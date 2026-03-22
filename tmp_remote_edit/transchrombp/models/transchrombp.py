"""Top-level TransChromBP model.

Architecture:
1) Convolutional feature extraction from one-hot sequence.
2) Optional Transformer encoder for long-range dependency modeling.
3) Optional ChromBPNet-style bias decomposition and additive fusion.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

import torch
import torch.nn.functional as F
from torch import Tensor, nn

from .bias_branch import ChromBPNetBiasBranch, ResidualDilatedConvBlock, center_crop_1d
from .genos_adapter import GenosCountProj, GenosGatedAdapter, GenosSummaryFiLM
from .transformer_encoder import SequenceTransformerEncoder


def _to_channel_first(seq_onehot: Tensor) -> Tensor:
    if seq_onehot.dim() != 3:
        raise ValueError(
            f"Expected sequence tensor with 3 dims [B, L, 4] or [B, 4, L], got {tuple(seq_onehot.shape)}"
        )
    if seq_onehot.size(1) == 4:
        return seq_onehot
    if seq_onehot.size(2) == 4:
        return seq_onehot.transpose(1, 2)
    raise ValueError(f"Expected a 4-channel axis in sequence tensor, got {tuple(seq_onehot.shape)}")


class ConvStem(nn.Module):
    """Convolutional feature extractor before transformer."""

    def __init__(
        self,
        input_channels: int,
        hidden_channels: int,
        stem_kernel_size: int,
        n_conv_layers: int,
        conv_kernel_size: int,
        dropout: float,
    ) -> None:
        super().__init__()
        if stem_kernel_size % 2 == 0 or conv_kernel_size % 2 == 0:
            raise ValueError("stem_kernel_size and conv_kernel_size must be odd")
        if n_conv_layers <= 0:
            raise ValueError(f"n_conv_layers must be positive, got {n_conv_layers}")

        layers = [
            nn.Conv1d(
                in_channels=input_channels,
                out_channels=hidden_channels,
                kernel_size=stem_kernel_size,
                padding=stem_kernel_size // 2,
            ),
            nn.GELU(),
            nn.Dropout(dropout),
        ]
        for _ in range(n_conv_layers - 1):
            layers.extend(
                [
                    nn.Conv1d(
                        in_channels=hidden_channels,
                        out_channels=hidden_channels,
                        kernel_size=conv_kernel_size,
                        padding=conv_kernel_size // 2,
                    ),
                    nn.GELU(),
                    nn.Dropout(dropout),
                ]
            )

        self.net = nn.Sequential(*layers)

    def forward(self, x: Tensor) -> Tensor:
        return self.net(x)


class LocalDilatedTower(nn.Module):
    """Optional BPNet-style local context tower before the transformer."""

    def __init__(
        self,
        channels: int,
        kernel_size: int,
        n_dil_layers: int,
        dilation_cycle_length: int,
        dropout: float,
    ) -> None:
        super().__init__()
        if n_dil_layers < 0:
            raise ValueError(f"n_dil_layers must be non-negative, got {n_dil_layers}")
        if dilation_cycle_length <= 0:
            raise ValueError(f"dilation_cycle_length must be positive, got {dilation_cycle_length}")

        if n_dil_layers == 0:
            self.net = nn.Identity()
            return

        blocks = []
        for i in range(n_dil_layers):
            dilation = 2 ** (i % dilation_cycle_length)
            blocks.append(
                ResidualDilatedConvBlock(
                    channels=channels,
                    kernel_size=kernel_size,
                    dilation=dilation,
                    dropout=dropout,
                )
            )
        self.net = nn.Sequential(*blocks)

    def forward(self, x: Tensor) -> Tensor:
        return self.net(x)


def _validate_heads_config(heads_cfg: Dict[str, Any]) -> None:
    count_head_type = str(heads_cfg.get("count_head", "linear")).lower()
    if count_head_type != "linear":
        raise ValueError(f"Unsupported heads.count_head={count_head_type!r}; expected 'linear'")


def _load_checkpoint(path: Path) -> Dict[str, Any]:
    try:
        return torch.load(path, map_location="cpu", weights_only=False)
    except TypeError:
        return torch.load(path, map_location="cpu")


def _inverse_softplus(value: float, eps: float = 1e-6) -> float:
    value = max(float(value), eps)
    return float(value + torch.log(-torch.expm1(torch.tensor(-value))).item())


@dataclass
class TransChromBPOutput:
    """Structured output container."""

    profile_logits_full: Tensor
    logcount_full: Tensor
    profile_logits_debiased: Tensor
    logcount_debiased: Tensor
    profile_bias: Tensor
    count_bias: Tensor


class TransChromBP(nn.Module):
    """TransChromBP model: conv -> optional transformer + optional bias decomposition."""

    def __init__(
        self,
        input_channels: int = 4,
        d_model: int = 256,
        n_heads: int = 8,
        n_layers: int = 6,
        ff_mult: int = 4,
        dropout: float = 0.1,
        output_len: int = 1000,
        stem_kernel_size: int = 21,
        conv_kernel_size: int = 7,
        n_conv_layers: int = 2,
        signal_dil_kernel_size: int = 3,
        signal_n_dil_layers: int = 0,
        signal_dilation_cycle_length: int = 8,
        use_transformer: bool = True,
        use_rope: bool = True,
        rope_theta: float = 10000.0,
        transformer_max_len: int = 32768,
        use_sdpa: bool = True,
        use_bias_decomposition: bool = True,
        bias_hidden_channels: int = 128,
        bias_kernel_size: int = 21,
        bias_n_dil_layers: int = 4,
        learnable_scales: bool = True,
        profile_fusion: str = "add",
        profile_scale_init: float = 1.0,
        count_scale_init: float = 1.0,
        positive_scales: bool = True,
        count_fusion: str = "logsumexp",
        bias_pretrained_path: str = "",
        freeze_bias_core: bool = False,
        profile_bias_stop_gradient: bool = False,
        bias_profile_pool_factor: int = 0,
    ) -> None:
        super().__init__()
        self.output_len = output_len
        self.use_transformer = bool(use_transformer)
        self.use_bias_decomposition = bool(use_bias_decomposition)
        self.learnable_scales = bool(learnable_scales)
        self.positive_scales = bool(positive_scales)
        self.profile_fusion = str(profile_fusion).lower()
        self.count_fusion = str(count_fusion).lower()
        if self.profile_fusion != "add":
            raise ValueError(f"Unsupported profile_fusion={profile_fusion!r}; expected 'add'")
        if self.count_fusion not in {"add", "logsumexp"}:
            raise ValueError(f"Unsupported count_fusion={count_fusion!r}; expected 'add' or 'logsumexp'")
        self.profile_bias_stop_gradient = bool(profile_bias_stop_gradient)

        self.conv_stem = ConvStem(
            input_channels=input_channels,
            hidden_channels=d_model,
            stem_kernel_size=stem_kernel_size,
            n_conv_layers=n_conv_layers,
            conv_kernel_size=conv_kernel_size,
            dropout=dropout,
        )
        self.local_tower = LocalDilatedTower(
            channels=d_model,
            kernel_size=signal_dil_kernel_size,
            n_dil_layers=signal_n_dil_layers,
            dilation_cycle_length=signal_dilation_cycle_length,
            dropout=dropout,
        )

        if self.use_transformer:
            self.transformer = SequenceTransformerEncoder(
                d_model=d_model,
                n_heads=n_heads,
                n_layers=n_layers,
                ff_mult=ff_mult,
                dropout=dropout,
                activation="gelu",
                use_positional_encoding=use_rope,
                max_len=transformer_max_len,
                rope_theta=rope_theta,
                use_sdpa=use_sdpa,
            )
        else:
            self.transformer = nn.Identity()

        self.profile_signal_head = nn.Linear(d_model, 1)
        self.count_signal_head = nn.Sequential(
            nn.LayerNorm(d_model),
            nn.Linear(d_model, d_model // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_model // 2, 1),
        )

        if self.use_bias_decomposition:
            self.bias_branch: Optional[ChromBPNetBiasBranch] = ChromBPNetBiasBranch(
                input_channels=input_channels,
                hidden_channels=bias_hidden_channels,
                kernel_size=bias_kernel_size,
                n_dil_layers=bias_n_dil_layers,
                output_len=output_len,
                dropout=dropout,
                profile_pool_factor=bias_profile_pool_factor,
            )
            if bias_pretrained_path:
                self.load_bias_branch_weights(bias_pretrained_path, freeze=freeze_bias_core)
            elif freeze_bias_core:
                raise ValueError("freeze_bias_core=True requires bias_branch.pretrained_path to be set")
        else:
            self.bias_branch = None
            if bias_pretrained_path:
                raise ValueError("bias_branch.pretrained_path cannot be used when bias decomposition is disabled")

        # ── Optional Genos gated fusion adapter ──
        self.genos_adapter: Optional[GenosGatedAdapter] = None
        # ── Optional cached Genos summary adapters ──
        self.genos_film: Optional[GenosSummaryFiLM] = None
        self.genos_count_proj: Optional[GenosCountProj] = None

        profile_scale_value = float(profile_scale_init)
        count_scale_value = float(count_scale_init)
        if self.positive_scales and self.learnable_scales:
            profile_scale = torch.tensor(_inverse_softplus(profile_scale_value), dtype=torch.float32)
            count_scale = torch.tensor(_inverse_softplus(count_scale_value), dtype=torch.float32)
        else:
            profile_scale = torch.tensor(profile_scale_value, dtype=torch.float32)
            count_scale = torch.tensor(count_scale_value, dtype=torch.float32)
        if learnable_scales:
            self.profile_scale = nn.Parameter(profile_scale)
            self.count_scale = nn.Parameter(count_scale)
        else:
            self.register_buffer("profile_scale", profile_scale)
            self.register_buffer("count_scale", count_scale)

    @staticmethod
    def _extract_bias_state_dict(payload: Dict[str, Any]) -> Dict[str, Tensor]:
        state: Dict[str, Any] = payload
        if "bias_branch_state" in state and isinstance(state["bias_branch_state"], dict):
            state = state["bias_branch_state"]
        elif "model_state" in state and isinstance(state["model_state"], dict):
            state = state["model_state"]

        if any(key.startswith("bias_branch.") for key in state):
            return {
                key.split("bias_branch.", 1)[1]: value
                for key, value in state.items()
                if key.startswith("bias_branch.")
            }
        return state

    def load_bias_branch_weights(self, checkpoint_path: str, freeze: bool = False, strict: bool = True) -> None:
        """Load bias-branch weights from either a raw state dict or a saved training checkpoint."""
        if self.bias_branch is None:
            raise RuntimeError("Cannot load bias weights when bias decomposition is disabled")

        path = Path(checkpoint_path)
        if not path.exists():
            raise FileNotFoundError(f"Bias checkpoint not found: {checkpoint_path}")

        payload = _load_checkpoint(path)
        if not isinstance(payload, dict):
            raise ValueError(f"Unsupported bias checkpoint format at {checkpoint_path}")

        bias_state = self._extract_bias_state_dict(payload)
        missing, unexpected = self.bias_branch.load_state_dict(bias_state, strict=strict)
        if strict and (missing or unexpected):
            raise RuntimeError(
                "Failed to load bias branch strictly: "
                f"missing={missing} unexpected={unexpected}"
            )
        if freeze:
            self.bias_branch.freeze_core()

    def _resolve_scale(self, scale: Tensor, reference: Tensor) -> Tensor:
        scale = scale.to(device=reference.device, dtype=reference.dtype)
        if self.positive_scales:
            if self.learnable_scales:
                scale = F.softplus(scale)
            else:
                scale = torch.clamp_min(scale, 0.0)
        return scale

    def forward(
        self,
        seq_onehot: Tensor,
        bias_profile: Optional[Tensor] = None,
        bias_count: Optional[Tensor] = None,
        genos_feat: Optional[Tensor] = None,
        genos_summary: Optional[Tensor] = None,
    ) -> TransChromBPOutput:
        """Forward pass.

        Args:
            seq_onehot: [B, L, 4] or [B, 4, L]
            bias_profile: optional external profile bias logits [B, output_len]
            bias_count: optional external count bias [B, 1]
            genos_feat: optional Genos hidden-state features [B, L, genos_dim]
            genos_summary: optional cached Genos global summary [B, genos_dim]
        """
        x = _to_channel_first(seq_onehot)

        conv_feat = self.conv_stem(x)  # [B, D, L]
        local_feat = self.local_tower(conv_feat)  # [B, D, L]
        tokens = local_feat.transpose(1, 2)  # [B, L, D] batch-first

        if genos_feat is not None and self.genos_adapter is not None:
            tokens = self.genos_adapter(tokens, genos_feat)

        encoded = self.transformer(tokens)

        # P1: late FiLM modulation after transformer
        if genos_summary is not None and self.genos_film is not None:
            encoded = self.genos_film(encoded, genos_summary)

        profile_signal = self.profile_signal_head(encoded).squeeze(-1)  # [B, L]
        profile_signal = center_crop_1d(profile_signal, self.output_len)

        pooled = encoded.mean(dim=1)

        # P2: inject genos into count head hidden layer
        if genos_summary is not None and self.genos_count_proj is not None:
            # Manually run count head layers to inject at hidden layer
            count_layers = list(self.count_signal_head)
            # layers: [LN, Linear(d_model, d_model//2), GELU, Dropout, Linear(d_model//2, 1)]
            h = count_layers[0](pooled)   # LayerNorm
            h = count_layers[1](h)        # Linear -> d_model//2
            h = count_layers[2](h)        # GELU
            h = self.genos_count_proj(h, genos_summary)  # additive injection
            h = count_layers[3](h)        # Dropout
            count_signal = count_layers[4](h)  # Linear -> 1
        else:
            count_signal = self.count_signal_head(pooled)  # [B, 1]

        if (bias_profile is None) != (bias_count is None):
            raise ValueError("bias_profile and bias_count must be both provided or both None")

        if bias_profile is not None and bias_count is not None:
            # Case 1: Use provided external bias
            profile_bias, count_bias = bias_profile, bias_count
        elif self.use_bias_decomposition:
            # Case 2: Compute bias using internal branch
            if self.bias_branch is None:
                raise RuntimeError("Bias decomposition is enabled but bias_branch is not initialized")
            profile_bias, count_bias = self.bias_branch(seq_onehot)
        else:
            # Case 3: No bias correction
            profile_bias = torch.zeros_like(profile_signal)
            count_bias = torch.zeros_like(count_signal)

        if profile_bias.shape != profile_signal.shape:
            raise ValueError(
                f"profile_bias shape {tuple(profile_bias.shape)} does not match profile_signal shape {tuple(profile_signal.shape)}"
            )
        if count_bias.shape != count_signal.shape:
            raise ValueError(
                f"count_bias shape {tuple(count_bias.shape)} does not match count_signal shape {tuple(count_signal.shape)}"
            )

        profile_scale = self._resolve_scale(self.profile_scale, profile_signal)
        count_scale = self._resolve_scale(self.count_scale, count_signal)

        if self.profile_bias_stop_gradient:
            profile_full = profile_signal + profile_scale * profile_bias.detach()
        else:
            profile_full = profile_signal + profile_scale * profile_bias
        if self.count_fusion == "logsumexp":
            log_bias_scale = torch.where(
                count_scale > 0,
                torch.log(count_scale),
                torch.full_like(count_scale, -torch.inf),
            )
            count_full = torch.logsumexp(
                torch.stack([count_signal, count_bias + log_bias_scale.expand_as(count_bias)], dim=0),
                dim=0,
            )
        else:
            count_full = count_signal + count_scale * count_bias

        return TransChromBPOutput(
            profile_logits_full=profile_full,
            logcount_full=count_full,
            profile_logits_debiased=profile_signal,
            logcount_debiased=count_signal,
            profile_bias=profile_bias,
            count_bias=count_bias,
        )

    def freeze_bias_core(self) -> None:
        """Freeze bias branch parameters (no-op if bias branch disabled)."""
        if self.bias_branch is None:
            return
        self.bias_branch.freeze_core()

    def unfreeze_bias_core(self) -> None:
        """Unfreeze bias branch parameters (no-op if bias branch disabled)."""
        if self.bias_branch is None:
            return
        self.bias_branch.unfreeze_core()


def build_transchrombp_from_config(config: Dict[str, Any]) -> TransChromBP:
    """Build model from nested config dict loaded from YAML/JSON.

    Expected structure:
      sequence_encoder: enabled, d_model, n_heads, n_layers, ff_mult, dropout,
                        use_rope, rope_theta, max_len, use_sdpa
      conv_stem: stem_kernel_size, conv_kernel_size, n_conv_layers
      local_tower: kernel_size, n_dil_layers, dilation_cycle_length
      bias_branch: enabled, hidden_channels, kernel_size, n_dil_layers, freeze_bias_core
      fusion: profile_scale_init, count_scale_init, learnable_scales
      heads: profile_output_len
    """
    seq_cfg = config.get("sequence_encoder", {})
    conv_cfg = config.get("conv_stem", {})
    local_cfg = config.get("local_tower", {})
    bias_cfg = config.get("bias_branch", {})
    fusion_cfg = config.get("fusion", {})
    heads_cfg = config.get("heads", {})
    _validate_heads_config(heads_cfg)

    genos_cfg = config.get("genos_branch", {})

    model = TransChromBP(
        input_channels=int(config.get("input_channels", 4)),
        d_model=int(seq_cfg.get("d_model", 256)),
        n_heads=int(seq_cfg.get("n_heads", 8)),
        n_layers=int(seq_cfg.get("n_layers", 6)),
        ff_mult=int(seq_cfg.get("ff_mult", 4)),
        dropout=float(seq_cfg.get("dropout", 0.1)),
        use_transformer=bool(seq_cfg.get("enabled", True)),
        use_rope=bool(seq_cfg.get("use_rope", True)),
        rope_theta=float(seq_cfg.get("rope_theta", 10000.0)),
        transformer_max_len=int(seq_cfg.get("max_len", 32768)),
        use_sdpa=bool(seq_cfg.get("use_sdpa", True)),
        use_bias_decomposition=bool(bias_cfg.get("enabled", True)),
        output_len=int(heads_cfg.get("profile_output_len", 1000)),
        stem_kernel_size=int(conv_cfg.get("stem_kernel_size", 21)),
        conv_kernel_size=int(conv_cfg.get("conv_kernel_size", 7)),
        n_conv_layers=int(conv_cfg.get("n_conv_layers", 2)),
        signal_dil_kernel_size=int(local_cfg.get("kernel_size", 3)),
        signal_n_dil_layers=int(local_cfg.get("n_dil_layers", 0)),
        signal_dilation_cycle_length=int(local_cfg.get("dilation_cycle_length", 8)),
        bias_hidden_channels=int(bias_cfg.get("hidden_channels", 128)),
        bias_kernel_size=int(bias_cfg.get("kernel_size", 21)),
        bias_n_dil_layers=int(bias_cfg.get("n_dil_layers", 4)),
        learnable_scales=bool(fusion_cfg.get("learnable_scales", True)),
        profile_fusion=str(fusion_cfg.get("profile_fusion", "add")),
        profile_scale_init=float(fusion_cfg.get("profile_scale_init", 1.0)),
        count_scale_init=float(fusion_cfg.get("count_scale_init", 1.0)),
        positive_scales=bool(fusion_cfg.get("positive_scales", True)),
        count_fusion=str(fusion_cfg.get("count_fusion", "logsumexp")),
        bias_pretrained_path=str(bias_cfg.get("pretrained_path", "")),
        freeze_bias_core=bool(bias_cfg.get("freeze_bias_core", False)),
        profile_bias_stop_gradient=bool(fusion_cfg.get("profile_bias_stop_gradient", False)),
        bias_profile_pool_factor=int(bias_cfg.get("profile_pool_factor", 0)),
    )

    # ── Attach optional Genos gated adapter (online mode) ──
    if genos_cfg.get("enabled", False):
        d_model = int(seq_cfg.get("d_model", 256))
        fusion_mode = str(genos_cfg.get("fusion_mode", "gate"))

        if fusion_mode == "gate":
            model.genos_adapter = GenosGatedAdapter(
                d_model=d_model,
                genos_hidden_size=int(genos_cfg.get("hidden_size", 1024)),
                gate_bias_init=float(genos_cfg.get("gate_bias_init", -2.0)),
                pool_mode=str(genos_cfg.get("pool_mode", "none")),
            )
        else:
            raise ValueError(
                f"Unsupported genos_branch.fusion_mode={fusion_mode!r}; expected 'gate'"
            )

    # ── Attach optional cached Genos summary adapters ──
    genos_cached_cfg = config.get("genos_cached", {})
    if genos_cached_cfg.get("enabled", False):
        d_model = int(seq_cfg.get("d_model", 256))
        genos_dim = int(genos_cached_cfg.get("genos_dim", 1024))
        cached_fusion = str(genos_cached_cfg.get("fusion_mode", "late_film"))

        if cached_fusion in {"late_film", "late_film_and_count"}:
            model.genos_film = GenosSummaryFiLM(
                d_model=d_model,
                genos_dim=genos_dim,
            )
        if cached_fusion in {"count_only", "late_film_and_count"}:
            model.genos_count_proj = GenosCountProj(
                d_count_hidden=d_model // 2,
                genos_dim=genos_dim,
            )
        if cached_fusion not in {"late_film", "count_only", "late_film_and_count"}:
            raise ValueError(
                "Unsupported genos_cached.fusion_mode="
                f"{cached_fusion!r}; expected one of "
                "'late_film', 'count_only', 'late_film_and_count'"
            )

    return model


def build_bias_branch_from_config(config: Dict[str, Any]) -> ChromBPNetBiasBranch:
    """Build just the bias branch so it can be pretrained separately from the main model."""
    seq_cfg = config.get("sequence_encoder", {})
    bias_cfg = config.get("bias_branch", {})
    heads_cfg = config.get("heads", {})
    _validate_heads_config(heads_cfg)
    return ChromBPNetBiasBranch(
        input_channels=int(config.get("input_channels", 4)),
        hidden_channels=int(bias_cfg.get("hidden_channels", 128)),
        kernel_size=int(bias_cfg.get("kernel_size", 21)),
        n_dil_layers=int(bias_cfg.get("n_dil_layers", 4)),
        output_len=int(heads_cfg.get("profile_output_len", 1000)),
        dropout=float(seq_cfg.get("dropout", 0.1)),
        profile_pool_factor=int(bias_cfg.get("profile_pool_factor", 0)),
    )
