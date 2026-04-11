"""Hierarchical encoder-decoder TransChromBP variant."""

from __future__ import annotations

import math
from typing import Any, Dict, Sequence

import torch
import torch.nn.functional as F
from torch import Tensor, nn

from .bias_branch import ChromBPNetBiasBranch, center_crop_1d
from .transformer_encoder import SequenceTransformerEncoder
from .transchrombp import (
    TransChromBPOutput,
    _inverse_softplus,
    _load_checkpoint,
    _to_channel_first,
    _validate_heads_config,
)


class DownsampleStage(nn.Module):
    """Stride-2 conv stage with a small residual refinement block."""

    def __init__(self, in_channels: int, out_channels: int, kernel_size: int, dropout: float) -> None:
        super().__init__()
        if kernel_size % 2 == 0:
            raise ValueError("kernel_size must be odd for same-length padding")

        self.downsample = nn.Sequential(
            nn.Conv1d(
                in_channels=in_channels,
                out_channels=out_channels,
                kernel_size=kernel_size,
                stride=2,
                padding=kernel_size // 2,
            ),
            nn.GELU(),
        )
        self.refine = nn.Sequential(
            nn.Conv1d(out_channels, out_channels, kernel_size=3, padding=1),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Conv1d(out_channels, out_channels, kernel_size=3, padding=1),
        )
        self.norm = nn.BatchNorm1d(out_channels)
        self.act = nn.GELU()

    def forward(self, x: Tensor) -> Tensor:
        x = self.downsample(x)
        residual = x
        x = self.refine(x)
        return self.act(self.norm(x + residual))


class UpsampleStage(nn.Module):
    """Upsample to a skip resolution, then fuse and refine."""

    def __init__(self, in_channels: int, skip_channels: int, out_channels: int, dropout: float) -> None:
        super().__init__()
        self.input_proj = nn.Conv1d(in_channels, out_channels, kernel_size=1)
        self.skip_proj = nn.Conv1d(skip_channels, out_channels, kernel_size=1)
        self.refine = nn.Sequential(
            nn.GELU(),
            nn.Conv1d(out_channels, out_channels, kernel_size=3, padding=1),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Conv1d(out_channels, out_channels, kernel_size=3, padding=1),
        )
        self.norm = nn.BatchNorm1d(out_channels)
        self.act = nn.GELU()

    def forward(self, x: Tensor, skip: Tensor) -> Tensor:
        x = F.interpolate(x, size=skip.size(-1), mode="linear", align_corners=False)
        fused = self.input_proj(x) + self.skip_proj(skip)
        refined = self.refine(fused)
        return self.act(self.norm(refined + fused))


class HierarchicalTransChromBP(nn.Module):
    """Minimal hierarchical encoder-decoder for corrected-B E2 experiments."""

    def __init__(
        self,
        input_channels: int = 4,
        d_model: int = 256,
        stem_channels: int = 128,
        downsample_channels: Sequence[int] = (128, 256),
        n_heads: int = 8,
        n_layers: int = 6,
        ff_mult: int = 4,
        dropout: float = 0.1,
        output_len: int = 1000,
        stem_kernel_size: int = 21,
        stage_kernel_size: int = 7,
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
        count_pool_mode: str = "full",
    ) -> None:
        super().__init__()
        if len(tuple(downsample_channels)) != 2:
            raise ValueError(f"downsample_channels must contain exactly 2 stages, got {downsample_channels}")
        if stem_kernel_size % 2 == 0 or stage_kernel_size % 2 == 0:
            raise ValueError("stem_kernel_size and stage_kernel_size must be odd")

        stage1_channels, stage2_channels = [int(channels) for channels in downsample_channels]
        if stage2_channels != d_model:
            raise ValueError(
                f"Final downsample stage must match d_model so transformer dims align, got {stage2_channels} != {d_model}"
            )

        self.output_len = int(output_len)
        self.total_stride = 4
        self.use_transformer = bool(use_transformer)
        self.use_bias_decomposition = bool(use_bias_decomposition)
        self.learnable_scales = bool(learnable_scales)
        self.positive_scales = bool(positive_scales)
        self.profile_fusion = str(profile_fusion).lower()
        self.count_fusion = str(count_fusion).lower()
        self.count_pool_mode = str(count_pool_mode).lower()
        self.profile_bias_stop_gradient = bool(profile_bias_stop_gradient)
        if self.profile_fusion != "add":
            raise ValueError(f"Unsupported profile_fusion={profile_fusion!r}; expected 'add'")
        if self.count_fusion not in {"add", "logsumexp"}:
            raise ValueError(f"Unsupported count_fusion={count_fusion!r}; expected 'add' or 'logsumexp'")
        if self.count_pool_mode not in {"full", "center", "attention"}:
            raise ValueError(
                f"Unsupported count_pool_mode={count_pool_mode!r}; expected 'full', 'center', or 'attention'"
            )

        self.stem = nn.Sequential(
            nn.Conv1d(
                in_channels=input_channels,
                out_channels=stem_channels,
                kernel_size=stem_kernel_size,
                padding=stem_kernel_size // 2,
            ),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Conv1d(stem_channels, stem_channels, kernel_size=stage_kernel_size, padding=stage_kernel_size // 2),
            nn.GELU(),
            nn.Dropout(dropout),
        )
        self.downsample1 = DownsampleStage(
            in_channels=stem_channels,
            out_channels=stage1_channels,
            kernel_size=stage_kernel_size,
            dropout=dropout,
        )
        self.downsample2 = DownsampleStage(
            in_channels=stage1_channels,
            out_channels=stage2_channels,
            kernel_size=stage_kernel_size,
            dropout=dropout,
        )

        if self.use_transformer:
            compressed_max_len = max(1, int(math.ceil(transformer_max_len / self.total_stride)))
            self.transformer = SequenceTransformerEncoder(
                d_model=d_model,
                n_heads=n_heads,
                n_layers=n_layers,
                ff_mult=ff_mult,
                dropout=dropout,
                activation="gelu",
                use_positional_encoding=use_rope,
                max_len=compressed_max_len,
                rope_theta=rope_theta,
                use_sdpa=use_sdpa,
            )
        else:
            self.transformer = nn.Identity()

        self.upsample1 = UpsampleStage(
            in_channels=d_model,
            skip_channels=stage1_channels,
            out_channels=stage1_channels,
            dropout=dropout,
        )
        self.upsample2 = UpsampleStage(
            in_channels=stage1_channels,
            skip_channels=stem_channels,
            out_channels=stem_channels,
            dropout=dropout,
        )
        self.profile_head = nn.Sequential(
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Conv1d(stem_channels, 1, kernel_size=1),
        )

        self.count_signal_head = nn.Sequential(
            nn.LayerNorm(d_model),
            nn.Linear(d_model, d_model // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_model // 2, 1),
        )
        if self.count_pool_mode == "attention":
            self.count_pool_proj: nn.Module | None = nn.Linear(d_model, 1)
        else:
            self.count_pool_proj = None

        if self.use_bias_decomposition:
            self.bias_branch: ChromBPNetBiasBranch | None = ChromBPNetBiasBranch(
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

        profile_scale_value = float(profile_scale_init)
        count_scale_value = float(count_scale_init)
        if self.positive_scales and self.learnable_scales:
            profile_scale = torch.tensor(_inverse_softplus(profile_scale_value), dtype=torch.float32)
            count_scale = torch.tensor(_inverse_softplus(count_scale_value), dtype=torch.float32)
        else:
            profile_scale = torch.tensor(profile_scale_value, dtype=torch.float32)
            count_scale = torch.tensor(count_scale_value, dtype=torch.float32)
        if self.learnable_scales:
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
        if self.bias_branch is None:
            raise RuntimeError("Cannot load bias weights when bias decomposition is disabled")

        payload = _load_checkpoint(checkpoint_path)
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

    def _pool_for_count(self, encoded: Tensor) -> Tensor:
        if encoded.dim() != 3:
            raise ValueError(f"Expected encoded shape [B, L, D], got {tuple(encoded.shape)}")

        if self.count_pool_mode == "full":
            return encoded.mean(dim=1)

        if self.count_pool_mode == "center":
            seq_len = encoded.size(1)
            center_tokens = max(1, int(math.ceil(self.output_len / self.total_stride)))
            if center_tokens > seq_len:
                raise ValueError(
                    f"count center pool requires {center_tokens} compressed tokens but only {seq_len} are available"
                )
            start = (seq_len - center_tokens) // 2
            end = start + center_tokens
            return encoded[:, start:end, :].mean(dim=1)

        if self.count_pool_mode == "attention":
            if self.count_pool_proj is None:
                raise RuntimeError("count_pool_mode='attention' requires count_pool_proj to be initialized")
            attn_logits = self.count_pool_proj(encoded).squeeze(-1)
            attn_weights = torch.softmax(attn_logits, dim=1)
            return torch.sum(encoded * attn_weights.unsqueeze(-1), dim=1)

        raise RuntimeError(f"Unhandled count_pool_mode={self.count_pool_mode!r}")

    def forward(
        self,
        seq_onehot: Tensor,
        bias_profile: Tensor | None = None,
        bias_count: Tensor | None = None,
        **extra_inputs: Tensor,
    ) -> TransChromBPOutput:
        if extra_inputs:
            unsupported_keys = ", ".join(sorted(extra_inputs))
            raise ValueError(
                "HierarchicalTransChromBP does not support extra forward inputs: "
                f"{unsupported_keys}"
            )

        x = _to_channel_first(seq_onehot)

        skip0 = self.stem(x)
        skip1 = self.downsample1(skip0)
        bottleneck = self.downsample2(skip1)
        encoded = self.transformer(bottleneck.transpose(1, 2))

        decoded = self.upsample1(encoded.transpose(1, 2), skip1)
        decoded = self.upsample2(decoded, skip0)
        profile_signal = self.profile_head(decoded).squeeze(1)
        profile_signal = center_crop_1d(profile_signal, self.output_len)

        pooled = self._pool_for_count(encoded)
        count_signal = self.count_signal_head(pooled)

        if (bias_profile is None) != (bias_count is None):
            raise ValueError("bias_profile and bias_count must be both provided or both None")

        if bias_profile is not None and bias_count is not None:
            profile_bias, count_bias = bias_profile, bias_count
        elif self.use_bias_decomposition:
            if self.bias_branch is None:
                raise RuntimeError("Bias decomposition is enabled but bias_branch is not initialized")
            profile_bias, count_bias = self.bias_branch(seq_onehot)
        else:
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


def _raise_on_unsupported_hierarchical_config(config: Dict[str, Any]) -> None:
    unsupported_sections: list[str] = []
    for section_name in ("genos_branch", "caduceus_branch", "genos_cached", "foundation_model"):
        section_cfg = config.get(section_name, {})
        if isinstance(section_cfg, dict) and bool(section_cfg.get("enabled", False)):
            unsupported_sections.append(section_name)

    profile_decoder_cfg = config.get("profile_decoder", {})
    if isinstance(profile_decoder_cfg, dict) and profile_decoder_cfg:
        unsupported_sections.append("profile_decoder")

    if unsupported_sections:
        raise ValueError(
            "HierarchicalTransChromBP does not support optional config sections: "
            + ", ".join(sorted(set(unsupported_sections)))
        )


def build_hierarchical_transchrombp_from_config(config: Dict[str, Any]) -> HierarchicalTransChromBP:
    """Build the hierarchical encoder-decoder variant from a model config."""

    _raise_on_unsupported_hierarchical_config(config)

    architecture_cfg = config.get("architecture", {})
    seq_cfg = config.get("sequence_encoder", {})
    conv_cfg = config.get("conv_stem", {})
    bias_cfg = config.get("bias_branch", {})
    fusion_cfg = config.get("fusion", {})
    heads_cfg = config.get("heads", {})
    _validate_heads_config(heads_cfg)

    d_model = int(seq_cfg.get("d_model", 256))
    stem_channels = int(architecture_cfg.get("stem_channels", max(64, d_model // 2)))
    downsample_channels = architecture_cfg.get("downsample_channels", [stem_channels, d_model])

    return HierarchicalTransChromBP(
        input_channels=int(config.get("input_channels", 4)),
        d_model=d_model,
        stem_channels=stem_channels,
        downsample_channels=[int(channels) for channels in downsample_channels],
        n_heads=int(seq_cfg.get("n_heads", 8)),
        n_layers=int(seq_cfg.get("n_layers", 6)),
        ff_mult=int(seq_cfg.get("ff_mult", 4)),
        dropout=float(seq_cfg.get("dropout", 0.1)),
        output_len=int(heads_cfg.get("profile_output_len", 1000)),
        stem_kernel_size=int(conv_cfg.get("stem_kernel_size", 21)),
        stage_kernel_size=int(conv_cfg.get("conv_kernel_size", 7)),
        use_transformer=bool(seq_cfg.get("enabled", True)),
        use_rope=bool(seq_cfg.get("use_rope", True)),
        rope_theta=float(seq_cfg.get("rope_theta", 10000.0)),
        transformer_max_len=int(seq_cfg.get("max_len", 32768)),
        use_sdpa=bool(seq_cfg.get("use_sdpa", True)),
        use_bias_decomposition=bool(bias_cfg.get("enabled", True)),
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
        count_pool_mode=str(heads_cfg.get("count_pool_mode", "full")),
    )
