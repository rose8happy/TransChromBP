# ChatGPT 咨询包 01：当前最佳模型代码与配置

## 这份文件是什么
这份文件的目标是让你在不访问本地仓库的情况下，也能理解我们当前默认最好模型的实现细节、训练语义和关键开关。

## 先给结论
- 当前 paper-facing 默认模型：`corrected B = Transformer + center pool + sg=true + deb2`
- 代码主体目前在：`vendor/transchrombp/transchrombp/`
- 最关键的设计不是“单纯换成 Transformer”，而是：
  1. `conv stem + local dilated tower + transformer encoder`
  2. 独立 `bias branch`
  3. `full` 与 `debiased` 双输出
  4. count head 使用 `center pool`
  5. 训练时显式加 `debiased_profile_weight = 2.0`

## 需要你特别注意的 5 个行为语义
1. `profile_logits_full = profile_signal + scale * bias`，而 `profile_logits_debiased = profile_signal`；我们会同时评估这两者。
2. `profile_bias_stop_gradient=true` 时，profile 融合路径会对 bias 分支做 `detach()`。
3. `count_fusion=logsumexp`，不是简单相加。
4. 当前默认最好 readout 不是 attention pool，而是 `count_pool_mode: center`。
5. clean-matrix 里真正抑制 bias reliance 的主因更像 `debiased_profile_weight=2.0`，而不是 stop-gradient 本身。

## 当前默认最好模型的最小重建清单
- 模型配置：`vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool.yaml`
- 代表性训练配置：`vendor/transchrombp/transchrombp/configs/train/train_v2fix_cpool_s1234_6000_single.yaml`
- 主模型代码：`vendor/transchrombp/transchrombp/models/transchrombp.py`
- bias branch：`vendor/transchrombp/transchrombp/models/bias_branch.py`
- transformer encoder：`vendor/transchrombp/transchrombp/models/transformer_encoder.py`

## 模型配置（paper-facing corrected B）
```yaml
# Teacher v2 scaffold with count head center pooling.
# This is the corrected B baseline: keep local_tower, bias profile pooling,
# and stop-gradient semantics from teacher_v2, and only switch count pooling
# from full-sequence mean to center-window mean.

model_name: transchrombp_teacher_v2_center_pool

sequence_encoder:
  type: transformer
  enabled: true
  d_model: 256
  n_heads: 8
  n_layers: 6
  ff_mult: 4
  dropout: 0.1
  use_rope: true
  rope_theta: 10000.0
  max_len: 2114
  use_sdpa: true

conv_stem:
  stem_kernel_size: 21
  conv_kernel_size: 7
  n_conv_layers: 2

local_tower:
  kernel_size: 3
  n_dil_layers: 8
  dilation_cycle_length: 8

bias_branch:
  type: chrombpnet_bias
  enabled: true
  hidden_channels: 128
  kernel_size: 21
  n_dil_layers: 4
  pretrained_path: ""
  freeze_bias_core: false
  profile_pool_factor: 32

fusion:
  profile_fusion: add
  count_fusion: logsumexp
  profile_scale_init: 1.0
  count_scale_init: 1.0
  learnable_scales: true
  positive_scales: true
  profile_bias_stop_gradient: true

heads:
  profile_output_len: 1000
  count_head: linear
  count_pool_mode: center
```

## 代表性训练配置（B=center pool, deb2）
```yaml
# V2fix experiment B on 6000 single GPU.
# Keep seed/global-batch semantics close to the 6002 F_s1234 run:
# bs/gpu=16 + grad_accum=2 => effective batch ~= 32.
seed: 1234
max_epochs: 40

trainer:
  backend: nccl
  precision: bf16
  grad_accum_steps: 2
  clip_grad_norm: 1.0
  log_every_steps: 20
  validate_every_epochs: 1
  checkpoint_every_epochs: 1
  find_unused_parameters: false
  compile: false
  strict_max_len_check: true
  best_metric: peak.profile_target_jsd_full_mean
  best_metric_mode: min
  early_stop_patience: 10
  early_stop_min_delta: 0.0

optimizer:
  name: adamw
  learning_rate: 5.0e-4
  weight_decay: 0.01
  betas: [0.9, 0.95]

loss:
  profile_weight: 1.0
  count_weight: 0.1
  debiased_profile_weight: 2.0
  debiased_count_weight: 0.0

schedule:
  name: cosine
  warmup_steps: 0
  warmup_ratio: 0.06
  min_lr_ratio: 0.1

data:
  source: chrombpnet_bigwig
  config_path: configs/data/data_tutorial_canonical_v1.yaml
  max_seq_len: 2114
  input_len: 2114
  output_len: 1000
  supervised_bp: 1000
  profile_bin_size: 1
  batch_size_per_gpu: 16
  num_workers: 2
  pin_memory: true
  prefetch_factor: 2
  persistent_workers: true
  region_source: both
  train_region_source: both
  val_region_source: both
  nonpeak_ratio: 0.1
  max_jitter: 500
  peak_max_jitter: 500
  nonpeak_max_jitter: 0
  train_revcomp: true
  val_revcomp: false
  revcomp_prob: 0.5
  track_total_count_target: 0.0

logging:
  output_dir: /data1/zhoujiazhen/bylw_atac/TransChromBP/outputs
  run_name: v2fix_cpool_s1234_6000_single
```

## 如何理解这份配置
- `sequence_encoder.enabled=true`：启用 Transformer backbone。
- `local_tower.n_dil_layers=8`：保留 BPNet 风格的局部扩张卷积塔。
- `bias_branch.profile_pool_factor=32`：把 bias profile 压成低频，避免 bias branch 直接表示细粒度 motif 形状。
- `fusion.profile_bias_stop_gradient=true`：对 profile bias 融合路径阻断梯度。
- `loss.debiased_profile_weight=2.0`：显式要求 signal branch 自己也把 profile 学好，这是当前最关键的“bias-safe”设计。
- `heads.count_pool_mode=center`：count 只看中间 1000bp 的 encoded token 平均，这个设计比 attention pool 更稳。

## 训练/诊断语义摘要
- profile loss：multinomial NLL
- count loss：`MSE(log1p(total_count))`
- total loss：`profile + 0.1 * count + 2.0 * debiased_profile`
- 关键诊断：
  - `effective_profile_scale`
  - `profile_bias_rms_over_signal_rms`
  - `profile_full_debiased_jsd`
  - `count_full_debiased_abs`
- best checkpoint 默认按 `peak.profile_target_jsd_full_mean` 选

## 主模型源码：transchrombp.py
```python
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
    count_pool_mode = str(heads_cfg.get("count_pool_mode", "full")).lower()
    if count_pool_mode not in {"full", "center", "attention"}:
        raise ValueError(
            f"Unsupported heads.count_pool_mode={count_pool_mode!r}; "
            "expected 'full', 'center', or 'attention'"
        )


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
        count_pool_mode: str = "full",
    ) -> None:
        super().__init__()
        self.output_len = output_len
        self.use_transformer = bool(use_transformer)
        self.use_bias_decomposition = bool(use_bias_decomposition)
        self.learnable_scales = bool(learnable_scales)
        self.positive_scales = bool(positive_scales)
        self.profile_fusion = str(profile_fusion).lower()
        self.count_fusion = str(count_fusion).lower()
        self.count_pool_mode = str(count_pool_mode).lower()
        if self.profile_fusion != "add":
            raise ValueError(f"Unsupported profile_fusion={profile_fusion!r}; expected 'add'")
        if self.count_fusion not in {"add", "logsumexp"}:
            raise ValueError(f"Unsupported count_fusion={count_fusion!r}; expected 'add' or 'logsumexp'")
        if self.count_pool_mode not in {"full", "center", "attention"}:
            raise ValueError(
                f"Unsupported count_pool_mode={count_pool_mode!r}; "
                "expected 'full', 'center', or 'attention'"
            )
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
        self.count_pool_proj: Optional[nn.Linear]
        if self.count_pool_mode == "attention":
            self.count_pool_proj = nn.Linear(d_model, 1)
        else:
            self.count_pool_proj = None

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

    def _pool_for_count(self, encoded: Tensor) -> Tensor:
        if encoded.dim() != 3:
            raise ValueError(f"Expected encoded shape [B, L, D], got {tuple(encoded.shape)}")

        if self.count_pool_mode == "full":
            return encoded.mean(dim=1)

        if self.count_pool_mode == "center":
            seq_len = encoded.size(1)
            if self.output_len > seq_len:
                raise ValueError(
                    f"count center pool requires output_len <= seq_len, got {self.output_len} > {seq_len}"
                )
            start = (seq_len - self.output_len) // 2
            end = start + self.output_len
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

        pooled = self._pool_for_count(encoded)

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
        count_pool_mode=str(heads_cfg.get("count_pool_mode", "full")),
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
```

## Bias Branch 源码：bias_branch.py
```python
"""ChromBPNet-style bias branch for TransChromBP."""

from __future__ import annotations

from typing import Tuple

from torch import Tensor, nn


def center_crop_1d(x: Tensor, target_len: int) -> Tensor:
    """Center-crop last dimension to target length."""
    if x.dim() < 1:
        raise ValueError("center_crop_1d expects tensor with at least 1 dimension")
    if target_len <= 0:
        raise ValueError(f"target_len must be positive, got {target_len}")

    seq_len = x.size(-1)
    if seq_len < target_len:
        raise ValueError(f"Cannot crop from length {seq_len} to larger target {target_len}")
    if seq_len == target_len:
        return x

    start = (seq_len - target_len) // 2
    end = start + target_len
    return x[..., start:end]


def _to_channel_first(seq_onehot: Tensor) -> Tensor:
    """Convert sequence tensor to [B, C, L] where C=4."""
    if seq_onehot.dim() != 3:
        raise ValueError(
            f"Expected sequence tensor with 3 dims [B, L, 4] or [B, 4, L], got {tuple(seq_onehot.shape)}"
        )

    if seq_onehot.size(1) == 4:
        return seq_onehot
    if seq_onehot.size(2) == 4:
        return seq_onehot.transpose(1, 2)

    raise ValueError(
        f"Unable to infer channel axis in sequence tensor with shape {tuple(seq_onehot.shape)}; expected a 4-channel axis"
    )


class ResidualDilatedConvBlock(nn.Module):
    """Simple residual dilated conv block."""

    def __init__(self, channels: int, kernel_size: int, dilation: int, dropout: float) -> None:
        super().__init__()
        if kernel_size % 2 == 0:
            raise ValueError("kernel_size must be odd for same-length padding")

        padding = (kernel_size // 2) * dilation
        self.block = nn.Sequential(
            nn.Conv1d(channels, channels, kernel_size=kernel_size, padding=padding, dilation=dilation),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Conv1d(channels, channels, kernel_size=1),
        )
        self.norm = nn.BatchNorm1d(channels)
        self.act = nn.GELU()

    def forward(self, x: Tensor) -> Tensor:
        residual = x
        x = self.block(x)
        x = self.norm(x + residual)
        return self.act(x)


class ChromBPNetBiasBranch(nn.Module):
    """Bias branch that predicts profile/count bias logits from sequence."""

    def __init__(
        self,
        input_channels: int = 4,
        hidden_channels: int = 128,
        kernel_size: int = 21,
        n_dil_layers: int = 4,
        output_len: int = 1000,
        dropout: float = 0.1,
        profile_pool_factor: int = 0,
    ) -> None:
        super().__init__()
        if kernel_size % 2 == 0:
            raise ValueError("kernel_size must be odd for same-length padding")
        if n_dil_layers <= 0:
            raise ValueError(f"n_dil_layers must be positive, got {n_dil_layers}")

        self.output_len = output_len
        self.stem = nn.Sequential(
            nn.Conv1d(
                in_channels=input_channels,
                out_channels=hidden_channels,
                kernel_size=kernel_size,
                padding=kernel_size // 2,
            ),
            nn.GELU(),
            nn.Dropout(dropout),
        )

        blocks = []
        for i in range(n_dil_layers):
            dilation = 2**i
            blocks.append(
                ResidualDilatedConvBlock(
                    channels=hidden_channels,
                    kernel_size=3,
                    dilation=dilation,
                    dropout=dropout,
                )
            )
        self.dilated_tower = nn.Sequential(*blocks)

        self.profile_head = nn.Conv1d(hidden_channels, 1, kernel_size=1)
        self.count_pool = nn.AdaptiveAvgPool1d(1)
        self.count_head = nn.Linear(hidden_channels, 1)

        # Resolution bottleneck: limits profile output to low-frequency patterns.
        # When profile_pool_factor > 0, the profile is downsampled then upsampled,
        # preventing the bias branch from representing fine-grained motif shapes.
        self.profile_pool_factor = int(profile_pool_factor)
        if self.profile_pool_factor > 0:
            bottleneck_len = max(1, output_len // self.profile_pool_factor)
            self.profile_bottleneck = nn.Sequential(
                nn.AdaptiveAvgPool1d(bottleneck_len),
                nn.Upsample(size=output_len, mode="linear", align_corners=False),
            )
        else:
            self.profile_bottleneck = None

    def forward(self, seq_onehot: Tensor) -> Tuple[Tensor, Tensor]:
        """Forward pass.

        Args:
            seq_onehot: [B, L, 4] or [B, 4, L]

        Returns:
            profile_bias_logits: [B, output_len]
            count_bias: [B, 1]
        """
        x = _to_channel_first(seq_onehot)
        x = self.stem(x)
        x = self.dilated_tower(x)

        profile_bias_logits = self.profile_head(x).squeeze(1)
        profile_bias_logits = center_crop_1d(profile_bias_logits, self.output_len)

        if self.profile_bottleneck is not None:
            profile_bias_logits = self.profile_bottleneck(
                profile_bias_logits.unsqueeze(1)
            ).squeeze(1)

        pooled = self.count_pool(x).squeeze(-1)
        count_bias = self.count_head(pooled)

        return profile_bias_logits, count_bias

    def freeze_core(self) -> None:
        """Freeze all bias branch parameters."""
        for p in self.parameters():
            p.requires_grad = False

    def unfreeze_core(self) -> None:
        """Unfreeze all bias branch parameters."""
        for p in self.parameters():
            p.requires_grad = True
```

## Transformer Encoder 源码：transformer_encoder.py
```python
"""Transformer encoder modules with Rotary Positional Embeddings (RoPE).

Design notes:
- Inputs are batch-first: [B, L, D].
- Attention is non-causal by default (genomic context is bidirectional).
- RoPE cache dynamically expands when sequence length exceeds cached size.
- Optional SDPA path uses torch.scaled_dot_product_attention for speed.
"""

from __future__ import annotations

from typing import Optional

import torch
import torch.nn.functional as F
from torch import Tensor, nn


class RotaryEmbedding(nn.Module):
    """Rotary positional embedding with dynamic cache growth."""

    def __init__(self, dim: int, max_seq_len: int = 32768, theta: float = 10000.0) -> None:
        super().__init__()
        if dim <= 0:
            raise ValueError(f"dim must be positive, got {dim}")
        if dim % 2 != 0:
            raise ValueError(f"RoPE dim must be even, got {dim}")
        if max_seq_len <= 0:
            raise ValueError(f"max_seq_len must be positive, got {max_seq_len}")
        if theta <= 0:
            raise ValueError(f"theta must be positive, got {theta}")

        self.dim = dim
        self.theta = float(theta)
        inv_freq = 1.0 / (self.theta ** (torch.arange(0, dim, 2, dtype=torch.float32) / dim))
        self.register_buffer("inv_freq", inv_freq, persistent=False)

        self.max_seq_len_cached = 0
        self.register_buffer("cos_cached", torch.empty(0), persistent=False)
        self.register_buffer("sin_cached", torch.empty(0), persistent=False)
        self._build_cache(max_seq_len)

    def _build_cache(self, seq_len: int) -> None:
        t = torch.arange(seq_len, device=self.inv_freq.device, dtype=self.inv_freq.dtype)
        freqs = torch.einsum("i,j->ij", t, self.inv_freq)
        emb = torch.cat((freqs, freqs), dim=-1)
        self.cos_cached = emb.cos()[None, None, :, :]
        self.sin_cached = emb.sin()[None, None, :, :]
        self.max_seq_len_cached = seq_len

    def forward(self, seq_len: int, device: torch.device, dtype: torch.dtype) -> tuple[Tensor, Tensor]:
        """Return cos/sin with shape [1, 1, seq_len, dim]."""
        if seq_len <= 0:
            raise ValueError(f"seq_len must be positive, got {seq_len}")
        if seq_len > self.max_seq_len_cached:
            self._build_cache(seq_len)

        cos = self.cos_cached[:, :, :seq_len, :].to(device=device, dtype=dtype)
        sin = self.sin_cached[:, :, :seq_len, :].to(device=device, dtype=dtype)
        return cos, sin


def rotate_half(x: Tensor) -> Tensor:
    """Rotate last-dim halves: [x1, x2] -> [-x2, x1]."""
    x1, x2 = x.chunk(2, dim=-1)
    return torch.cat((-x2, x1), dim=-1)


def apply_rotary_pos_emb(q: Tensor, k: Tensor, cos: Tensor, sin: Tensor) -> tuple[Tensor, Tensor]:
    """Apply RoPE to query and key tensors.

    Args:
        q, k: [B, H, L, Dh]
        cos, sin: [1, 1, L, Dh]
    """
    q_out = (q * cos) + (rotate_half(q) * sin)
    k_out = (k * cos) + (rotate_half(k) * sin)
    return q_out, k_out


class RoPEMultiheadAttention(nn.Module):
    """Batch-first multi-head self-attention with RoPE.

    Notes:
    - Non-causal by default.
    - key_padding_mask uses PyTorch convention: True means "ignore this token".
    """

    def __init__(
        self,
        d_model: int,
        n_heads: int,
        dropout: float = 0.1,
        rotary_emb: Optional[RotaryEmbedding] = None,
        use_sdpa: bool = True,
    ) -> None:
        super().__init__()
        if d_model <= 0:
            raise ValueError(f"d_model must be positive, got {d_model}")
        if n_heads <= 0:
            raise ValueError(f"n_heads must be positive, got {n_heads}")
        if d_model % n_heads != 0:
            raise ValueError(f"d_model ({d_model}) must be divisible by n_heads ({n_heads})")

        self.d_model = d_model
        self.n_heads = n_heads
        self.head_dim = d_model // n_heads
        self.scale = self.head_dim**-0.5

        self.q_proj = nn.Linear(d_model, d_model)
        self.k_proj = nn.Linear(d_model, d_model)
        self.v_proj = nn.Linear(d_model, d_model)
        self.out_proj = nn.Linear(d_model, d_model)

        self.dropout_p = float(dropout)
        self.dropout = nn.Dropout(dropout)
        self.rotary_emb = rotary_emb
        self.use_sdpa = bool(use_sdpa) and hasattr(F, "scaled_dot_product_attention")

    def _shape_qkv(self, x: Tensor) -> tuple[Tensor, Tensor, Tensor]:
        bsz, seq_len, _ = x.shape
        q = self.q_proj(x).view(bsz, seq_len, self.n_heads, self.head_dim).transpose(1, 2)
        k = self.k_proj(x).view(bsz, seq_len, self.n_heads, self.head_dim).transpose(1, 2)
        v = self.v_proj(x).view(bsz, seq_len, self.n_heads, self.head_dim).transpose(1, 2)
        return q, k, v

    def forward(self, x: Tensor, key_padding_mask: Optional[Tensor] = None) -> Tensor:
        """Forward self-attention.

        Args:
            x: [B, L, D]
            key_padding_mask: optional [B, L], True means token is padding/ignored.
        """
        if x.dim() != 3:
            raise ValueError(f"Expected x with shape [B, L, D], got {tuple(x.shape)}")

        bsz, seq_len, _ = x.shape
        if key_padding_mask is not None:
            if key_padding_mask.shape != (bsz, seq_len):
                raise ValueError(
                    "key_padding_mask shape mismatch: "
                    f"expected {(bsz, seq_len)}, got {tuple(key_padding_mask.shape)}"
                )
            key_padding_mask = key_padding_mask.to(dtype=torch.bool)
            if torch.any(torch.all(key_padding_mask, dim=1)):
                raise ValueError("Found sequence with all tokens masked in key_padding_mask")

        q, k, v = self._shape_qkv(x)

        if self.rotary_emb is not None:
            cos, sin = self.rotary_emb(seq_len=seq_len, device=q.device, dtype=q.dtype)
            q, k = apply_rotary_pos_emb(q, k, cos, sin)

        if self.use_sdpa:
            attn_mask = None
            if key_padding_mask is not None:
                # SDPA bool mask uses True = keep, False = mask.
                keep_mask = (~key_padding_mask).unsqueeze(1).unsqueeze(2)  # [B,1,1,L]
                attn_mask = keep_mask

            out = F.scaled_dot_product_attention(
                q,
                k,
                v,
                attn_mask=attn_mask,
                dropout_p=self.dropout_p if self.training else 0.0,
                is_causal=False,
            )
        else:
            attn_weights = torch.matmul(q, k.transpose(-2, -1)) * self.scale
            if key_padding_mask is not None:
                mask = key_padding_mask.unsqueeze(1).unsqueeze(2)  # [B,1,1,L]
                attn_weights = attn_weights.masked_fill(mask, float("-inf"))
            attn_weights = F.softmax(attn_weights, dim=-1)
            attn_weights = self.dropout(attn_weights)
            out = torch.matmul(attn_weights, v)

        out = out.transpose(1, 2).contiguous().view(bsz, seq_len, self.d_model)
        return self.out_proj(out)


class RoPETransformerLayer(nn.Module):
    """Pre-Norm transformer encoder layer with RoPE attention."""

    def __init__(
        self,
        d_model: int,
        n_heads: int,
        dim_feedforward: int,
        dropout: float,
        activation: str,
        rotary_emb: Optional[RotaryEmbedding] = None,
        use_sdpa: bool = True,
    ) -> None:
        super().__init__()
        self.norm1 = nn.LayerNorm(d_model)
        self.attn = RoPEMultiheadAttention(
            d_model=d_model,
            n_heads=n_heads,
            dropout=dropout,
            rotary_emb=rotary_emb,
            use_sdpa=use_sdpa,
        )
        self.dropout1 = nn.Dropout(dropout)

        self.norm2 = nn.LayerNorm(d_model)
        self.linear1 = nn.Linear(d_model, dim_feedforward)
        self.dropout = nn.Dropout(dropout)
        self.linear2 = nn.Linear(dim_feedforward, d_model)
        self.dropout2 = nn.Dropout(dropout)

        if activation == "relu":
            self.act = nn.ReLU()
        elif activation == "gelu":
            self.act = nn.GELU()
        else:
            raise ValueError(f"Unsupported activation: {activation}")

    def forward(self, x: Tensor, key_padding_mask: Optional[Tensor] = None) -> Tensor:
        x_norm = self.norm1(x)
        attn_out = self.attn(x_norm, key_padding_mask=key_padding_mask)
        x = x + self.dropout1(attn_out)

        x_norm = self.norm2(x)
        ff_out = self.linear2(self.dropout(self.act(self.linear1(x_norm))))
        x = x + self.dropout2(ff_out)
        return x


class SequenceTransformerEncoder(nn.Module):
    """Transformer encoder with optional RoPE and optional SDPA acceleration."""

    def __init__(
        self,
        d_model: int,
        n_heads: int,
        n_layers: int,
        ff_mult: int = 4,
        dropout: float = 0.1,
        activation: str = "gelu",
        use_positional_encoding: bool = True,
        max_len: int = 32768,
        rope_theta: float = 10000.0,
        use_sdpa: bool = True,
    ) -> None:
        super().__init__()
        if d_model <= 0:
            raise ValueError(f"d_model must be positive, got {d_model}")
        if n_heads <= 0:
            raise ValueError(f"n_heads must be positive, got {n_heads}")
        if d_model % n_heads != 0:
            raise ValueError(f"d_model ({d_model}) must be divisible by n_heads ({n_heads})")
        if n_layers <= 0:
            raise ValueError(f"n_layers must be positive, got {n_layers}")
        if ff_mult <= 0:
            raise ValueError(f"ff_mult must be positive, got {ff_mult}")

        head_dim = d_model // n_heads
        self.rotary_emb = (
            RotaryEmbedding(dim=head_dim, max_seq_len=max_len, theta=rope_theta)
            if use_positional_encoding
            else None
        )

        self.layers = nn.ModuleList(
            [
                RoPETransformerLayer(
                    d_model=d_model,
                    n_heads=n_heads,
                    dim_feedforward=d_model * ff_mult,
                    dropout=dropout,
                    activation=activation,
                    rotary_emb=self.rotary_emb,
                    use_sdpa=use_sdpa,
                )
                for _ in range(n_layers)
            ]
        )
        self.final_norm = nn.LayerNorm(d_model)

    def forward(self, x: Tensor, key_padding_mask: Optional[Tensor] = None) -> Tensor:
        if x.dim() != 3:
            raise ValueError(f"Expected x with shape [B, L, D], got {tuple(x.shape)}")

        for layer in self.layers:
            x = layer(x, key_padding_mask=key_padding_mask)
        return self.final_norm(x)
```

## 给 ChatGPT 的两个阅读提示
1. 这不是一个“普通 backbone 替换”项目。真正的核心是：Transformer 被接入了一个带独立 bias branch 的 factorized 框架，因此模型是否“安全”不能只看最终 full output。
2. 如果你在判断方法贡献，请把 “模型结构收益” 和 “full/debiased 诊断框架” 分开看；后者现在反而是更稳的贡献点。
