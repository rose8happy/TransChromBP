# AlphaGenome-like Factor Ladder Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a three-step AlphaGenome-like factor ladder to the TransChromBP archive repo so we can run `E1` long-context-only, `E2` full-model hierarchical encoder-decoder, and `E3` model-teacher distillation with consistent configs, launchers, tests, and docs.

**Architecture:** Keep `E1` on the existing `corrected B` codepath and change only the context length. Introduce `E2` as a separate hierarchical model variant that still returns `TransChromBPOutput` and preserves the bias-safe/count-head contract, then layer `E3` on top by exporting `profile16/logcount` teacher caches from an `E2` checkpoint without changing the training loop semantics.

**Tech Stack:** Python, PyTorch, YAML configs, pytest, bash launchers, existing `transchrombp.training.train_ddp` and `transchrombp.data.ChromBPNetBigWigDataset`

**Execution Location:** Implement and validate code in the remote runtime repos on `6000` (`/data1/zhoujiazhen/bylw_atac/TransChromBP`) and `6002` (`/home/zhengwei/bylw_atac/TransChromBP`). The local archive repo stays canonical for specs/plans/reports only. When following file paths below, map local archive `vendor/transchrombp/transchrombp/...` paths to remote runtime `src/transchrombp/...` or remote `configs/...` / `scripts/...` peers before editing.

---

## File Structure

**Create:**

- `vendor/transchrombp/transchrombp/configs/data/data_tutorial_canonical_v1_longctx4096.yaml`
- `vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_longctx4096.yaml`
- `vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_hierdec4096.yaml`
- `vendor/transchrombp/transchrombp/configs/train/train_tutorial_teacher_v2_longctx4096_short10.yaml`
- `vendor/transchrombp/transchrombp/configs/train/train_tutorial_teacher_v2_hierdec4096_short10.yaml`
- `vendor/transchrombp/transchrombp/configs/train/train_tutorial_teacher_v2_hierdec4096_teacher30.yaml`
- `vendor/transchrombp/transchrombp/configs/train/train_tutorial_teacher_v2_hierdec4096_distill_short10.yaml`
- `vendor/transchrombp/transchrombp/models/hierarchical_transchrombp.py`
- `vendor/transchrombp/transchrombp/evaluation/model_teacher_cache_export.py`
- `vendor/transchrombp/transchrombp/scripts/run_teacher_v2_longctx4096_probe.sh`
- `vendor/transchrombp/transchrombp/scripts/run_teacher_v2_hierdec4096_probe.sh`
- `vendor/transchrombp/transchrombp/scripts/run_teacher_v2_hierdec4096_distill.sh`
- `tests/test_long_context_probe_launcher.py`
- `tests/test_hierarchical_transchrombp.py`
- `tests/test_model_teacher_cache_export.py`
- `tests/test_factor_ladder_launchers.py`
- `tests/test_factor_ladder_docs.py`

**Modify:**

- `vendor/transchrombp/transchrombp/models/transchrombp.py`
- `vendor/transchrombp/transchrombp/models/__init__.py`
- `TRACKING.md`
- `docs/experiments/registry.md`

**Responsibility Split:**

- Config files describe `E1/E2/E3` experiments without hiding any logic in ad-hoc CLI flags.
- `hierarchical_transchrombp.py` owns the new `E2`/`E3` model family so `transchrombp.py` remains the baseline path.
- `model_teacher_cache_export.py` owns `E2 teacher -> record-aligned profile16/logcount cache` export, leaving `train_ddp.py` unchanged.
- Launcher scripts own runtime orchestration only; they do not embed model logic.
- Root `tests/` hold regression checks for configs, builder dispatch, export format, and launcher wiring.

### Task 1: Scaffold `E1` Long-Context Baseline

**Files:**

- Create: `vendor/transchrombp/transchrombp/configs/data/data_tutorial_canonical_v1_longctx4096.yaml`
- Create: `vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_longctx4096.yaml`
- Create: `vendor/transchrombp/transchrombp/configs/train/train_tutorial_teacher_v2_longctx4096_short10.yaml`
- Create: `vendor/transchrombp/transchrombp/scripts/run_teacher_v2_longctx4096_probe.sh`
- Test: `tests/test_long_context_probe_launcher.py`

- [ ] **Step 1: Write the failing config/launcher smoke test**

```python
from __future__ import annotations

import subprocess
from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_long_context_probe_scaffolding_exists_and_pins_4096() -> None:
    data_cfg_path = (
        REPO_ROOT / "vendor" / "transchrombp" / "transchrombp" / "configs" / "data"
        / "data_tutorial_canonical_v1_longctx4096.yaml"
    )
    model_cfg_path = (
        REPO_ROOT / "vendor" / "transchrombp" / "transchrombp" / "configs" / "model"
        / "transchrombp_teacher_v2_center_pool_longctx4096.yaml"
    )
    train_cfg_path = (
        REPO_ROOT / "vendor" / "transchrombp" / "transchrombp" / "configs" / "train"
        / "train_tutorial_teacher_v2_longctx4096_short10.yaml"
    )
    launcher_path = (
        REPO_ROOT / "vendor" / "transchrombp" / "transchrombp" / "scripts"
        / "run_teacher_v2_longctx4096_probe.sh"
    )

    data_cfg = yaml.safe_load(data_cfg_path.read_text(encoding="utf-8"))
    model_cfg = yaml.safe_load(model_cfg_path.read_text(encoding="utf-8"))
    train_cfg = yaml.safe_load(train_cfg_path.read_text(encoding="utf-8"))

    assert data_cfg["window"]["input_len"] == 4096
    assert model_cfg["sequence_encoder"]["max_len"] == 4096
    assert train_cfg["data"]["input_len"] == 4096
    assert train_cfg["data"]["output_len"] == 1000
    assert train_cfg["data"]["supervised_bp"] == 1000
    assert train_cfg["loss"]["distill_profile_weight"] == 0.0
    assert train_cfg["data"]["teacher_target_names"] == []

    launcher_text = launcher_path.read_text(encoding="utf-8")
    assert "train_tutorial_teacher_v2_longctx4096_short10.yaml" in launcher_text
    assert "transchrombp_teacher_v2_center_pool_longctx4096.yaml" in launcher_text
    assert "data_tutorial_canonical_v1_longctx4096.yaml" in launcher_text

    result = subprocess.run(["bash", "-n", str(launcher_path)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
```

- [ ] **Step 2: Run the test to verify it fails**

Run:

```bash
pytest tests/test_long_context_probe_launcher.py -q
```

Expected:

```text
E   FileNotFoundError: ...data_tutorial_canonical_v1_longctx4096.yaml
```

- [ ] **Step 3: Add the `4096` data/model/train configs**

```yaml
# vendor/transchrombp/transchrombp/configs/data/data_tutorial_canonical_v1_longctx4096.yaml
genome_fasta: /data1/zhoujiazhen/bylw_atac/hg38.fa
chrom_sizes: /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes
folds_json: /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/folds.json

input:
  peaks_bed: /data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/preprocessing/tutorial_canonical_v1/step4_filtering/filtered.peaks.bed
  nonpeaks_bed: /data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/preprocessing/tutorial_canonical_v1/step4_filtering/filtered.nonpeaks.bed
  bigwig: /data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/preprocessing/tutorial_canonical_v1/step2_bigwig/merged_unstranded.bw

window:
  input_len: 4096
  output_len: 1000
  max_jitter: 500

sampling:
  nonpeak_ratio: 0.1
```

```yaml
# vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_longctx4096.yaml
model_name: transchrombp_teacher_v2_center_pool_longctx4096

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
  max_len: 4096
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

```yaml
# vendor/transchrombp/transchrombp/configs/train/train_tutorial_teacher_v2_longctx4096_short10.yaml
seed: 42
max_epochs: 10

trainer:
  backend: nccl
  precision: bf16
  grad_accum_steps: 1
  clip_grad_norm: 1.0
  log_every_steps: 20
  validate_every_epochs: 1
  checkpoint_every_epochs: 1
  find_unused_parameters: false
  compile: false
  strict_max_len_check: true
  best_metric: peak.profile_target_jsd_full_mean
  best_metric_mode: min
  early_stop_patience: 3
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
  distill_profile_weight: 0.0
  distill_count_weight: 0.0
  distill_rank_weight: 0.0

schedule:
  name: cosine
  warmup_steps: 0
  warmup_ratio: 0.06
  min_lr_ratio: 0.1

data:
  source: chrombpnet_bigwig
  config_path: configs/data/data_tutorial_canonical_v1_longctx4096.yaml
  max_seq_len: 4096
  input_len: 4096
  output_len: 1000
  supervised_bp: 1000
  profile_bin_size: 1
  batch_size_per_gpu: 8
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
  teacher_cache_dir: ""
  teacher_target_names: []

logging:
  output_dir: /data1/zhoujiazhen/bylw_atac/TransChromBP/outputs
  run_name: teacher_v2_center_pool_longctx4096_short10
```

- [ ] **Step 4: Add the `E1` launcher**

```bash
#!/usr/bin/env bash
set -euo pipefail

if [ -z "${BASH_VERSION:-}" ]; then
  echo "[error] Please run this script with bash." >&2
  exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TRAIN_GPU_IDS="${TRAIN_GPU_IDS:-${GPU_IDS:-${1:-0}}}"
NPROC_PER_NODE="${NPROC_PER_NODE:-1}"
OUTPUT_BASE="${OUTPUT_BASE:-${ROOT_DIR}/outputs}"
TRAIN_CONFIG="${TRAIN_CONFIG:-${ROOT_DIR}/configs/train/train_tutorial_teacher_v2_longctx4096_short10.yaml}"
MODEL_CONFIG="${MODEL_CONFIG:-${ROOT_DIR}/configs/model/transchrombp_teacher_v2_center_pool_longctx4096.yaml}"
DATA_CONFIG="${DATA_CONFIG:-${ROOT_DIR}/configs/data/data_tutorial_canonical_v1_longctx4096.yaml}"
RUN_NAME="${RUN_NAME:-teacher_v2_center_pool_longctx4096_short10_s42}"
BATCH_SIZE_PER_GPU="${BATCH_SIZE_PER_GPU:-8}"

require_file() {
  local path="$1"
  local label="$2"
  if [ ! -f "${path}" ]; then
    echo "[error] Missing ${label}: ${path}" >&2
    exit 1
  fi
}

require_file "${TRAIN_CONFIG}" "train config"
require_file "${MODEL_CONFIG}" "model config"
require_file "${DATA_CONFIG}" "data config"

export CUDA_VISIBLE_DEVICES="${TRAIN_GPU_IDS}"
export PYTHONPATH="${ROOT_DIR}:${PYTHONPATH:-}"

python -m transchrombp.training.train_ddp \
  --train-config "${TRAIN_CONFIG}" \
  --model-config "${MODEL_CONFIG}" \
  --data-config "${DATA_CONFIG}" \
  --run-name "${RUN_NAME}" \
  --output-dir "${OUTPUT_BASE}" \
  --batch-size-per-gpu "${BATCH_SIZE_PER_GPU}"
```

- [ ] **Step 5: Re-run the smoke test**

Run:

```bash
pytest tests/test_long_context_probe_launcher.py -q
```

Expected:

```text
1 passed
```

- [ ] **Step 6: Commit**

```bash
git add \
  vendor/transchrombp/transchrombp/configs/data/data_tutorial_canonical_v1_longctx4096.yaml \
  vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_longctx4096.yaml \
  vendor/transchrombp/transchrombp/configs/train/train_tutorial_teacher_v2_longctx4096_short10.yaml \
  vendor/transchrombp/transchrombp/scripts/run_teacher_v2_longctx4096_probe.sh \
  tests/test_long_context_probe_launcher.py
git commit -m "feat: scaffold long-context baseline probe"
```

### Task 2: Add the Hierarchical `E2` Model Variant

**Files:**

- Create: `vendor/transchrombp/transchrombp/models/hierarchical_transchrombp.py`
- Modify: `vendor/transchrombp/transchrombp/models/transchrombp.py`
- Modify: `vendor/transchrombp/transchrombp/models/__init__.py`
- Create: `tests/test_hierarchical_transchrombp.py`

- [ ] **Step 1: Write the failing builder/shape test**

```python
from __future__ import annotations

import sys
from pathlib import Path

import torch
import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "vendor" / "transchrombp"))

from transchrombp.models import TransChromBPOutput, build_transchrombp_from_config


def test_hierarchical_variant_builds_and_preserves_output_contract() -> None:
    cfg_path = (
        REPO_ROOT / "vendor" / "transchrombp" / "transchrombp" / "configs" / "model"
        / "transchrombp_teacher_v2_hierdec4096.yaml"
    )
    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    model = build_transchrombp_from_config(cfg)
    x = torch.randn(2, 4096, 4)

    outputs = model(x)

    assert isinstance(outputs, TransChromBPOutput)
    assert outputs.profile_logits_full.shape == (2, 1000)
    assert outputs.profile_logits_debiased.shape == (2, 1000)
    assert outputs.logcount_full.shape == (2, 1)
    assert outputs.logcount_debiased.shape == (2, 1)
    assert getattr(model, "output_len", None) == 1000
    assert hasattr(model, "bias_branch")
```

- [ ] **Step 2: Run the test to verify it fails**

Run:

```bash
PYTHONPATH=vendor/transchrombp pytest tests/test_hierarchical_transchrombp.py -q
```

Expected:

```text
E   FileNotFoundError: ...transchrombp_teacher_v2_hierdec4096.yaml
```

- [ ] **Step 3: Add the new hierarchical model file**

```python
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict

import torch
import torch.nn.functional as F
from torch import Tensor, nn

from .bias_branch import ChromBPNetBiasBranch, center_crop_1d
from .transchrombp import TransChromBPOutput, _inverse_softplus
from .transformer_encoder import SequenceTransformerEncoder


class DownsampleStage(nn.Module):
    def __init__(self, in_channels: int, out_channels: int, dropout: float) -> None:
        super().__init__()
        self.block = nn.Sequential(
            nn.Conv1d(in_channels, out_channels, kernel_size=5, stride=2, padding=2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Conv1d(out_channels, out_channels, kernel_size=3, padding=1),
            nn.GELU(),
        )

    def forward(self, x: Tensor) -> Tensor:
        return self.block(x)


class UpsampleStage(nn.Module):
    def __init__(self, in_channels: int, skip_channels: int, out_channels: int, dropout: float) -> None:
        super().__init__()
        self.fuse = nn.Sequential(
            nn.Conv1d(in_channels + skip_channels, out_channels, kernel_size=3, padding=1),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Conv1d(out_channels, out_channels, kernel_size=3, padding=1),
            nn.GELU(),
        )

    def forward(self, x: Tensor, skip: Tensor) -> Tensor:
        x = F.interpolate(x, size=skip.size(-1), mode="linear", align_corners=False)
        return self.fuse(torch.cat([x, skip], dim=1))


class HierarchicalTransChromBP(nn.Module):
    def __init__(
        self,
        d_model: int = 256,
        output_len: int = 1000,
        dropout: float = 0.1,
        count_pool_mode: str = "center",
        profile_fusion: str = "add",
        count_fusion: str = "logsumexp",
        profile_scale_init: float = 1.0,
        count_scale_init: float = 1.0,
        learnable_scales: bool = True,
        positive_scales: bool = True,
        profile_bias_stop_gradient: bool = True,
        bias_hidden_channels: int = 128,
        bias_kernel_size: int = 21,
        bias_n_dil_layers: int = 4,
        bias_profile_pool_factor: int = 32,
        n_heads: int = 8,
        n_layers: int = 6,
        ff_mult: int = 4,
        transformer_max_len: int = 4096,
        rope_theta: float = 10000.0,
        use_sdpa: bool = True,
    ) -> None:
        super().__init__()
        self.output_len = output_len
        self.count_pool_mode = count_pool_mode
        self.profile_fusion = profile_fusion
        self.count_fusion = count_fusion
        self.learnable_scales = learnable_scales
        self.positive_scales = positive_scales
        self.profile_bias_stop_gradient = profile_bias_stop_gradient

        self.stem = nn.Conv1d(4, d_model, kernel_size=21, padding=10)
        self.enc1 = DownsampleStage(d_model, d_model, dropout)
        self.enc2 = DownsampleStage(d_model, d_model, dropout)
        self.transformer = SequenceTransformerEncoder(
            d_model=d_model,
            n_heads=n_heads,
            n_layers=n_layers,
            ff_mult=ff_mult,
            dropout=dropout,
            activation="gelu",
            use_positional_encoding=True,
            max_len=transformer_max_len // 4,
            rope_theta=rope_theta,
            use_sdpa=use_sdpa,
        )
        self.dec1 = UpsampleStage(d_model, d_model, d_model, dropout)
        self.dec2 = UpsampleStage(d_model, d_model, d_model, dropout)
        self.profile_signal_head = nn.Conv1d(d_model, 1, kernel_size=1)
        self.count_signal_head = nn.Sequential(
            nn.LayerNorm(d_model),
            nn.Linear(d_model, d_model // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_model // 2, 1),
        )
        self.bias_branch = ChromBPNetBiasBranch(
            input_channels=4,
            hidden_channels=bias_hidden_channels,
            kernel_size=bias_kernel_size,
            n_dil_layers=bias_n_dil_layers,
            output_len=output_len,
            dropout=dropout,
            profile_pool_factor=bias_profile_pool_factor,
        )

        if positive_scales and learnable_scales:
            self.profile_scale = nn.Parameter(torch.tensor(_inverse_softplus(profile_scale_init)))
            self.count_scale = nn.Parameter(torch.tensor(_inverse_softplus(count_scale_init)))
        else:
            self.profile_scale = nn.Parameter(torch.tensor(profile_scale_init), requires_grad=learnable_scales)
            self.count_scale = nn.Parameter(torch.tensor(count_scale_init), requires_grad=learnable_scales)

    def _resolve_scale(self, scale: Tensor, ref: Tensor) -> Tensor:
        scale = scale.to(device=ref.device, dtype=ref.dtype)
        if self.positive_scales:
            scale = F.softplus(scale) if self.learnable_scales else torch.clamp_min(scale, 0.0)
        return scale

    def forward(self, seq_onehot: Tensor, bias_profile: Tensor | None = None, bias_count: Tensor | None = None, **_: Any) -> TransChromBPOutput:
        x = seq_onehot.transpose(1, 2) if seq_onehot.size(-1) == 4 else seq_onehot
        s0 = F.gelu(self.stem(x))
        s1 = self.enc1(s0)
        s2 = self.enc2(s1)
        tokens = self.transformer(s2.transpose(1, 2)).transpose(1, 2)
        d1 = self.dec1(tokens, s1)
        d2 = self.dec2(d1, s0)
        profile_signal = center_crop_1d(self.profile_signal_head(d2).squeeze(1), self.output_len)
        pooled = center_crop_1d(d2, self.output_len).mean(dim=-1)
        count_signal = self.count_signal_head(pooled)

        if bias_profile is None or bias_count is None:
            bias_profile, bias_count = self.bias_branch(x)
        if self.profile_bias_stop_gradient:
            bias_profile = bias_profile.detach()

        profile_scale = self._resolve_scale(self.profile_scale, profile_signal)
        count_scale = self._resolve_scale(self.count_scale, count_signal)
        profile_logits_full = profile_signal + (profile_scale * bias_profile)
        logcount_full = torch.logsumexp(torch.stack([count_signal, count_scale * bias_count], dim=0), dim=0)

        return TransChromBPOutput(
            profile_logits_full=profile_logits_full,
            logcount_full=logcount_full,
            profile_logits_debiased=profile_signal,
            logcount_debiased=count_signal,
            profile_bias=bias_profile,
            count_bias=bias_count,
        )


def build_hierarchical_transchrombp_from_config(config: Dict[str, Any]) -> HierarchicalTransChromBP:
    seq_cfg = config.get("sequence_encoder", {})
    bias_cfg = config.get("bias_branch", {})
    fusion_cfg = config.get("fusion", {})
    heads_cfg = config.get("heads", {})
    return HierarchicalTransChromBP(
        d_model=int(seq_cfg.get("d_model", 256)),
        output_len=int(heads_cfg.get("profile_output_len", 1000)),
        dropout=float(seq_cfg.get("dropout", 0.1)),
        count_pool_mode=str(heads_cfg.get("count_pool_mode", "center")),
        profile_fusion=str(fusion_cfg.get("profile_fusion", "add")),
        count_fusion=str(fusion_cfg.get("count_fusion", "logsumexp")),
        profile_scale_init=float(fusion_cfg.get("profile_scale_init", 1.0)),
        count_scale_init=float(fusion_cfg.get("count_scale_init", 1.0)),
        learnable_scales=bool(fusion_cfg.get("learnable_scales", True)),
        positive_scales=bool(fusion_cfg.get("positive_scales", True)),
        profile_bias_stop_gradient=bool(fusion_cfg.get("profile_bias_stop_gradient", True)),
        bias_hidden_channels=int(bias_cfg.get("hidden_channels", 128)),
        bias_kernel_size=int(bias_cfg.get("kernel_size", 21)),
        bias_n_dil_layers=int(bias_cfg.get("n_dil_layers", 4)),
        bias_profile_pool_factor=int(bias_cfg.get("profile_pool_factor", 32)),
        n_heads=int(seq_cfg.get("n_heads", 8)),
        n_layers=int(seq_cfg.get("n_layers", 6)),
        ff_mult=int(seq_cfg.get("ff_mult", 4)),
        transformer_max_len=int(seq_cfg.get("max_len", 4096)),
        rope_theta=float(seq_cfg.get("rope_theta", 10000.0)),
        use_sdpa=bool(seq_cfg.get("use_sdpa", True)),
    )
```

- [ ] **Step 4: Dispatch the new architecture from the existing builder**

```python
# vendor/transchrombp/transchrombp/models/transchrombp.py
def build_transchrombp_from_config(config: Dict[str, Any]) -> TransChromBP:
    arch_cfg = config.get("architecture", {})
    variant = str(arch_cfg.get("variant", "baseline_v1")).lower()
    if variant == "hierarchical_encoder_decoder_v1":
        from .hierarchical_transchrombp import build_hierarchical_transchrombp_from_config

        return build_hierarchical_transchrombp_from_config(config)

    seq_cfg = config.get("sequence_encoder", {})
    ...
```

```python
# vendor/transchrombp/transchrombp/models/__init__.py
from .hierarchical_transchrombp import HierarchicalTransChromBP, build_hierarchical_transchrombp_from_config

__all__ = [
    ...
    "HierarchicalTransChromBP",
    "build_hierarchical_transchrombp_from_config",
]
```

```yaml
# vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_hierdec4096.yaml
model_name: transchrombp_teacher_v2_hierdec4096

architecture:
  variant: hierarchical_encoder_decoder_v1

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
  max_len: 4096
  use_sdpa: true

bias_branch:
  enabled: true
  hidden_channels: 128
  kernel_size: 21
  n_dil_layers: 4
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

- [ ] **Step 5: Re-run the builder/shape test**

Run:

```bash
PYTHONPATH=vendor/transchrombp pytest tests/test_hierarchical_transchrombp.py -q
```

Expected:

```text
1 passed
```

- [ ] **Step 6: Commit**

```bash
git add \
  vendor/transchrombp/transchrombp/models/hierarchical_transchrombp.py \
  vendor/transchrombp/transchrombp/models/transchrombp.py \
  vendor/transchrombp/transchrombp/models/__init__.py \
  vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_hierdec4096.yaml \
  tests/test_hierarchical_transchrombp.py
git commit -m "feat: add hierarchical encoder-decoder variant"
```

### Task 3: Wire the `E2` Launcher and Training Configs

**Files:**

- Create: `vendor/transchrombp/transchrombp/configs/train/train_tutorial_teacher_v2_hierdec4096_short10.yaml`
- Create: `vendor/transchrombp/transchrombp/configs/train/train_tutorial_teacher_v2_hierdec4096_teacher30.yaml`
- Create: `vendor/transchrombp/transchrombp/scripts/run_teacher_v2_hierdec4096_probe.sh`
- Create: `tests/test_factor_ladder_launchers.py`

- [ ] **Step 1: Write the failing launcher/config smoke test**

```python
from __future__ import annotations

import subprocess
from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_hierdec_probe_launcher_targets_hierdec4096_assets() -> None:
    short_cfg_path = (
        REPO_ROOT / "vendor" / "transchrombp" / "transchrombp" / "configs" / "train"
        / "train_tutorial_teacher_v2_hierdec4096_short10.yaml"
    )
    teacher_cfg_path = (
        REPO_ROOT / "vendor" / "transchrombp" / "transchrombp" / "configs" / "train"
        / "train_tutorial_teacher_v2_hierdec4096_teacher30.yaml"
    )
    launcher_path = (
        REPO_ROOT / "vendor" / "transchrombp" / "transchrombp" / "scripts"
        / "run_teacher_v2_hierdec4096_probe.sh"
    )

    short_cfg = yaml.safe_load(short_cfg_path.read_text(encoding="utf-8"))
    teacher_cfg = yaml.safe_load(teacher_cfg_path.read_text(encoding="utf-8"))
    launcher_text = launcher_path.read_text(encoding="utf-8")

    assert short_cfg["data"]["input_len"] == 4096
    assert short_cfg["logging"]["run_name"] == "teacher_v2_hierdec4096_short10"
    assert teacher_cfg["max_epochs"] == 30
    assert "transchrombp_teacher_v2_hierdec4096.yaml" in launcher_text
    assert "train_tutorial_teacher_v2_hierdec4096_short10.yaml" in launcher_text
    assert "data_tutorial_canonical_v1_longctx4096.yaml" in launcher_text

    result = subprocess.run(["bash", "-n", str(launcher_path)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
```

- [ ] **Step 2: Run the test to verify it fails**

Run:

```bash
pytest tests/test_factor_ladder_launchers.py::test_hierdec_probe_launcher_targets_hierdec4096_assets -q
```

Expected:

```text
E   FileNotFoundError: ...train_tutorial_teacher_v2_hierdec4096_short10.yaml
```

- [ ] **Step 3: Add the `E2` training configs**

```yaml
# vendor/transchrombp/transchrombp/configs/train/train_tutorial_teacher_v2_hierdec4096_short10.yaml
seed: 42
max_epochs: 10

trainer:
  backend: nccl
  precision: bf16
  grad_accum_steps: 1
  clip_grad_norm: 1.0
  log_every_steps: 20
  validate_every_epochs: 1
  checkpoint_every_epochs: 1
  find_unused_parameters: false
  compile: false
  strict_max_len_check: true
  best_metric: peak.profile_target_jsd_full_mean
  best_metric_mode: min
  early_stop_patience: 3
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
  distill_profile_weight: 0.0
  distill_count_weight: 0.0
  distill_rank_weight: 0.0

schedule:
  name: cosine
  warmup_steps: 0
  warmup_ratio: 0.06
  min_lr_ratio: 0.1

data:
  source: chrombpnet_bigwig
  config_path: configs/data/data_tutorial_canonical_v1_longctx4096.yaml
  max_seq_len: 4096
  input_len: 4096
  output_len: 1000
  supervised_bp: 1000
  profile_bin_size: 1
  batch_size_per_gpu: 8
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
  teacher_cache_dir: ""
  teacher_target_names: []

logging:
  output_dir: /data1/zhoujiazhen/bylw_atac/TransChromBP/outputs
  run_name: teacher_v2_hierdec4096_short10
```

```yaml
# vendor/transchrombp/transchrombp/configs/train/train_tutorial_teacher_v2_hierdec4096_teacher30.yaml
seed: 42
max_epochs: 30

trainer:
  backend: nccl
  precision: bf16
  grad_accum_steps: 1
  clip_grad_norm: 1.0
  log_every_steps: 20
  validate_every_epochs: 1
  checkpoint_every_epochs: 1
  find_unused_parameters: false
  compile: false
  strict_max_len_check: true
  best_metric: peak.profile_target_jsd_full_mean
  best_metric_mode: min
  early_stop_patience: 6
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
  distill_profile_weight: 0.0
  distill_count_weight: 0.0
  distill_rank_weight: 0.0

schedule:
  name: cosine
  warmup_steps: 0
  warmup_ratio: 0.06
  min_lr_ratio: 0.1

data:
  source: chrombpnet_bigwig
  config_path: configs/data/data_tutorial_canonical_v1_longctx4096.yaml
  max_seq_len: 4096
  input_len: 4096
  output_len: 1000
  supervised_bp: 1000
  profile_bin_size: 1
  batch_size_per_gpu: 8
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
  teacher_cache_dir: ""
  teacher_target_names: []

logging:
  output_dir: /data1/zhoujiazhen/bylw_atac/TransChromBP/outputs
  run_name: teacher_v2_hierdec4096_teacher30
```

- [ ] **Step 4: Add the `E2` probe launcher**

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TRAIN_GPU_IDS="${TRAIN_GPU_IDS:-${GPU_IDS:-${1:-0}}}"
NPROC_PER_NODE="${NPROC_PER_NODE:-1}"
OUTPUT_BASE="${OUTPUT_BASE:-${ROOT_DIR}/outputs}"
TRAIN_CONFIG="${TRAIN_CONFIG:-${ROOT_DIR}/configs/train/train_tutorial_teacher_v2_hierdec4096_short10.yaml}"
MODEL_CONFIG="${MODEL_CONFIG:-${ROOT_DIR}/configs/model/transchrombp_teacher_v2_hierdec4096.yaml}"
DATA_CONFIG="${DATA_CONFIG:-${ROOT_DIR}/configs/data/data_tutorial_canonical_v1_longctx4096.yaml}"
RUN_NAME="${RUN_NAME:-teacher_v2_hierdec4096_short10_s42}"
BATCH_SIZE_PER_GPU="${BATCH_SIZE_PER_GPU:-8}"

require_file() {
  local path="$1"
  local label="$2"
  if [ ! -f "${path}" ]; then
    echo "[error] Missing ${label}: ${path}" >&2
    exit 1
  fi
}

require_file "${TRAIN_CONFIG}" "train config"
require_file "${MODEL_CONFIG}" "model config"
require_file "${DATA_CONFIG}" "data config"

export CUDA_VISIBLE_DEVICES="${TRAIN_GPU_IDS}"
export PYTHONPATH="${ROOT_DIR}:${PYTHONPATH:-}"

python -m transchrombp.training.train_ddp \
  --train-config "${TRAIN_CONFIG}" \
  --model-config "${MODEL_CONFIG}" \
  --data-config "${DATA_CONFIG}" \
  --run-name "${RUN_NAME}" \
  --output-dir "${OUTPUT_BASE}" \
  --batch-size-per-gpu "${BATCH_SIZE_PER_GPU}"
```

- [ ] **Step 5: Re-run the launcher test**

Run:

```bash
pytest tests/test_factor_ladder_launchers.py::test_hierdec_probe_launcher_targets_hierdec4096_assets -q
```

Expected:

```text
1 passed
```

- [ ] **Step 6: Commit**

```bash
git add \
  vendor/transchrombp/transchrombp/configs/train/train_tutorial_teacher_v2_hierdec4096_short10.yaml \
  vendor/transchrombp/transchrombp/configs/train/train_tutorial_teacher_v2_hierdec4096_teacher30.yaml \
  vendor/transchrombp/transchrombp/scripts/run_teacher_v2_hierdec4096_probe.sh \
  tests/test_factor_ladder_launchers.py
git commit -m "feat: wire hierdec4096 probe configs"
```

### Task 4: Export Record-Aligned `E2` Teacher Caches

**Files:**

- Create: `vendor/transchrombp/transchrombp/evaluation/model_teacher_cache_export.py`
- Create: `tests/test_model_teacher_cache_export.py`

- [ ] **Step 1: Write the failing teacher-cache export test**

```python
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import torch


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "vendor" / "transchrombp"))

from transchrombp.evaluation.model_teacher_cache_export import (
    export_arrays,
    pool_profile_probs_to_bins,
    write_teacher_manifest,
)


def test_model_teacher_export_helpers_write_profile16_logcount_manifest(tmp_path: Path) -> None:
    probs = torch.full((2, 1000), 1 / 1000, dtype=torch.float32)
    pooled = pool_profile_probs_to_bins(probs, num_bins=16)
    assert pooled.shape == (2, 16)
    assert torch.allclose(pooled.sum(dim=-1), torch.ones(2))

    export_arrays(tmp_path, "valid", pooled.numpy(), np.ones((2, 1), dtype=np.float32))
    write_teacher_manifest(
        tmp_path,
        split="valid",
        n_records=2,
        record_sha1="abc",
        model_checkpoint="/tmp/model.pt",
        targets=["profile16", "logcount"],
    )

    manifest = json.loads((tmp_path / "teacher_manifest_valid.json").read_text(encoding="utf-8"))
    assert manifest["targets"] == ["profile16", "logcount"]
    assert manifest["record_sha1"] == "abc"
    assert (tmp_path / "valid_profile16.f32.npy").is_file()
    assert (tmp_path / "valid_logcount.f32.npy").is_file()
```

- [ ] **Step 2: Run the test to verify it fails**

Run:

```bash
PYTHONPATH=vendor/transchrombp pytest tests/test_model_teacher_cache_export.py -q
```

Expected:

```text
E   ModuleNotFoundError: No module named 'transchrombp.evaluation.model_teacher_cache_export'
```

- [ ] **Step 3: Add the model-teacher export module**

```python
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import torch
import yaml
from torch.utils.data import DataLoader

from transchrombp.data import ChromBPNetBigWigDataset, compute_record_sha1
from transchrombp.models import build_transchrombp_from_config
from transchrombp.training.train_ddp import apply_external_data_config_defaults, resolve_data_config_path


def pool_profile_probs_to_bins(profile_probs: torch.Tensor, num_bins: int = 16) -> torch.Tensor:
    if profile_probs.size(-1) % num_bins != 0:
        raise ValueError(f"profile length {profile_probs.size(-1)} must be divisible by {num_bins}")
    bin_width = profile_probs.size(-1) // num_bins
    pooled = profile_probs.reshape(profile_probs.size(0), num_bins, bin_width).sum(dim=-1)
    return pooled / pooled.sum(dim=-1, keepdim=True)


def export_arrays(output_dir: Path, split: str, profile16: np.ndarray, logcount: np.ndarray) -> None:
    np.save(output_dir / f"{split}_profile16.f32.npy", profile16.astype(np.float32, copy=False))
    np.save(output_dir / f"{split}_logcount.f32.npy", logcount.astype(np.float32, copy=False))


def write_teacher_manifest(
    output_dir: Path,
    *,
    split: str,
    n_records: int,
    record_sha1: str,
    model_checkpoint: str,
    targets: Sequence[str],
) -> None:
    manifest = {
        "split": split,
        "n_records": int(n_records),
        "record_sha1": str(record_sha1),
        "model_checkpoint": str(model_checkpoint),
        "targets": list(targets),
    }
    (output_dir / f"teacher_manifest_{split}.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def export_model_teacher_cache(
    *,
    checkpoint_path: str,
    model_config: str,
    train_config: str,
    data_config: str,
    output_dir: str,
    splits: Sequence[str] = ("train", "valid"),
    batch_size: int = 8,
) -> None:
    model_cfg = yaml.safe_load(Path(model_config).read_text(encoding="utf-8"))
    train_cfg = yaml.safe_load(Path(train_config).read_text(encoding="utf-8"))
    data_cfg_path = resolve_data_config_path(train_cfg, data_config)
    data_source_cfg = yaml.safe_load(Path(data_cfg_path).read_text(encoding="utf-8"))
    apply_external_data_config_defaults(train_cfg, data_source_cfg)

    model = build_transchrombp_from_config(model_cfg)
    payload = torch.load(checkpoint_path, map_location="cpu")
    model.load_state_dict(payload["model_state"] if "model_state" in payload else payload)
    model.eval()

    out_root = Path(output_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    for split in splits:
        ds = ChromBPNetBigWigDataset(
            genome_fasta=str(train_cfg["data"]["genome_fasta"]),
            bigwig_path=str(train_cfg["data"]["bigwig"]),
            peaks_bed=str(train_cfg["data"]["peaks_bed"]),
            nonpeaks_bed=str(train_cfg["data"]["nonpeaks_bed"]),
            folds_json=str(train_cfg["data"]["folds_json"]),
            split=split,
            input_len=int(train_cfg["data"]["input_len"]),
            supervised_bp=int(train_cfg["data"]["supervised_bp"]),
            profile_bin_size=int(train_cfg["data"]["profile_bin_size"]),
            nonpeak_ratio=float(train_cfg["data"]["nonpeak_ratio"]),
        )
        loader = DataLoader(ds, batch_size=batch_size, shuffle=False)
        profile_rows = []
        count_rows = []
        with torch.no_grad():
            for batch in loader:
                outputs = model(batch["seq"].float())
                probs = torch.softmax(outputs.profile_logits_debiased, dim=-1)
                profile_rows.append(pool_profile_probs_to_bins(probs).cpu())
                count_rows.append(outputs.logcount_debiased.cpu())
        profile16 = torch.cat(profile_rows, dim=0).numpy()
        logcount = torch.cat(count_rows, dim=0).numpy()
        export_arrays(out_root, split, profile16, logcount)
        write_teacher_manifest(
            out_root,
            split=split,
            n_records=len(ds),
            record_sha1=compute_record_sha1(ds.records),
            model_checkpoint=checkpoint_path,
            targets=["profile16", "logcount"],
        )


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint-path", required=True)
    parser.add_argument("--model-config", required=True)
    parser.add_argument("--train-config", required=True)
    parser.add_argument("--data-config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--splits", nargs="+", default=["train", "valid"])
    parser.add_argument("--batch-size", type=int, default=8)
    return parser.parse_args(list(argv) if argv is not None else None)
```

- [ ] **Step 4: Re-run the export helper test**

Run:

```bash
PYTHONPATH=vendor/transchrombp pytest tests/test_model_teacher_cache_export.py -q
```

Expected:

```text
1 passed
```

- [ ] **Step 5: Commit**

```bash
git add \
  vendor/transchrombp/transchrombp/evaluation/model_teacher_cache_export.py \
  tests/test_model_teacher_cache_export.py
git commit -m "feat: export model teacher caches"
```

### Task 5: Wire `E3` Distillation Around the Exported Teacher Cache

**Files:**

- Create: `vendor/transchrombp/transchrombp/configs/train/train_tutorial_teacher_v2_hierdec4096_distill_short10.yaml`
- Create: `vendor/transchrombp/transchrombp/scripts/run_teacher_v2_hierdec4096_distill.sh`
- Modify: `tests/test_factor_ladder_launchers.py`

- [ ] **Step 1: Extend the launcher test with a failing distill case**

```python
def test_hierdec_distill_launcher_exports_teacher_cache_then_trains() -> None:
    launcher_path = (
        REPO_ROOT / "vendor" / "transchrombp" / "transchrombp" / "scripts"
        / "run_teacher_v2_hierdec4096_distill.sh"
    )
    train_cfg_path = (
        REPO_ROOT / "vendor" / "transchrombp" / "transchrombp" / "configs" / "train"
        / "train_tutorial_teacher_v2_hierdec4096_distill_short10.yaml"
    )

    train_cfg = yaml.safe_load(train_cfg_path.read_text(encoding="utf-8"))
    launcher_text = launcher_path.read_text(encoding="utf-8")

    assert train_cfg["loss"]["distill_profile_weight"] == 0.5
    assert train_cfg["loss"]["distill_count_weight"] == 0.05
    assert train_cfg["data"]["teacher_target_names"] == ["profile16", "logcount"]
    assert "transchrombp.evaluation.model_teacher_cache_export" in launcher_text
    assert "teacher_cache_dir" in launcher_text

    result = subprocess.run(["bash", "-n", str(launcher_path)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
```

- [ ] **Step 2: Run the test to verify it fails**

Run:

```bash
pytest tests/test_factor_ladder_launchers.py::test_hierdec_distill_launcher_exports_teacher_cache_then_trains -q
```

Expected:

```text
E   FileNotFoundError: ...train_tutorial_teacher_v2_hierdec4096_distill_short10.yaml
```

- [ ] **Step 3: Add the distill config**

```yaml
# vendor/transchrombp/transchrombp/configs/train/train_tutorial_teacher_v2_hierdec4096_distill_short10.yaml
seed: 42
max_epochs: 10

trainer:
  backend: nccl
  precision: bf16
  grad_accum_steps: 1
  clip_grad_norm: 1.0
  log_every_steps: 20
  validate_every_epochs: 1
  checkpoint_every_epochs: 1
  find_unused_parameters: false
  compile: false
  strict_max_len_check: true
  best_metric: peak.profile_target_jsd_full_mean
  best_metric_mode: min
  early_stop_patience: 3
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
  distill_profile_weight: 0.5
  distill_count_weight: 0.05
  distill_rank_weight: 0.0

schedule:
  name: cosine
  warmup_steps: 0
  warmup_ratio: 0.06
  min_lr_ratio: 0.1

data:
  source: chrombpnet_bigwig
  config_path: configs/data/data_tutorial_canonical_v1_longctx4096.yaml
  max_seq_len: 4096
  input_len: 4096
  output_len: 1000
  supervised_bp: 1000
  profile_bin_size: 1
  batch_size_per_gpu: 8
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
  teacher_cache_dir: outputs/teacher_cache/hierdec4096_teacher30
  teacher_target_names: ["profile16", "logcount"]

logging:
  output_dir: /data1/zhoujiazhen/bylw_atac/TransChromBP/outputs
  run_name: teacher_v2_hierdec4096_distill_short10
```

- [ ] **Step 4: Add the distill launcher**

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TRAIN_GPU_IDS="${TRAIN_GPU_IDS:-${GPU_IDS:-${1:-0}}}"
OUTPUT_BASE="${OUTPUT_BASE:-${ROOT_DIR}/outputs}"
MODEL_CONFIG="${MODEL_CONFIG:-${ROOT_DIR}/configs/model/transchrombp_teacher_v2_hierdec4096.yaml}"
TEACHER_TRAIN_CONFIG="${TEACHER_TRAIN_CONFIG:-${ROOT_DIR}/configs/train/train_tutorial_teacher_v2_hierdec4096_teacher30.yaml}"
STUDENT_TRAIN_CONFIG="${STUDENT_TRAIN_CONFIG:-${ROOT_DIR}/configs/train/train_tutorial_teacher_v2_hierdec4096_distill_short10.yaml}"
DATA_CONFIG="${DATA_CONFIG:-${ROOT_DIR}/configs/data/data_tutorial_canonical_v1_longctx4096.yaml}"
TEACHER_CKPT="${TEACHER_CKPT:?set TEACHER_CKPT to the E2 teacher checkpoint}"
TEACHER_CACHE_DIR="${TEACHER_CACHE_DIR:-${OUTPUT_BASE}/teacher_cache/hierdec4096_teacher30}"
RUN_NAME="${RUN_NAME:-teacher_v2_hierdec4096_distill_short10_s42}"

require_file() {
  local path="$1"
  local label="$2"
  if [ ! -f "${path}" ]; then
    echo "[error] Missing ${label}: ${path}" >&2
    exit 1
  fi
}

require_file "${MODEL_CONFIG}" "model config"
require_file "${TEACHER_TRAIN_CONFIG}" "teacher train config"
require_file "${STUDENT_TRAIN_CONFIG}" "student train config"
require_file "${DATA_CONFIG}" "data config"
require_file "${TEACHER_CKPT}" "teacher checkpoint"

mkdir -p "${TEACHER_CACHE_DIR}"
export CUDA_VISIBLE_DEVICES="${TRAIN_GPU_IDS}"
export PYTHONPATH="${ROOT_DIR}:${PYTHONPATH:-}"

if [ ! -f "${TEACHER_CACHE_DIR}/teacher_manifest_valid.json" ]; then
  python -m transchrombp.evaluation.model_teacher_cache_export \
    --checkpoint-path "${TEACHER_CKPT}" \
    --model-config "${MODEL_CONFIG}" \
    --train-config "${TEACHER_TRAIN_CONFIG}" \
    --data-config "${DATA_CONFIG}" \
    --output-dir "${TEACHER_CACHE_DIR}" \
    --splits train valid
fi

python -m transchrombp.training.train_ddp \
  --train-config "${STUDENT_TRAIN_CONFIG}" \
  --model-config "${MODEL_CONFIG}" \
  --data-config "${DATA_CONFIG}" \
  --teacher-cache-dir "${TEACHER_CACHE_DIR}" \
  --run-name "${RUN_NAME}" \
  --output-dir "${OUTPUT_BASE}" \
  --batch-size-per-gpu 8
```

- [ ] **Step 5: Re-run the distill launcher test**

Run:

```bash
pytest tests/test_factor_ladder_launchers.py::test_hierdec_distill_launcher_exports_teacher_cache_then_trains -q
```

Expected:

```text
1 passed
```

- [ ] **Step 6: Commit**

```bash
git add \
  vendor/transchrombp/transchrombp/configs/train/train_tutorial_teacher_v2_hierdec4096_distill_short10.yaml \
  vendor/transchrombp/transchrombp/scripts/run_teacher_v2_hierdec4096_distill.sh \
  tests/test_factor_ladder_launchers.py
git commit -m "feat: add hierdec4096 distill wiring"
```

### Task 6: Register the New Family in Live Docs Before Any Run Launch

**Files:**

- Modify: `TRACKING.md`
- Modify: `docs/experiments/registry.md`

- [ ] **Step 1: Write the failing docs regression check**

```python
from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_factor_ladder_docs_are_registered() -> None:
    tracking = (REPO_ROOT / "TRACKING.md").read_text(encoding="utf-8")
    registry = (REPO_ROOT / "docs" / "experiments" / "registry.md").read_text(encoding="utf-8")

    assert "AlphaGenome-like factor ladder" in tracking
    assert "`alphagenome_factor_ladder`" in registry
```

- [ ] **Step 2: Run the docs regression check to verify it fails**

Run:

```bash
pytest tests/test_factor_ladder_docs.py -q
```

Expected:

```text
E   AssertionError: assert 'AlphaGenome-like factor ladder' in ...
```

- [ ] **Step 3: Add the docs regression file**

```python
# tests/test_factor_ladder_docs.py
from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_factor_ladder_docs_are_registered() -> None:
    tracking = (REPO_ROOT / "TRACKING.md").read_text(encoding="utf-8")
    registry = (REPO_ROOT / "docs" / "experiments" / "registry.md").read_text(encoding="utf-8")

    assert "AlphaGenome-like factor ladder" in tracking
    assert "`alphagenome_factor_ladder`" in registry
```

- [ ] **Step 4: Update `TRACKING.md` and `registry.md`**

```markdown
| AlphaGenome-like factor ladder（E1/E2/E3） | 待处理 | `E1=4096 corrected B`、`E2=hierdec4096`、`E3=model-teacher distill` 的实现计划已写定，尚未起跑。 | 先完成代码与 launcher，再按 stop-rule 依次启动 `Phase 0` smoke。 | `docs/superpowers/specs/2026-04-11-alphagenome-like-factor-ladder-design.md`、`docs/superpowers/plans/2026-04-11-alphagenome-like-factor-ladder.md` |
```

```markdown
| `alphagenome_factor_ladder` | `planned` | `master` | `n/a` | `n/a` | `docs/superpowers/specs/2026-04-11-alphagenome-like-factor-ladder-design.md` | 仅允许按 `E1 -> E2 -> E3` 与 stop-rule 启动 | `n/a` | 新 family 尚未起跑；必须先完成 4096 long-context、hierdec 和 model-teacher cache wiring。 |
```

- [ ] **Step 5: Run the docs regression check**

Run:

```bash
pytest tests/test_factor_ladder_docs.py -q
```

Expected:

```text
1 passed
```

- [ ] **Step 6: Commit**

```bash
git add TRACKING.md docs/experiments/registry.md tests/test_factor_ladder_docs.py
git commit -m "docs: register alphagenome factor ladder family"
```

## Self-Review

- Spec coverage: `E1` long-context config/launcher is covered in Task 1; `E2` hierarchical model and probe wiring are covered in Tasks 2-3; `E3` model-teacher distillation is covered in Tasks 4-5; live-doc registration before run launch is covered in Task 6.
- Placeholder scan: no `TBD`, `TODO`, or “similar to Task N” shortcuts remain.
- Type consistency: all distill tasks use the existing `profile16` / `logcount` interface, and all `E2` paths reference the same `transchrombp_teacher_v2_hierdec4096.yaml` model config.
