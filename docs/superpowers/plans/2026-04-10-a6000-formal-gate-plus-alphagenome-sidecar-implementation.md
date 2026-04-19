# A6000 Formal Gate Plus AlphaGenome Sidecar Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement a clean `corrected B + multiscale/local-skip decoder v2` 30-epoch A6000 formal gate while preparing and launching the `AlphaGenome 12-20 loci` sidecar as an independent, non-GPU blocking track.

**Architecture:** Treat the work as two coupled but independent deliverables. First add the missing local decoder code path, configs, tests, and launch assets inside the canonical repo so the A6000 gate is scientifically clean and reproducible; then add a deterministic `AlphaGenome v2` matched-panel builder and committed CSV so the sidecar is reproducible too. After the local assets are verified, deploy the exact runtime files into the 6000 `src/transchrombp` workspace, launch the GPU gate and the API sidecar as separate background jobs, update live docs immediately, and finally close out each track with its own report and archive move.

**Tech Stack:** Python, PyTorch, YAML, Bash, `torchrun`, `pyBigWig`, AlphaGenome SDK/API, SSH/SCP, Markdown

---

## File Map

- Create: `vendor/transchrombp/transchrombp/models/profile_decoder.py`
  Purpose: house the new `MultiScaleLocalSkipDecoderV2` module used only by the profile head.
- Modify: `vendor/transchrombp/transchrombp/models/transchrombp.py`
  Purpose: wire the optional profile decoder into `TransChromBP` while keeping bias/count semantics unchanged.
- Modify: `vendor/transchrombp/transchrombp/models/__init__.py`
  Purpose: export the new decoder class cleanly.
- Create: `vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_msdls_v2.yaml`
  Purpose: define the clean `corrected B + multiscale/local-skip decoder v2` model config.
- Create: `tests/test_multiscale_local_skip_v2.py`
  Purpose: enforce shape, config, and bias-safe contract checks for the new decoder path.
- Create: `vendor/transchrombp/transchrombp/scripts/run_msdls_v2_gate.sh`
  Purpose: foreground launcher for the A6000 30-epoch formal gate that writes a runtime train config with `seed=42`.
- Create: `scripts/alphagenome_pilot/build_matched_panel_v2.py`
  Purpose: build a deterministic `12-20 loci` matched panel from tutorial canonical test chromosomes.
- Create: `tests/test_alphagenome_panel_builder.py`
  Purpose: unit-test the deterministic panel-selection logic without requiring remote files.
- Create: `scripts/alphagenome_pilot/regions_k562_tutorial_matched_panel_v2.csv`
  Purpose: committed `AlphaGenome v2` sidecar panel artifact; it also serves as the local comparator table because it carries the matched loci metadata and `observed_total` for the expanded panel.
- Create: `scripts/deploy_a6000_dual_track_to_6000.sh`
  Purpose: sync the exact local code/config/panel files into the 6000 runtime repo’s `src/transchrombp`, `configs`, and `scripts` paths.
- Modify: `TRACKING.md`
  Purpose: record live launch state, then final closeout state.
- Modify: `docs/plan/dual_machine_experiment_charter_20260409.md`
  Purpose: record both active tracks and their post-closeout next-step boundaries.
- Create: `reports/closeout/multiscale_local_skip_v2_30ep_gate_closeout_20260410.md`
  Purpose: hard verdict for the A6000 30-epoch formal gate.
- Create: `reports/closeout/alphagenome_matched_raw_track_slice_v2_closeout_20260410.md`
  Purpose: hard verdict for the `12-20 loci` AlphaGenome sidecar.
- Modify: `TRACKING_archive.md`
  Purpose: archive whichever track reaches terminal state first after closeout.

---

### Task 1: Add the local multiscale/local-skip decoder path

**Files:**
- Create: `vendor/transchrombp/transchrombp/models/profile_decoder.py`
- Modify: `vendor/transchrombp/transchrombp/models/transchrombp.py`
- Modify: `vendor/transchrombp/transchrombp/models/__init__.py`
- Test: `tests/test_multiscale_local_skip_v2.py`

- [ ] **Step 1: Write the failing decoder-path tests**

Create `tests/test_multiscale_local_skip_v2.py` with this test module:

```python
from __future__ import annotations

from copy import deepcopy

import pytest
import torch

from transchrombp.models.transchrombp import build_transchrombp_from_config


def _base_config() -> dict:
    return {
        "sequence_encoder": {
            "enabled": True,
            "d_model": 64,
            "n_heads": 4,
            "n_layers": 2,
            "ff_mult": 2,
            "dropout": 0.0,
            "use_rope": True,
            "rope_theta": 10000.0,
            "max_len": 2114,
            "use_sdpa": True,
        },
        "conv_stem": {
            "stem_kernel_size": 21,
            "conv_kernel_size": 7,
            "n_conv_layers": 2,
        },
        "local_tower": {
            "kernel_size": 3,
            "n_dil_layers": 4,
            "dilation_cycle_length": 4,
        },
        "bias_branch": {
            "enabled": True,
            "hidden_channels": 32,
            "kernel_size": 21,
            "n_dil_layers": 2,
            "profile_pool_factor": 32,
        },
        "fusion": {
            "profile_fusion": "add",
            "count_fusion": "logsumexp",
            "profile_scale_init": 1.0,
            "count_scale_init": 1.0,
            "learnable_scales": True,
            "positive_scales": True,
            "profile_bias_stop_gradient": True,
        },
        "heads": {
            "profile_output_len": 1000,
            "count_head": "linear",
            "count_pool_mode": "center",
        },
        "profile_decoder": {
            "mode": "multiscale_local_skip_v2",
            "hidden_channels": 64,
            "decoder_channels": [64, 48, 32],
            "dropout": 0.0,
            "upsample_mode": "linear",
        },
    }


def test_multiscale_local_skip_v2_builds_and_preserves_output_shapes():
    model = build_transchrombp_from_config(_base_config())
    seq = torch.randn(2, 2114, 4)
    out = model(seq)
    assert out.profile_logits_full.shape == (2, 1000)
    assert out.profile_logits_debiased.shape == (2, 1000)
    assert out.logcount_full.shape == (2, 1)
    assert out.logcount_debiased.shape == (2, 1)


def test_multiscale_local_skip_v2_keeps_count_pool_mode_center():
    model = build_transchrombp_from_config(_base_config())
    assert model.count_pool_mode == "center"
    assert model.profile_bias_stop_gradient is True


def test_unknown_profile_decoder_mode_raises():
    cfg = deepcopy(_base_config())
    cfg["profile_decoder"]["mode"] = "does_not_exist"
    with pytest.raises(ValueError, match="profile_decoder.mode"):
        build_transchrombp_from_config(cfg)
```

- [ ] **Step 2: Run the new tests and verify they fail**

Run:

```bash
PYTHONPATH=vendor/transchrombp:${PYTHONPATH:-} python -m pytest tests/test_multiscale_local_skip_v2.py -q
```

Expected:

- FAIL because `build_transchrombp_from_config` does not yet understand `profile_decoder.mode`, or because the decoder module does not exist.

- [ ] **Step 3: Implement the decoder module and wire it into `TransChromBP`**

Create `vendor/transchrombp/transchrombp/models/profile_decoder.py` with this implementation:

```python
from __future__ import annotations

from typing import Sequence

import torch
import torch.nn.functional as F
from torch import Tensor, nn

from .bias_branch import center_crop_1d


class ConvRefineBlock(nn.Module):
    def __init__(self, channels: int, dropout: float) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv1d(channels, channels, kernel_size=3, padding=1),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Conv1d(channels, channels, kernel_size=3, padding=1),
        )

    def forward(self, x: Tensor) -> Tensor:
        return x + self.net(x)


class MultiScaleLocalSkipDecoderV2(nn.Module):
    def __init__(
        self,
        d_model: int,
        output_len: int,
        hidden_channels: int,
        decoder_channels: Sequence[int],
        dropout: float = 0.1,
        upsample_mode: str = "linear",
    ) -> None:
        super().__init__()
        if not decoder_channels:
            raise ValueError("decoder_channels must be non-empty")
        if upsample_mode != "linear":
            raise ValueError(f"Unsupported upsample_mode={upsample_mode!r}")
        self.output_len = int(output_len)
        self.upsample_mode = upsample_mode
        self.input_proj = nn.Conv1d(d_model, hidden_channels, kernel_size=1)
        self.skip_proj = nn.Conv1d(d_model, hidden_channels, kernel_size=1)
        self.down_blocks = nn.ModuleList()
        self.up_blocks = nn.ModuleList()
        in_channels = hidden_channels
        for out_channels in decoder_channels:
            self.down_blocks.append(
                nn.Sequential(
                    nn.Conv1d(in_channels, out_channels, kernel_size=3, stride=2, padding=1),
                    nn.GELU(),
                    nn.Dropout(dropout),
                )
            )
            in_channels = out_channels
        for out_channels in reversed(decoder_channels[:-1] + [hidden_channels]):
            self.up_blocks.append(
                nn.Sequential(
                    nn.Conv1d(in_channels, out_channels, kernel_size=3, padding=1),
                    nn.GELU(),
                    nn.Dropout(dropout),
                    ConvRefineBlock(out_channels, dropout),
                )
            )
            in_channels = out_channels
        self.out_proj = nn.Conv1d(hidden_channels, 1, kernel_size=1)

    def forward(self, encoded_tokens: Tensor, local_features_cf: Tensor) -> Tensor:
        x = encoded_tokens.transpose(1, 2)
        x = self.input_proj(x)
        skip = self.skip_proj(local_features_cf)
        pyramids: list[Tensor] = []
        for block in self.down_blocks:
            x = block(x)
            pyramids.append(x)
        for idx, block in enumerate(self.up_blocks):
            target_len = skip.size(-1) if idx == len(self.up_blocks) - 1 else pyramids[-idx - 2].size(-1)
            x = F.interpolate(x, size=target_len, mode=self.upsample_mode, align_corners=False)
            if idx < len(self.up_blocks) - 1:
                x = x + pyramids[-idx - 2]
            else:
                x = x + skip
            x = block(x)
        profile = self.out_proj(x).squeeze(1)
        return center_crop_1d(profile, self.output_len)
```

Modify `vendor/transchrombp/transchrombp/models/transchrombp.py` in three places:

```python
from .profile_decoder import MultiScaleLocalSkipDecoderV2
```

```python
        profile_decoder_cfg: Optional[Dict[str, Any]] = None,
```

```python
        self.profile_decoder_mode = "linear"
        self.profile_decoder: Optional[nn.Module] = None
        if profile_decoder_cfg:
            self.profile_decoder_mode = str(profile_decoder_cfg.get("mode", "linear")).lower()
            if self.profile_decoder_mode == "multiscale_local_skip_v2":
                self.profile_decoder = MultiScaleLocalSkipDecoderV2(
                    d_model=d_model,
                    output_len=output_len,
                    hidden_channels=int(profile_decoder_cfg.get("hidden_channels", d_model)),
                    decoder_channels=list(profile_decoder_cfg.get("decoder_channels", [d_model, d_model // 2, d_model // 2])),
                    dropout=float(profile_decoder_cfg.get("dropout", dropout)),
                    upsample_mode=str(profile_decoder_cfg.get("upsample_mode", "linear")),
                )
            elif self.profile_decoder_mode != "linear":
                raise ValueError(
                    f"Unsupported profile_decoder.mode={self.profile_decoder_mode!r}; expected 'linear' or 'multiscale_local_skip_v2'"
                )
```

```python
        if self.profile_decoder is None:
            profile_signal = self.profile_signal_head(encoded).squeeze(-1)
            profile_signal = center_crop_1d(profile_signal, self.output_len)
        else:
            profile_signal = self.profile_decoder(encoded, local_feat)
```

Also pass the config through `build_transchrombp_from_config`:

```python
    profile_decoder_cfg = config.get("profile_decoder", {})
```

```python
        profile_decoder_cfg=profile_decoder_cfg,
```

Modify `vendor/transchrombp/transchrombp/models/__init__.py`:

```python
from .profile_decoder import MultiScaleLocalSkipDecoderV2
```

```python
    "MultiScaleLocalSkipDecoderV2",
```

- [ ] **Step 4: Run the decoder tests again and verify they pass**

Run:

```bash
PYTHONPATH=vendor/transchrombp:${PYTHONPATH:-} python -m pytest tests/test_multiscale_local_skip_v2.py -q
```

Expected:

- PASS with `3 passed`.

- [ ] **Step 5: Commit the decoder-path implementation**

Run:

```bash
git add \
  vendor/transchrombp/transchrombp/models/profile_decoder.py \
  vendor/transchrombp/transchrombp/models/transchrombp.py \
  vendor/transchrombp/transchrombp/models/__init__.py \
  tests/test_multiscale_local_skip_v2.py
git commit -m "feat: add multiscale local skip decoder v2"
```

### Task 2: Add the A6000 formal-gate config and launcher

**Files:**
- Create: `vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_msdls_v2.yaml`
- Create: `vendor/transchrombp/transchrombp/scripts/run_msdls_v2_gate.sh`

- [ ] **Step 1: Add the new model YAML**

Create `vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_msdls_v2.yaml`:

```yaml
# corrected B + multiscale/local-skip decoder v2
model_name: transchrombp_teacher_v2_center_pool_msdls_v2

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

profile_decoder:
  mode: multiscale_local_skip_v2
  hidden_channels: 256
  decoder_channels: [256, 192, 128]
  dropout: 0.1
  upsample_mode: linear
```

- [ ] **Step 2: Add the foreground gate launcher**

Create `vendor/transchrombp/transchrombp/scripts/run_msdls_v2_gate.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TRAIN_GPU_IDS="${TRAIN_GPU_IDS:-0,1}"
NPROC_PER_NODE="${NPROC_PER_NODE:-2}"
MASTER_ADDR="${MASTER_ADDR:-127.0.0.1}"
MASTER_PORT_BASE="${MASTER_PORT_BASE:-29920}"
VENV_DIR="${VENV_DIR:-/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b}"
OUTPUT_BASE="${OUTPUT_BASE:-${ROOT_DIR}/outputs}"
TRAIN_CONFIG="${TRAIN_CONFIG:-${ROOT_DIR}/configs/train/train_tutorial_teacher_v2_main.yaml}"
DATA_CONFIG="${DATA_CONFIG:-${ROOT_DIR}/configs/data/data_tutorial_canonical_v1.yaml}"
MODEL_CONFIG="${MODEL_CONFIG:-${ROOT_DIR}/configs/model/transchrombp_teacher_v2_center_pool_msdls_v2.yaml}"
RUN_NAME="${RUN_NAME:-teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1}"
SEED="${SEED:-42}"
BATCH_SIZE_PER_GPU="${BATCH_SIZE_PER_GPU:-12}"
NUM_WORKERS="${NUM_WORKERS:-4}"
PYTHON_BIN="${PYTHON_BIN:-python}"
DRY_RUN="${DRY_RUN:-0}"

export PATH="${VENV_DIR}/bin:${PATH:-}"
export PYTHONPATH="${ROOT_DIR}/src:${PYTHONPATH:-}"

RUNTIME_DIR="${OUTPUT_BASE}/runtime/msdls_v2_gate"
mkdir -p "${RUNTIME_DIR}"
RUNTIME_TRAIN="${RUNTIME_DIR}/train_${RUN_NAME}.yaml"

"${PYTHON_BIN}" - "${TRAIN_CONFIG}" "${RUNTIME_TRAIN}" "${SEED}" "${RUN_NAME}" "${BATCH_SIZE_PER_GPU}" "${NUM_WORKERS}" <<'PY'
import sys
from pathlib import Path
import yaml

src, dst, seed, run_name, batch_size, num_workers = sys.argv[1:7]
cfg = yaml.safe_load(Path(src).read_text(encoding="utf-8"))
cfg["seed"] = int(seed)
cfg.setdefault("data", {})["batch_size_per_gpu"] = int(batch_size)
cfg["data"]["num_workers"] = int(num_workers)
cfg.setdefault("logging", {})["run_name"] = run_name
Path(dst).write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
PY

echo "RUN_NAME=${RUN_NAME}"
echo "TRAIN_GPU_IDS=${TRAIN_GPU_IDS}"
echo "MASTER_PORT_BASE=${MASTER_PORT_BASE}"
echo "MODEL_CONFIG=${MODEL_CONFIG}"
echo "RUNTIME_TRAIN=${RUNTIME_TRAIN}"

if [[ "${DRY_RUN}" == "1" ]]; then
  exit 0
fi

export CUDA_VISIBLE_DEVICES="${TRAIN_GPU_IDS}"
torchrun \
  --nnodes=1 \
  --nproc_per_node="${NPROC_PER_NODE}" \
  --master_addr="${MASTER_ADDR}" \
  --master_port="${MASTER_PORT_BASE}" \
  -m transchrombp.training.train_ddp \
  --train-config "${RUNTIME_TRAIN}" \
  --model-config "${MODEL_CONFIG}" \
  --data-config "${DATA_CONFIG}" \
  --run-name "${RUN_NAME}" \
  --output-dir "${OUTPUT_BASE}" \
  --batch-size-per-gpu "${BATCH_SIZE_PER_GPU}"
```

- [ ] **Step 3: Run syntax and dry-run checks**

Run:

```bash
python -m py_compile vendor/transchrombp/transchrombp/models/profile_decoder.py \
  vendor/transchrombp/transchrombp/models/transchrombp.py
bash -n vendor/transchrombp/transchrombp/scripts/run_msdls_v2_gate.sh
DRY_RUN=1 bash vendor/transchrombp/transchrombp/scripts/run_msdls_v2_gate.sh
```

Expected:

- `py_compile` exits `0`
- `bash -n` exits `0`
- `DRY_RUN=1` prints `RUN_NAME`, `TRAIN_GPU_IDS`, `MASTER_PORT_BASE`, `MODEL_CONFIG`, `RUNTIME_TRAIN`

- [ ] **Step 4: Commit the gate config and launcher**

Run:

```bash
git add \
  vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_msdls_v2.yaml \
  vendor/transchrombp/transchrombp/scripts/run_msdls_v2_gate.sh
git commit -m "feat: add msdls v2 gate launcher"
```

### Task 3: Add the deterministic AlphaGenome v2 panel builder

**Files:**
- Create: `scripts/alphagenome_pilot/build_matched_panel_v2.py`
- Create: `tests/test_alphagenome_panel_builder.py`
- Create: `scripts/alphagenome_pilot/regions_k562_tutorial_matched_panel_v2.csv`

- [ ] **Step 1: Write the failing panel-selection test**

Create `tests/test_alphagenome_panel_builder.py`:

```python
from __future__ import annotations

from scripts.alphagenome_pilot.build_matched_panel_v2 import choose_panel_rows


def test_choose_panel_rows_returns_unique_quantile_anchors():
    rows = [
        {"label": f"peak_{i}", "source": "peak", "observed_total": float(i)}
        for i in range(1, 21)
    ] + [
        {"label": f"nonpeak_{i}", "source": "nonpeak", "observed_total": float(i)}
        for i in range(1, 21)
    ]
    panel = choose_panel_rows(rows, quantiles=[0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9, 0.98])
    assert len(panel) == 16
    assert len({row["label"] for row in panel}) == 16
    assert {row["source"] for row in panel} == {"peak", "nonpeak"}
```

- [ ] **Step 2: Run the test and verify it fails**

Run:

```bash
python -m pytest tests/test_alphagenome_panel_builder.py -q
```

Expected:

- FAIL because the builder script and `choose_panel_rows` do not exist yet.

- [ ] **Step 3: Implement the panel builder**

Create `scripts/alphagenome_pilot/build_matched_panel_v2.py`:

```python
#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pyBigWig
import yaml

from vendor.transchrombp.transchrombp.data.real_data import (
    filter_records_by_chroms,
    load_bigwig_chrom_sizes,
    load_fold_chroms,
    load_regions_from_bed,
)

PANEL_QUANTILES = [0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9, 0.98]
WINDOW_BP = 1000


def _window_total(bigwig: pyBigWig.pyBigWig, chrom: str, center: int) -> float:
    start = max(0, center - WINDOW_BP // 2)
    end = center + WINDOW_BP // 2
    values = np.nan_to_num(bigwig.values(chrom, start, end, numpy=True), nan=0.0)
    return float(values.sum())


def choose_panel_rows(rows: list[dict], quantiles: list[float]) -> list[dict]:
    panel: list[dict] = []
    used_labels: set[str] = set()
    for source in ("peak", "nonpeak"):
        source_rows = [row for row in rows if row["source"] == source]
        source_rows = sorted(source_rows, key=lambda row: row["observed_total"])
        totals = np.asarray([row["observed_total"] for row in source_rows], dtype=np.float64)
        for q in quantiles:
            target = float(np.quantile(totals, q))
            idx = min(range(len(source_rows)), key=lambda i: abs(source_rows[i]["observed_total"] - target))
            while source_rows[idx]["label"] in used_labels:
                idx = min(
                    (i for i in range(len(source_rows)) if source_rows[i]["label"] not in used_labels),
                    key=lambda i: abs(source_rows[i]["observed_total"] - target),
                )
                break
            row = dict(source_rows[idx])
            row["quantile"] = q
            panel.append(row)
            used_labels.add(row["label"])
    return panel


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-config", required=True)
    parser.add_argument("--split", default="test")
    parser.add_argument("--output-csv", required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(Path(args.data_config).read_text(encoding="utf-8"))
    split_chroms = load_fold_chroms(cfg["folds_json"], args.split)
    bigwig_path = cfg["input"]["bigwig"]
    peaks_bed = cfg["input"]["peaks_bed"]
    nonpeaks_bed = cfg["input"]["nonpeaks_bed"]
    chrom_sizes = load_bigwig_chrom_sizes(bigwig_path)
    peaks = filter_records_by_chroms(load_regions_from_bed(peaks_bed, split_chroms, source="peak"), chrom_sizes.keys(), peaks_bed, "peak")
    nonpeaks = filter_records_by_chroms(load_regions_from_bed(nonpeaks_bed, split_chroms, source="nonpeak"), chrom_sizes.keys(), nonpeaks_bed, "nonpeak")

    rows: list[dict] = []
    with pyBigWig.open(bigwig_path) as bw:
        for idx, record in enumerate(peaks + nonpeaks):
            rows.append(
                {
                    "label": f"{record.source}_{idx}",
                    "source": record.source,
                    "chrom": record.chrom,
                    "center": int(record.center),
                    "observed_total": _window_total(bw, record.chrom, record.center),
                    "notes": "tutorial canonical test split quantile anchor",
                }
            )

    panel = choose_panel_rows(rows, PANEL_QUANTILES)
    out = Path(args.output_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["label", "source", "quantile", "chrom", "center", "observed_total", "notes"])
        writer.writeheader()
        writer.writerows(panel)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run the unit test and generate the committed CSV**

Run:

```bash
python -m pytest tests/test_alphagenome_panel_builder.py -q
python scripts/alphagenome_pilot/build_matched_panel_v2.py \
  --data-config vendor/transchrombp/transchrombp/configs/data/data_tutorial_canonical_v1.yaml \
  --split test \
  --output-csv scripts/alphagenome_pilot/regions_k562_tutorial_matched_panel_v2.csv
wc -l scripts/alphagenome_pilot/regions_k562_tutorial_matched_panel_v2.csv
sed -n '1,10p' scripts/alphagenome_pilot/regions_k562_tutorial_matched_panel_v2.csv
```

Expected:

- `pytest` PASS
- `wc -l` should print `17` (header + 16 loci)
- CSV preview should show both `peak` and `nonpeak` rows with `quantile` column

- [ ] **Step 5: Commit the sidecar-panel assets**

Run:

```bash
git add \
  scripts/alphagenome_pilot/build_matched_panel_v2.py \
  scripts/alphagenome_pilot/regions_k562_tutorial_matched_panel_v2.csv \
  tests/test_alphagenome_panel_builder.py
git commit -m "feat: add alphagenome matched panel v2"
```

### Task 4: Deploy the dual-track assets to 6000 and launch both jobs

**Files:**
- Create: `scripts/deploy_a6000_dual_track_to_6000.sh`
- Modify: `TRACKING.md`
- Modify: `docs/plan/dual_machine_experiment_charter_20260409.md`

- [ ] **Step 1: Add the deployment helper**

Create `scripts/deploy_a6000_dual_track_to_6000.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REMOTE="zhoujiazhen@127.0.0.1"
PORT="6000"
REMOTE_ROOT="/data1/zhoujiazhen/bylw_atac/TransChromBP"

ssh -p "${PORT}" "${REMOTE}" \
  "mkdir -p \
    '${REMOTE_ROOT}/src/transchrombp/models' \
    '${REMOTE_ROOT}/configs/model' \
    '${REMOTE_ROOT}/scripts/alphagenome_pilot'"

scp -P "${PORT}" \
  "${REPO_ROOT}/vendor/transchrombp/transchrombp/models/profile_decoder.py" \
  "${REPO_ROOT}/vendor/transchrombp/transchrombp/models/transchrombp.py" \
  "${REPO_ROOT}/vendor/transchrombp/transchrombp/models/__init__.py" \
  "${REMOTE}:${REMOTE_ROOT}/src/transchrombp/models/"

scp -P "${PORT}" \
  "${REPO_ROOT}/vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_msdls_v2.yaml" \
  "${REMOTE}:${REMOTE_ROOT}/configs/model/"

scp -P "${PORT}" \
  "${REPO_ROOT}/vendor/transchrombp/transchrombp/scripts/run_msdls_v2_gate.sh" \
  "${REMOTE}:${REMOTE_ROOT}/scripts/"

scp -P "${PORT}" \
  "${REPO_ROOT}/scripts/alphagenome_pilot/run_alphagenome_pilot.py" \
  "${REPO_ROOT}/scripts/alphagenome_pilot/merge_locus_totals.py" \
  "${REPO_ROOT}/scripts/alphagenome_pilot/build_matched_panel_v2.py" \
  "${REPO_ROOT}/scripts/alphagenome_pilot/regions_k562_tutorial_matched_panel_v2.csv" \
  "${REMOTE}:${REMOTE_ROOT}/scripts/alphagenome_pilot/"
```

- [ ] **Step 2: Deploy and verify remote syntax**

Run:

```bash
bash scripts/deploy_a6000_dual_track_to_6000.sh
ssh -p 6000 zhoujiazhen@127.0.0.1 \
  "cd /data1/zhoujiazhen/bylw_atac/TransChromBP && \
   /data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b/bin/python -m py_compile \
     src/transchrombp/models/profile_decoder.py \
     src/transchrombp/models/transchrombp.py && \
   bash -n scripts/run_msdls_v2_gate.sh && \
   /data1/zhoujiazhen/bylw_atac/.mamba/envs/alphagenome/bin/python -m py_compile \
     scripts/alphagenome_pilot/build_matched_panel_v2.py \
     scripts/alphagenome_pilot/run_alphagenome_pilot.py \
     scripts/alphagenome_pilot/merge_locus_totals.py"
```

Expected:

- All syntax checks exit `0`.

- [ ] **Step 3: Launch the A6000 formal gate**

Run:

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 \
  "cd /data1/zhoujiazhen/bylw_atac/TransChromBP && \
   nohup bash -lc '
     TRAIN_GPU_IDS=0,1 \
     NPROC_PER_NODE=2 \
     MASTER_PORT_BASE=29920 \
     RUN_NAME=teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1 \
     bash scripts/run_msdls_v2_gate.sh
   ' > logs/teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1.log 2>&1 & \
   echo \$!"
```

Expected:

- Prints one PID
- `nvidia-smi` should show both A6000 cards engaged within the first few minutes

- [ ] **Step 4: Launch the AlphaGenome sidecar**

Run:

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 \
  "cd /data1/zhoujiazhen/bylw_atac/TransChromBP && \
   nohup bash -lc '
     /data1/zhoujiazhen/bylw_atac/.mamba/envs/alphagenome/bin/python \
       scripts/alphagenome_pilot/run_alphagenome_pilot.py \
       --regions-csv scripts/alphagenome_pilot/regions_k562_tutorial_matched_panel_v2.csv \
       --output-dir outputs/alphagenome_pilot/alphagenome_matched_raw_track_slice_v2_20260410 && \
     /data1/zhoujiazhen/bylw_atac/.mamba/envs/alphagenome/bin/python \
       scripts/alphagenome_pilot/merge_locus_totals.py \
       --alpha-summary outputs/alphagenome_pilot/alphagenome_matched_raw_track_slice_v2_20260410/summary.csv \
       --local-totals scripts/alphagenome_pilot/regions_k562_tutorial_matched_panel_v2.csv \
       --output-csv outputs/alphagenome_pilot/alphagenome_matched_raw_track_slice_v2_20260410/merged_locus_totals.csv
   ' > logs/alphagenome_matched_raw_track_slice_v2_20260410.log 2>&1 & \
   echo \$!"
```

Expected:

- Prints one PID
- This sidecar is allowed to finish quickly and does not need A6000 utilization
- `merged_locus_totals.csv` will join the AlphaGenome summary against the committed panel metadata and `observed_total`, so the expanded `12-20 loci` run does not depend on the old 4-locus totals table.

- [ ] **Step 5: Update the live docs immediately after both launches are verified**

Edit `TRACKING.md` and `docs/plan/dual_machine_experiment_charter_20260409.md` so they record:

```md
- `6000` active GPU run: `teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1`
- `AlphaGenome` active sidecar: `alphagenome_matched_raw_track_slice_v2_20260410`
- 两条线互不等待
- A6000 gate 只和 `corrected B` 比
- AlphaGenome 仍是窄 external coordinate
```

Use these monitoring commands in the doc and user update:

```bash
tail -f /data1/zhoujiazhen/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1.log
tail -f /data1/zhoujiazhen/bylw_atac/TransChromBP/logs/alphagenome_matched_raw_track_slice_v2_20260410.log
```

Write the A6000 ETA as an initial rough window of `+5h` to `+8h` from launch, and tighten it after the first completed validation epoch appears in the log.

- [ ] **Step 6: Commit the launch-state doc update**

Run:

```bash
git add \
  scripts/deploy_a6000_dual_track_to_6000.sh \
  TRACKING.md \
  docs/plan/dual_machine_experiment_charter_20260409.md
git commit -m "docs: track dual a6000 launch state"
```

### Task 5: Close out both tracks and archive terminal items

**Files:**
- Create: `reports/closeout/multiscale_local_skip_v2_30ep_gate_closeout_20260410.md`
- Create: `reports/closeout/alphagenome_matched_raw_track_slice_v2_closeout_20260410.md`
- Modify: `TRACKING.md`
- Modify: `TRACKING_archive.md`
- Modify: `docs/plan/dual_machine_experiment_charter_20260409.md`

- [ ] **Step 1: Re-check the finished A6000 run and classify the gate**

Run:

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 \
  "cd /data1/zhoujiazhen/bylw_atac/TransChromBP && \
   tail -n 80 logs/teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1.log && \
   echo '---META---' && \
   sed -n '1,80p' outputs/logs/teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1/run_meta.json"
```

Expected:

- You can read the best epoch, best validation metrics, and terminal stop condition.

- [ ] **Step 2: Write the A6000 closeout report**

Create `reports/closeout/multiscale_local_skip_v2_30ep_gate_closeout_20260410.md` with this skeleton:

```md
# Multiscale Local-Skip Decoder v2 30-Epoch Gate Closeout (2026-04-10)

## Verdict
- `pass` / `fail`

## Run Identity
- run name
- machine / GPU
- log path
- launch / finish time

## Baseline Comparator
- `corrected B`
- comparator metric lines

## Gate Evidence
- peak profile JSD delta
- peak count_r delta
- full/debiased gap behavior
- whether training was stable

## Decision
- if `pass`: only second-seed or fuller budget allowed
- if `fail`: stop the v2 family; do not branch into `v2a/v2b/v2c`
```

- [ ] **Step 3: Re-check the finished AlphaGenome sidecar**

Run:

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 \
  "cd /data1/zhoujiazhen/bylw_atac/TransChromBP && \
   tail -n 60 logs/alphagenome_matched_raw_track_slice_v2_20260410.log && \
   echo '---FILES---' && \
   ls outputs/alphagenome_pilot/alphagenome_matched_raw_track_slice_v2_20260410"
```

Expected:

- You can confirm the larger sidecar still produced `summary`, `metadata`, `profiles`, `run_meta`, and `merged`.

- [ ] **Step 4: Write the AlphaGenome v2 closeout report**

Create `reports/closeout/alphagenome_matched_raw_track_slice_v2_closeout_20260410.md` with this skeleton:

```md
# AlphaGenome Matched Raw-Track Slice v2 Closeout (2026-04-10)

## Verdict
- `pass` / `fail`

## Run Identity
- run name
- panel size
- log path
- output dir

## Evidence
- number of loci completed
- number of loci with usable ATAC tracks after filtering
- whether merged table formed

## Decision
- if `pass`: sidecar remains an external coordinate only
- if `fail`: stop the v2 panel and document the blocker
```

- [ ] **Step 5: Update live docs and archive terminal items**

Apply these rules:

```md
- A6000 run remains in `TRACKING.md` only if a follow-up seed/budget is actually pending
- AlphaGenome sidecar row leaves `TRACKING.md` once its closeout is terminal and the top summary row already carries the conclusion
- `TRACKING_archive.md` gets one completed-item row per terminal track
```

- [ ] **Step 6: Commit the closeout docs**

Run:

```bash
git add \
  reports/closeout/multiscale_local_skip_v2_30ep_gate_closeout_20260410.md \
  reports/closeout/alphagenome_matched_raw_track_slice_v2_closeout_20260410.md \
  TRACKING.md \
  TRACKING_archive.md \
  docs/plan/dual_machine_experiment_charter_20260409.md
git commit -m "docs: close out dual a6000 experiment tracks"
```
