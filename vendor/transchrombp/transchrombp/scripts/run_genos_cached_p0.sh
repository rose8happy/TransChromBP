#!/usr/bin/env bash
# Genos cached-fusion P0 launcher
# 单独运行 matched-control baseline_short10；不加载 Genos cache
# 默认支持双卡 DDP matched recipe：
#   - nproc_per_node=2
#   - batch_size_per_gpu=10
#   - global_batch=20

set -euo pipefail

if [ -z "${BASH_VERSION:-}" ]; then
    echo "[error] Please run this script with bash." >&2
    exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GPU_IDS="${GPU_IDS:-${1:-${GPU_ID:-0,1}}}"
NPROC_PER_NODE="${NPROC_PER_NODE:-2}"
MASTER_ADDR="${MASTER_ADDR:-127.0.0.1}"
MASTER_PORT="${MASTER_PORT:-29731}"
VENV_DIR="${VENV_DIR:-/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b}"
OUTPUT_BASE="${OUTPUT_BASE:-${ROOT_DIR}/outputs}"
TRAIN_CONFIG="${TRAIN_CONFIG:-${ROOT_DIR}/configs/train/train_genos_cached_short10.yaml}"
MODEL_CONFIG="${MODEL_CONFIG:-${ROOT_DIR}/configs/model/v2fix_baseline.yaml}"
DATA_CONFIG="${DATA_CONFIG:-${ROOT_DIR}/configs/data/data_tutorial_canonical_v1.yaml}"
PYTHON_BIN="${PYTHON_BIN:-python}"
RUN_NAME="${RUN_NAME:-genos_cached_P0_baseline_$(date +%Y%m%d_%H%M%S)}"
if [ -z "${BATCH_SIZE_PER_GPU:-}" ]; then
    if [ "${NPROC_PER_NODE}" -gt 1 ]; then
        BATCH_SIZE_PER_GPU=10
    else
        BATCH_SIZE_PER_GPU=20
    fi
fi
if [ -z "${NUM_WORKERS:-}" ]; then
    if [ "${NPROC_PER_NODE}" -gt 1 ]; then
        NUM_WORKERS=4
    else
        NUM_WORKERS=2
    fi
fi
if [ -z "${DDP_FIND_UNUSED_PARAMETERS:-}" ]; then
    DDP_FIND_UNUSED_PARAMETERS=false
fi

export CUDA_VISIBLE_DEVICES="${GPU_IDS}"
export PATH="${VENV_DIR}/bin:${PATH:-}"
export PYTHONPATH="${ROOT_DIR}/src:${PYTHONPATH:-}"
export PYTHONUNBUFFERED=1

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

mkdir -p "${OUTPUT_BASE}"

GLOBAL_BATCH=$((BATCH_SIZE_PER_GPU * NPROC_PER_NODE))
RUNTIME_DIR="${OUTPUT_BASE}/runtime/${RUN_NAME}"
mkdir -p "${RUNTIME_DIR}"
EFFECTIVE_TRAIN_CONFIG="${RUNTIME_DIR}/train_config_runtime.yaml"
"${PYTHON_BIN}" - "${TRAIN_CONFIG}" "${EFFECTIVE_TRAIN_CONFIG}" "${BATCH_SIZE_PER_GPU}" "${NUM_WORKERS}" "${DDP_FIND_UNUSED_PARAMETERS}" <<'PY'
import sys
from pathlib import Path
import yaml

src, dst, batch_size, num_workers, find_unused = sys.argv[1:6]
with open(src, "r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
cfg.setdefault("data", {})["batch_size_per_gpu"] = int(batch_size)
cfg["data"]["num_workers"] = int(num_workers)
cfg.setdefault("trainer", {})["find_unused_parameters"] = str(find_unused).lower() == "true"
Path(dst).write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
PY

echo "=============================================="
echo "Genos Cached P0 Baseline"
echo "ROOT_DIR: ${ROOT_DIR}"
echo "GPU_IDS: ${GPU_IDS}"
echo "NPROC_PER_NODE: ${NPROC_PER_NODE}"
echo "MASTER_ADDR: ${MASTER_ADDR}"
echo "MASTER_PORT: ${MASTER_PORT}"
echo "RUN_NAME: ${RUN_NAME}"
echo "TRAIN_CONFIG: ${TRAIN_CONFIG}"
echo "EFFECTIVE_TRAIN_CONFIG: ${EFFECTIVE_TRAIN_CONFIG}"
echo "MODEL_CONFIG: ${MODEL_CONFIG}"
echo "DATA_CONFIG: ${DATA_CONFIG}"
echo "BATCH_SIZE_PER_GPU: ${BATCH_SIZE_PER_GPU}"
echo "NUM_WORKERS: ${NUM_WORKERS}"
echo "DDP_FIND_UNUSED_PARAMETERS: ${DDP_FIND_UNUSED_PARAMETERS}"
echo "GLOBAL_BATCH: ${GLOBAL_BATCH}"
echo "=============================================="

if [ "${NPROC_PER_NODE}" -gt 1 ]; then
    if ! command -v torchrun >/dev/null 2>&1; then
        echo "[error] torchrun not found in PATH" >&2
        exit 1
    fi
    torchrun \
        --nnodes=1 \
        --nproc_per_node="${NPROC_PER_NODE}" \
        --master_addr="${MASTER_ADDR}" \
        --master_port="${MASTER_PORT}" \
        -m transchrombp.training.train_ddp \
        --train-config "${EFFECTIVE_TRAIN_CONFIG}" \
        --model-config "${MODEL_CONFIG}" \
        --data-config "${DATA_CONFIG}" \
        --run-name "${RUN_NAME}" \
        --output-dir "${OUTPUT_BASE}" \
        --batch-size-per-gpu "${BATCH_SIZE_PER_GPU}"
else
    "${PYTHON_BIN}" -m transchrombp.training.train_ddp \
        --train-config "${EFFECTIVE_TRAIN_CONFIG}" \
        --model-config "${MODEL_CONFIG}" \
        --data-config "${DATA_CONFIG}" \
        --run-name "${RUN_NAME}" \
        --output-dir "${OUTPUT_BASE}" \
        --batch-size-per-gpu "${BATCH_SIZE_PER_GPU}"
fi
