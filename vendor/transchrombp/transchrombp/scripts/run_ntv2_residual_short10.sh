#!/usr/bin/env bash
# NT v2 cached residual-head pilot launcher.
# Sequence:
#   1. build dataset-aligned NT v2 cache (train + valid)
#   2. optionally run unified valid probe against baseline checkpoint
#   3. run corrected-B + NT v2 residual-head short10 training

set -euo pipefail

if [ -z "${BASH_VERSION:-}" ]; then
    echo "[error] Please run this script with bash." >&2
    exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TRAIN_GPU_IDS="${TRAIN_GPU_IDS:-${GPU_IDS:-${1:-${GPU_ID:-0,1}}}}"
CACHE_GPU_IDS="${CACHE_GPU_IDS:-${TRAIN_GPU_IDS}}"
CACHE_GPU_ID="${CACHE_GPU_ID:-${CACHE_GPU_IDS%%,*}}"
NPROC_PER_NODE="${NPROC_PER_NODE:-2}"
CACHE_NPROC_PER_NODE="${CACHE_NPROC_PER_NODE:-${NPROC_PER_NODE}}"
MASTER_ADDR="${MASTER_ADDR:-127.0.0.1}"
MASTER_PORT_BASE="${MASTER_PORT_BASE:-29820}"
VENV_DIR="${VENV_DIR:-/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b}"
OUTPUT_BASE="${OUTPUT_BASE:-${ROOT_DIR}/outputs}"
FOUNDATION_CACHE_DIR="${FOUNDATION_CACHE_DIR:-${OUTPUT_BASE}/foundation_cache/ntv2_tutorial_canonical_v1}"
NT_MODEL_DIR="${NT_MODEL_DIR:-/data1/zhoujiazhen/bylw_atac/foundation_models/nucleotide_transformer/nt-v2-500m-multi-species}"
TRAIN_CONFIG="${TRAIN_CONFIG:-${ROOT_DIR}/configs/train/train_tutorial_foundation_short10.yaml}"
DATA_CONFIG="${DATA_CONFIG:-${ROOT_DIR}/configs/data/data_tutorial_canonical_v1.yaml}"
MODEL_CONFIG="${MODEL_CONFIG:-${ROOT_DIR}/configs/model/transchrombp_teacher_v2_center_pool_ntv2_residual.yaml}"
BASELINE_CHECKPOINT="${BASELINE_CHECKPOINT:-}"
BASELINE_MODEL_CONFIG="${BASELINE_MODEL_CONFIG:-${ROOT_DIR}/configs/model/transchrombp_teacher_v2_center_pool.yaml}"
PYTHON_BIN="${PYTHON_BIN:-python}"

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

export PATH="${VENV_DIR}/bin:${PATH:-}"
if [ -d "${ROOT_DIR}/src/transchrombp" ]; then
    PACKAGE_IMPORT_ROOT="${ROOT_DIR}/src"
elif [ -d "${ROOT_DIR}/../transchrombp" ]; then
    PACKAGE_IMPORT_ROOT="$(cd "${ROOT_DIR}/.." && pwd)"
else
    echo "[error] Could not locate transchrombp package import root from ${ROOT_DIR}" >&2
    exit 1
fi
export PYTHONPATH="${PACKAGE_IMPORT_ROOT}:${PYTHONPATH:-}"

require_file() {
    local path="$1"
    local label="$2"
    if [ ! -f "${path}" ]; then
        echo "[error] Missing ${label}: ${path}" >&2
        exit 1
    fi
}

require_dir() {
    local path="$1"
    local label="$2"
    if [ ! -d "${path}" ]; then
        echo "[error] Missing ${label}: ${path}" >&2
        exit 1
    fi
}

count_csv_items() {
    local csv="$1"
    if [ -z "${csv}" ]; then
        echo 0
        return 0
    fi
    awk -F',' '{print NF}' <<<"${csv}"
}

launch_train() {
    local master_port="$1"
    shift
    export CUDA_VISIBLE_DEVICES="${TRAIN_GPU_IDS}"
    if [ "${NPROC_PER_NODE}" -gt 1 ]; then
        torchrun \
            --nnodes=1 \
            --nproc_per_node="${NPROC_PER_NODE}" \
            --master_addr="${MASTER_ADDR}" \
            --master_port="${master_port}" \
            -m transchrombp.training.train_ddp \
            "$@"
    else
        "${PYTHON_BIN}" -m transchrombp.training.train_ddp "$@"
    fi
}

require_file "${TRAIN_CONFIG}" "train config"
require_file "${DATA_CONFIG}" "data config"
require_file "${MODEL_CONFIG}" "model config"
require_dir "${NT_MODEL_DIR}" "NT v2 model directory"

TRAIN_GPU_COUNT="$(count_csv_items "${TRAIN_GPU_IDS}")"
if [ "${NPROC_PER_NODE}" -gt "${TRAIN_GPU_COUNT}" ]; then
    echo "[error] NPROC_PER_NODE=${NPROC_PER_NODE} exceeds visible TRAIN_GPU_IDS=${TRAIN_GPU_IDS}" >&2
    exit 1
fi
CACHE_GPU_COUNT="$(count_csv_items "${CACHE_GPU_IDS}")"
if [ "${CACHE_NPROC_PER_NODE}" -gt "${CACHE_GPU_COUNT}" ]; then
    echo "[error] CACHE_NPROC_PER_NODE=${CACHE_NPROC_PER_NODE} exceeds visible CACHE_GPU_IDS=${CACHE_GPU_IDS}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_BASE}" "${FOUNDATION_CACHE_DIR}"

RUNTIME_DIR="${OUTPUT_BASE}/runtime/ntv2_residual_short10"
mkdir -p "${RUNTIME_DIR}"
EFFECTIVE_TRAIN_CONFIG="${RUNTIME_DIR}/train_tutorial_foundation_short10_runtime.yaml"
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
echo "NT v2 Residual Short10"
echo "ROOT_DIR: ${ROOT_DIR}"
echo "CACHE_GPU_IDS: ${CACHE_GPU_IDS}"
echo "CACHE_GPU_ID: ${CACHE_GPU_ID} (used only when CACHE_NPROC_PER_NODE=1)"
echo "CACHE_NPROC_PER_NODE: ${CACHE_NPROC_PER_NODE}"
echo "TRAIN_GPU_IDS: ${TRAIN_GPU_IDS} (preferred dual-GPU training)"
echo "NPROC_PER_NODE: ${NPROC_PER_NODE}"
echo "MASTER_PORT_BASE: ${MASTER_PORT_BASE}"
echo "Cache dir: ${FOUNDATION_CACHE_DIR}"
echo "NT model dir: ${NT_MODEL_DIR}"
echo "Train config: ${TRAIN_CONFIG}"
echo "Effective train config: ${EFFECTIVE_TRAIN_CONFIG}"
echo "Model config: ${MODEL_CONFIG}"
echo "Data config: ${DATA_CONFIG}"
echo "BATCH_SIZE_PER_GPU: ${BATCH_SIZE_PER_GPU}"
echo "NUM_WORKERS: ${NUM_WORKERS}"
echo "GLOBAL_BATCH: $((BATCH_SIZE_PER_GPU * NPROC_PER_NODE))"
echo "=============================================="

if [ ! -f "${FOUNDATION_CACHE_DIR}/manifest_train.json" ] || [ ! -f "${FOUNDATION_CACHE_DIR}/manifest_valid.json" ]; then
    echo ""
    if [ "${CACHE_NPROC_PER_NODE}" -gt 1 ]; then
        echo "[Step 0] Building dataset-aligned NT v2 cache on CACHE_GPU_IDS=${CACHE_GPU_IDS} with ${CACHE_NPROC_PER_NODE} workers..."
        CUDA_VISIBLE_DEVICES="${CACHE_GPU_IDS}" \
        torchrun \
            --nnodes=1 \
            --nproc_per_node="${CACHE_NPROC_PER_NODE}" \
            --master_addr="${MASTER_ADDR}" \
            --master_port="$((MASTER_PORT_BASE - 1))" \
            "${ROOT_DIR}/scripts/build_foundation_cache.py" \
            --data_config "${DATA_CONFIG}" \
            --train_config "${EFFECTIVE_TRAIN_CONFIG}" \
            --output_dir "${FOUNDATION_CACHE_DIR}" \
            --backend nt_v2 \
            --model_dir "${NT_MODEL_DIR}" \
            --splits train valid \
            --layers 7,14 \
            --feature_types global_mean bins4_mean \
            --record_splits valid \
            --batch_size 8 \
            --dtype float16
    else
        echo "[Step 0] Building dataset-aligned NT v2 cache on single GPU ${CACHE_GPU_ID}..."
        CUDA_VISIBLE_DEVICES="${CACHE_GPU_ID}" "${PYTHON_BIN}" "${ROOT_DIR}/scripts/build_foundation_cache.py" \
            --data_config "${DATA_CONFIG}" \
            --train_config "${EFFECTIVE_TRAIN_CONFIG}" \
            --output_dir "${FOUNDATION_CACHE_DIR}" \
            --backend nt_v2 \
            --model_dir "${NT_MODEL_DIR}" \
            --splits train valid \
            --layers 7,14 \
            --feature_types global_mean bins4_mean \
            --record_splits valid \
            --batch_size 8 \
            --dtype float16
    fi
else
    echo "[Step 0] Cache already exists, skipping."
fi

if [ -n "${BASELINE_CHECKPOINT}" ] && [ -f "${BASELINE_CHECKPOINT}" ] && [ -f "${FOUNDATION_CACHE_DIR}/records_valid.jsonl" ]; then
    echo ""
    echo "[Step 1] Running unified valid probe..."
    "${PYTHON_BIN}" "${ROOT_DIR}/../../scripts/foundation_model_probe.py" \
        --records-jsonl "${FOUNDATION_CACHE_DIR}/records_valid.jsonl" \
        --cache-dir "${FOUNDATION_CACHE_DIR}" \
        --split valid \
        --baseline-checkpoint "${BASELINE_CHECKPOINT}" \
        --baseline-model-config "${BASELINE_MODEL_CONFIG}" \
        --transchrombp-root "${ROOT_DIR}" \
        --output-dir "${OUTPUT_BASE}/foundation_probes/ntv2_tutorial_valid"
else
    echo "[Step 1] Probe skipped. Set BASELINE_CHECKPOINT and ensure records_valid.jsonl exists."
fi

echo ""
echo "[Step 2] Launching residual-head short10 training on TRAIN_GPU_IDS=${TRAIN_GPU_IDS}..."
RUN_NAME="${RUN_NAME:-ntv2_residual_short10_s42}"
launch_train "${MASTER_PORT_BASE}" \
    --train-config "${EFFECTIVE_TRAIN_CONFIG}" \
    --model-config "${MODEL_CONFIG}" \
    --data-config "${DATA_CONFIG}" \
    --run-name "${RUN_NAME}" \
    --output-dir "${OUTPUT_BASE}" \
    --foundation-cache-dir "${FOUNDATION_CACHE_DIR}" \
    --batch-size-per-gpu "${BATCH_SIZE_PER_GPU}"

echo ""
echo "=============================================="
echo "NT v2 residual short10 completed."
echo "Results in: ${OUTPUT_BASE}/checkpoints, ${OUTPUT_BASE}/logs"
echo "=============================================="
