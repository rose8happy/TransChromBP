#!/usr/bin/env bash
# Corrected-B + multiscale decoder probe launcher.
# Default target is the tutorial-canonical full-budget s42 run.
# Override TRAIN_CONFIG to a short-budget config when you only want a structure smoke.

set -euo pipefail

if [ -z "${BASH_VERSION:-}" ]; then
    echo "[error] Please run this script with bash." >&2
    exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TRAIN_GPU_IDS="${TRAIN_GPU_IDS:-${GPU_IDS:-${1:-${GPU_ID:-0,1}}}}"
NPROC_PER_NODE="${NPROC_PER_NODE:-2}"
MASTER_ADDR="${MASTER_ADDR:-127.0.0.1}"
MASTER_PORT_BASE="${MASTER_PORT_BASE:-29920}"
VENV_DIR="${VENV_DIR:-/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b}"
OUTPUT_BASE="${OUTPUT_BASE:-${ROOT_DIR}/outputs}"
TRAIN_CONFIG="${TRAIN_CONFIG:-${ROOT_DIR}/configs/train/train_tutorial_teacher_v2_main.yaml}"
DATA_CONFIG="${DATA_CONFIG:-${ROOT_DIR}/configs/data/data_tutorial_canonical_v1.yaml}"
MODEL_CONFIG="${MODEL_CONFIG:-${ROOT_DIR}/configs/model/transchrombp_teacher_v2_center_pool_msdec_v1.yaml}"
PYTHON_BIN="${PYTHON_BIN:-python}"
TRAIN_SEED="${TRAIN_SEED:-42}"
RUN_NAME="${RUN_NAME:-teacher_v2_center_pool_msdec_v1_s${TRAIN_SEED}}"

if [ -z "${BATCH_SIZE_PER_GPU:-}" ]; then
    if [ "${NPROC_PER_NODE}" -gt 1 ]; then
        BATCH_SIZE_PER_GPU=12
    else
        BATCH_SIZE_PER_GPU=24
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

TRAIN_GPU_COUNT="$(count_csv_items "${TRAIN_GPU_IDS}")"
if [ "${NPROC_PER_NODE}" -gt "${TRAIN_GPU_COUNT}" ]; then
    echo "[error] NPROC_PER_NODE=${NPROC_PER_NODE} exceeds visible TRAIN_GPU_IDS=${TRAIN_GPU_IDS}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_BASE}"

RUNTIME_DIR="${OUTPUT_BASE}/runtime/multiscale_decoder_probe"
mkdir -p "${RUNTIME_DIR}"
TRAIN_CONFIG_BASENAME="$(basename "${TRAIN_CONFIG}" .yaml)"
EFFECTIVE_TRAIN_CONFIG="${RUNTIME_DIR}/${TRAIN_CONFIG_BASENAME}_${RUN_NAME}_runtime.yaml"
"${PYTHON_BIN}" - "${TRAIN_CONFIG}" "${EFFECTIVE_TRAIN_CONFIG}" "${BATCH_SIZE_PER_GPU}" "${NUM_WORKERS}" "${DDP_FIND_UNUSED_PARAMETERS}" "${TRAIN_SEED}" "${RUN_NAME}" <<'PY'
import sys
from pathlib import Path
import yaml

src, dst, batch_size, num_workers, find_unused, seed, run_name = sys.argv[1:8]
with open(src, "r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
cfg["seed"] = int(seed)
cfg.setdefault("data", {})["batch_size_per_gpu"] = int(batch_size)
cfg["data"]["num_workers"] = int(num_workers)
cfg.setdefault("trainer", {})["find_unused_parameters"] = str(find_unused).lower() == "true"
cfg.setdefault("logging", {})["run_name"] = run_name
Path(dst).write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
PY

echo "=============================================="
echo "Corrected-B + Multiscale Decoder Probe"
echo "ROOT_DIR: ${ROOT_DIR}"
echo "TRAIN_GPU_IDS: ${TRAIN_GPU_IDS}"
echo "NPROC_PER_NODE: ${NPROC_PER_NODE}"
echo "MASTER_PORT_BASE: ${MASTER_PORT_BASE}"
echo "Train config: ${TRAIN_CONFIG}"
echo "Effective train config: ${EFFECTIVE_TRAIN_CONFIG}"
echo "Model config: ${MODEL_CONFIG}"
echo "Data config: ${DATA_CONFIG}"
echo "RUN_NAME: ${RUN_NAME}"
echo "TRAIN_SEED: ${TRAIN_SEED}"
echo "BATCH_SIZE_PER_GPU: ${BATCH_SIZE_PER_GPU}"
echo "NUM_WORKERS: ${NUM_WORKERS}"
echo "GLOBAL_BATCH: $((BATCH_SIZE_PER_GPU * NPROC_PER_NODE))"
echo "DDP_FIND_UNUSED_PARAMETERS: ${DDP_FIND_UNUSED_PARAMETERS}"
echo "=============================================="

echo ""
echo "[Step 1] Launching multiscale decoder probe on TRAIN_GPU_IDS=${TRAIN_GPU_IDS}..."
launch_train "${MASTER_PORT_BASE}" \
    --train-config "${EFFECTIVE_TRAIN_CONFIG}" \
    --model-config "${MODEL_CONFIG}" \
    --data-config "${DATA_CONFIG}" \
    --run-name "${RUN_NAME}" \
    --output-dir "${OUTPUT_BASE}" \
    --batch-size-per-gpu "${BATCH_SIZE_PER_GPU}"

echo ""
echo "=============================================="
echo "Multiscale decoder probe completed."
echo "Results in: ${OUTPUT_BASE}/checkpoints, ${OUTPUT_BASE}/logs"
echo "=============================================="
