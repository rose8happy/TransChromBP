#!/usr/bin/env bash
# NT v2 teacher-distill launcher.
# Sequence:
#   1. build dataset-aligned NT v2 bins16 cache
#   2. export record-aligned teacher cache (profile16 + logcount)
#   3. train corrected-B student with distillation losses

set -euo pipefail

if [ -z "${BASH_VERSION:-}" ]; then
    echo "[error] Please run this script with bash." >&2
    exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TRAIN_GPU_IDS="${TRAIN_GPU_IDS:-${GPU_IDS:-${1:-${GPU_ID:-0,1}}}}"
CACHE_GPU_IDS="${CACHE_GPU_IDS:-${TRAIN_GPU_IDS}}"
NPROC_PER_NODE="${NPROC_PER_NODE:-2}"
CACHE_NPROC_PER_NODE="${CACHE_NPROC_PER_NODE:-${NPROC_PER_NODE}}"
MASTER_ADDR="${MASTER_ADDR:-127.0.0.1}"
MASTER_PORT_BASE="${MASTER_PORT_BASE:-29980}"
VENV_DIR="${VENV_DIR:-/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b}"
OUTPUT_BASE="${OUTPUT_BASE:-${ROOT_DIR}/outputs}"
FOUNDATION_CACHE_DIR="${FOUNDATION_CACHE_DIR:-${OUTPUT_BASE}/foundation_cache/ntv2_tutorial_canonical_v1_bins16}"
TEACHER_CACHE_DIR="${TEACHER_CACHE_DIR:-${OUTPUT_BASE}/teacher_cache/ntv2_tutorial_teacher_v1}"
NT_MODEL_DIR="${NT_MODEL_DIR:-/data1/zhoujiazhen/bylw_atac/foundation_models/nucleotide_transformer/nt-v2-500m-multi-species}"
TRAIN_CONFIG="${TRAIN_CONFIG:-${ROOT_DIR}/configs/train/train_tutorial_teacher_v2_ntv2_distill_short10.yaml}"
FULL_TRAIN_CONFIG_DEFAULT="${ROOT_DIR}/configs/train/train_tutorial_teacher_v2_ntv2_distill_full.yaml"
DATA_CONFIG="${DATA_CONFIG:-${ROOT_DIR}/configs/data/data_tutorial_canonical_v1.yaml}"
MODEL_CONFIG="${MODEL_CONFIG:-${ROOT_DIR}/configs/model/transchrombp_teacher_v2_center_pool_ntv2_distill.yaml}"
FEATURE_NAME="${FEATURE_NAME:-layer_07__bins16_mean}"
CACHE_SPLITS="${CACHE_SPLITS:-train valid test}"
PREDICT_SPLITS="${PREDICT_SPLITS:-train valid test}"
PROFILE_PROBE_BINS="${PROFILE_PROBE_BINS:-16}"
RIDGE_ALPHA="${RIDGE_ALPHA:-1.0}"
PYTHON_BIN="${PYTHON_BIN:-python}"
TRAIN_SEED="${TRAIN_SEED:-42}"
RUN_NAME="${RUN_NAME:-ntv2_teacher_distill_s${TRAIN_SEED}}"

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

read -r -a CACHE_SPLIT_ARR <<<"${CACHE_SPLITS}"
read -r -a PREDICT_SPLIT_ARR <<<"${PREDICT_SPLITS}"

require_file "${TRAIN_CONFIG}" "train config"
require_file "${FULL_TRAIN_CONFIG_DEFAULT}" "full train config"
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

mkdir -p "${OUTPUT_BASE}" "${FOUNDATION_CACHE_DIR}" "${TEACHER_CACHE_DIR}"

RUNTIME_DIR="${OUTPUT_BASE}/runtime/ntv2_teacher_distill"
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
echo "NT v2 Teacher Distill"
echo "ROOT_DIR: ${ROOT_DIR}"
echo "CACHE_GPU_IDS: ${CACHE_GPU_IDS}"
echo "TRAIN_GPU_IDS: ${TRAIN_GPU_IDS}"
echo "NPROC_PER_NODE: ${NPROC_PER_NODE}"
echo "CACHE_NPROC_PER_NODE: ${CACHE_NPROC_PER_NODE}"
echo "MASTER_PORT_BASE: ${MASTER_PORT_BASE}"
echo "Train config: ${TRAIN_CONFIG}"
echo "Full train config default: ${FULL_TRAIN_CONFIG_DEFAULT}"
echo "Effective train config: ${EFFECTIVE_TRAIN_CONFIG}"
echo "Model config: ${MODEL_CONFIG}"
echo "Data config: ${DATA_CONFIG}"
echo "RUN_NAME: ${RUN_NAME}"
echo "TRAIN_SEED: ${TRAIN_SEED}"
echo "FOUNDATION_CACHE_DIR: ${FOUNDATION_CACHE_DIR}"
echo "TEACHER_CACHE_DIR: ${TEACHER_CACHE_DIR}"
echo "FEATURE_NAME: ${FEATURE_NAME}"
echo "CACHE_SPLITS: ${CACHE_SPLITS}"
echo "PREDICT_SPLITS: ${PREDICT_SPLITS}"
echo "PROFILE_PROBE_BINS: ${PROFILE_PROBE_BINS}"
echo "RIDGE_ALPHA: ${RIDGE_ALPHA}"
echo "BATCH_SIZE_PER_GPU: ${BATCH_SIZE_PER_GPU}"
echo "NUM_WORKERS: ${NUM_WORKERS}"
echo "GLOBAL_BATCH: $((BATCH_SIZE_PER_GPU * NPROC_PER_NODE))"
echo "DDP_FIND_UNUSED_PARAMETERS: ${DDP_FIND_UNUSED_PARAMETERS}"
echo "=============================================="

missing_cache=0
for split in "${CACHE_SPLIT_ARR[@]}"; do
    if [ ! -f "${FOUNDATION_CACHE_DIR}/manifest_${split}.json" ]; then
        missing_cache=1
        break
    fi
done

if [ "${missing_cache}" -eq 1 ]; then
    echo ""
    echo "[Step 0] Building NT v2 bins16 cache..."
    if [ "${CACHE_NPROC_PER_NODE}" -gt 1 ]; then
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
            --splits "${CACHE_SPLIT_ARR[@]}" \
            --layers 7 \
            --feature_types bins16_mean \
            --record_splits "${PREDICT_SPLIT_ARR[@]}" \
            --profile_probe_bins "${PROFILE_PROBE_BINS}" \
            --batch_size 8 \
            --dtype float16
    else
        CUDA_VISIBLE_DEVICES="${CACHE_GPU_IDS%%,*}" "${PYTHON_BIN}" "${ROOT_DIR}/scripts/build_foundation_cache.py" \
            --data_config "${DATA_CONFIG}" \
            --train_config "${EFFECTIVE_TRAIN_CONFIG}" \
            --output_dir "${FOUNDATION_CACHE_DIR}" \
            --backend nt_v2 \
            --model_dir "${NT_MODEL_DIR}" \
            --splits "${CACHE_SPLIT_ARR[@]}" \
            --layers 7 \
            --feature_types bins16_mean \
            --record_splits "${PREDICT_SPLIT_ARR[@]}" \
            --profile_probe_bins "${PROFILE_PROBE_BINS}" \
            --batch_size 8 \
            --dtype float16
    fi
else
    echo "[Step 0] Foundation cache already exists, skipping."
fi

missing_teacher=0
for split in "${PREDICT_SPLIT_ARR[@]}"; do
    if [ ! -f "${TEACHER_CACHE_DIR}/teacher_manifest_${split}.json" ]; then
        missing_teacher=1
        break
    fi
done

if [ "${missing_teacher}" -eq 1 ]; then
    echo ""
    echo "[Step 1] Exporting teacher cache via teacher_cache_export..."
    "${PYTHON_BIN}" -m transchrombp.evaluation.teacher_cache_export \
        --cache-dir "${FOUNDATION_CACHE_DIR}" \
        --output-dir "${TEACHER_CACHE_DIR}" \
        --feature-name "${FEATURE_NAME}" \
        --train-split train \
        --predict-splits "${PREDICT_SPLIT_ARR[@]}" \
        --alpha "${RIDGE_ALPHA}"
else
    echo "[Step 1] Teacher cache already exists, skipping."
fi

echo ""
echo "[Step 2] Launching corrected-B student distill training on TRAIN_GPU_IDS=${TRAIN_GPU_IDS}..."
launch_train "${MASTER_PORT_BASE}" \
    --train-config "${EFFECTIVE_TRAIN_CONFIG}" \
    --model-config "${MODEL_CONFIG}" \
    --data-config "${DATA_CONFIG}" \
    --run-name "${RUN_NAME}" \
    --output-dir "${OUTPUT_BASE}" \
    --teacher-cache-dir "${TEACHER_CACHE_DIR}" \
    --batch-size-per-gpu "${BATCH_SIZE_PER_GPU}"

echo ""
echo "=============================================="
echo "NT v2 teacher distill completed."
echo "Results in: ${OUTPUT_BASE}/checkpoints, ${OUTPUT_BASE}/logs"
echo "=============================================="
