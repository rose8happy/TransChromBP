#!/usr/bin/env bash
# Genos Cached Pilot Launcher
# 串行运行 P0 (baseline_short10) → P1 (global_late_film) → P2 (global_count_only)

set -euo pipefail

if [ -z "${BASH_VERSION:-}" ]; then
    echo "[error] Please run this script with bash." >&2
    exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PACKAGE_PARENT="$(cd "${ROOT_DIR}/.." && pwd)"
GPU_ID="${1:-${GPU_ID:-0}}"
VENV_DIR="${VENV_DIR:-/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b}"
OUTPUT_BASE="${OUTPUT_BASE:-${ROOT_DIR}/outputs}"
GENOS_CACHE_DIR="${GENOS_CACHE_DIR:-${ROOT_DIR}/genos_cache}"
GENOS_MODEL_PATH="${GENOS_MODEL_PATH:-/data1/zhoujiazhen/bylw_atac/foundation_models/genos/Genos-1.2B}"
TRAIN_CONFIG="${TRAIN_CONFIG:-${ROOT_DIR}/configs/train/train_genos_cached_short10.yaml}"
DATA_CONFIG="${DATA_CONFIG:-${ROOT_DIR}/configs/data/data_tutorial_canonical_v1.yaml}"
PYTHON_BIN="${PYTHON_BIN:-python}"

export CUDA_VISIBLE_DEVICES="${GPU_ID}"
export PATH="${VENV_DIR}/bin:${PATH:-}"
export PYTHONPATH="${PACKAGE_PARENT}:${ROOT_DIR}/src:${PYTHONPATH:-}"

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

find_g0_best() {
    if [ -n "${G0_BEST_PT:-}" ] && [ -f "${G0_BEST_PT}" ]; then
        printf '%s\n' "${G0_BEST_PT}"
        return 0
    fi

    local candidates=(
        "${OUTPUT_BASE}/checkpoints/genos_20260321_baseline_s42/best.pt"
        "${OUTPUT_BASE}/checkpoints/genos_cached_P0_baseline/best.pt"
        "${OUTPUT_BASE}/genos_pilot_G0_baseline/best.pt"
    )

    local path
    for path in "${candidates[@]}"; do
        if [ -f "${path}" ]; then
            printf '%s\n' "${path}"
            return 0
        fi
    done
    return 1
}

G0_BEST_PT="$(find_g0_best || true)"

echo "=============================================="
echo "Genos Cached Pilot"
echo "ROOT_DIR: ${ROOT_DIR}"
echo "GPU: ${GPU_ID}"
echo "Cache dir: ${GENOS_CACHE_DIR}"
echo "Train config: ${TRAIN_CONFIG}"
echo "Data config: ${DATA_CONFIG}"
echo "=============================================="

require_file "${TRAIN_CONFIG}" "train config"
require_file "${DATA_CONFIG}" "data config"

mkdir -p "${OUTPUT_BASE}" "${GENOS_CACHE_DIR}"

# Step 0: Build cache if missing
if [ ! -f "${GENOS_CACHE_DIR}/manifest_train.json" ] || [ ! -f "${GENOS_CACHE_DIR}/manifest_valid.json" ]; then
    require_dir "${GENOS_MODEL_PATH}" "Genos model directory"
    echo ""
    echo "[Step 0] Building Genos summary cache..."
    "${PYTHON_BIN}" "${ROOT_DIR}/scripts/build_genos_summary_cache.py" \
        --data_config "${DATA_CONFIG}" \
        --train_config "${TRAIN_CONFIG}" \
        --genos_model_path "${GENOS_MODEL_PATH}" \
        --output_dir "${GENOS_CACHE_DIR}" \
        --splits train valid \
        --features global_mean bins4_mean \
        --batch_size 8 \
        --bins4_budget_gib 4.0
else
    echo "[Step 0] Cache already exists, skipping."
fi

# Step 1: Run probes
if [ -n "${G0_BEST_PT}" ]; then
    echo ""
    echo "[Step 1] Running Genos summary probes with ${G0_BEST_PT}..."
    "${PYTHON_BIN}" "${ROOT_DIR}/scripts/run_genos_summary_probe.py" \
        --g0_checkpoint "${G0_BEST_PT}" \
        --model_config "${ROOT_DIR}/configs/model/v2fix_baseline.yaml" \
        --train_config "${TRAIN_CONFIG}" \
        --data_config "${DATA_CONFIG}" \
        --genos_cache_dir "${GENOS_CACHE_DIR}" \
        --output_dir "${OUTPUT_BASE}/genos_cached_probes" \
        --batch_size 64
else
    echo "[Step 1] SKIP: G0 best checkpoint not found. Set G0_BEST_PT to override."
fi

# Step 2: P0 baseline_short10
echo ""
echo "[Step 2] P0: baseline_short10"
RUN_NAME="genos_cached_P0_baseline"
"${PYTHON_BIN}" -m transchrombp.training.train_ddp \
    --train-config "${TRAIN_CONFIG}" \
    --model-config "${ROOT_DIR}/configs/model/v2fix_baseline.yaml" \
    --data-config "${DATA_CONFIG}" \
    --run-name "${RUN_NAME}" \
    --output-dir "${OUTPUT_BASE}"

# Step 3: P1 global_late_film
echo ""
echo "[Step 3] P1: global_late_film"
RUN_NAME="genos_cached_P1_late_film"
"${PYTHON_BIN}" -m transchrombp.training.train_ddp \
    --train-config "${TRAIN_CONFIG}" \
    --model-config "${ROOT_DIR}/configs/model/v2fix_genos_global_late_film.yaml" \
    --data-config "${DATA_CONFIG}" \
    --run-name "${RUN_NAME}" \
    --output-dir "${OUTPUT_BASE}" \
    --genos-cache-dir "${GENOS_CACHE_DIR}"

# Step 4: P2 global_count_only
echo ""
echo "[Step 4] P2: global_count_only"
RUN_NAME="genos_cached_P2_count_only"
"${PYTHON_BIN}" -m transchrombp.training.train_ddp \
    --train-config "${TRAIN_CONFIG}" \
    --model-config "${ROOT_DIR}/configs/model/v2fix_genos_global_count.yaml" \
    --data-config "${DATA_CONFIG}" \
    --run-name "${RUN_NAME}" \
    --output-dir "${OUTPUT_BASE}" \
    --genos-cache-dir "${GENOS_CACHE_DIR}"

echo ""
echo "=============================================="
echo "All pilot runs completed."
echo "Results in: ${OUTPUT_BASE}/checkpoints, ${OUTPUT_BASE}/logs"
echo "=============================================="
