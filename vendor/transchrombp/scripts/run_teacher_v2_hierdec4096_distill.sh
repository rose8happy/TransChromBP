#!/usr/bin/env bash
# Export the hierdec4096 teacher cache on demand, then train the short10 distill student.

set -euo pipefail

if [ -z "${BASH_VERSION:-}" ]; then
    echo "[error] Please run this script with bash." >&2
    exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_CONFIG="${DATA_CONFIG:-${ROOT_DIR}/configs/data/data_tutorial_canonical_v1_longctx4096.yaml}"
MODEL_CONFIG="${MODEL_CONFIG:-${ROOT_DIR}/configs/model/transchrombp_teacher_v2_hierdec4096.yaml}"
TEACHER_TRAIN_CONFIG="${TEACHER_TRAIN_CONFIG:-${ROOT_DIR}/configs/train/train_tutorial_teacher_v2_hierdec4096_teacher30.yaml}"
TRAIN_CONFIG="${TRAIN_CONFIG:-${ROOT_DIR}/configs/train/train_tutorial_teacher_v2_hierdec4096_distill_short10.yaml}"
OUTPUT_BASE="${OUTPUT_BASE:-${ROOT_DIR}/outputs}"
TEACHER_CACHE_DIR="${TEACHER_CACHE_DIR:-${OUTPUT_BASE}/teacher_cache/teacher_v2_hierdec4096_teacher30}"
RUN_NAME="${RUN_NAME:-teacher_v2_hierdec4096_distill_short10}"
TEACHER_CKPT="${TEACHER_CKPT:-}"
VENV_DIR="${VENV_DIR:-/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b}"

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

require_file "${DATA_CONFIG}" "data config"
require_file "${MODEL_CONFIG}" "model config"
require_file "${TEACHER_TRAIN_CONFIG}" "teacher train config"
require_file "${TRAIN_CONFIG}" "distill train config"

mkdir -p "${OUTPUT_BASE}" "${TEACHER_CACHE_DIR}"

missing_teacher_cache=0
for split in train valid; do
    if [ ! -f "${TEACHER_CACHE_DIR}/teacher_manifest_${split}.json" ]; then
        missing_teacher_cache=1
        break
    fi
done

if [ "${missing_teacher_cache}" -eq 1 ]; then
    if [ -z "${TEACHER_CKPT}" ]; then
        echo "[error] TEACHER_CKPT must point to an existing teacher checkpoint." >&2
        exit 1
    fi
    require_file "${TEACHER_CKPT}" "teacher checkpoint"
    python -m transchrombp.evaluation.model_teacher_cache_export \
        --checkpoint "${TEACHER_CKPT}" \
        --data-config "${DATA_CONFIG}" \
        --output-dir "${TEACHER_CACHE_DIR}" \
        --splits train valid
fi

python -m transchrombp.training.train_ddp \
    --model-config "${MODEL_CONFIG}" \
    --train-config "${TRAIN_CONFIG}" \
    --data-config "${DATA_CONFIG}" \
    --output-dir "${OUTPUT_BASE}" \
    --run-name "${RUN_NAME}" \
    --teacher-cache-dir "${TEACHER_CACHE_DIR}" \
    "$@"
