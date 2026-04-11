#!/usr/bin/env bash
# Launch the hierdec4096 E2 probe on longctx4096 inputs with the short10 budget.

set -euo pipefail

if [ -z "${BASH_VERSION:-}" ]; then
    echo "[error] Please run this script with bash." >&2
    exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_CONFIG="${DATA_CONFIG:-${ROOT_DIR}/configs/data/data_tutorial_canonical_v1_longctx4096.yaml}"
MODEL_CONFIG="${MODEL_CONFIG:-${ROOT_DIR}/configs/model/transchrombp_teacher_v2_hierdec4096.yaml}"
TRAIN_CONFIG="${TRAIN_CONFIG:-${ROOT_DIR}/configs/train/train_tutorial_teacher_v2_hierdec4096_short10.yaml}"
OUTPUT_DIR="${OUTPUT_DIR:-${ROOT_DIR}/outputs}"
RUN_NAME="${RUN_NAME:-teacher_v2_hierdec4096_short10}"
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
require_file "${TRAIN_CONFIG}" "train config"

exec python -m transchrombp.training.train_ddp \
    --model-config "${MODEL_CONFIG}" \
    --train-config "${TRAIN_CONFIG}" \
    --data-config "${DATA_CONFIG}" \
    --output-dir "${OUTPUT_DIR}" \
    --run-name "${RUN_NAME}" \
    "$@"
