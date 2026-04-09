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

GPU_IDS_CLEAN="${TRAIN_GPU_IDS//[[:space:]]/}"
IFS=',' read -r -a GPU_ID_LIST <<< "${GPU_IDS_CLEAN}"
if [ "${#GPU_ID_LIST[@]}" -lt "${NPROC_PER_NODE}" ]; then
    echo "[error] TRAIN_GPU_IDS (${TRAIN_GPU_IDS}) has fewer GPUs than NPROC_PER_NODE (${NPROC_PER_NODE})" >&2
    exit 1
fi

require_file "${TRAIN_CONFIG}" "train config"

RUNTIME_DIR="${OUTPUT_BASE}/runtime/msdls_v2_gate"
mkdir -p "${RUNTIME_DIR}"
RUNTIME_TRAIN="${RUNTIME_DIR}/train_${RUN_NAME}.yaml"

"${PYTHON_BIN}" - "${TRAIN_CONFIG}" "${RUNTIME_TRAIN}" "${SEED}" "${RUN_NAME}" "${BATCH_SIZE_PER_GPU}" "${NUM_WORKERS}" "${DATA_CONFIG}" "${OUTPUT_BASE}" <<'PY'
import sys
from pathlib import Path
import yaml

src, dst, seed, run_name, batch_size, num_workers, data_config, output_base = sys.argv[1:9]
cfg = yaml.safe_load(Path(src).read_text(encoding="utf-8"))
cfg["seed"] = int(seed)
cfg.setdefault("data", {})["batch_size_per_gpu"] = int(batch_size)
cfg["data"]["num_workers"] = int(num_workers)
cfg["data"]["config_path"] = data_config
cfg.setdefault("logging", {})["run_name"] = run_name
cfg["logging"]["output_dir"] = output_base
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

require_file "${DATA_CONFIG}" "data config"
require_file "${MODEL_CONFIG}" "model config"

if [ "${NPROC_PER_NODE}" -gt 1 ]; then
    if ! command -v torchrun >/dev/null 2>&1; then
        echo "[error] torchrun not found in PATH" >&2
        exit 1
    fi
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
