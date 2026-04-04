#!/usr/bin/env bash
# ============================================================================
# Caduceus-PS tutorial matched A/B launcher
# ============================================================================
# A: corrected-B baseline
# B: corrected-B + frozen online Caduceus-PS token fusion
#
# Default settings keep the matched tutorial strict-compare semantics:
#   - 2-rank DDP
#   - batch_size_per_gpu=16
#   - global batch=32
#   - peak_max_jitter=500
#   - train_revcomp=true
#   - debiased_profile_weight=2.0
#   - count_pool_mode=center
#   - select ckpt by peak.profile_target_jsd_full_mean
# ============================================================================
set -euo pipefail

if [ -z "${BASH_VERSION:-}" ]; then
    echo "[error] Please run this script with bash." >&2
    exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GPU_IDS="${GPU_IDS:-0,1}"
NPROC_PER_NODE="${NPROC_PER_NODE:-2}"
MASTER_ADDR="${MASTER_ADDR:-127.0.0.1}"
MASTER_PORT_BASE="${MASTER_PORT_BASE:-29840}"
VENV_DIR="${VENV_DIR:-/data1/zhoujiazhen/bylw_atac/.venvs/caduceus-ps}"
OUTPUT_DIR="${OUTPUT_DIR:-/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs}"
LOG_DIR="${LOG_DIR:-/data1/zhoujiazhen/bylw_atac/logs}"
BIAS_CKPT="${BIAS_CKPT:-${OUTPUT_DIR}/checkpoints/ablation_tf_20260318_shared_bias/best.pt}"
CADUCEUS_MODEL_PATH="${CADUCEUS_MODEL_PATH:-/data1/zhoujiazhen/bylw_atac/foundation_models/caduceus/caduceus-ps_seqlen-131k_d_model-256_n_layer-16}"
CADUCEUS_LOCAL_FILES_ONLY="${CADUCEUS_LOCAL_FILES_ONLY:-true}"
BATCH_SIZE_PER_GPU="${BATCH_SIZE_PER_GPU:-16}"
NUM_WORKERS="${NUM_WORKERS:-2}"
SEED="${SEED:-42}"
SMOKE_STEPS="${SMOKE_STEPS:-2}"
DRY_RUN="${DRY_RUN:-0}"
FORCE_RERUN="${FORCE_RERUN:-0}"
RUN_LABEL="${RUN_LABEL:-$(date +%Y%m%d_%H%M%S)}"
ARMS=()
SMOKE=false

MODEL_A="${ROOT_DIR}/configs/model/transchrombp_teacher_v2_center_pool.yaml"
MODEL_B="${ROOT_DIR}/configs/model/transchrombp_teacher_v2_center_pool_caduceus_ps.yaml"
TRAIN_CFG="${ROOT_DIR}/configs/train/train_tutorial_corrected_b_strict_compare_6000.yaml"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --arms) IFS=',' read -ra ARMS <<< "$2"; shift 2 ;;
        --gpu-ids) GPU_IDS="$2"; shift 2 ;;
        --seed) SEED="$2"; shift 2 ;;
        --smoke) SMOKE=true; shift ;;
        --run-label) RUN_LABEL="$2"; shift 2 ;;
        *) echo "[error] Unknown option: $1" >&2; exit 1 ;;
    esac
done

if [[ ${#ARMS[@]} -eq 0 ]]; then
    ARMS=(A B)
fi

if [[ ! -d "${VENV_DIR}" ]]; then
    echo "[error] Missing Caduceus env: ${VENV_DIR}" >&2
    exit 1
fi
if [[ ! -x "${VENV_DIR}/bin/python" ]]; then
    echo "[error] Missing python in Caduceus env: ${VENV_DIR}/bin/python" >&2
    exit 1
fi
if [[ ! -x "${VENV_DIR}/bin/torchrun" ]]; then
    echo "[error] Missing torchrun in Caduceus env: ${VENV_DIR}/bin/torchrun" >&2
    exit 1
fi
if [[ ! -f "${BIAS_CKPT}" ]]; then
    echo "[error] Missing bias checkpoint: ${BIAS_CKPT}" >&2
    exit 1
fi
if [[ "${CADUCEUS_LOCAL_FILES_ONLY}" == "true" ]] && [[ ! -d "${CADUCEUS_MODEL_PATH}" ]]; then
    echo "[error] Missing local Caduceus model directory: ${CADUCEUS_MODEL_PATH}" >&2
    exit 1
fi

export CUDA_VISIBLE_DEVICES="${GPU_IDS}"
export PATH="${VENV_DIR}/bin:${PATH:-}"
export PYTHONPATH="${ROOT_DIR}/src:${PYTHONPATH:-}"

mkdir -p "${LOG_DIR}" "${OUTPUT_DIR}/runtime/caduceus_tutorial_ab_${RUN_LABEL}"
RUNTIME_DIR="${OUTPUT_DIR}/runtime/caduceus_tutorial_ab_${RUN_LABEL}"
PYTHON_BIN="${VENV_DIR}/bin/python"
TORCHRUN_BIN="${VENV_DIR}/bin/torchrun"
cd "${ROOT_DIR}"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$1"; }

generate_configs() {
    local arm="$1"
    local run_name="$2"
    local model_src="$3"
    local out_model="${RUNTIME_DIR}/model_${run_name}.yaml"
    local out_train="${RUNTIME_DIR}/train_${run_name}.yaml"

    "${PYTHON_BIN}" - "$model_src" "$out_model" "$BIAS_CKPT" "$CADUCEUS_MODEL_PATH" "$CADUCEUS_LOCAL_FILES_ONLY" <<'PY'
import sys
from pathlib import Path
import yaml

src, dst, bias_ckpt, caduceus_model_path, local_only = sys.argv[1:6]
with open(src, "r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
cfg.setdefault("bias_branch", {})["pretrained_path"] = bias_ckpt
cfg["bias_branch"]["freeze_bias_core"] = True
if "caduceus_branch" in cfg:
    cfg["caduceus_branch"]["model_path"] = caduceus_model_path
    cfg["caduceus_branch"]["local_files_only"] = str(local_only).lower() == "true"
Path(dst).write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
PY

    "${PYTHON_BIN}" - "$TRAIN_CFG" "$out_train" "$SEED" "$run_name" "$ROOT_DIR" "$BATCH_SIZE_PER_GPU" "$NUM_WORKERS" <<'PY'
import sys
from pathlib import Path
import yaml

src, dst, seed, run_name, root, batch_size, num_workers = sys.argv[1:8]
with open(src, "r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
cfg["seed"] = int(seed)
cfg.setdefault("logging", {})["run_name"] = run_name
cfg.setdefault("data", {})["batch_size_per_gpu"] = int(batch_size)
cfg["data"]["num_workers"] = int(num_workers)
data_cfg = cfg["data"].get("config_path", "")
if data_cfg and not Path(data_cfg).is_absolute():
    cfg["data"]["config_path"] = str((Path(root) / data_cfg).resolve())
Path(dst).write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
PY
}

launch_train() {
    local run_name="$1"
    local model_cfg="$2"
    local train_cfg="$3"
    local master_port="$4"
    local log_path="${LOG_DIR}/${run_name}.log"

    if [[ "${FORCE_RERUN}" != "1" ]] && [[ -f "${OUTPUT_DIR}/checkpoints/${run_name}/best.pt" ]]; then
        log "[skip] ${run_name} already has best.pt"
        return 0
    fi

    log "[run] ${run_name} gpu_ids=${GPU_IDS} log=${log_path}"
    if [[ "${DRY_RUN}" == "1" ]]; then
        log "[dry-run] torchrun --nproc_per_node=${NPROC_PER_NODE} -m transchrombp.training.train_ddp ..."
        return 0
    fi

    local dry_flag=()
    if ${SMOKE}; then
        dry_flag=(--dry-run-steps "${SMOKE_STEPS}")
    fi

    "${TORCHRUN_BIN}" \
        --nnodes=1 \
        --nproc_per_node="${NPROC_PER_NODE}" \
        --master_addr="${MASTER_ADDR}" \
        --master_port="${master_port}" \
        -m transchrombp.training.train_ddp \
        --model-config "${model_cfg}" \
        --train-config "${train_cfg}" \
        --output-dir "${OUTPUT_DIR}" \
        --run-name "${run_name}" \
        --batch-size-per-gpu "${BATCH_SIZE_PER_GPU}" \
        "${dry_flag[@]}" \
        >> "${log_path}" 2>&1

    log "[done] ${run_name}"
}

log "========== Caduceus tutorial matched A/B =========="
log "[config] RUN_LABEL=${RUN_LABEL} ARMS=${ARMS[*]} GPU_IDS=${GPU_IDS} SEED=${SEED} SMOKE=${SMOKE}"
log "[config] BIAS_CKPT=${BIAS_CKPT}"
log "[config] CADUCEUS_MODEL_PATH=${CADUCEUS_MODEL_PATH} local_only=${CADUCEUS_LOCAL_FILES_ONLY}"
log "[config] BATCH_SIZE_PER_GPU=${BATCH_SIZE_PER_GPU} NPROC_PER_NODE=${NPROC_PER_NODE}"

for arm in "${ARMS[@]}"; do
    case "${arm^^}" in
        A)
            run_name="caduceus_tutorial_A_corrected_b_s${SEED}_${RUN_LABEL}"
            generate_configs "A" "${run_name}" "${MODEL_A}"
            launch_train "${run_name}" "${RUNTIME_DIR}/model_${run_name}.yaml" "${RUNTIME_DIR}/train_${run_name}.yaml" "${MASTER_PORT_BASE}"
            ;;
        B)
            run_name="caduceus_tutorial_B_caduceus_ps_s${SEED}_${RUN_LABEL}"
            generate_configs "B" "${run_name}" "${MODEL_B}"
            launch_train "${run_name}" "${RUNTIME_DIR}/model_${run_name}.yaml" "${RUNTIME_DIR}/train_${run_name}.yaml" "$((MASTER_PORT_BASE + 1))"
            ;;
        *)
            echo "[error] Unknown arm: ${arm}" >&2
            exit 1
            ;;
    esac
done

log "========== Caduceus tutorial A/B complete =========="
