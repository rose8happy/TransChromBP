#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

case "$ROOT_DIR" in
  /data1/zhoujiazhen/bylw_atac/TransChromBP)
    HOST_TAG="6000"
    ENV_DIR="/data1/zhoujiazhen/bylw_atac/.mamba/envs/transchrombp"
    OUTPUT_DIR_DEFAULT="/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs"
    LOG_DIR_DEFAULT="/data1/zhoujiazhen/bylw_atac/logs"
    DATA_CONFIG_REL="configs/data/data_tutorial_canonical_v1.yaml"
    BIAS_CKPT_DEFAULT="/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/checkpoints/ablation_tf_20260318_shared_bias/best.pt"
    ;;
  /home/zhengwei/bylw_atac/TransChromBP)
    HOST_TAG="6002"
    ENV_DIR="/home/zhengwei/bylw_atac/.mamba/envs/transchrombp"
    OUTPUT_DIR_DEFAULT="/home/zhengwei/bylw_atac/TransChromBP/outputs"
    LOG_DIR_DEFAULT="/home/zhengwei/bylw_atac/logs"
    DATA_CONFIG_REL="configs/data/data_tutorial_canonical_v1_6002.yaml"
    BIAS_CKPT_DEFAULT="/home/zhengwei/bylw_atac/TransChromBP/outputs/checkpoints/ablation_tf_20260318_shared_bias/best.pt"
    ;;
  *)
    echo "[error] Unsupported ROOT_DIR=$ROOT_DIR" >&2
    exit 1
    ;;
esac

VARIANT="${1:-}"
if [[ -z "$VARIANT" ]]; then
  echo "Usage: bash scripts/launch_profile_shortcut_probe.sh <tf_sg0_deb2|tf_sg1_deb0|notf_sg0_deb2|notf_sg1_deb0|tf_center_sg1_deb2>" >&2
  exit 1
fi

case "$VARIANT" in
  tf_sg0_deb2)
    MODEL_REL="configs/model/transchrombp_teacher_v2_nosg.yaml"
    DEBIASED_PROFILE_WEIGHT="2.0"
    ;;
  tf_sg1_deb0)
    MODEL_REL="configs/model/transchrombp_teacher_v2.yaml"
    DEBIASED_PROFILE_WEIGHT="0.0"
    ;;
  notf_sg0_deb2)
    MODEL_REL="configs/model/transchrombp_teacher_v2_noTF_nosg.yaml"
    DEBIASED_PROFILE_WEIGHT="2.0"
    ;;
  notf_sg1_deb0)
    MODEL_REL="configs/model/transchrombp_teacher_v2_noTF.yaml"
    DEBIASED_PROFILE_WEIGHT="0.0"
    ;;
  tf_center_sg1_deb2)
    MODEL_REL="configs/model/transchrombp_teacher_v2_center_pool.yaml"
    DEBIASED_PROFILE_WEIGHT="2.0"
    ;;
  *)
    echo "[error] Unsupported variant: $VARIANT" >&2
    exit 1
    ;;
esac

SEED="${SEED:-42}"
RUN_TAG="${RUN_TAG:-profile_shortcut_$(date +%Y%m%d)}"
MAX_EPOCHS="${MAX_EPOCHS:-30}"
BATCH_SIZE_PER_GPU="${BATCH_SIZE_PER_GPU:-16}"
GRAD_ACCUM_STEPS="${GRAD_ACCUM_STEPS:-2}"
NPROC_PER_NODE="${NPROC_PER_NODE:-1}"
MASTER_PORT="${MASTER_PORT:-29591}"
OUTPUT_DIR="${OUTPUT_DIR:-$OUTPUT_DIR_DEFAULT}"
LOG_DIR="${LOG_DIR:-$LOG_DIR_DEFAULT}"
BIAS_CKPT="${BIAS_CKPT:-$BIAS_CKPT_DEFAULT}"
RUN_NAME="${RUN_NAME:-${RUN_TAG}_${VARIANT}_s${SEED}_${HOST_TAG}}"

BASE_TRAIN_TEMPLATE="${BASE_TRAIN_TEMPLATE:-$ROOT_DIR/configs/train/train_ablation_v2_main.yaml}"
MODEL_SRC="${ROOT_DIR}/${MODEL_REL}"
RUNTIME_DIR="${OUTPUT_DIR}/runtime/${RUN_TAG}"
RUNTIME_CONFIG_DIR="${RUNTIME_DIR}/configs"
MODEL_RT="${RUNTIME_CONFIG_DIR}/model_${RUN_NAME}.yaml"
TRAIN_RT="${RUNTIME_CONFIG_DIR}/train_${RUN_NAME}.yaml"
LOG_FILE="${LOG_FILE:-${LOG_DIR}/${RUN_NAME}.log}"

mkdir -p "$RUNTIME_CONFIG_DIR" "$LOG_DIR" "$OUTPUT_DIR"

if [[ ! -f "$MODEL_SRC" ]]; then
  echo "[error] Missing model config: $MODEL_SRC" >&2
  exit 1
fi
if [[ ! -f "$BASE_TRAIN_TEMPLATE" ]]; then
  echo "[error] Missing train template: $BASE_TRAIN_TEMPLATE" >&2
  exit 1
fi
if [[ ! -f "$BIAS_CKPT" ]]; then
  echo "[error] Missing bias checkpoint: $BIAS_CKPT" >&2
  exit 1
fi

python3 - "$MODEL_SRC" "$MODEL_RT" "$BIAS_CKPT" <<'PY'
import sys, yaml
from pathlib import Path

src, dst, bias_ckpt = sys.argv[1:]
cfg = yaml.safe_load(Path(src).read_text())
cfg.setdefault("bias_branch", {})["pretrained_path"] = bias_ckpt
cfg["bias_branch"]["freeze_bias_core"] = True
Path(dst).write_text(yaml.safe_dump(cfg, sort_keys=False))
PY

python3 - "$BASE_TRAIN_TEMPLATE" "$TRAIN_RT" "$SEED" "$RUN_NAME" "$OUTPUT_DIR" "$ROOT_DIR" "$DATA_CONFIG_REL" "$BATCH_SIZE_PER_GPU" "$GRAD_ACCUM_STEPS" "$MAX_EPOCHS" "$DEBIASED_PROFILE_WEIGHT" <<'PY'
import sys, yaml
from pathlib import Path

(
    src,
    dst,
    seed,
    run_name,
    output_dir,
    root_dir,
    data_cfg_rel,
    batch_size,
    grad_accum,
    max_epochs,
    debiased_profile_weight,
) = sys.argv[1:]

cfg = yaml.safe_load(Path(src).read_text())
cfg["seed"] = int(seed)
cfg["max_epochs"] = int(max_epochs)
cfg.setdefault("trainer", {})["grad_accum_steps"] = int(grad_accum)
cfg.setdefault("loss", {})["debiased_profile_weight"] = float(debiased_profile_weight)
cfg.setdefault("data", {})["batch_size_per_gpu"] = int(batch_size)
cfg["data"]["config_path"] = str((Path(root_dir) / data_cfg_rel).resolve())
cfg.setdefault("logging", {})["run_name"] = run_name
cfg["logging"]["output_dir"] = output_dir
Path(dst).write_text(yaml.safe_dump(cfg, sort_keys=False))
PY

export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="$ROOT_DIR/src:${PYTHONPATH:-}"

cd "$ROOT_DIR"

echo "[launch] host=$HOST_TAG variant=$VARIANT run_name=$RUN_NAME"
echo "[launch] model_config=$MODEL_RT"
echo "[launch] train_config=$TRAIN_RT"
echo "[launch] log_file=$LOG_FILE"
echo "[launch] batch_size_per_gpu=$BATCH_SIZE_PER_GPU grad_accum_steps=$GRAD_ACCUM_STEPS max_epochs=$MAX_EPOCHS"

torchrun \
  --standalone \
  --nproc_per_node="$NPROC_PER_NODE" \
  --master_port="$MASTER_PORT" \
  -m transchrombp.training.train_ddp \
  --model-config "$MODEL_RT" \
  --train-config "$TRAIN_RT" \
  --run-name "$RUN_NAME" \
  --output-dir "$OUTPUT_DIR" \
  2>&1 | tee "$LOG_FILE"
