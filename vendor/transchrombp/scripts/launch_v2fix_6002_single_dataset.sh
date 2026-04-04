#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="/home/zhengwei/bylw_atac/TransChromBP"
ENV_DIR="/home/zhengwei/bylw_atac/.mamba/envs/transchrombp"
OUTPUT_DIR="${OUTPUT_DIR:-$ROOT_DIR/outputs}"
LOG_DIR="${LOG_DIR:-/home/zhengwei/bylw_atac/logs}"

MODEL_CONFIG="${MODEL_CONFIG:-${1:-$ROOT_DIR/configs/model/v2fix_baseline.yaml}}"
if [[ $# -gt 0 ]]; then shift; fi
TRAIN_CONFIG="${TRAIN_CONFIG:-${1:-$ROOT_DIR/configs/train/train_gm12878_v2_6002_smoke.yaml}}"
if [[ $# -gt 0 ]]; then shift; fi
RUN_NAME="${RUN_NAME:-$(basename "${TRAIN_CONFIG%.yaml}")}"

NPROC_PER_NODE="${NPROC_PER_NODE:-1}"
MASTER_PORT="${MASTER_PORT:-29571}"
LOG_FILE="${LOG_FILE:-$LOG_DIR/${RUN_NAME}.log}"

mkdir -p "$LOG_DIR" "$OUTPUT_DIR"
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="$ROOT_DIR/src:${PYTHONPATH:-}"

cd "$ROOT_DIR"

echo "[launch] run_name=$RUN_NAME"
echo "[launch] model_config=$MODEL_CONFIG"
echo "[launch] train_config=$TRAIN_CONFIG"
echo "[launch] log_file=$LOG_FILE"

torchrun \
  --standalone \
  --nproc_per_node="$NPROC_PER_NODE" \
  --master_port="$MASTER_PORT" \
  -m transchrombp.training.train_ddp \
  --model-config "$MODEL_CONFIG" \
  --train-config "$TRAIN_CONFIG" \
  --run-name "$RUN_NAME" \
  --output-dir "$OUTPUT_DIR" \
  "$@" 2>&1 | tee "$LOG_FILE"
