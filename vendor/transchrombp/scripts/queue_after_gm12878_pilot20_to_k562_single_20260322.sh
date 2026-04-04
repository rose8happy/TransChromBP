#!/usr/bin/env bash
set -euo pipefail

export TZ=Asia/Shanghai
CURRENT_RUN="gm12878_v2_6002_pilot20_20260322_191746"
NOT_BEFORE="2026-03-23 06:00:00"
ROOT="/home/zhengwei/bylw_atac/TransChromBP"
LOG_DIR="/home/zhengwei/bylw_atac/logs"
MODEL_CONFIG="$ROOT/configs/model/v2fix_baseline.yaml"
TRAIN_CONFIG="$ROOT/configs/train/train_k562_v2_6002_single.yaml"
MASTER_PORT=29574

mkdir -p "$LOG_DIR"
NOT_BEFORE_EPOCH=$(date -d "$NOT_BEFORE" +%s)

echo "[queue] $(date '+%F %T %Z') wait current_run=$CURRENT_RUN not_before=$NOT_BEFORE"
while pgrep -f "$CURRENT_RUN" >/dev/null 2>&1; do
  echo "[queue] $(date '+%F %T %Z') current run still active"
  sleep 120
done

while [ "$(date +%s)" -lt "$NOT_BEFORE_EPOCH" ]; do
  remain=$((NOT_BEFORE_EPOCH - $(date +%s)))
  echo "[queue] $(date '+%F %T %Z') current run done; hold until $NOT_BEFORE remain=${remain}s"
  sleep 120
done

RUN_TS=$(date +%Y%m%d_%H%M%S)
RUN_NAME="k562_v2_6002_single_${RUN_TS}"
TRAIN_LOG="$LOG_DIR/${RUN_NAME}.log"
LAUNCH_LOG="$LOG_DIR/${RUN_NAME}.launch.log"

echo "[queue] $(date '+%F %T %Z') launching $RUN_NAME"
nohup env \
  RUN_NAME="$RUN_NAME" \
  LOG_FILE="$TRAIN_LOG" \
  MASTER_PORT="$MASTER_PORT" \
  bash "$ROOT/scripts/launch_v2fix_6002_single_dataset.sh" \
  "$MODEL_CONFIG" \
  "$TRAIN_CONFIG" \
  > "$LAUNCH_LOG" 2>&1 < /dev/null &
CHILD=$!
echo "[queue] child_pid=$CHILD"
echo "[queue] train_log=$TRAIN_LOG"
echo "[queue] launch_log=$LAUNCH_LOG"
