#!/usr/bin/env bash
# ============================================================================
# Transformer Ablation Experiment — V2-full vs V2-noTF
# ============================================================================
# Design:
#   1. One shared bias pretrain (Transformer not used during bias_pretrain)
#   2. 6 main runs: 2 models × 3 seeds
#   3. Checkpoint selection: peak.profile_target_jsd_debiased_mean (min)
#   4. No test-set peeking during training
#
# Usage:
#   bash scripts/run_ablation_transformer.sh
#
# Environment variables (optional overrides):
#   NPROC_PER_NODE   — GPUs per node (default: 2)
#   DRY_RUN          — set to "1" to only print commands without running
#   ABLATION_TAG     — resume or override experiment tag (default: ablation_tf_YYYYMMDD)
#   WAIT_FOR_IDLE_GPUS — wait for GPUs to become idle before each launch (default: 1)
# ============================================================================
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_DIR="/data1/zhoujiazhen/bylw_atac/.mamba/envs/transchrombp"
LOG_DIR="/data1/zhoujiazhen/bylw_atac/logs"
OUTPUT_DIR="$ROOT_DIR/outputs"

export PATH="/usr/bin:/bin:${ENV_DIR}/bin:${PATH:-}"
export LD_LIBRARY_PATH="${ENV_DIR}/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="$ROOT_DIR/src:${PYTHONPATH:-}"

# Make all relative paths resolve from the repository root, even under nohup/cron.
cd "$ROOT_DIR"

NPROC_PER_NODE="${NPROC_PER_NODE:-2}"
DRY_RUN="${DRY_RUN:-0}"
MASTER_PORT_BASE="${MASTER_PORT_BASE:-29530}"
WAIT_FOR_IDLE_GPUS="${WAIT_FOR_IDLE_GPUS:-1}"
GPU_IDLE_CHECK_INTERVAL="${GPU_IDLE_CHECK_INTERVAL:-60}"
GPU_IDLE_USED_MEM_THRESHOLD_MIB="${GPU_IDLE_USED_MEM_THRESHOLD_MIB:-1024}"
RUN_RETRIES_ON_OOM="${RUN_RETRIES_ON_OOM:-3}"

# ---------- Configs ----------
MODEL_V2_FULL="$ROOT_DIR/configs/model/transchrombp_teacher_v2.yaml"
MODEL_V2_NOTF="$ROOT_DIR/configs/model/transchrombp_teacher_v2_noTF.yaml"
TRAIN_BIAS="$ROOT_DIR/configs/train/train_bias_tutorial_teacher_v2.yaml"
TRAIN_MAIN="$ROOT_DIR/configs/train/train_ablation_v2_main.yaml"

# ---------- Seeds ----------
SEEDS=(42 1234 2024)

# ---------- Naming ----------
ABLATION_TAG="${ABLATION_TAG:-ablation_tf_$(date +%Y%m%d)}"
BIAS_RUN_NAME="${ABLATION_TAG}_shared_bias"

# ---------- Runtime config dir ----------
RUNTIME_DIR="$OUTPUT_DIR/runtime/$ABLATION_TAG"
CONFIG_DIR="$RUNTIME_DIR/configs"
META_DIR="$RUNTIME_DIR/meta"
mkdir -p "$LOG_DIR" "$CONFIG_DIR" "$META_DIR"

# ---------- Validate prerequisites ----------
for f in "$MODEL_V2_FULL" "$MODEL_V2_NOTF" "$TRAIN_BIAS" "$TRAIN_MAIN"; do
  if [[ ! -f "$f" ]]; then
    echo "[error] missing config: $f" >&2
    exit 1
  fi
done

# Verify train_ddp.py has validation metrics support
if ! grep -q 'profile_target_jsd_debiased_mean' "$ROOT_DIR/src/transchrombp/training/train_ddp.py"; then
  echo "[error] train_ddp.py does not have validation metric alignment code." >&2
  echo "[error] Deploy the updated train_ddp.py first." >&2
  exit 1
fi

log() { printf '[%(%F %T)T] %s\n' -1 "$1"; }

count_compute_pids() {
  nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null | awk 'NF{count++} END{print count+0}'
}

gpu_snapshot() {
  nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv,noheader,nounits 2>/dev/null
}

all_gpus_idle() {
  local mem_ok=1
  while IFS=',' read -r gpu_idx gpu_mem gpu_util; do
    gpu_mem="${gpu_mem// /}"
    [[ -z "$gpu_mem" ]] && continue
    if (( gpu_mem > GPU_IDLE_USED_MEM_THRESHOLD_MIB )); then
      mem_ok=0
      break
    fi
  done < <(gpu_snapshot)

  [[ "$mem_ok" -eq 1 ]] || return 1
  [[ "$(count_compute_pids)" -eq 0 ]]
}

wait_for_idle_gpus() {
  [[ "$WAIT_FOR_IDLE_GPUS" == "1" ]] || return 0

  local warned=0
  until all_gpus_idle; do
    if [[ "$warned" -eq 0 ]]; then
      log "[wait] GPUs busy; waiting for compute slots to clear (mem_threshold=${GPU_IDLE_USED_MEM_THRESHOLD_MIB} MiB)"
      warned=1
    fi
    log "[wait] gpu_snapshot=$(gpu_snapshot | paste -sd ';' -) compute_pids=$(count_compute_pids)"
    sleep "$GPU_IDLE_CHECK_INTERVAL"
  done

  if [[ "$warned" -eq 1 ]]; then
    log "[wait] GPUs look idle; resuming launch"
  fi
}

run_best_ckpt_path() {
  local run_name="$1"
  printf '%s/checkpoints/%s/best.pt\n' "$OUTPUT_DIR" "$run_name"
}

run_meta_path() {
  local run_name="$1"
  printf '%s/logs/%s/run_meta.json\n' "$OUTPUT_DIR" "$run_name"
}

run_is_complete() {
  local run_name="$1"
  [[ -f "$(run_best_ckpt_path "$run_name")" && -f "$(run_meta_path "$run_name")" ]]
}

run_is_running() {
  local run_name="$1"
  pgrep -f -- "run-name ${run_name}" >/dev/null 2>&1
}

wait_for_existing_run_completion() {
  local run_name="$1"
  local waited=0

  while run_is_running "$run_name"; do
    if [[ "$waited" -eq 0 ]]; then
      log "[wait] $run_name is already running; waiting for it to finish"
      waited=1
    fi
    sleep "$GPU_IDLE_CHECK_INTERVAL"
  done

  if [[ "$waited" -eq 0 ]]; then
    return 2
  fi

  if run_is_complete "$run_name"; then
    log "[skip] $run_name finished while waiting"
    return 0
  fi

  log "[error] $run_name stopped without a completion marker (best.pt + run_meta.json)"
  return 1
}

run_training() {
  local run_name="$1"
  local model_config="$2"
  local train_config="$3"
  local master_port="$4"
  local stage_log="$LOG_DIR/${run_name}.log"
  local attempt=1
  local max_attempts=$((RUN_RETRIES_ON_OOM + 1))

  log "[run] $run_name port=$master_port log=$stage_log"

  if [[ "$DRY_RUN" == "1" ]]; then
    log "[dry-run] would run: $run_name"
    return 0
  fi

  while (( attempt <= max_attempts )); do
    wait_for_idle_gpus
    printf '\n[%(%F %T)T] ===== launcher attempt %d/%d for %s =====\n' -1 "$attempt" "$max_attempts" "$run_name" >> "$stage_log"

    # seed and run_name are baked into the train_config YAML
    if torchrun \
      --standalone \
      --nproc_per_node="$NPROC_PER_NODE" \
      --master_port="$master_port" \
      -m transchrombp.training.train_ddp \
      --model-config "$model_config" \
      --train-config "$train_config" \
      --run-name "$run_name" \
      >> "$stage_log" 2>&1; then
      log "[done] $run_name"
      return 0
    fi

    if (( attempt >= max_attempts )); then
      log "[error] $run_name failed after $attempt attempt(s)"
      return 1
    fi

    if grep -q 'CUDA out of memory' "$stage_log"; then
      log "[warn] $run_name hit CUDA OOM on attempt $attempt/$max_attempts; will wait and retry"
      attempt=$((attempt + 1))
      sleep "$GPU_IDLE_CHECK_INTERVAL"
      continue
    fi

    log "[error] $run_name failed with a non-OOM error; not retrying"
    return 1
  done
}

# ============================================================================
# Phase 1: Shared bias pretrain
# ============================================================================
log "========== Phase 1: Shared bias pretrain =========="

BIAS_BEST_CKPT="$OUTPUT_DIR/checkpoints/$BIAS_RUN_NAME/best.pt"

if run_is_complete "$BIAS_RUN_NAME"; then
  log "[skip] bias pretrain already complete: $BIAS_BEST_CKPT"
elif wait_for_existing_run_completion "$BIAS_RUN_NAME"; then
  log "[skip] bias pretrain already completed by another launcher: $BIAS_BEST_CKPT"
else
  wait_status=$?
  if [[ "$wait_status" -eq 1 ]]; then
    exit 1
  fi

  # Generate bias train config with fixed seed and run_name
  BIAS_TRAIN_RT="$CONFIG_DIR/train_bias.yaml"
  python3 - "$TRAIN_BIAS" "$BIAS_TRAIN_RT" "$BIAS_RUN_NAME" "$ROOT_DIR" <<'PY'
import sys, yaml
from pathlib import Path
base_path, out_path, run_name, root_dir = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
with open(base_path) as f:
    cfg = yaml.safe_load(f)
cfg["seed"] = 1234
cfg.setdefault("logging", {})["run_name"] = run_name
data_cfg = cfg.get("data", {}).get("config_path")
if data_cfg and not Path(data_cfg).is_absolute():
    cfg["data"]["config_path"] = str((Path(root_dir) / data_cfg).resolve())
Path(out_path).write_text(yaml.safe_dump(cfg, sort_keys=False))
PY

  # Use V2-full model config for bias (identical bias branch for both models)
  run_training "$BIAS_RUN_NAME" "$MODEL_V2_FULL" "$BIAS_TRAIN_RT" "$MASTER_PORT_BASE"

  if [[ "$DRY_RUN" != "1" ]] && [[ ! -f "$BIAS_BEST_CKPT" ]]; then
    echo "[error] bias pretrain failed — no best.pt at $BIAS_BEST_CKPT" >&2
    exit 1
  fi
fi

log "[bias] best_ckpt=$BIAS_BEST_CKPT"

# ============================================================================
# Phase 2: Main training — 2 models × 3 seeds
# ============================================================================
log "========== Phase 2: Main training (2 models × 3 seeds) =========="

declare -A MODEL_CONFIGS
MODEL_CONFIGS[full]="$MODEL_V2_FULL"
MODEL_CONFIGS[noTF]="$MODEL_V2_NOTF"

port_offset=1
for model_tag in full noTF; do
  base_model_config="${MODEL_CONFIGS[$model_tag]}"

  for seed in "${SEEDS[@]}"; do
    run_name="${ABLATION_TAG}_${model_tag}_s${seed}"
    main_model_config="$CONFIG_DIR/model_main_${model_tag}_s${seed}.yaml"
    main_train_config="$CONFIG_DIR/train_main_${model_tag}_s${seed}.yaml"
    master_port=$((MASTER_PORT_BASE + port_offset))
    port_offset=$((port_offset + 1))

    # Check if already completed
    if run_is_complete "$run_name"; then
      log "[skip] $run_name already complete"
      continue
    elif wait_for_existing_run_completion "$run_name"; then
      continue
    else
      wait_status=$?
      if [[ "$wait_status" -eq 1 ]]; then
        exit 1
      fi
    fi

    # Generate runtime model config with pretrained bias path
    python3 - "$base_model_config" "$main_model_config" "$BIAS_BEST_CKPT" <<'PY'
import sys, yaml
from pathlib import Path

base_cfg_path, out_path, bias_ckpt = sys.argv[1], sys.argv[2], sys.argv[3]
with open(base_cfg_path) as f:
    cfg = yaml.safe_load(f)

cfg.setdefault("bias_branch", {})["pretrained_path"] = bias_ckpt
cfg["bias_branch"]["freeze_bias_core"] = True
cfg.setdefault("fusion", {})["learnable_scales"] = True

Path(out_path).write_text(yaml.safe_dump(cfg, sort_keys=False))
PY

    # Generate runtime train config with correct seed and run_name
    python3 - "$TRAIN_MAIN" "$main_train_config" "$seed" "$run_name" "$ROOT_DIR" <<'PY'
import sys, yaml
from pathlib import Path

base_path, out_path, seed, run_name, root_dir = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
with open(base_path) as f:
    cfg = yaml.safe_load(f)

cfg["seed"] = int(seed)
cfg.setdefault("logging", {})["run_name"] = run_name
data_cfg = cfg.get("data", {}).get("config_path")
if data_cfg and not Path(data_cfg).is_absolute():
    cfg["data"]["config_path"] = str((Path(root_dir) / data_cfg).resolve())

Path(out_path).write_text(yaml.safe_dump(cfg, sort_keys=False))
PY

    # Copy configs to meta dir for reproducibility
    cp "$main_model_config" "$META_DIR/"
    cp "$main_train_config" "$META_DIR/"

    run_training "$run_name" "$main_model_config" "$main_train_config" "$master_port"
  done
done

# ============================================================================
# Phase 3: Summary
# ============================================================================
log "========== Ablation complete =========="
log "Bias checkpoint: $BIAS_BEST_CKPT"
log "Run configs saved to: $META_DIR/"

echo ""
echo "Next steps:"
echo "  1. Evaluate:  bash scripts/evaluate_ablation.sh $ABLATION_TAG"
echo "  2. Analyze:   python scripts/analyze_ablation.py $ABLATION_TAG"
