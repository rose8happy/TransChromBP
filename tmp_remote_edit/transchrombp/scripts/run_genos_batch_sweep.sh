#!/usr/bin/env bash
# ============================================================================
# Genos Phase 1 Quick Batch Sweep
# ============================================================================
# 目的：
#   在 corrected smoke 通过后，用 G1 短程 dry-run 快速比较不同 batch size 的
#   显存与耗时，选出正式 pilot 统一使用的 batch。
#
# 默认策略：
#   - 单卡运行，不启用 DDP
#   - 仅测试 G1（逐位置 Genos 注入）
#   - 每个 batch 运行固定 dry-run steps
#
# Usage:
#   bash scripts/run_genos_batch_sweep.sh
#   GPU=1 BATCH_SIZES="8 12 16 20" DRY_RUN_STEPS=50 bash scripts/run_genos_batch_sweep.sh
#
# Environment variables:
#   GPU                default 0
#   GROUP              default G1
#   SEED               default 42
#   BATCH_SIZES        default "8 12 16 20"
#   DRY_RUN_STEPS      default 50
#   POLL_INTERVAL      default 1
#   STOP_ON_FAILURE    default 1
#   DRY_RUN            default 0
#   BIAS_CKPT          optional override
# ============================================================================
set -euo pipefail

if [ -z "${BASH_VERSION:-}" ]; then
    echo "[error] Please run this script with bash, e.g. 'bash scripts/run_genos_batch_sweep.sh'" >&2
    exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_DIR="/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b"
LOG_DIR="/data1/zhoujiazhen/bylw_atac/logs"
OUTPUT_DIR="$ROOT_DIR/outputs"

export PATH="${VENV_DIR}/bin:${PATH:-}"
export PYTHONPATH="$ROOT_DIR/src:${PYTHONPATH:-}"

cd "$ROOT_DIR"

GPU="${GPU:-0}"
GROUP="${GROUP:-G1}"
SEED="${SEED:-42}"
BATCH_SIZES="${BATCH_SIZES:-8 12 16 20}"
DRY_RUN_STEPS="${DRY_RUN_STEPS:-50}"
POLL_INTERVAL="${POLL_INTERVAL:-1}"
STOP_ON_FAILURE="${STOP_ON_FAILURE:-1}"
DRY_RUN="${DRY_RUN:-0}"
RUN_STAMP="$(date +%Y%m%d_%H%M%S)"

case "${GROUP^^}" in
    G0) MODEL_SRC="$ROOT_DIR/configs/model/v2fix_baseline.yaml" ;;
    G1) MODEL_SRC="$ROOT_DIR/configs/model/v2fix_genos_gate.yaml" ;;
    G2) MODEL_SRC="$ROOT_DIR/configs/model/v2fix_genos_mean.yaml" ;;
    *)
        echo "[error] Unsupported GROUP=$GROUP (expected G0, G1, or G2)" >&2
        exit 1
        ;;
esac

TRAIN_CFG="$ROOT_DIR/configs/train/train_v2fix_genos_profile_select.yaml"
if [[ -z "${BIAS_CKPT:-}" ]]; then
    BIAS_CKPT="$OUTPUT_DIR/checkpoints/ablation_tf_20260318_shared_bias/best.pt"
fi
if [[ ! -f "$BIAS_CKPT" ]]; then
    echo "[error] Bias checkpoint not found: $BIAS_CKPT" >&2
    exit 1
fi

RUNTIME_DIR="$OUTPUT_DIR/runtime/genos_batch_sweep_${RUN_STAMP}"
CONFIG_DIR="$RUNTIME_DIR/configs"
mkdir -p "$LOG_DIR" "$CONFIG_DIR" "$OUTPUT_DIR/logs"

MAIN_LOG="$LOG_DIR/genos_batch_sweep_${RUN_STAMP}.log"
SUMMARY_TSV="$OUTPUT_DIR/logs/genos_batch_sweep_${RUN_STAMP}.tsv"

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$1" | tee -a "$MAIN_LOG"
}

generate_configs() {
    local batch_size="$1"
    local run_name="$2"
    local out_model="$CONFIG_DIR/model_${run_name}.yaml"
    local out_train="$CONFIG_DIR/train_${run_name}.yaml"

    python3 - "$MODEL_SRC" "$out_model" "$BIAS_CKPT" <<'PY'
import sys, yaml
from pathlib import Path
src, dst, bias = sys.argv[1], sys.argv[2], sys.argv[3]
with open(src) as f:
    cfg = yaml.safe_load(f)
cfg.setdefault("bias_branch", {})["pretrained_path"] = bias
cfg["bias_branch"]["freeze_bias_core"] = True
cfg.setdefault("fusion", {})["learnable_scales"] = True
if int(cfg.get("local_tower", {}).get("n_dil_layers", 0)) <= 0:
    raise SystemExit(f"Expected local_tower to be enabled in {src}")
if int(cfg["bias_branch"].get("profile_pool_factor", 0)) <= 0:
    raise SystemExit(f"Expected bias_branch.profile_pool_factor > 0 in {src}")
if not bool(cfg.get("fusion", {}).get("profile_bias_stop_gradient", False)):
    raise SystemExit(f"Expected fusion.profile_bias_stop_gradient=true in {src}")
Path(dst).write_text(yaml.safe_dump(cfg, sort_keys=False))
PY

    python3 - "$TRAIN_CFG" "$out_train" "$SEED" "$run_name" "$ROOT_DIR" "$batch_size" <<'PY'
import sys, yaml
from pathlib import Path
src, dst, seed, name, root, batch = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]
with open(src) as f:
    cfg = yaml.safe_load(f)
cfg["seed"] = int(seed)
cfg.setdefault("logging", {})["run_name"] = name
cfg.setdefault("data", {})["batch_size_per_gpu"] = int(batch)
dp = cfg.get("data", {}).get("config_path")
if dp and not Path(dp).is_absolute():
    cfg["data"]["config_path"] = str((Path(root) / dp).resolve())
Path(dst).write_text(yaml.safe_dump(cfg, sort_keys=False))
PY

    echo "$out_model" "$out_train"
}

sample_gpu() {
    local gpu="$1"
    local util mem
    util="$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits -i "$gpu" | tr -d ' ')"
    mem="$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$gpu" | tr -d ' ')"
    printf '%s,%s,%s\n' "$(date +%s)" "$mem" "$util"
}

printf 'run_name\tgroup\tbatch_size\tdry_run_steps\texit_code\twall_seconds\tpeak_mem_mib\tstage_log\tmonitor_log\n' > "$SUMMARY_TSV"

log "========== Genos Quick Batch Sweep =========="
log "[config] GROUP=$GROUP GPU=$GPU SEED=$SEED BATCH_SIZES=$BATCH_SIZES DRY_RUN_STEPS=$DRY_RUN_STEPS POLL_INTERVAL=$POLL_INTERVAL STOP_ON_FAILURE=$STOP_ON_FAILURE"
log "[config] MODEL_SRC=$MODEL_SRC"
log "[config] BIAS_CKPT=$BIAS_CKPT"
log "[config] MAIN_LOG=$MAIN_LOG"
log "[config] SUMMARY_TSV=$SUMMARY_TSV"

for batch_size in $BATCH_SIZES; do
    run_name="genos_bsweep_${RUN_STAMP}_${GROUP,,}_bs${batch_size}_s${SEED}"
    read -r model_cfg train_cfg < <(generate_configs "$batch_size" "$run_name")
    stage_log="$LOG_DIR/${run_name}.log"
    monitor_log="$LOG_DIR/${run_name}.gpu.csv"

    log "[run] $run_name batch=$batch_size stage_log=$stage_log"

    if [[ "$DRY_RUN" == "1" ]]; then
        log "[dry-run] would run batch=$batch_size model=$model_cfg train=$train_cfg"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$run_name" "$GROUP" "$batch_size" "$DRY_RUN_STEPS" "0" "0" "0" "$stage_log" "$monitor_log" >> "$SUMMARY_TSV"
        continue
    fi

    start_ts="$(date +%s)"
    CUDA_VISIBLE_DEVICES="$GPU" \
    python -m transchrombp.training.train_ddp \
        --model-config "$model_cfg" \
        --train-config "$train_cfg" \
        --run-name "$run_name" \
        --batch-size-per-gpu "$batch_size" \
        --dry-run-steps "$DRY_RUN_STEPS" \
        >> "$stage_log" 2>&1 &
    child_pid=$!

    : > "$monitor_log"
    while kill -0 "$child_pid" 2>/dev/null; do
        sample_gpu "$GPU" >> "$monitor_log"
        sleep "$POLL_INTERVAL"
    done
    if wait "$child_pid"; then
        exit_code=0
    else
        exit_code=$?
    fi
    end_ts="$(date +%s)"
    wall_seconds=$((end_ts - start_ts))
    peak_mem="$(awk -F, 'BEGIN{m=0} NF>=2 && $2+0>m {m=$2+0} END{print m}' "$monitor_log")"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$run_name" "$GROUP" "$batch_size" "$DRY_RUN_STEPS" "$exit_code" "$wall_seconds" "$peak_mem" "$stage_log" "$monitor_log" >> "$SUMMARY_TSV"

    if [[ "$exit_code" -eq 0 ]]; then
        log "[done] $run_name batch=$batch_size wall=${wall_seconds}s peak_mem=${peak_mem}MiB"
    else
        log "[fail] $run_name batch=$batch_size exit_code=$exit_code wall=${wall_seconds}s peak_mem=${peak_mem}MiB"
        if [[ "$STOP_ON_FAILURE" == "1" ]]; then
            log "[stop] STOP_ON_FAILURE=1 -> abort remaining sweep"
            break
        fi
    fi
done

log "========== Sweep Summary =========="
cat "$SUMMARY_TSV" | tee -a "$MAIN_LOG"
log "========== Genos Quick Batch Sweep Complete =========="
