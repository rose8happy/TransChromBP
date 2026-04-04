#!/usr/bin/env bash
# ============================================================================
# Genos Phase 1 Minimal Pilot — G0 / G1 / G2
# ============================================================================
# 用 genos-1.2b venv 在 A6000 上运行。
# Phase 1 默认不启用 2-rank DDP；双卡通过“每卡一条单卡 run”并行利用，
# 这样可以保持 global batch 恒定，避免把实验口径和调度策略混在一起。
#
# Usage:
#   # Smoke test (默认 2 步，验证链路；run_name 与正式 pilot 分离)
#   bash scripts/run_genos_pilot.sh --smoke
#
#   # G0 baseline only
#   bash scripts/run_genos_pilot.sh --groups G0
#
#   # G1 + G2 并行 (分别占一张 A6000；当前正式 pilot 默认每卡 batch=20)
#   bash scripts/run_genos_pilot.sh --groups G1 --gpu 0 &
#   bash scripts/run_genos_pilot.sh --groups G2 --gpu 1 &
#   wait
#
#   # 全部串行
#   bash scripts/run_genos_pilot.sh --groups G0,G1,G2
#
# Environment variables:
#   DRY_RUN            — set to "1" to only print commands
#   BIAS_CKPT          — path to shared bias pretrain checkpoint
#   BATCH_SIZE_PER_GPU — default 20 (selected by 2026-03-21 G1 quick sweep);
#                         if later testing DDP separately, lower this to keep
#                         global batch matched; if rerunning smoke, override
#                         explicitly rather than relying on this default
#   SMOKE_STEPS        — default 2
#   FORCE_RERUN        — set to "1" to ignore existing best/meta outputs
# ============================================================================
set -euo pipefail

if [ -z "${BASH_VERSION:-}" ]; then
    echo "[error] Please run this script with bash, e.g. 'bash scripts/run_genos_pilot.sh --groups G0'" >&2
    exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_DIR="/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b"
LOG_DIR="/data1/zhoujiazhen/bylw_atac/logs"
OUTPUT_DIR="$ROOT_DIR/outputs"

export PATH="${VENV_DIR}/bin:${PATH:-}"
export PYTHONPATH="$ROOT_DIR/src:${PYTHONPATH:-}"

cd "$ROOT_DIR"

DRY_RUN="${DRY_RUN:-0}"
RUN_DATE="$(date +%Y%m%d)"
SEED="${SEED:-42}"
GPU="${GPU:-0}"
BATCH_SIZE_PER_GPU="${BATCH_SIZE_PER_GPU:-20}"
SMOKE_STEPS="${SMOKE_STEPS:-2}"
FORCE_RERUN="${FORCE_RERUN:-0}"
SMOKE=false
RUN_GROUPS=()

# ---------- Configs ----------
MODEL_G0="$ROOT_DIR/configs/model/v2fix_baseline.yaml"
MODEL_G1="$ROOT_DIR/configs/model/v2fix_genos_gate.yaml"
MODEL_G2="$ROOT_DIR/configs/model/v2fix_genos_mean.yaml"
TRAIN_CFG="$ROOT_DIR/configs/train/train_v2fix_genos_profile_select.yaml"

# ---------- Argument parsing ----------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --smoke)     SMOKE=true; shift ;;
        --groups)    IFS=',' read -ra RUN_GROUPS <<< "$2"; shift 2 ;;
        --gpu)       GPU="$2"; shift 2 ;;
        --seed)      SEED="$2"; shift 2 ;;
        *)           echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

if [[ ${#RUN_GROUPS[@]} -eq 0 ]]; then
    RUN_GROUPS=(G0 G1 G2)
fi

# ---------- Bias checkpoint ----------
if [[ -z "${BIAS_CKPT:-}" ]]; then
    BIAS_CKPT="$OUTPUT_DIR/checkpoints/ablation_tf_20260318_shared_bias/best.pt"
fi
if [[ ! -f "$BIAS_CKPT" ]]; then
    echo "[error] Bias checkpoint not found: $BIAS_CKPT" >&2
    exit 1
fi

# ---------- Runtime config generation ----------
RUN_LABEL="$RUN_DATE"
if $SMOKE; then
    RUN_LABEL="smoke_${RUN_DATE}"
fi

RUNTIME_DIR="$OUTPUT_DIR/runtime/genos_pilot_${RUN_LABEL}"
CONFIG_DIR="$RUNTIME_DIR/configs"
mkdir -p "$LOG_DIR" "$CONFIG_DIR"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$1"; }

generate_configs() {
    local group="$1"
    local run_name="$2"
    local model_src="$3"
    local out_model="$CONFIG_DIR/model_${run_name}.yaml"
    local out_train="$CONFIG_DIR/train_${run_name}.yaml"

    python3 - "$model_src" "$out_model" "$BIAS_CKPT" <<'PY'
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

    python3 - "$TRAIN_CFG" "$out_train" "$SEED" "$run_name" "$ROOT_DIR" "$BATCH_SIZE_PER_GPU" <<'PY'
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

run_experiment() {
    local run_name="$1"
    local model_cfg="$2"
    local train_cfg="$3"
    local stage_log="$LOG_DIR/${run_name}.log"
    local dry_flag=""

    if $SMOKE; then
        dry_flag="--dry-run-steps ${SMOKE_STEPS}"
    fi

    # Check if already complete
    local best_pt="$OUTPUT_DIR/checkpoints/$run_name/best.pt"
    local meta_json="$OUTPUT_DIR/logs/$run_name/run_meta.json"
    if [[ "$FORCE_RERUN" != "1" ]] && [[ -f "$best_pt" && -f "$meta_json" ]] && ! $SMOKE; then
        log "[skip] $run_name already complete"
        return 0
    fi

    log "[run] $run_name gpu=$GPU log=$stage_log"

    if [[ "$DRY_RUN" == "1" ]]; then
        log "[dry-run] would run: CUDA_VISIBLE_DEVICES=$GPU python -m transchrombp.training.train_ddp ..."
        return 0
    fi

    CUDA_VISIBLE_DEVICES="$GPU" \
    python -m transchrombp.training.train_ddp \
        --model-config "$model_cfg" \
        --train-config "$train_cfg" \
        --run-name "$run_name" \
        --batch-size-per-gpu "$BATCH_SIZE_PER_GPU" \
        $dry_flag \
        >> "$stage_log" 2>&1

    log "[done] $run_name"
}

# ---------- Main ----------
log "========== Genos Phase 1 Pilot =========="
log "[config] GROUPS=${RUN_GROUPS[*]} GPU=$GPU SEED=$SEED SMOKE=$SMOKE"
log "[config] RUN_LABEL=$RUN_LABEL BATCH_SIZE_PER_GPU=$BATCH_SIZE_PER_GPU SMOKE_STEPS=$SMOKE_STEPS FORCE_RERUN=$FORCE_RERUN"
log "[config] BIAS_CKPT=$BIAS_CKPT"
log "[config] VENV=$VENV_DIR"
log "[config] RUNTIME_DIR=$RUNTIME_DIR"

for group in "${RUN_GROUPS[@]}"; do
    case "${group^^}" in
        G0)
            run_name="genos_${RUN_LABEL}_baseline_s${SEED}"
            model_src="$MODEL_G0"
            ;;
        G1)
            run_name="genos_${RUN_LABEL}_gate_s${SEED}"
            model_src="$MODEL_G1"
            ;;
        G2)
            run_name="genos_${RUN_LABEL}_mean_s${SEED}"
            model_src="$MODEL_G2"
            ;;
        *)
            echo "[error] Unknown group: $group (expected G0, G1, G2)" >&2
            exit 1
            ;;
    esac

    read -r runtime_model runtime_train < <(generate_configs "$group" "$run_name" "$model_src")
    run_experiment "$run_name" "$runtime_model" "$runtime_train"
done

log "========== Genos pilot complete =========="
