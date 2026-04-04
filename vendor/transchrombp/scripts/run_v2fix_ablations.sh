#!/usr/bin/env bash
# ============================================================================
# V2 Code Improvement Ablation Experiments (v2, 2026-03-20)
# ============================================================================
# Design:
#   Shared bias pretrain from existing ablation_tf experiment.
#   6 experiment groups × 3 seeds each:
#     A: freeze fix only (baseline with fix)
#     B: center count pooling
#     C: profile_scale_init=0.1 (optional, 2nd batch)
#     D: bias GroupNorm (optional, needs own bias pretrain)
#     F: attention pooling for count head (NEW)
#     G: profile refinement conv (NEW)
#
# Usage:
#   # 第一批: A/B/F/G (共享已有 bias checkpoint, 40 epoch)
#   bash scripts/run_v2fix_ablations.sh --batch 1
#
#   # 只跑第一批中的 A/B，且只跑 seed=42（推荐首次启动）
#   bash scripts/run_v2fix_ablations.sh --batch 1 --groups A,B --seeds 42
#
#   # 第二批: C/D (可选)
#   bash scripts/run_v2fix_ablations.sh --batch 2
#
#   # 只跑 soup 评估 (零训练成本)
#   bash scripts/run_v2fix_ablations.sh --soup-only
#
# Environment variables:
#   NPROC_PER_NODE   — GPUs per node (default: 2)
#   DRY_RUN          — set to "1" to only print commands
#   BIAS_CKPT        — path to shared bias pretrain checkpoint (required for batch 1)
#   ABLATION_DATE    — override date tag (default: YYYYMMDD)
#   TRAIN_STRATEGY   — scientific | throughput_smoke | legacy_loss_total
# ============================================================================
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_DIR="/data1/zhoujiazhen/bylw_atac/.mamba/envs/transchrombp"
LOG_DIR="/data1/zhoujiazhen/bylw_atac/logs"
OUTPUT_DIR="$ROOT_DIR/outputs"

export PATH="${ENV_DIR}/bin:/usr/bin:/bin:${PATH:-}"
export LD_LIBRARY_PATH="${ENV_DIR}/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="$ROOT_DIR/src:${PYTHONPATH:-}"

cd "$ROOT_DIR"

NPROC_PER_NODE="${NPROC_PER_NODE:-2}"
DRY_RUN="${DRY_RUN:-0}"
MASTER_PORT_BASE="${MASTER_PORT_BASE:-29550}"
ABLATION_DATE="${ABLATION_DATE:-$(date +%Y%m%d)}"
TRAIN_STRATEGY="${TRAIN_STRATEGY:-scientific}"
SOUP_TOPK_METRIC="${SOUP_TOPK_METRIC:-peak.profile_target_jsd_full_mean}"

SEEDS=(42 1234 2024)
GROUP_FILTER=()

# ---------- Configs ----------
# Model configs for each experiment
MODEL_A="$ROOT_DIR/configs/model/v2fix_baseline.yaml"
MODEL_B="$ROOT_DIR/configs/model/v2fix_center_pool.yaml"
MODEL_C="$ROOT_DIR/configs/model/v2fix_pscale01.yaml"
MODEL_D="$ROOT_DIR/configs/model/v2fix_groupnorm.yaml"
MODEL_F="$ROOT_DIR/configs/model/v2fix_attn_pool.yaml"           # [NEW] 实验 F
MODEL_G="$ROOT_DIR/configs/model/v2fix_profile_refine.yaml"      # [NEW] 实验 G

# Train config strategy
case "$TRAIN_STRATEGY" in
    scientific)
        TRAIN_MAIN="$ROOT_DIR/configs/train/train_ablation_v2_main_profile_select.yaml"
        ;;
    throughput_smoke)
        TRAIN_MAIN="$ROOT_DIR/configs/train/train_ablation_v2_main_throughput_smoke.yaml"
        ;;
    legacy_loss_total)
        TRAIN_MAIN="$ROOT_DIR/configs/train/train_ablation_v2_main.yaml"
        ;;
    *)
        echo "[error] Unsupported TRAIN_STRATEGY=$TRAIN_STRATEGY" >&2
        echo "[error] Expected one of: scientific, throughput_smoke, legacy_loss_total" >&2
        exit 1
        ;;
esac
EVAL_TRAIN_CONFIG="${EVAL_TRAIN_CONFIG:-$TRAIN_MAIN}"

# Bias pretrain (for experiment D only)
TRAIN_BIAS="$ROOT_DIR/configs/train/train_bias_tutorial_teacher_v2.yaml"

# ---------- Runtime dir ----------
TAG="v2fix_${ABLATION_DATE}"
RUNTIME_DIR="$OUTPUT_DIR/runtime/$TAG"
CONFIG_DIR="$RUNTIME_DIR/configs"
META_DIR="$RUNTIME_DIR/meta"
mkdir -p "$LOG_DIR" "$CONFIG_DIR" "$META_DIR"

log() { printf '[%(%F %T)T] %s\n' -1 "$1"; }

run_best_ckpt_path() {
    printf '%s/checkpoints/%s/best.pt\n' "$OUTPUT_DIR" "$1"
}

run_meta_path() {
    printf '%s/logs/%s/run_meta.json\n' "$OUTPUT_DIR" "$1"
}

run_is_complete() {
    [[ -f "$(run_best_ckpt_path "$1")" && -f "$(run_meta_path "$1")" ]]
}

metrics_support_metric() {
    local metrics_jsonl="$1"
    local metric_key="$2"
    python3 - "$metrics_jsonl" "$metric_key" <<'PY'
import json, sys
from pathlib import Path

metrics_path = Path(sys.argv[1])
metric_key = sys.argv[2]
parts = metric_key.split(".", 1)
if not metrics_path.exists():
    raise SystemExit(1)

with metrics_path.open() as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        rec = json.loads(line)
        val = rec.get("val", {})
        if len(parts) == 2:
            value = val.get(parts[0], {}).get(parts[1])
        else:
            value = val.get(parts[0])
        if value is not None:
            raise SystemExit(0)

raise SystemExit(1)
PY
}

pick_soup_metric() {
    local metrics_jsonl="$1"
    local preferred_metric="$2"
    local candidate
    local -a candidates=(
        "$preferred_metric"
        "peak.profile_target_jsd_debiased_mean"
        "peak.loss_profile"
    )
    local seen=" "

    for candidate in "${candidates[@]}"; do
        [[ -z "$candidate" ]] && continue
        if [[ "$seen" == *" $candidate "* ]]; then
            continue
        fi
        seen+="$candidate "
        if metrics_support_metric "$metrics_jsonl" "$candidate"; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

normalize_group_name() {
    local raw="${1,,}"
    case "$raw" in
        a|freeze)
            printf 'freeze\n'
            ;;
        b|cpool|center|center_pool)
            printf 'cpool\n'
            ;;
        c|pscale|pscale01)
            printf 'pscale01\n'
            ;;
        d|gn|groupnorm)
            printf 'gn\n'
            ;;
        f|attnpool|attention|attention_pool)
            printf 'attnpool\n'
            ;;
        g|profref|profile_refine|refine)
            printf 'profref\n'
            ;;
        *)
            return 1
            ;;
    esac
}

parse_seed_csv() {
    local csv="$1"
    local token
    local -a parsed=()
    IFS=',' read -r -a raw_tokens <<< "$csv"
    for token in "${raw_tokens[@]}"; do
        token="${token//[[:space:]]/}"
        [[ -z "$token" ]] && continue
        if [[ ! "$token" =~ ^[0-9]+$ ]]; then
            echo "[error] Invalid seed: $token" >&2
            exit 1
        fi
        parsed+=("$token")
    done
    if [[ ${#parsed[@]} -eq 0 ]]; then
        echo "[error] --seeds requires at least one integer seed" >&2
        exit 1
    fi
    SEEDS=("${parsed[@]}")
}

parse_group_csv() {
    local csv="$1"
    local token normalized
    local -a parsed=()
    IFS=',' read -r -a raw_tokens <<< "$csv"
    for token in "${raw_tokens[@]}"; do
        token="${token//[[:space:]]/}"
        [[ -z "$token" ]] && continue
        if ! normalized="$(normalize_group_name "$token")"; then
            echo "[error] Invalid group selector: $token" >&2
            echo "[error] Supported groups: A,B,C,D,F,G or freeze,cpool,pscale01,gn,attnpool,profref" >&2
            exit 1
        fi
        if [[ " ${parsed[*]} " != *" $normalized "* ]]; then
            parsed+=("$normalized")
        fi
    done
    if [[ ${#parsed[@]} -eq 0 ]]; then
        echo "[error] --groups requires at least one valid group" >&2
        exit 1
    fi
    GROUP_FILTER=("${parsed[@]}")
}

group_is_selected() {
    local group="$1"
    local selected
    if [[ ${#GROUP_FILTER[@]} -eq 0 ]]; then
        return 0
    fi
    for selected in "${GROUP_FILTER[@]}"; do
        if [[ "$selected" == "$group" ]]; then
            return 0
        fi
    done
    return 1
}

generate_runtime_configs() {
    local model_config="$1"
    local run_name="$2"
    local seed="$3"
    local bias_ckpt="$4"
    local out_model="$CONFIG_DIR/model_${run_name}.yaml"
    local out_train="$CONFIG_DIR/train_${run_name}.yaml"

    # Inject bias pretrain path + freeze into model config
    python3 - "$model_config" "$out_model" "$bias_ckpt" <<'PY'
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

    # Inject seed + run_name into train config
    python3 - "$TRAIN_MAIN" "$out_train" "$seed" "$run_name" "$ROOT_DIR" <<'PY'
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

    echo "$out_model" "$out_train"
}

run_training() {
    local run_name="$1"
    local model_config="$2"
    local train_config="$3"
    local master_port="$4"
    local stage_log="$LOG_DIR/${run_name}.log"

    if run_is_complete "$run_name"; then
        log "[skip] $run_name already complete"
        return 0
    fi

    log "[run] $run_name port=$master_port log=$stage_log"

    if [[ "$DRY_RUN" == "1" ]]; then
        log "[dry-run] would run: $run_name"
        return 0
    fi

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

    log "[error] $run_name failed"
    return 1
}

run_experiment_group() {
    local exp_tag="$1"     # e.g. "freeze", "cpool", "pscale01"
    local model_config="$2"
    local bias_ckpt="$3"
    local port_start="$4"

    log "========== Experiment: $exp_tag =========="
    local port_offset=0

    for seed in "${SEEDS[@]}"; do
        local run_name="${TAG}_${exp_tag}_s${seed}"
        local master_port=$((port_start + port_offset))
        port_offset=$((port_offset + 1))

        read -r runtime_model runtime_train < <(
            generate_runtime_configs "$model_config" "$run_name" "$seed" "$bias_ckpt"
        )

        cp "$runtime_model" "$META_DIR/"
        cp "$runtime_train" "$META_DIR/"

        run_training "$run_name" "$runtime_model" "$runtime_train" "$master_port"
    done
}

# ============================================================================
# Soup evaluation (实验 E)
# ============================================================================
run_soup_eval() {
    log "========== Experiment E: Checkpoint Soup =========="

    # 对每个已完成的 baseline run 做 soup
    for seed in "${SEEDS[@]}"; do
        local baseline_run="ablation_tf_20260318_full_s${seed}"
        local ckpt_dir="$OUTPUT_DIR/checkpoints/$baseline_run"
        local metrics_file="$OUTPUT_DIR/logs/$baseline_run/epoch_metrics.jsonl"
        local eval_dir="$OUTPUT_DIR/eval"
        mkdir -p "$eval_dir"

        if [[ ! -d "$ckpt_dir" ]]; then
            log "[skip] $baseline_run: no checkpoint dir"
            continue
        fi

        # Soup: last 5
        local soup_last5="$ckpt_dir/soup_last5.pt"
        if [[ ! -f "$soup_last5" ]]; then
            log "[soup] $baseline_run: creating soup_last5"
            python3 scripts/checkpoint_soup.py \
                --ckpt-dir "$ckpt_dir" \
                --mode last_k --k 5 \
                --output "$soup_last5" || log "[warn] soup_last5 failed for $baseline_run"
        fi

        # Soup: top 3 by profile-facing metric, with fallback for older runs
        local soup_top3="$ckpt_dir/soup_top3.pt"
        if [[ ! -f "$soup_top3" ]] && [[ -f "$metrics_file" ]]; then
            local selected_metric=""
            if selected_metric="$(pick_soup_metric "$metrics_file" "$SOUP_TOPK_METRIC")"; then
                if [[ "$selected_metric" == "$SOUP_TOPK_METRIC" ]]; then
                    log "[soup] $baseline_run: creating soup_top3 metric=$selected_metric"
                else
                    log "[soup] $baseline_run: fallback metric=$selected_metric (preferred=$SOUP_TOPK_METRIC missing)"
                fi
                python3 scripts/checkpoint_soup.py \
                    --ckpt-dir "$ckpt_dir" \
                    --mode top_k_jsd --k 3 \
                    --metrics-jsonl "$metrics_file" \
                    --metric-key "$selected_metric" \
                    --output "$soup_top3" || log "[warn] soup_top3 failed for $baseline_run"
            else
                log "[skip] $baseline_run: metrics file lacks all supported soup metrics"
            fi
        fi

        # Evaluate all variants
        local model_cfg="$ROOT_DIR/configs/model/v2fix_baseline.yaml"
        for variant in best soup_last5 soup_top3; do
            local ckpt_path="$ckpt_dir/${variant}.pt"
            local eval_json="$eval_dir/${baseline_run}_${variant}.json"

            if [[ ! -f "$ckpt_path" ]]; then
                log "[skip] $ckpt_path not found"
                continue
            fi
            if [[ -f "$eval_json" ]]; then
                log "[skip] $eval_json already exists"
                continue
            fi

            log "[eval] $baseline_run $variant"
            if [[ "$DRY_RUN" == "1" ]]; then
                log "[dry-run] would evaluate $ckpt_path"
                continue
            fi
            python3 scripts/evaluate_checkpoint.py \
                --model-config "$model_cfg" \
                --train-config "$EVAL_TRAIN_CONFIG" \
                --split valid \
                --checkpoint "$ckpt_path" \
                --output-json "$eval_json" || log "[warn] eval failed for $ckpt_path"
        done
    done
}

# ============================================================================
# Parse arguments
# ============================================================================
BATCH=""
SOUP_ONLY=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --batch)
            BATCH="$2"; shift 2 ;;
        --groups)
            parse_group_csv "$2"; shift 2 ;;
        --seeds)
            parse_seed_csv "$2"; shift 2 ;;
        --soup-only)
            SOUP_ONLY=true; shift ;;
        *)
            echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# ============================================================================
# Main
# ============================================================================

if $SOUP_ONLY; then
    log "[config] SOUP_ONLY=true SEEDS=${SEEDS[*]} GROUPS=${GROUP_FILTER[*]:-all}"
    run_soup_eval
    log "Soup evaluation complete."
    exit 0
fi

if [[ "$BATCH" == "1" || -z "$BATCH" ]]; then
    # 第一批: A/B/F/G — 需要 BIAS_CKPT, 40 epoch
    log "[config] TRAIN_STRATEGY=$TRAIN_STRATEGY TRAIN_MAIN=$TRAIN_MAIN SEEDS=${SEEDS[*]} GROUPS=${GROUP_FILTER[*]:-all}"
    if [[ -z "${BIAS_CKPT:-}" ]]; then
        # 尝试自动发现
        BIAS_CKPT="$OUTPUT_DIR/checkpoints/ablation_tf_20260318_shared_bias/best.pt"
        if [[ ! -f "$BIAS_CKPT" ]]; then
            echo "[error] BIAS_CKPT not set and auto-discovery failed." >&2
            echo "[error] Set BIAS_CKPT=/path/to/shared_bias/best.pt" >&2
            exit 1
        fi
        log "[auto] Using bias checkpoint: $BIAS_CKPT"
    fi

    # Validate shared prerequisite
    for f in "$TRAIN_MAIN"; do
        if [[ ! -f "$f" ]]; then
            echo "[error] Missing config: $f" >&2
            exit 1
        fi
    done

    # Run experiments in priority order: B (center pool) -> A (freeze fix) -> F -> G
    batch1_selected=0
    for spec in \
        "cpool:$MODEL_B:$MASTER_PORT_BASE" \
        "freeze:$MODEL_A:$((MASTER_PORT_BASE + 3))" \
        "attnpool:$MODEL_F:$((MASTER_PORT_BASE + 6))" \
        "profref:$MODEL_G:$((MASTER_PORT_BASE + 9))"; do
        IFS=':' read -r exp_tag model_cfg port_start <<< "$spec"
        if ! group_is_selected "$exp_tag"; then
            log "[skip] batch1 group $exp_tag filtered out"
            continue
        fi
        batch1_selected=1
        if [[ ! -f "$model_cfg" ]]; then
            echo "[error] Missing config: $model_cfg" >&2
            exit 1
        fi
        run_experiment_group "$exp_tag" "$model_cfg" "$BIAS_CKPT" "$port_start"
    done

    if [[ "$batch1_selected" != "1" ]]; then
        log "[skip] No selected groups for batch 1"
    fi

    log "========== Batch 1 complete (A/B/F/G) =========="
fi

if [[ "$BATCH" == "2" ]]; then
    # 第二批 (可选): C (pscale) + D (GN, 需先训练 bias)
    log "[config] TRAIN_STRATEGY=$TRAIN_STRATEGY TRAIN_MAIN=$TRAIN_MAIN SEEDS=${SEEDS[*]} GROUPS=${GROUP_FILTER[*]:-all}"
    if [[ -z "${BIAS_CKPT:-}" ]]; then
        BIAS_CKPT="$OUTPUT_DIR/checkpoints/ablation_tf_20260318_shared_bias/best.pt"
        if [[ ! -f "$BIAS_CKPT" ]]; then
            echo "[error] BIAS_CKPT not set and auto-discovery failed." >&2
            exit 1
        fi
        log "[auto] Using bias checkpoint: $BIAS_CKPT"
    fi

    batch2_selected=0

    # C: profile_scale_init=0.1
    if group_is_selected "pscale01"; then
        batch2_selected=1
        if [[ ! -f "$MODEL_C" ]]; then
            echo "[error] Missing config: $MODEL_C" >&2
            exit 1
        fi
        run_experiment_group "pscale01" "$MODEL_C" "$BIAS_CKPT" "$MASTER_PORT_BASE"
    fi

    # D: GroupNorm — 需要先训练 GN bias branch
    if group_is_selected "gn"; then
        batch2_selected=1
        if [[ ! -f "$MODEL_D" ]]; then
            echo "[error] Missing config: $MODEL_D" >&2
            exit 1
        fi
        GN_BIAS_RUN="${TAG}_gn_shared_bias"
        GN_BIAS_CKPT="$(run_best_ckpt_path "$GN_BIAS_RUN")"

        if run_is_complete "$GN_BIAS_RUN"; then
            log "[skip] GN bias pretrain already complete: $GN_BIAS_CKPT"
        else
            log "========== Phase D.1: GroupNorm bias pretrain =========="

            BIAS_TRAIN_RT="$CONFIG_DIR/train_gn_bias.yaml"
            python3 - "$TRAIN_BIAS" "$BIAS_TRAIN_RT" "$GN_BIAS_RUN" "$ROOT_DIR" <<'PY'
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

            run_training "$GN_BIAS_RUN" "$MODEL_D" "$BIAS_TRAIN_RT" "$MASTER_PORT_BASE"

            if [[ "$DRY_RUN" != "1" ]] && [[ ! -f "$GN_BIAS_CKPT" ]]; then
                echo "[error] GN bias pretrain failed — no best.pt" >&2
                exit 1
            fi
        fi

        log "========== Phase D.2: GroupNorm main training =========="
        run_experiment_group "gn" "$MODEL_D" "$GN_BIAS_CKPT" $((MASTER_PORT_BASE + 4))
    fi

    if [[ "$batch2_selected" != "1" ]]; then
        log "[skip] No selected groups for batch 2"
    fi

    log "========== Batch 2 complete (C/D) =========="
fi

# Print summary
log "========== All requested experiments complete =========="
log "Runtime configs: $META_DIR/"
echo ""
echo "Next steps:"
echo "  1. Evaluate soup:   bash scripts/run_v2fix_ablations.sh --soup-only"
echo "  2. Compare results: python scripts/compare_v2fix_results.py $TAG"
