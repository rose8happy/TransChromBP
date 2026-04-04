#!/usr/bin/env bash
# ============================================================================
# Evaluate all ablation checkpoints on held-out test set (chr1)
# ============================================================================
# Usage:
#   bash scripts/evaluate_ablation.sh ablation_tf_20260318
#
# Runs evaluate_checkpoint.py for each of the 6 best.pt checkpoints,
# producing standardized test JSON files.
# ============================================================================
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_DIR="/data1/zhoujiazhen/bylw_atac/.mamba/envs/transchrombp"
OUTPUT_DIR="$ROOT_DIR/outputs"
METRICS_DIR="$OUTPUT_DIR/metrics"

export PATH="/usr/bin:/bin:${ENV_DIR}/bin:${PATH:-}"
export LD_LIBRARY_PATH="${ENV_DIR}/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="$ROOT_DIR/src:${PYTHONPATH:-}"
cd "$ROOT_DIR"

ABLATION_TAG="${1:?Usage: $0 <ablation_tag> e.g. ablation_tf_20260318}"
NONPEAK_RATIO="${NONPEAK_RATIO:-0.1}"
DEVICE="${DEVICE:-cuda:0}"

SEEDS=(42 1234 2024)
MODELS=(full noTF)

mkdir -p "$METRICS_DIR"

log() { printf '[%(%F %T)T] %s\n' -1 "$1"; }

for model_tag in "${MODELS[@]}"; do
  for seed in "${SEEDS[@]}"; do
    run_name="${ABLATION_TAG}_${model_tag}_s${seed}"
    ckpt="$OUTPUT_DIR/checkpoints/$run_name/best.pt"
    out_json="$METRICS_DIR/${run_name}_test.json"

    if [[ ! -f "$ckpt" ]]; then
      log "[warn] checkpoint missing: $ckpt — skipping"
      continue
    fi

    if [[ -f "$out_json" ]]; then
      log "[skip] already evaluated: $out_json"
      continue
    fi

    log "[eval] $run_name → $out_json"
    python -u -m transchrombp.evaluation.evaluate_checkpoint \
      --checkpoint "$ckpt" \
      --nonpeak-ratio "$NONPEAK_RATIO" \
      --device "$DEVICE" \
      --output "$out_json"
    log "[done] $run_name"
  done
done

log "========== All evaluations complete =========="
echo "Results in: $METRICS_DIR/${ABLATION_TAG}_*_test.json"
echo "Run analysis: python scripts/analyze_ablation.py $ABLATION_TAG"
