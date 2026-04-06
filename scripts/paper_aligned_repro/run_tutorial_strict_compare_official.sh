#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

function usage() {
  cat <<'USAGE'
Usage:
  bash scripts/paper_aligned_repro/run_tutorial_strict_compare_official.sh \
    --mode fidelity|controlled \
    --work-root /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare \
    --seed 42 \
    --folds "fold_0" \
    --gpus 0,1

Optional shared-region overrides:
  --run-suffix "_L3"
  --shared-peaks /path/to/shared_filtered_peaks.bed
  --shared-nonpeaks /path/to/shared_filtered_nonpeaks.bed
  --shared-bigwig /path/to/shared_unstranded.bw
  --official-root /path/to/official/chrombpnet

Modes:
  fidelity
    Official-faithful anchor on tutorial data. Uses a larger batch and default-ish early stop.

  controlled
    Fair-comparison official arm on tutorial data. Locks global batch to 32 across 2 GPUs
    and saves every epoch checkpoint for external best-epoch selection.
USAGE
}

MODE="controlled"
WORK_ROOT="/data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare"
SEED="42"
FOLDS="fold_0"
GPUS="0,1"
MAX_PARALLEL="2"
NAME=""
MULTI_GPU_TRAIN="0"
RUN_SUFFIX=""
SHARED_PEAKS=""
SHARED_NONPEAKS=""
SHARED_BIGWIG=""
OFFICIAL_ROOT="${CHROMBPNET_OFFICIAL_ROOT:-}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    --work-root) WORK_ROOT="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --folds) FOLDS="$2"; shift 2 ;;
    --gpus) GPUS="$2"; shift 2 ;;
    --max-parallel) MAX_PARALLEL="$2"; shift 2 ;;
    --name) NAME="$2"; shift 2 ;;
    --run-suffix) RUN_SUFFIX="$2"; shift 2 ;;
    --shared-peaks) SHARED_PEAKS="$2"; shift 2 ;;
    --shared-nonpeaks) SHARED_NONPEAKS="$2"; shift 2 ;;
    --shared-bigwig) SHARED_BIGWIG="$2"; shift 2 ;;
    --official-root) OFFICIAL_ROOT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

case "${MODE}" in
  fidelity)
    BATCH_SIZE="64"
    EPOCHS="50"
    EARLY_STOP="5"
    SAVE_ALL_CHECKPOINTS="0"
    MULTI_GPU_TRAIN="0"
    RUN_NAME_DEFAULT="tutorial_official_fidelity_s${SEED}"
    ;;
  controlled)
    BATCH_SIZE="32"
    EPOCHS="50"
    EARLY_STOP="50"
    SAVE_ALL_CHECKPOINTS="1"
    MULTI_GPU_TRAIN="1"
    MAX_PARALLEL="1"
    RUN_NAME_DEFAULT="tutorial_official_controlled_s${SEED}"
    ;;
  *)
    echo "ERROR: --mode must be fidelity or controlled" >&2
    exit 1
    ;;
esac

if [[ -n "${SHARED_PEAKS}" && -z "${SHARED_NONPEAKS}" ]] || [[ -z "${SHARED_PEAKS}" && -n "${SHARED_NONPEAKS}" ]]; then
  echo "ERROR: --shared-peaks and --shared-nonpeaks must be provided together" >&2
  exit 1
fi

if [[ -n "${SHARED_PEAKS}" && ! -f "${SHARED_PEAKS}" ]]; then
  echo "ERROR: shared peaks not found: ${SHARED_PEAKS}" >&2
  exit 1
fi

if [[ -n "${SHARED_NONPEAKS}" && ! -f "${SHARED_NONPEAKS}" ]]; then
  echo "ERROR: shared nonpeaks not found: ${SHARED_NONPEAKS}" >&2
  exit 1
fi

if [[ -n "${SHARED_BIGWIG}" && ! -f "${SHARED_BIGWIG}" ]]; then
  echo "ERROR: shared bigwig not found: ${SHARED_BIGWIG}" >&2
  exit 1
fi

if [[ -z "${OFFICIAL_ROOT}" ]]; then
  echo "ERROR: missing official ChromBPNet root; pass --official-root or set CHROMBPNET_OFFICIAL_ROOT" >&2
  exit 1
fi

if [[ -z "${NAME}" ]]; then
  NAME="${RUN_NAME_DEFAULT}${RUN_SUFFIX}"
fi

PEAKS_PATH="/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/overlap.bed.gz"
if [[ -n "${SHARED_PEAKS}" ]]; then
  PEAKS_PATH="${SHARED_PEAKS}"
fi

FOLD_DIR="${WORK_ROOT}/folds"
mkdir -p "${WORK_ROOT}" "${WORK_ROOT}/logs"

if [[ ! -f "${FOLD_DIR}/fold_0.json" ]]; then
  python3 "${SCRIPT_DIR}/generate_fold_set.py" \
    --chrom-sizes /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes \
    --output-dir "${FOLD_DIR}"
fi

CMD=(
  bash "${SCRIPT_DIR}/run_paper_aligned_fast_1seed.sh"
  --name "${NAME}"
  --genome /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.fa
  --chrom-sizes /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes
  --bam /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged.bam
  --peaks "${PEAKS_PATH}"
  --blacklist /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/blacklist.bed.gz
  --fold-dir "${FOLD_DIR}"
  --work-root "${WORK_ROOT}"
  --seed "${SEED}"
  --gpus "${GPUS}"
  --max-parallel "${MAX_PARALLEL}"
  --batch-size "${BATCH_SIZE}"
  --predict-batch-size 512
  --epochs "${EPOCHS}"
  --early-stop "${EARLY_STOP}"
  --folds "${FOLDS}"
)

if [[ "${SAVE_ALL_CHECKPOINTS}" == "1" ]]; then
  CMD+=(--save-all-checkpoints)
fi
if [[ "${MULTI_GPU_TRAIN:-0}" == "1" ]]; then
  CMD+=(--multi-gpu-train)
fi
if [[ -n "${SHARED_NONPEAKS}" ]]; then
  CMD+=(--external-nonpeaks "${SHARED_NONPEAKS}" --use-input-peaks-for-eval)
fi
if [[ -n "${SHARED_BIGWIG}" ]]; then
  CMD+=(--eval-bigwig "${SHARED_BIGWIG}")
fi
CMD+=(--official-root "${OFFICIAL_ROOT}")

echo "[strict-compare] mode=${MODE}"
echo "[strict-compare] name=${NAME}"
echo "[strict-compare] seed=${SEED}"
echo "[strict-compare] folds=${FOLDS}"
echo "[strict-compare] gpus=${GPUS}"
echo "[strict-compare] batch_size=${BATCH_SIZE}"
echo "[strict-compare] epochs=${EPOCHS}"
echo "[strict-compare] early_stop=${EARLY_STOP}"
echo "[strict-compare] max_parallel=${MAX_PARALLEL}"
echo "[strict-compare] multi_gpu_train=${MULTI_GPU_TRAIN:-0}"
echo "[strict-compare] work_root=${WORK_ROOT}"
echo "[strict-compare] shared_peaks=${SHARED_PEAKS:-<none>}"
echo "[strict-compare] shared_nonpeaks=${SHARED_NONPEAKS:-<none>}"
echo "[strict-compare] shared_bigwig=${SHARED_BIGWIG:-<none>}"
echo "[strict-compare] official_root=${OFFICIAL_ROOT}"

"${CMD[@]}"
