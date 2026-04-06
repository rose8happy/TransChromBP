#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

function usage() {
  cat <<'USAGE'
Usage:
  bash scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh \
    --name K562_ATAC \
    --genome /path/to/hg38.fa \
    --chrom-sizes /path/to/hg38.chrom.sizes \
    --bam /path/to/merged.bam \
    --peaks /path/to/overlap.bed.gz \
    --blacklist /path/to/blacklist.bed.gz \
    --fold-dir /path/to/folds \
    --work-root /path/to/chrombpnet_paper_repro \
    --seed 1234 \
    --gpus 0,1

Required:
  --name
  --genome
  --chrom-sizes
  --bam
  --peaks
  --blacklist
  --fold-dir
  --work-root

Optional:
  --data-type ATAC|DNASE     (default: ATAC)
  --seed INT                 (default: 1234)
  --gpus GPU_IDS             (default: 0,1)
  --max-parallel INT         (default: number of gpus)
  --bias-threshold-factor X  (default: 0.5)
  --batch-size INT           (default: 64)
  --predict-batch-size INT   (default: 512)
  --epochs INT               (default: 50)
  --early-stop INT           (default: 5)
  --external-nonpeaks BED    Reuse a precomputed nonpeaks BED and skip `prep nonpeaks`
  --eval-bigwig BIGWIG       Use this bigWig for evaluation instead of the run-generated one
  --use-input-peaks-for-eval Evaluate with the input `--peaks` file instead of `filtered.peaks.bed`
  --save-all-checkpoints     Also save one checkpoint per epoch for external best-epoch selection
  --multi-gpu-train          Train bias/chrombpnet with MirroredStrategy across all GPUs in --gpus
  --official-root PATH       Official ChromBPNet repo root (default: \$CHROMBPNET_OFFICIAL_ROOT)
  --folds "fold_0 fold_1"    (default: all fold_*.json in fold-dir)
USAGE
}

function require_cmd() {
  local cmd="$1"
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "ERROR: command not found: ${cmd}" >&2
    exit 1
  fi
}

NAME=""
GENOME=""
CHROM_SIZES=""
BAM=""
PEAKS=""
BLACKLIST=""
FOLD_DIR=""
WORK_ROOT=""
DATA_TYPE="ATAC"
SEED="1234"
GPUS="0,1"
MAX_PARALLEL=""
BIAS_THRESHOLD_FACTOR="0.5"
BATCH_SIZE="64"
PREDICT_BATCH_SIZE="512"
EPOCHS="50"
EARLY_STOP="5"
FOLDS=""
SAVE_ALL_CHECKPOINTS="0"
MULTI_GPU_TRAIN="0"
EXTERNAL_NONPEAKS=""
EVAL_BIGWIG=""
USE_INPUT_PEAKS_FOR_EVAL="0"
OFFICIAL_ROOT="${CHROMBPNET_OFFICIAL_ROOT:-}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --name) NAME="$2"; shift 2 ;;
    --genome) GENOME="$2"; shift 2 ;;
    --chrom-sizes) CHROM_SIZES="$2"; shift 2 ;;
    --bam) BAM="$2"; shift 2 ;;
    --peaks) PEAKS="$2"; shift 2 ;;
    --blacklist) BLACKLIST="$2"; shift 2 ;;
    --fold-dir) FOLD_DIR="$2"; shift 2 ;;
    --work-root) WORK_ROOT="$2"; shift 2 ;;
    --data-type) DATA_TYPE="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --gpus) GPUS="$2"; shift 2 ;;
    --max-parallel) MAX_PARALLEL="$2"; shift 2 ;;
    --bias-threshold-factor) BIAS_THRESHOLD_FACTOR="$2"; shift 2 ;;
    --batch-size) BATCH_SIZE="$2"; shift 2 ;;
    --predict-batch-size) PREDICT_BATCH_SIZE="$2"; shift 2 ;;
    --epochs) EPOCHS="$2"; shift 2 ;;
    --early-stop) EARLY_STOP="$2"; shift 2 ;;
    --external-nonpeaks) EXTERNAL_NONPEAKS="$2"; shift 2 ;;
    --eval-bigwig) EVAL_BIGWIG="$2"; shift 2 ;;
    --use-input-peaks-for-eval) USE_INPUT_PEAKS_FOR_EVAL="1"; shift 1 ;;
    --save-all-checkpoints) SAVE_ALL_CHECKPOINTS="1"; shift 1 ;;
    --multi-gpu-train) MULTI_GPU_TRAIN="1"; shift 1 ;;
    --official-root) OFFICIAL_ROOT="$2"; shift 2 ;;
    --folds) FOLDS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "${NAME}" || -z "${GENOME}" || -z "${CHROM_SIZES}" || -z "${BAM}" || -z "${PEAKS}" || -z "${BLACKLIST}" || -z "${FOLD_DIR}" || -z "${WORK_ROOT}" ]]; then
  echo "ERROR: missing required arguments" >&2
  usage
  exit 1
fi

if [[ -n "${EXTERNAL_NONPEAKS}" && ! -f "${EXTERNAL_NONPEAKS}" ]]; then
  echo "ERROR: external nonpeaks not found: ${EXTERNAL_NONPEAKS}" >&2
  exit 1
fi

if [[ -n "${EVAL_BIGWIG}" && ! -f "${EVAL_BIGWIG}" ]]; then
  echo "ERROR: evaluation bigwig not found: ${EVAL_BIGWIG}" >&2
  exit 1
fi

if [[ -z "${OFFICIAL_ROOT}" ]]; then
  echo "ERROR: missing official ChromBPNet root; pass --official-root or set CHROMBPNET_OFFICIAL_ROOT" >&2
  exit 1
fi

OFFICIAL_PREDICT="${OFFICIAL_ROOT}/chrombpnet/training/predict.py"
if [[ ! -f "${OFFICIAL_PREDICT}" ]]; then
  echo "ERROR: official predict.py not found: ${OFFICIAL_PREDICT}" >&2
  exit 1
fi

require_cmd chrombpnet
require_cmd python3

export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

RUN_ROOT="${WORK_ROOT}/runs/${NAME}"
LOG_DIR="${RUN_ROOT}/logs_fast_seed_${SEED}"
mkdir -p "${LOG_DIR}"

IFS=',' read -r -a GPU_ARRAY <<< "${GPUS}"
if [[ ${#GPU_ARRAY[@]} -eq 0 ]]; then
  echo "ERROR: parsed empty GPU list from --gpus='${GPUS}'" >&2
  exit 1
fi

if [[ -z "${MAX_PARALLEL}" ]]; then
  if [[ "${MULTI_GPU_TRAIN}" == "1" ]]; then
    MAX_PARALLEL="1"
  else
    MAX_PARALLEL="${#GPU_ARRAY[@]}"
  fi
fi

if [[ "${MAX_PARALLEL}" -lt 1 ]]; then
  echo "ERROR: --max-parallel must be >=1" >&2
  exit 1
fi

if [[ "${MULTI_GPU_TRAIN}" == "1" && "${MAX_PARALLEL}" -ne 1 ]]; then
  echo "WARN: --multi-gpu-train requires --max-parallel=1, overriding from ${MAX_PARALLEL}" >&2
  MAX_PARALLEL="1"
fi

declare -a FOLD_FILES
if [[ -n "${FOLDS}" ]]; then
  read -r -a FOLD_KEYS <<< "${FOLDS}"
  for fold_key in "${FOLD_KEYS[@]}"; do
    fold_json="${FOLD_DIR}/${fold_key}.json"
    if [[ ! -f "${fold_json}" ]]; then
      echo "ERROR: fold json not found: ${fold_json}" >&2
      exit 1
    fi
    FOLD_FILES+=("${fold_json}")
  done
else
  mapfile -t FOLD_FILES < <(find "${FOLD_DIR}" -maxdepth 1 -type f -name 'fold_*.json' | sort -V)
fi

if [[ ${#FOLD_FILES[@]} -eq 0 ]]; then
  echo "ERROR: no fold files found" >&2
  exit 1
fi

echo "[INFO] run root: ${RUN_ROOT}"
echo "[INFO] mode: fast learning reproduction (single seed, no modisco)"
echo "[INFO] seed: ${SEED}"
echo "[INFO] gpus: ${GPUS}"
echo "[INFO] max parallel folds: ${MAX_PARALLEL}"
echo "[INFO] total folds: ${#FOLD_FILES[@]}"
echo "[INFO] save all checkpoints: ${SAVE_ALL_CHECKPOINTS}"
echo "[INFO] multi gpu train: ${MULTI_GPU_TRAIN}"
echo "[INFO] official root: ${OFFICIAL_ROOT}"
echo "[INFO] external nonpeaks: ${EXTERNAL_NONPEAKS:-<none>}"
echo "[INFO] eval bigwig override: ${EVAL_BIGWIG:-<none>}"
echo "[INFO] use input peaks for eval: ${USE_INPUT_PEAKS_FOR_EVAL}"

run_fold() {
  local fold_json="$1"
  local gpu="$2"
  local fold_key
  fold_key="$(basename "${fold_json}" .json)"
  local train_visible_devices="${gpu}"
  local predict_gpu="${gpu}"
  local train_multi_gpu="0"
  if [[ "${MULTI_GPU_TRAIN}" == "1" ]]; then
    train_visible_devices="${GPUS}"
    predict_gpu="${GPU_ARRAY[0]}"
    train_multi_gpu="1"
  fi

  local fold_dir="${RUN_ROOT}/${fold_key}"
  local bg_dir="${fold_dir}/background"
  mkdir -p "${bg_dir}"

  local nonpeak_prefix="${bg_dir}/nonpeaks"
  local nonpeaks=""

  echo "[INFO][${fold_key}] gpu=${gpu} train_visible_devices=${train_visible_devices} multi_gpu=${train_multi_gpu}: start"

  if [[ -n "${EXTERNAL_NONPEAKS}" ]]; then
    nonpeaks="${EXTERNAL_NONPEAKS}"
    echo "[INFO][${fold_key}] prep nonpeaks: skip (using external nonpeaks)"
  else
    nonpeaks="${nonpeak_prefix}_negatives.bed"
    if [[ ! -f "${nonpeaks}" ]]; then
      if [[ -d "${nonpeak_prefix}_auxiliary" ]]; then
        rm -rf "${nonpeak_prefix}_auxiliary"
      fi
      echo "[INFO][${fold_key}] prep nonpeaks"
      chrombpnet prep nonpeaks \
        -g "${GENOME}" \
        -o "${nonpeak_prefix}" \
        -p "${PEAKS}" \
        -c "${CHROM_SIZES}" \
        -fl "${fold_json}" \
        -il 2114 \
        -st 1000 \
        -br "${BLACKLIST}"
    else
      echo "[INFO][${fold_key}] prep nonpeaks: skip"
    fi
  fi

  local seed_dir="${fold_dir}/seed_${SEED}"
  mkdir -p "${seed_dir}/logs"

  local bias_dir="${seed_dir}/bias"
  local chrom_dir="${seed_dir}/chrombpnet"
  local bias_model="${bias_dir}/models/bias.h5"
  local bias_metrics="${bias_dir}/evaluation/bias_metrics.json"
  local chrom_model="${chrom_dir}/models/chrombpnet.h5"
  local chrom_metrics="${chrom_dir}/evaluation/chrombpnet_metrics.json"

  if [[ ! -f "${bias_model}" ]]; then
    if [[ -d "${bias_dir}" ]]; then
      rm -rf "${bias_dir}"
    fi
    echo "[INFO][${fold_key}] bias train seed=${SEED}"
    CUDA_VISIBLE_DEVICES="${train_visible_devices}" CHROMBPNET_MULTI_GPU="${train_multi_gpu}" chrombpnet bias train \
      -g "${GENOME}" \
      -c "${CHROM_SIZES}" \
      -ibam "${BAM}" \
      -o "${bias_dir}" \
      -d "${DATA_TYPE}" \
      -p "${PEAKS}" \
      -n "${nonpeaks}" \
      -fl "${fold_json}" \
      -b "${BIAS_THRESHOLD_FACTOR}" \
      -s "${SEED}" \
      -bs "${BATCH_SIZE}" \
      -e "${EPOCHS}" \
      -es "${EARLY_STOP}"
  else
    echo "[INFO][${fold_key}] bias train: skip"
  fi

  if [[ ! -f "${bias_metrics}" ]]; then
    local bias_bigwig="${bias_dir}/auxiliary/data_unstranded.bw"
    if [[ ! -f "${bias_bigwig}" ]]; then
      echo "ERROR: bias bigwig missing: ${bias_bigwig}" >&2
      return 1
    fi
    echo "[INFO][${fold_key}] bias predict metrics"
    CUDA_VISIBLE_DEVICES="${predict_gpu}" CHROMBPNET_MULTI_GPU=0 python3 "${OFFICIAL_PREDICT}" \
      -g "${GENOME}" \
      -b "${bias_bigwig}" \
      -p "${PEAKS}" \
      -n "${nonpeaks}" \
      -o "${bias_dir}/evaluation/bias" \
      -fl "${fold_json}" \
      -m "${bias_model}" \
      -bs "${PREDICT_BATCH_SIZE}" \
      -il 2114 \
      -ol 1000 \
      -s "${SEED}"
  else
    echo "[INFO][${fold_key}] bias metrics: skip"
  fi

  if [[ ! -f "${chrom_model}" ]]; then
    if [[ -d "${chrom_dir}" ]]; then
      rm -rf "${chrom_dir}"
    fi
    echo "[INFO][${fold_key}] chrombpnet train seed=${SEED}"
    local chrom_cmd=(
      chrombpnet train
      -g "${GENOME}"
      -c "${CHROM_SIZES}"
      -ibam "${BAM}"
      -o "${chrom_dir}"
      -d "${DATA_TYPE}"
      -p "${PEAKS}"
      -n "${nonpeaks}"
      -fl "${fold_json}"
      -b "${bias_model}"
      -s "${SEED}"
      -bs "${BATCH_SIZE}"
      -e "${EPOCHS}"
      -es "${EARLY_STOP}"
    )
    if [[ "${SAVE_ALL_CHECKPOINTS}" == "1" ]]; then
      chrom_cmd+=(--save-all-checkpoints)
    fi
    CUDA_VISIBLE_DEVICES="${train_visible_devices}" CHROMBPNET_MULTI_GPU="${train_multi_gpu}" "${chrom_cmd[@]}"
  else
    echo "[INFO][${fold_key}] chrombpnet train: skip"
  fi

  if [[ ! -f "${chrom_metrics}" ]]; then
    local chrom_bigwig="${chrom_dir}/auxiliary/data_unstranded.bw"
    if [[ ! -f "${chrom_bigwig}" ]]; then
      echo "ERROR: chrombpnet bigwig missing: ${chrom_bigwig}" >&2
      return 1
    fi
    local eval_bigwig="${chrom_bigwig}"
    if [[ -n "${EVAL_BIGWIG}" ]]; then
      eval_bigwig="${EVAL_BIGWIG}"
    fi
    local eval_peaks="${chrom_dir}/auxiliary/filtered.peaks.bed"
    local eval_nonpeaks="None"
    if [[ -n "${EXTERNAL_NONPEAKS}" ]]; then
      eval_nonpeaks="${nonpeaks}"
    fi
    if [[ "${USE_INPUT_PEAKS_FOR_EVAL}" == "1" ]]; then
      eval_peaks="${PEAKS}"
    elif [[ ! -f "${eval_peaks}" ]]; then
      echo "ERROR: filtered peaks missing: ${eval_peaks}" >&2
      return 1
    fi
    echo "[INFO][${fold_key}] chrombpnet predict metrics"
    CUDA_VISIBLE_DEVICES="${predict_gpu}" CHROMBPNET_MULTI_GPU=0 python3 "${OFFICIAL_PREDICT}" \
      -g "${GENOME}" \
      -b "${eval_bigwig}" \
      -p "${eval_peaks}" \
      -n "${eval_nonpeaks}" \
      -o "${chrom_dir}/evaluation/chrombpnet" \
      -fl "${fold_json}" \
      -m "${chrom_model}" \
      -bs "${PREDICT_BATCH_SIZE}" \
      -il 2114 \
      -ol 1000 \
      -s "${SEED}"
  else
    echo "[INFO][${fold_key}] chrombpnet metrics: skip"
  fi

  echo "[INFO][${fold_key}] gpu=${gpu}: done"
}

for idx in "${!FOLD_FILES[@]}"; do
  fold_json="${FOLD_FILES[$idx]}"
  fold_key="$(basename "${fold_json}" .json)"
  gpu="${GPU_ARRAY[$((idx % ${#GPU_ARRAY[@]}))]}"
  fold_log="${LOG_DIR}/${fold_key}_seed_${SEED}_gpu${gpu}.log"

  while [[ "$(jobs -rp | wc -l)" -ge "${MAX_PARALLEL}" ]]; do
    wait -n
  done

  echo "[INFO] launch ${fold_key} on gpu ${gpu} -> ${fold_log}"
  run_fold "${fold_json}" "${gpu}" > "${fold_log}" 2>&1 &
done

wait

echo "[INFO] summarize metrics"
python3 "${SCRIPT_DIR}/summarize_metrics.py" \
  --run-root "${RUN_ROOT}" \
  --output-dir "${RUN_ROOT}/summary_fast_seed_${SEED}"

echo "[INFO] finished fast run: ${RUN_ROOT}"
