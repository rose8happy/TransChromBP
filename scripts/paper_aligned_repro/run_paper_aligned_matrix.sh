#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

function usage() {
  cat <<'USAGE'
Usage:
  bash scripts/paper_aligned_repro/run_paper_aligned_matrix.sh \
    --name K562_ATAC \
    --genome /path/to/hg38.fa \
    --chrom-sizes /path/to/hg38.chrom.sizes \
    --bam /path/to/merged.bam \
    --peaks /path/to/overlap.bed.gz \
    --blacklist /path/to/blacklist.bed.gz \
    --fold-dir /path/to/folds \
    --work-root /path/to/chrombpnet_paper_repro \
    --seeds "1234 2345 3456" \
    --gpu 0,1

Required:
  --name
  --genome
  --chrom-sizes
  --bam
  --peaks
  --blacklist
  --fold-dir
  --work-root
  --seeds

Optional:
  --data-type ATAC|DNASE    (default: ATAC)
  --gpu GPU_IDS             (default: 0)
  --bias-threshold-factor X (default: 0.5)
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
SEEDS=""
DATA_TYPE="ATAC"
GPU="0"
BIAS_THRESHOLD_FACTOR="0.5"

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
    --seeds) SEEDS="$2"; shift 2 ;;
    --data-type) DATA_TYPE="$2"; shift 2 ;;
    --gpu) GPU="$2"; shift 2 ;;
    --bias-threshold-factor) BIAS_THRESHOLD_FACTOR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "${NAME}" || -z "${GENOME}" || -z "${CHROM_SIZES}" || -z "${BAM}" || -z "${PEAKS}" || -z "${BLACKLIST}" || -z "${FOLD_DIR}" || -z "${WORK_ROOT}" || -z "${SEEDS}" ]]; then
  echo "ERROR: missing required arguments" >&2
  usage
  exit 1
fi

require_cmd chrombpnet
require_cmd python3

RUN_ROOT="${WORK_ROOT}/runs/${NAME}"
LOG_DIR="${RUN_ROOT}/logs"
mkdir -p "${LOG_DIR}"

echo "[INFO] run root: ${RUN_ROOT}"
echo "[INFO] data type: ${DATA_TYPE}"
echo "[INFO] seeds: ${SEEDS}"
echo "[INFO] gpu: ${GPU}"

read -r -a SEED_ARRAY <<< "${SEEDS}"
if [[ ${#SEED_ARRAY[@]} -eq 0 ]]; then
  echo "ERROR: seeds parsed empty from '${SEEDS}'" >&2
  exit 1
fi

mapfile -t FOLD_FILES < <(find "${FOLD_DIR}" -maxdepth 1 -type f -name 'fold_*.json' | sort -V)
if [[ ${#FOLD_FILES[@]} -eq 0 ]]; then
  echo "ERROR: no fold_*.json found in ${FOLD_DIR}" >&2
  exit 1
fi

for fold_json in "${FOLD_FILES[@]}"; do
  fold_key="$(basename "${fold_json}" .json)"
  fold_dir="${RUN_ROOT}/${fold_key}"
  bg_dir="${fold_dir}/background"
  mkdir -p "${bg_dir}"

  nonpeak_prefix="${bg_dir}/nonpeaks"
  nonpeaks="${nonpeak_prefix}_negatives.bed"

  if [[ ! -f "${nonpeaks}" ]]; then
    if [[ -d "${nonpeak_prefix}_auxiliary" ]]; then
      rm -rf "${nonpeak_prefix}_auxiliary"
    fi
    echo "[INFO] prep nonpeaks for ${fold_key}"
    chrombpnet prep nonpeaks \
      -g "${GENOME}" \
      -o "${nonpeak_prefix}" \
      -p "${PEAKS}" \
      -c "${CHROM_SIZES}" \
      -fl "${fold_json}" \
      -il 2114 \
      -st 1000 \
      -br "${BLACKLIST}" \
      > "${fold_dir}/step3_prep_nonpeaks.log" 2>&1
  else
    echo "[INFO] prep nonpeaks: skip ${fold_key} (exists ${nonpeaks})"
  fi

  for seed in "${SEED_ARRAY[@]}"; do
    seed_dir="${fold_dir}/seed_${seed}"
    mkdir -p "${seed_dir}/logs"

    bias_dir="${seed_dir}/bias"
    chrom_dir="${seed_dir}/chrombpnet"
    bias_metrics="${bias_dir}/evaluation/bias_metrics.json"
    chrom_metrics="${chrom_dir}/evaluation/chrombpnet_metrics.json"

    if [[ ! -f "${bias_metrics}" ]]; then
      if [[ -d "${bias_dir}" ]]; then
        rm -rf "${bias_dir}"
      fi
      echo "[INFO] bias pipeline ${fold_key} seed=${seed}"
      CUDA_VISIBLE_DEVICES="${GPU}" chrombpnet bias pipeline \
        -g "${GENOME}" \
        -c "${CHROM_SIZES}" \
        -ibam "${BAM}" \
        -o "${bias_dir}" \
        -d "${DATA_TYPE}" \
        -p "${PEAKS}" \
        -n "${nonpeaks}" \
        -fl "${fold_json}" \
        -b "${BIAS_THRESHOLD_FACTOR}" \
        -s "${seed}" \
        > "${seed_dir}/logs/step4_bias_pipeline.log" 2>&1
    else
      echo "[INFO] bias pipeline: skip ${fold_key} seed=${seed} (exists ${bias_metrics})"
    fi

    if [[ ! -f "${chrom_metrics}" ]]; then
      if [[ -d "${chrom_dir}" ]]; then
        rm -rf "${chrom_dir}"
      fi
      bias_model="${bias_dir}/models/bias.h5"
      if [[ ! -f "${bias_model}" ]]; then
        echo "ERROR: expected bias model missing: ${bias_model}" >&2
        exit 1
      fi

      echo "[INFO] chrombpnet pipeline ${fold_key} seed=${seed}"
      CUDA_VISIBLE_DEVICES="${GPU}" chrombpnet pipeline \
        -g "${GENOME}" \
        -c "${CHROM_SIZES}" \
        -ibam "${BAM}" \
        -o "${chrom_dir}" \
        -d "${DATA_TYPE}" \
        -p "${PEAKS}" \
        -n "${nonpeaks}" \
        -fl "${fold_json}" \
        -b "${bias_model}" \
        -s "${seed}" \
        > "${seed_dir}/logs/step6_chrombpnet_pipeline.log" 2>&1
    else
      echo "[INFO] chrombpnet pipeline: skip ${fold_key} seed=${seed} (exists ${chrom_metrics})"
    fi
  done
done

echo "[INFO] summarize metrics"
python3 "${SCRIPT_DIR}/summarize_metrics.py" \
  --run-root "${RUN_ROOT}" \
  --output-dir "${RUN_ROOT}/summary"

echo "[INFO] finished: ${RUN_ROOT}"
