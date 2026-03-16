#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

WORK_ROOT="${1:-/data1/zhoujiazhen/bylw_atac/chrombpnet_paper_repro}"
SEEDS="${2:-1234 2345 3456}"
GPU="${3:-0,1}"

FOLD_DIR="${WORK_ROOT}/folds"
mkdir -p "${WORK_ROOT}" "${WORK_ROOT}/logs"

if [[ ! -f "${FOLD_DIR}/fold_0.json" ]]; then
  python3 "${SCRIPT_DIR}/generate_fold_set.py" \
    --chrom-sizes /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes \
    --output-dir "${FOLD_DIR}"
fi

nohup bash "${SCRIPT_DIR}/run_paper_aligned_matrix.sh" \
  --name K562_ATAC \
  --genome /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.fa \
  --chrom-sizes /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes \
  --bam /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged.bam \
  --peaks /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/overlap.bed.gz \
  --blacklist /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/blacklist.bed.gz \
  --fold-dir "${FOLD_DIR}" \
  --work-root "${WORK_ROOT}" \
  --seeds "${SEEDS}" \
  --gpu "${GPU}" \
  > "${WORK_ROOT}/logs/k562_matrix_$(date +%Y%m%d_%H%M%S).log" 2>&1 &

echo "Launched K562 paper-aligned matrix job."
echo "Work root: ${WORK_ROOT}"
echo "Monitor: ls -t ${WORK_ROOT}/logs/k562_matrix_*.log | head -n 1 | xargs tail -f"

