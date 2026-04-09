#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

select_best_epoch="${REPO_ROOT}/scripts/paper_aligned_repro/select_best_epoch.py"
run_fast_1seed="${REPO_ROOT}/scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh"
run_remote_dataset_prep="${REPO_ROOT}/scripts/run_remote_chrombpnet_dataset_prep.sh"
start_6000_dataset_prep="${REPO_ROOT}/scripts/start_6000_chrombpnet_dataset_prep.sh"
start_6002_dataset_prep="${REPO_ROOT}/scripts/start_6002_chrombpnet_dataset_prep.sh"
tutorial_step3="${REPO_ROOT}/workflows/tutorial/step3_get_background_regions.sh"

predict_entrypoint_pattern="${REPO_ROOT}/chrombpnet/training/predict.py|from chrombpnet.training import predict"
if rg -n "${predict_entrypoint_pattern}" "${select_best_epoch}" "${run_fast_1seed}"; then
  echo "ERROR: target files still reference the local chrombpnet predict entrypoint" >&2
  exit 1
fi

payload_package_suffix="egg-info"

for path in \
  "${REPO_ROOT}/chrombpnet" \
  "${REPO_ROOT}/chrombpnet.${payload_package_suffix}" \
  "${REPO_ROOT}/setup.py" \
  "${REPO_ROOT}/MANIFEST.in"
do
  if [[ -e "${path}" ]]; then
    echo "ERROR: expected payload/package path to be absent: ${path}" >&2
    exit 1
  fi
done

if ! python3 "${select_best_epoch}" --help 2>&1 | grep -q -- '--official-root'; then
  echo "ERROR: select_best_epoch.py --help does not expose --official-root" >&2
  exit 1
fi

if ! bash "${run_fast_1seed}" --help 2>&1 | grep -q -- '--official-root'; then
  echo "ERROR: run_paper_aligned_fast_1seed.sh --help does not expose --official-root" >&2
  exit 1
fi

if ! bash "${run_remote_dataset_prep}" --help 2>&1 | grep -q -- '--official-root'; then
  echo "ERROR: run_remote_chrombpnet_dataset_prep.sh --help does not expose --official-root" >&2
  exit 1
fi

if ! bash "${run_remote_dataset_prep}" --help 2>&1 | grep -q -- '--gc-helper-dir'; then
  echo "ERROR: run_remote_chrombpnet_dataset_prep.sh --help does not expose --gc-helper-dir" >&2
  exit 1
fi

if ! rg -n "CHROMBPNET_OFFICIAL_ROOT" "${run_fast_1seed}" >/dev/null; then
  echo "ERROR: run_paper_aligned_fast_1seed.sh does not mention CHROMBPNET_OFFICIAL_ROOT" >&2
  exit 1
fi

if rg -n 'REPO_ROOT/chrombpnet/helpers|../../chrombpnet/helpers' \
  "${start_6000_dataset_prep}" \
  "${start_6002_dataset_prep}" \
  "${tutorial_step3}"
then
  echo "ERROR: gc-helper scripts still reference retired local chrombpnet helpers" >&2
  exit 1
fi

if ! rg -n "CHROMBPNET_OFFICIAL_ROOT|MAKE_GC_MATCHED_NEGATIVES_SCRIPT" "${tutorial_step3}" >/dev/null; then
  echo "ERROR: tutorial step3 does not mention CHROMBPNET_OFFICIAL_ROOT or MAKE_GC_MATCHED_NEGATIVES_SCRIPT" >&2
  exit 1
fi

echo "official externalization guard passed"
