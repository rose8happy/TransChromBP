#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

select_best_epoch="${REPO_ROOT}/scripts/paper_aligned_repro/select_best_epoch.py"
run_fast_1seed="${REPO_ROOT}/scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh"

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

if ! rg -n "CHROMBPNET_OFFICIAL_ROOT" "${run_fast_1seed}" >/dev/null; then
  echo "ERROR: run_paper_aligned_fast_1seed.sh does not mention CHROMBPNET_OFFICIAL_ROOT" >&2
  exit 1
fi

echo "official externalization guard passed"
