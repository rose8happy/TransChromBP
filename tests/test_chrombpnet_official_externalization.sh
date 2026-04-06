#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

select_best_epoch="${REPO_ROOT}/scripts/paper_aligned_repro/select_best_epoch.py"
run_fast_1seed="${REPO_ROOT}/scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh"

if rg -n "REPO_ROOT/chrombpnet/training/predict.py|from chrombpnet.training import predict" "${select_best_epoch}" "${run_fast_1seed}"; then
  echo "ERROR: target files still reference the local chrombpnet predict entrypoint" >&2
  exit 1
fi

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
