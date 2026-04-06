#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

select_best_epoch="${REPO_ROOT}/scripts/paper_aligned_repro/select_best_epoch.py"
run_fast_1seed="${REPO_ROOT}/scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh"

local_predict_patterns=(
  'REPO_ROOT/chrombpnet/training/predict.py'
  '${REPO_ROOT}/chrombpnet/training/predict.py'
  'REPO_ROOT / "chrombpnet" / "training" / "predict.py"'
  'str(REPO_ROOT / "chrombpnet" / "training" / "predict.py")'
  'from chrombpnet.training import predict'
)

for pattern in "${local_predict_patterns[@]}"; do
  if rg -F -n "${pattern}" "${select_best_epoch}" "${run_fast_1seed}"; then
    echo "ERROR: target files still reference the local chrombpnet predict entrypoint: ${pattern}" >&2
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

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

fake_official_root="${tmpdir}/official"
mkdir -p "${fake_official_root}/chrombpnet/training" "${tmpdir}/data" "${tmpdir}/models" "${tmpdir}/out"

cat > "${fake_official_root}/chrombpnet/__init__.py" <<'PY'
from pathlib import Path

PACKAGE_FILE = str(Path(__file__).resolve())
PY

cat > "${fake_official_root}/chrombpnet/training/__init__.py" <<'PY'
PY

cat > "${fake_official_root}/chrombpnet/training/predict.py" <<'PY'
#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import chrombpnet


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-g")
    parser.add_argument("-b")
    parser.add_argument("-p")
    parser.add_argument("-n")
    parser.add_argument("-o", required=True)
    parser.add_argument("-fl")
    parser.add_argument("-m")
    parser.add_argument("-bs")
    parser.add_argument("-s")
    parser.add_argument("-il")
    parser.add_argument("-ol")
    parser.add_argument("--split")
    parser.add_argument("--metrics-only", action="store_true")
    args = parser.parse_args()

    output_prefix = Path(args.o)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "profile_metrics": {"peaks": {"median_jsd": 0.123}},
        "imported_chrombpnet_file": chrombpnet.PACKAGE_FILE,
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
    }
    output_prefix.with_name(f"{output_prefix.name}_metrics.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
PY

chmod +x "${fake_official_root}/chrombpnet/training/predict.py"

cat > "${tmpdir}/data/genome.fa" <<'EOF'
>chr1
ACGT
EOF
cat > "${tmpdir}/data/bigwig.bw" <<'EOF'
placeholder
EOF
cat > "${tmpdir}/data/fold.json" <<'EOF'
{"fold": 0}
EOF
cat > "${tmpdir}/data/peaks.bed" <<'EOF'
chr1	0	10
EOF
cat > "${tmpdir}/data/nonpeaks.bed" <<'EOF'
chr1	10	20
EOF
: > "${tmpdir}/models/chrombpnet.epoch_1.h5"

PYTHONPATH="${REPO_ROOT}${PYTHONPATH:+:${PYTHONPATH}}" \
CUDA_VISIBLE_DEVICES=9 \
python3 "${select_best_epoch}" \
  --model-glob "${tmpdir}/models/chrombpnet.epoch_*.h5" \
  --genome "${tmpdir}/data/genome.fa" \
  --bigwig "${tmpdir}/data/bigwig.bw" \
  --peaks "${tmpdir}/data/peaks.bed" \
  --nonpeaks "${tmpdir}/data/nonpeaks.bed" \
  --fold-json "${tmpdir}/data/fold.json" \
  --output-dir "${tmpdir}/out" \
  --official-root "${fake_official_root}"

python3 - "${tmpdir}/out/chrombpnet.epoch_1_metrics.json" "${fake_official_root}" <<'PY'
import json
import pathlib
import sys

payload = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
imported = payload["imported_chrombpnet_file"]
cuda_visible_devices = payload["cuda_visible_devices"]
if not imported.startswith(sys.argv[2]):
    raise SystemExit(f"import came from wrong package: {imported}")
if cuda_visible_devices != "9":
    raise SystemExit(f"CUDA_VISIBLE_DEVICES was not preserved: {cuda_visible_devices!r}")
PY

echo "official externalization guard passed"
