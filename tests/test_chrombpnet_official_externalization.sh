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

if ! bash "${REPO_ROOT}/scripts/paper_aligned_repro/run_tutorial_strict_compare_official.sh" --help 2>&1 | grep -q -- '--official-root'; then
  echo "ERROR: run_tutorial_strict_compare_official.sh --help does not expose --official-root" >&2
  exit 1
fi

if ! rg -n "CHROMBPNET_OFFICIAL_ROOT" "${run_fast_1seed}" >/dev/null; then
  echo "ERROR: run_paper_aligned_fast_1seed.sh does not mention CHROMBPNET_OFFICIAL_ROOT" >&2
  exit 1
fi

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

if ! rg -n "CHROMBPNET_OFFICIAL_ROOT|--official-root" "${REPO_ROOT}/scripts/paper_aligned_repro/run_tutorial_strict_compare_official.sh" >/dev/null; then
  echo "ERROR: run_tutorial_strict_compare_official.sh does not clearly handle official root" >&2
  exit 1
fi

missing_root="${tmpdir}/definitely_missing_official_root"
mkdir -p "${tmpdir}/missing_root_check"
if bash "${run_fast_1seed}" \
  --name missing_root_check \
  --genome "${tmpdir}/missing_root_check/genome.fa" \
  --chrom-sizes "${tmpdir}/missing_root_check/chrom.sizes" \
  --bam "${tmpdir}/missing_root_check/merged.bam" \
  --peaks "${tmpdir}/missing_root_check/peaks.bed" \
  --blacklist "${tmpdir}/missing_root_check/blacklist.bed" \
  --fold-dir "${tmpdir}/missing_root_check/folds" \
  --work-root "${tmpdir}/missing_root_check/work" \
  --official-root "${missing_root}" \
  2> "${tmpdir}/missing_root_check/stderr"; then
  echo "ERROR: run_paper_aligned_fast_1seed.sh unexpectedly succeeded with a missing official root" >&2
  exit 1
fi
if ! grep -q "official ChromBPNet root is not a directory" "${tmpdir}/missing_root_check/stderr"; then
  echo "ERROR: run_paper_aligned_fast_1seed.sh did not emit the expected directory validation error" >&2
  exit 1
fi
if grep -q "cd: no such file or directory" "${tmpdir}/missing_root_check/stderr"; then
  echo "ERROR: run_paper_aligned_fast_1seed.sh still emitted a raw cd error for missing official root" >&2
  exit 1
fi

write_fake_official_root() {
  local root="$1"
  mkdir -p "${root}/chrombpnet/training" "${root}/bin"

  cat > "${root}/chrombpnet/__init__.py" <<'PY'
from pathlib import Path

PACKAGE_FILE = str(Path(__file__).resolve())
PY

  cat > "${root}/chrombpnet/training/__init__.py" <<'PY'
PY

  cat > "${root}/chrombpnet/training/predict.py" <<'PY'
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
    split_metric = {
        "train": 0.789,
        "valid": 0.123,
        "test": 0.456,
    }.get(args.split or "", 0.654)
    bigwig_name = Path(args.b).name
    if args.split:
        median_jsd = split_metric
    else:
        median_jsd = 0.111 if "override" not in bigwig_name else 0.222

    payload = {
        "counts_metrics": {
            "peaks": {"pearsonr": 0.91, "spearmanr": 0.81, "mse": 0.19},
            "nonpeaks": {"pearsonr": 0.11},
        },
        "profile_metrics": {
            "peaks": {"mean_jsd": 0.321, "median_jsd": median_jsd, "median_norm_jsd": 0.222},
            "nonpeaks": {"median_jsd": 0.456},
        },
        "classification_metrics": {
            "peaks_vs_nonpeaks": {
                "logcounts": {"auroc": 0.73, "auprc": 0.64, "best_f1": 0.51}
            }
        },
        "imported_chrombpnet_file": chrombpnet.PACKAGE_FILE,
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
        "source_split": args.split or "",
        "source_bigwig": bigwig_name,
    }
    output_prefix.with_name(f"{output_prefix.name}_metrics.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
PY

  cat > "${root}/bin/chrombpnet" <<'PY'
#!/usr/bin/env python3

from __future__ import annotations

import os
import sys
from pathlib import Path


def get_arg(name: str, default: str = "") -> str:
    args = sys.argv[1:]
    for idx, value in enumerate(args):
        if value == name and idx + 1 < len(args):
            return args[idx + 1]
    return default


def main() -> None:
    if len(sys.argv) < 3:
        return
    subcommand = (sys.argv[1], sys.argv[2])
    if subcommand == ("prep", "nonpeaks"):
        prefix = Path(get_arg("-o"))
        prefix.parent.mkdir(parents=True, exist_ok=True)
        (Path(f"{prefix}_negatives.bed")).write_text("chr1\t0\t10\n", encoding="utf-8")
        aux = Path(f"{prefix}_auxiliary")
        aux.mkdir(parents=True, exist_ok=True)
        return
    if subcommand == ("bias", "train"):
        out_dir = Path(get_arg("-o"))
        (out_dir / "models").mkdir(parents=True, exist_ok=True)
        (out_dir / "auxiliary").mkdir(parents=True, exist_ok=True)
        (out_dir / "models" / "bias.h5").write_text("bias model\n", encoding="utf-8")
        (out_dir / "auxiliary" / "data_unstranded.bw").write_text("bw\n", encoding="utf-8")
        return
    if subcommand == ("train", "-g"):
        out_dir = Path(get_arg("-o"))
        (out_dir / "models").mkdir(parents=True, exist_ok=True)
        (out_dir / "auxiliary").mkdir(parents=True, exist_ok=True)
        (out_dir / "models" / "chrombpnet.h5").write_text("chrom model\n", encoding="utf-8")
        (out_dir / "auxiliary" / "data_unstranded.bw").write_text("bw\n", encoding="utf-8")
        return


if __name__ == "__main__":
    main()
PY

  chmod +x "${root}/chrombpnet/training/predict.py" "${root}/bin/chrombpnet"
}

fake_official_root_rel="${tmpdir}/official_rel"
write_fake_official_root "${fake_official_root_rel}"
mkdir -p "${tmpdir}/relwork/data" "${tmpdir}/relwork/models" "${tmpdir}/relwork/out"
cat > "${tmpdir}/relwork/data/genome.fa" <<'EOF'
>chr1
ACGT
EOF
cat > "${tmpdir}/relwork/data/bigwig.bw" <<'EOF'
placeholder
EOF
cat > "${tmpdir}/relwork/data/fold.json" <<'EOF'
{"fold": 0}
EOF
cat > "${tmpdir}/relwork/data/peaks.bed" <<'EOF'
chr1	0	10
EOF
cat > "${tmpdir}/relwork/data/nonpeaks.bed" <<'EOF'
chr1	10	20
EOF
: > "${tmpdir}/relwork/models/chrombpnet.epoch_1.h5"
cat > "${tmpdir}/relwork/out/chrombpnet.epoch_1_metrics.json" <<'EOF'
{
  "_official_predict_provenance": {
    "official_root": "/stale",
    "predict_py": "/stale"
  },
  "imported_chrombpnet_file": "stale",
  "profile_metrics": {
    "peaks": {
      "median_jsd": 9.999
    }
  }
}
EOF

(
  cd "${tmpdir}/relwork"
  PYTHONPATH="${REPO_ROOT}${PYTHONPATH:+:${PYTHONPATH}}" \
  CUDA_VISIBLE_DEVICES=9 \
  python3 "${select_best_epoch}" \
    --model-glob 'models/chrombpnet.epoch_*.h5' \
    --genome 'data/genome.fa' \
    --bigwig 'data/bigwig.bw' \
    --peaks 'data/peaks.bed' \
    --nonpeaks 'data/nonpeaks.bed' \
    --fold-json 'data/fold.json' \
    --output-dir 'out' \
    --official-root "${fake_official_root_rel}"
)

python3 - "${tmpdir}/relwork/out/chrombpnet.epoch_1_metrics.json" "${fake_official_root_rel}" "valid" 0.123 <<'PY'
import json
import pathlib
import sys

payload = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
imported = payload["imported_chrombpnet_file"]
cuda_visible_devices = payload["cuda_visible_devices"]
provenance = payload["_official_predict_provenance"]
if not imported.startswith(sys.argv[2]):
    raise SystemExit(f"import came from wrong package: {imported}")
if cuda_visible_devices != "9":
    raise SystemExit(f"CUDA_VISIBLE_DEVICES was not preserved: {cuda_visible_devices!r}")
if provenance["official_root"] != sys.argv[2]:
    raise SystemExit(f"selector provenance did not refresh: {provenance}")
if provenance["split"] != sys.argv[3]:
    raise SystemExit(f"selector provenance split mismatch: {provenance}")
if payload["profile_metrics"]["peaks"]["median_jsd"] != float(sys.argv[4]):
    raise SystemExit(f"selector metrics did not refresh: {payload['profile_metrics']['peaks']['median_jsd']}")
PY

(
  cd "${tmpdir}/relwork"
  PYTHONPATH="${REPO_ROOT}${PYTHONPATH:+:${PYTHONPATH}}" \
  CUDA_VISIBLE_DEVICES=9 \
  python3 "${select_best_epoch}" \
    --model-glob 'models/chrombpnet.epoch_*.h5' \
    --genome 'data/genome.fa' \
    --bigwig 'data/bigwig.bw' \
    --peaks 'data/peaks.bed' \
    --nonpeaks 'data/nonpeaks.bed' \
    --fold-json 'data/fold.json' \
    --output-dir 'out' \
    --split test \
    --official-root "${fake_official_root_rel}"
)

python3 - "${tmpdir}/relwork/out/chrombpnet.epoch_1_metrics.json" "${fake_official_root_rel}" "test" 0.456 <<'PY'
import json
import pathlib
import sys

payload = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
provenance = payload["_official_predict_provenance"]
if provenance["split"] != sys.argv[3]:
    raise SystemExit(f"selector did not refresh provenance after split change: {provenance}")
if payload["profile_metrics"]["peaks"]["median_jsd"] != float(sys.argv[4]):
    raise SystemExit(f"selector cache was not refreshed for changed split: {payload['profile_metrics']['peaks']['median_jsd']}")
PY

fake_official_root_run="${tmpdir}/official_run"
write_fake_official_root "${fake_official_root_run}"
mkdir -p "${tmpdir}/run/data" "${tmpdir}/run/folds" "${tmpdir}/run/work"
cat > "${tmpdir}/run/data/genome.fa" <<'EOF'
>chr1
ACGT
EOF
cat > "${tmpdir}/run/data/chrom.sizes" <<'EOF'
chr1	4
EOF
cat > "${tmpdir}/run/data/merged.bam" <<'EOF'
placeholder
EOF
cat > "${tmpdir}/run/data/peaks.bed" <<'EOF'
chr1	0	10
EOF
cat > "${tmpdir}/run/data/blacklist.bed" <<'EOF'
chr1	20	30
EOF
cat > "${tmpdir}/run/data/nonpeaks.bed" <<'EOF'
chr1	40	50
EOF
cat > "${tmpdir}/run/data/override_unstranded.bw" <<'EOF'
override
EOF
cat > "${tmpdir}/run/folds/fold_0.json" <<'EOF'
{"fold": 0}
EOF
mkdir -p "${tmpdir}/run/work/runs/test_official_externalization/fold_0/seed_7/bias/evaluation" \
  "${tmpdir}/run/work/runs/test_official_externalization/fold_0/seed_7/chrombpnet/evaluation"
cat > "${tmpdir}/run/work/runs/test_official_externalization/fold_0/seed_7/bias/evaluation/bias_metrics.json" <<'EOF'
{
  "_official_predict_provenance": {
    "official_root": "/stale",
    "predict_py": "/stale"
  },
  "imported_chrombpnet_file": "stale",
  "profile_metrics": {
    "peaks": {
      "median_jsd": 9.999
    }
  }
}
EOF
cat > "${tmpdir}/run/work/runs/test_official_externalization/fold_0/seed_7/chrombpnet/evaluation/chrombpnet_metrics.json" <<'EOF'
{
  "_official_predict_provenance": {
    "official_root": "/stale",
    "predict_py": "/stale"
  },
  "imported_chrombpnet_file": "stale",
  "profile_metrics": {
    "peaks": {
      "median_jsd": 9.999
    }
  }
}
EOF

PATH="${fake_official_root_run}/bin:${PATH}" \
CHROMBPNET_OFFICIAL_ROOT="${fake_official_root_run}" \
bash "${run_fast_1seed}" \
  --name test_official_externalization \
  --genome "${tmpdir}/run/data/genome.fa" \
  --chrom-sizes "${tmpdir}/run/data/chrom.sizes" \
  --bam "${tmpdir}/run/data/merged.bam" \
  --peaks "${tmpdir}/run/data/peaks.bed" \
  --blacklist "${tmpdir}/run/data/blacklist.bed" \
  --fold-dir "${tmpdir}/run/folds" \
  --work-root "${tmpdir}/run/work" \
  --seed 7 \
  --gpus 0 \
  --max-parallel 1 \
  --official-root "${fake_official_root_run}" \
  --use-input-peaks-for-eval

python3 - "${tmpdir}/run/work/runs/test_official_externalization/fold_0/seed_7/chrombpnet/evaluation/chrombpnet_metrics.json" "${fake_official_root_run}" "0.111" "data_unstranded.bw" <<'PY'
import json
import pathlib
import sys

payload = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
imported = payload["imported_chrombpnet_file"]
if not imported.startswith(sys.argv[2]):
    raise SystemExit(f"run script imported chrombpnet from wrong package: {imported}")
if payload["_official_predict_provenance"]["official_root"] != sys.argv[2]:
    raise SystemExit(f"run script provenance did not refresh: {payload['_official_predict_provenance']}")
if payload["profile_metrics"]["peaks"]["median_jsd"] != float(sys.argv[3]):
    raise SystemExit(f"run script metrics did not refresh: {payload['profile_metrics']['peaks']['median_jsd']}")
if payload["source_bigwig"] != sys.argv[4]:
    raise SystemExit(f"run script source bigwig mismatch: {payload['source_bigwig']}")
PY

python3 - "${tmpdir}/run/work/runs/test_official_externalization/fold_0/seed_7/bias/evaluation/bias_metrics.json" "${fake_official_root_run}" <<'PY'
import json
import pathlib
import sys

payload = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
if payload["imported_chrombpnet_file"] == "stale":
    raise SystemExit("bias metrics cache was not refreshed")
if payload["_official_predict_provenance"]["official_root"] != sys.argv[2]:
    raise SystemExit(f"bias metrics provenance did not refresh: {payload['_official_predict_provenance']}")
if payload["profile_metrics"]["peaks"]["median_jsd"] != 0.111:
    raise SystemExit(f"bias metrics did not refresh: {payload['profile_metrics']['peaks']['median_jsd']}")
PY

PATH="${fake_official_root_run}/bin:${PATH}" \
CHROMBPNET_OFFICIAL_ROOT="${fake_official_root_run}" \
bash "${run_fast_1seed}" \
  --name test_official_externalization \
  --genome "${tmpdir}/run/data/genome.fa" \
  --chrom-sizes "${tmpdir}/run/data/chrom.sizes" \
  --bam "${tmpdir}/run/data/merged.bam" \
  --peaks "${tmpdir}/run/data/peaks.bed" \
  --blacklist "${tmpdir}/run/data/blacklist.bed" \
  --fold-dir "${tmpdir}/run/folds" \
  --work-root "${tmpdir}/run/work" \
  --seed 7 \
  --gpus 0 \
  --max-parallel 1 \
  --official-root "${fake_official_root_run}" \
  --eval-bigwig "${tmpdir}/run/data/override_unstranded.bw" \
  --use-input-peaks-for-eval

python3 - "${tmpdir}/run/work/runs/test_official_externalization/fold_0/seed_7/chrombpnet/evaluation/chrombpnet_metrics.json" "${fake_official_root_run}" "0.222" "${tmpdir}/run/data/override_unstranded.bw" "override_unstranded.bw" <<'PY'
import json
import pathlib
import sys

payload = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
if not payload["imported_chrombpnet_file"].startswith(sys.argv[2]):
    raise SystemExit(f"run script imported chrombpnet from wrong package after eval_bigwig change: {payload['imported_chrombpnet_file']}")
if payload["_official_predict_provenance"]["effective_bigwig"] != str(pathlib.Path(sys.argv[4]).resolve()):
    raise SystemExit(f"run script provenance did not capture changed eval_bigwig: {payload['_official_predict_provenance']}")
if payload["profile_metrics"]["peaks"]["median_jsd"] != float(sys.argv[3]):
    raise SystemExit(f"run script cache was not refreshed for changed eval_bigwig: {payload['profile_metrics']['peaks']['median_jsd']}")
if payload["source_bigwig"] != sys.argv[5]:
    raise SystemExit(f"run script source bigwig did not change: {payload['source_bigwig']}")
PY

wrapper_tmp="${tmpdir}/wrapper"
mkdir -p "${wrapper_tmp}/work"
if bash "${REPO_ROOT}/scripts/paper_aligned_repro/run_tutorial_strict_compare_official.sh" \
  --mode controlled \
  --work-root "${wrapper_tmp}/work" \
  --seed 7 \
  --folds fold_0 \
  --gpus 0 2> "${wrapper_tmp}/stderr"; then
  echo "ERROR: wrapper preflight unexpectedly succeeded without --official-root" >&2
  exit 1
fi
if ! grep -q "pass --official-root or set CHROMBPNET_OFFICIAL_ROOT" "${wrapper_tmp}/stderr"; then
  echo "ERROR: wrapper preflight did not emit the expected official-root error" >&2
  exit 1
fi

echo "official externalization guard passed"
