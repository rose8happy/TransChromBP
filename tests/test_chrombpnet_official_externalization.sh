#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

select_best_epoch="${REPO_ROOT}/scripts/paper_aligned_repro/select_best_epoch.py"
run_fast_1seed="${REPO_ROOT}/scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh"
run_remote_dataset_prep="${REPO_ROOT}/scripts/run_remote_chrombpnet_dataset_prep.sh"
start_dataset_prep_6000="${REPO_ROOT}/scripts/start_6000_chrombpnet_dataset_prep.sh"
start_dataset_prep_6002="${REPO_ROOT}/scripts/start_6002_chrombpnet_dataset_prep.sh"
tutorial_step3="${REPO_ROOT}/workflows/tutorial/step3_get_background_regions.sh"
full_workflow_test="${REPO_ROOT}/tests/full_workflow.sh"

local_predict_patterns=(
  'REPO_ROOT''/chrombpnet/training/predict.py'
  '${REPO_ROOT}''/chrombpnet/training/predict.py'
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

if ! bash "${REPO_ROOT}/scripts/paper_aligned_repro/run_tutorial_strict_compare_official.sh" --help 2>&1 | grep -q -- '--official-root'; then
  echo "ERROR: run_tutorial_strict_compare_official.sh --help does not expose --official-root" >&2
  exit 1
fi

for flag in --official-root --gc-helper-dir; do
  if ! bash "${run_remote_dataset_prep}" --help 2>&1 | grep -q -- "${flag}"; then
    echo "ERROR: run_remote_chrombpnet_dataset_prep.sh --help does not expose ${flag}" >&2
    exit 1
  fi
done
if ! bash "${run_remote_dataset_prep}" --help 2>&1 | grep -q "pass exactly one of --official-root or --gc-helper-dir"; then
  echo "ERROR: run_remote_chrombpnet_dataset_prep.sh --help does not explain the helper source requirement" >&2
  exit 1
fi

if ! rg -n "CHROMBPNET_OFFICIAL_ROOT" "${run_fast_1seed}" >/dev/null; then
  echo "ERROR: run_paper_aligned_fast_1seed.sh does not mention CHROMBPNET_OFFICIAL_ROOT" >&2
  exit 1
fi

local_gc_helper_patterns=(
  'REPO_ROOT''/chrombpnet/helpers'
  '${REPO_ROOT}''/chrombpnet/helpers'
  '../../chrombpnet/helpers'
)

for pattern in "${local_gc_helper_patterns[@]}"; do
  if rg -F -n "${pattern}" "${start_dataset_prep_6000}" "${start_dataset_prep_6002}" "${tutorial_step3}"; then
    echo "ERROR: Task 2 targets still reference local GC helpers: ${pattern}" >&2
    exit 1
  fi
done

if ! rg -n "CHROMBPNET_OFFICIAL_ROOT" "${tutorial_step3}" >/dev/null; then
  echo "ERROR: step3_get_background_regions.sh does not mention CHROMBPNET_OFFICIAL_ROOT" >&2
  exit 1
fi

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

make_fake_remote_tools() {
  local fakebin="$1"

  mkdir -p "${fakebin}"

  cat > "${fakebin}/ssh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
log_file="${FAKE_SSH_LOG:?missing FAKE_SSH_LOG}"
{
  printf 'ssh'
  for arg in "$@"; do
    printf '\t%s' "$arg"
  done
  printf '\n'
} >> "${log_file}"
last_arg="${!#}"
if [[ "${last_arg}" == *'echo $!'* ]]; then
  printf '%s\n' "${FAKE_SSH_PID:-4242}"
fi
EOF

  cat > "${fakebin}/scp" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
log_file="${FAKE_SCP_LOG:?missing FAKE_SCP_LOG}"
{
  printf 'scp'
  for arg in "$@"; do
    printf '\t%s' "$arg"
  done
  printf '\n'
} >> "${log_file}"
EOF

  chmod +x "${fakebin}/ssh" "${fakebin}/scp"
}

assert_log_contains() {
  local needle="$1"
  local file="$2"
  if ! grep -F -q -- "${needle}" "${file}"; then
    echo "ERROR: missing log entry [${needle}] in ${file}" >&2
    exit 1
  fi
}

assert_log_lacks() {
  local needle="$1"
  local file="$2"
  if grep -F -q -- "${needle}" "${file}"; then
    echo "ERROR: unexpected log entry [${needle}] in ${file}" >&2
    exit 1
  fi
}

assert_log_empty() {
  local file="$1"
  if [[ -e "${file}" && -s "${file}" ]]; then
    echo "ERROR: expected no remote command log entries in ${file}" >&2
    exit 1
  fi
}

make_fake_full_workflow_tools() {
  local fakebin="$1"
  mkdir -p "${fakebin}"

  cat > "${fakebin}/step1_download_bams_and_peaks.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
data_dir="${1:?missing data_dir}"
mkdir -p "${data_dir}"
: > "${data_dir}/merged.bam"
cat > "${data_dir}/hg38.fa" <<'FA'
>chr1
ACGT
FA
cat > "${data_dir}/hg38.chrom.sizes" <<'EOF_SIZES'
chr1	1000
EOF_SIZES
: > "${data_dir}/blacklist.bed.gz"
: > "${data_dir}/overlap.bed.gz"
EOF

  cat > "${fakebin}/step2_make_bigwigs_from_bams.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
: > "${2}_unstranded.bw"
EOF

  cat > "${fakebin}/chrombpnet_make_splits" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
outdir=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o)
      outdir="$2"
      shift 2
      ;;
    *)
      shift
      ;;
  esac
done
mkdir -p "${outdir}"
cat > "${outdir}/fold_0.json" <<'JSON'
{"fold": 0}
JSON
EOF

  cat > "${fakebin}/chrombpnet_genomewide_gc" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
prefix=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o)
      prefix="$2"
      shift 2
      ;;
    *)
      shift
      ;;
  esac
done
: > "${prefix}.bed"
EOF

  cat > "${fakebin}/step3_get_background_regions.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "${CHROMBPNET_OFFICIAL_ROOT:-}" > "${FAKE_FULL_WORKFLOW_STEP3_ENV_FILE:?missing env capture file}"
mkdir -p "${7:?missing output dir}"
cat > "${7}/negatives_with_summit.bed" <<'BED'
chr1	100	200	.	.	.	.	.	.	1057
BED
EOF

  cat > "${fakebin}/step4_train_bias_model.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
mkdir -p "${7:?missing output dir}"
: > "${7}/bias.h5"
EOF

  cat > "${fakebin}/bedtools" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
case "${1:?missing subcommand}" in
  slop|intersect)
    printf 'chr1\t0\t10\n'
    ;;
  *)
    printf 'unexpected subcommand: %s\n' "$1" >&2
    exit 1
    ;;
esac
EOF

  cat > "${fakebin}/shuf" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
if [[ "${1:-}" == "-n" ]]; then
  shift 2
fi
cat
EOF

  cat > "${fakebin}/step5_interpret_bias_model.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
exit 0
EOF

  cat > "${fakebin}/step6_train_chrombpnet_model.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
mkdir -p "${7:?missing output dir}"
: > "${7}/chrombpnet_wo_bias.h5"
EOF

  cat > "${fakebin}/step7_interpret_chrombpnet_model.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
exit 0
EOF

  chmod +x \
    "${fakebin}/step1_download_bams_and_peaks.sh" \
    "${fakebin}/step2_make_bigwigs_from_bams.sh" \
    "${fakebin}/chrombpnet_make_splits" \
    "${fakebin}/chrombpnet_genomewide_gc" \
    "${fakebin}/step3_get_background_regions.sh" \
    "${fakebin}/step4_train_bias_model.sh" \
    "${fakebin}/bedtools" \
    "${fakebin}/shuf" \
    "${fakebin}/step5_interpret_bias_model.sh" \
    "${fakebin}/step6_train_chrombpnet_model.sh" \
    "${fakebin}/step7_interpret_chrombpnet_model.sh"
}

make_fake_run_remote_prep_fixture() {
  local base="$1"
  local root="${base}/root"
  local env_dir="${base}/env"
  local dataset_dir="${root}/chrombpnet_datasets/GM12878"
  local tutorial_dir="${root}/chrombpnet_tutorial/data"

  mkdir -p "${env_dir}/bin" "${dataset_dir}" "${tutorial_dir}"

  cat > "${env_dir}/bin/python" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
if [[ "${1:-}" == "-" ]]; then
  cat >/dev/null
  exit 0
fi
script_name="$(basename "${1:?missing script}")"
shift
case "${script_name}" in
  get_genomewide_gc_bins.py)
    prefix=""
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -o)
          prefix="$2"
          shift 2
          ;;
        *)
          shift
          ;;
      esac
    done
    cat > "${prefix}.bed" <<'BED'
chr1	0	10
BED
    ;;
  get_gc_content.py)
    prefix=""
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -op)
          prefix="$2"
          shift 2
          ;;
        *)
          shift
          ;;
      esac
    done
    cat > "${prefix}.bed" <<'BED'
chr1	30	40
BED
    ;;
  get_gc_matched_negatives.py)
    prefix=""
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -o)
          prefix="$2"
          shift 2
          ;;
        *)
          shift
          ;;
      esac
    done
    cat > "${prefix}.bed" <<'BED'
chr1	20	30
BED
    ;;
  *)
    printf 'unexpected fake python script: %s\n' "${script_name}" >&2
    exit 1
    ;;
esac
EOF

  cat > "${env_dir}/bin/samtools" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
subcommand="${1:?missing subcommand}"
shift
case "${subcommand}" in
  faidx)
    printf 'chr1\t4\t0\t4\t5\n' > "${1:?missing fasta}.fai"
    ;;
  quickcheck)
    exit 0
    ;;
  index)
    bam="${@: -1}"
    printf 'fake bai\n' > "${bam}.bai"
    ;;
  merge)
    out=""
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -f)
          out="$2"
          shift 2
          ;;
        -@)
          shift 2
          ;;
        *)
          cat "$1" > "${out}"
          break
          ;;
      esac
    done
    ;;
  *)
    printf 'unexpected fake samtools subcommand: %s\n' "${subcommand}" >&2
    exit 1
    ;;
esac
EOF

  cat > "${env_dir}/bin/bedtools" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
subcommand="${1:?missing subcommand}"
shift
case "${subcommand}" in
  bamtobed)
    printf 'chr1\t0\t10\t.\t0\t+\n'
    ;;
  genomecov)
    cat
    ;;
  slop)
    input=""
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -i)
          input="$2"
          shift 2
          ;;
        -g|-b)
          shift 2
          ;;
        *)
          shift
          ;;
      esac
    done
    cat "${input}"
    ;;
  sort|merge)
    cat
    ;;
  intersect)
    input=""
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -a)
          input="$2"
          shift 2
          ;;
        -b)
          shift 2
          ;;
        -v)
          shift
          ;;
        *)
          shift
          ;;
      esac
    done
    cat "${input}"
    ;;
  *)
    printf 'unexpected fake bedtools subcommand: %s\n' "${subcommand}" >&2
    exit 1
    ;;
esac
EOF

  cat > "${env_dir}/bin/bedGraphToBigWig" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
cp "${1:?missing bedgraph}" "${3:?missing bigwig}"
EOF

  cat > "${env_dir}/bin/gzip" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
if [[ "${1:-}" == "-dc" ]]; then
  cat "${2:?missing file}"
  exit 0
fi
printf 'unexpected fake gzip invocation\n' >&2
exit 1
EOF

  chmod +x \
    "${env_dir}/bin/python" \
    "${env_dir}/bin/samtools" \
    "${env_dir}/bin/bedtools" \
    "${env_dir}/bin/bedGraphToBigWig" \
    "${env_dir}/bin/gzip"

  cat > "${root}/hg38.fa" <<'EOF'
>chr1
ACGT
EOF
  cat > "${tutorial_dir}/hg38.chrom.sizes" <<'EOF'
chr1	1000
EOF
  cat > "${tutorial_dir}/blacklist.bed.gz" <<'EOF'
chr1	10	20
EOF
  cat > "${tutorial_dir}/folds.json" <<'EOF'
{"fold": 0}
EOF
  cat > "${dataset_dir}/overlap.bed" <<'EOF'
chr1	30	40
EOF
  cat > "${dataset_dir}/rep1.bam" <<'EOF'
fake bam
EOF
}

make_fake_official_helper_root() {
  local root="$1"
  mkdir -p "${root}/chrombpnet/helpers/make_gc_matched_negatives/get_genomewide_gc_buckets"
  : > "${root}/chrombpnet/helpers/make_gc_matched_negatives/get_gc_content.py"
  : > "${root}/chrombpnet/helpers/make_gc_matched_negatives/get_gc_matched_negatives.py"
  : > "${root}/chrombpnet/helpers/make_gc_matched_negatives/get_genomewide_gc_buckets/get_genomewide_gc_bins.py"
}

make_fake_staged_gc_helper_dir() {
  local root="$1"
  mkdir -p "${root}"
  : > "${root}/get_gc_content.py"
  : > "${root}/get_gc_matched_negatives.py"
  : > "${root}/get_genomewide_gc_bins.py"
}

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

prep_helper_tmp="${tmpdir}/dataset_prep_helper_checks"
mkdir -p "${prep_helper_tmp}/root" "${prep_helper_tmp}/env"

if bash "${run_remote_dataset_prep}" \
  --root "${prep_helper_tmp}/root" \
  --env-dir "${prep_helper_tmp}/env" \
  2> "${prep_helper_tmp}/missing_helper_source.stderr"; then
  echo "ERROR: run_remote_chrombpnet_dataset_prep.sh unexpectedly succeeded without helper source flags" >&2
  exit 1
fi
if ! grep -q "pass --official-root or --gc-helper-dir" "${prep_helper_tmp}/missing_helper_source.stderr"; then
  echo "ERROR: run_remote_chrombpnet_dataset_prep.sh did not emit the expected missing helper source error" >&2
  exit 1
fi

missing_gc_helper_dir="${prep_helper_tmp}/missing_gc_helper_dir"
if bash "${run_remote_dataset_prep}" \
  --root "${prep_helper_tmp}/root" \
  --env-dir "${prep_helper_tmp}/env" \
  --gc-helper-dir "${missing_gc_helper_dir}" \
  2> "${prep_helper_tmp}/missing_gc_helper_dir.stderr"; then
  echo "ERROR: run_remote_chrombpnet_dataset_prep.sh unexpectedly succeeded with a missing gc helper dir" >&2
  exit 1
fi
if ! grep -q "gc helper dir is not a directory" "${prep_helper_tmp}/missing_gc_helper_dir.stderr"; then
  echo "ERROR: run_remote_chrombpnet_dataset_prep.sh did not emit the expected gc helper dir validation error" >&2
  exit 1
fi

missing_official_helper_root="${prep_helper_tmp}/missing_official_root"
if bash "${run_remote_dataset_prep}" \
  --root "${prep_helper_tmp}/root" \
  --env-dir "${prep_helper_tmp}/env" \
  --official-root "${missing_official_helper_root}" \
  2> "${prep_helper_tmp}/missing_official_root.stderr"; then
  echo "ERROR: run_remote_chrombpnet_dataset_prep.sh unexpectedly succeeded with a missing official root" >&2
  exit 1
fi
if ! grep -q "official ChromBPNet root is not a directory" "${prep_helper_tmp}/missing_official_root.stderr"; then
  echo "ERROR: run_remote_chrombpnet_dataset_prep.sh did not emit the expected official root validation error" >&2
  exit 1
fi

mkdir -p "${prep_helper_tmp}/conflict_official/chrombpnet/helpers/make_gc_matched_negatives/get_genomewide_gc_buckets" \
  "${prep_helper_tmp}/conflict_gc_helper_dir"
if bash "${run_remote_dataset_prep}" \
  --root "${prep_helper_tmp}/root" \
  --env-dir "${prep_helper_tmp}/env" \
  --official-root "${prep_helper_tmp}/conflict_official" \
  --gc-helper-dir "${prep_helper_tmp}/conflict_gc_helper_dir" \
  2> "${prep_helper_tmp}/conflicting_helper_sources.stderr"; then
  echo "ERROR: run_remote_chrombpnet_dataset_prep.sh unexpectedly succeeded with conflicting helper source flags" >&2
  exit 1
fi
if ! grep -q "pass exactly one of --official-root or --gc-helper-dir" "${prep_helper_tmp}/conflicting_helper_sources.stderr"; then
  echo "ERROR: run_remote_chrombpnet_dataset_prep.sh did not emit the expected conflicting helper source error" >&2
  exit 1
fi

launcher_tmp="${tmpdir}/launcher_contracts"
mkdir -p "${launcher_tmp}/fakebin"
make_fake_remote_tools "${launcher_tmp}/fakebin"

FAKE_SSH_LOG="${launcher_tmp}/ssh.log" \
FAKE_SCP_LOG="${launcher_tmp}/scp.log" \
FAKE_SSH_PID=6000 \
PATH="${launcher_tmp}/fakebin:${PATH}" \
REMOTE_HOST="remote6000" \
REMOTE_PORT="6600" \
REMOTE_ROOT="/srv/remote6000" \
REMOTE_ENV="/envs/chrombpnet" \
REMOTE_PYTHON="/envs/chrombpnet/bin/python" \
CHROMBPNET_OFFICIAL_ROOT="/official/chrombpnet" \
DATASETS="GM12878" \
THREADS="8" \
NICE_LEVEL="3" \
RUN_TAG="unit_6000" \
bash "${start_dataset_prep_6000}" > "${launcher_tmp}/start6000.stdout"

assert_log_contains $'scp\t-P\t6600\t'"${REPO_ROOT}/scripts/run_remote_chrombpnet_dataset_prep.sh"$'\tremote6000:/srv/remote6000/.codex_jobs/chrombpnet_dataset_prep/unit_6000/' "${launcher_tmp}/scp.log"
assert_log_lacks "chrombpnet/helpers" "${launcher_tmp}/scp.log"
assert_log_contains "get_gc_content.py" "${launcher_tmp}/ssh.log"
assert_log_contains "get_gc_matched_negatives.py" "${launcher_tmp}/ssh.log"
assert_log_contains "get_genomewide_gc_bins.py" "${launcher_tmp}/ssh.log"
assert_log_contains "--official-root '/official/chrombpnet'" "${launcher_tmp}/ssh.log"
preflight_line="$(grep -n "get_gc_content.py" "${launcher_tmp}/ssh.log" | head -n 1 | cut -d: -f1)"
nohup_line="$(grep -n "nohup bash '/srv/remote6000/.codex_jobs/chrombpnet_dataset_prep/unit_6000/run_remote_chrombpnet_dataset_prep.sh'" "${launcher_tmp}/ssh.log" | head -n 1 | cut -d: -f1)"
if [[ -z "${preflight_line}" || -z "${nohup_line}" || "${preflight_line}" -ge "${nohup_line}" ]]; then
  echo "ERROR: start_6000 preflight did not run before background launch" >&2
  exit 1
fi

: > "${launcher_tmp}/ssh.log"
: > "${launcher_tmp}/scp.log"

FAKE_SSH_LOG="${launcher_tmp}/ssh.log" \
FAKE_SCP_LOG="${launcher_tmp}/scp.log" \
FAKE_SSH_PID=6002 \
PATH="${launcher_tmp}/fakebin:${PATH}" \
REMOTE_HOST="remote6002" \
REMOTE_PORT="6602" \
REMOTE_KEY="/tmp/fake_6002_key" \
REMOTE_ROOT="/srv/remote6002" \
REMOTE_ENV="/envs/transchrombp" \
REMOTE_PYTHON="/envs/transchrombp/bin/python" \
SOURCE_HOST="source6000" \
SOURCE_PORT="6610" \
SOURCE_OFFICIAL_ROOT="/official/source_root" \
DATASETS="K562" \
THREADS="6" \
NICE_LEVEL="4" \
RUN_TAG="unit_6002" \
bash "${start_dataset_prep_6002}" > "${launcher_tmp}/start6002.stdout"

assert_log_contains $'scp\t-P\t6610\tsource6000:/official/source_root/chrombpnet/helpers/make_gc_matched_negatives/get_gc_content.py' "${launcher_tmp}/scp.log"
assert_log_contains "source6000:/official/source_root/chrombpnet/helpers/make_gc_matched_negatives/get_gc_matched_negatives.py" "${launcher_tmp}/scp.log"
assert_log_contains "source6000:/official/source_root/chrombpnet/helpers/make_gc_matched_negatives/get_genomewide_gc_buckets/get_genomewide_gc_bins.py" "${launcher_tmp}/scp.log"
assert_log_contains $'scp\t-i\t/tmp/fake_6002_key\t-P\t6602\t'"${REPO_ROOT}/scripts/run_remote_chrombpnet_dataset_prep.sh" "${launcher_tmp}/scp.log"
assert_log_contains "get_gc_content.py" "${launcher_tmp}/scp.log"
assert_log_contains "get_gc_matched_negatives.py" "${launcher_tmp}/scp.log"
assert_log_contains "get_genomewide_gc_bins.py" "${launcher_tmp}/scp.log"
assert_log_contains "remote6002:/srv/remote6002/.codex_jobs/chrombpnet_dataset_prep/unit_6002/" "${launcher_tmp}/scp.log"
assert_log_contains "--gc-helper-dir '/srv/remote6002/.codex_jobs/chrombpnet_dataset_prep/unit_6002'" "${launcher_tmp}/ssh.log"
assert_log_lacks "${REPO_ROOT}/chrombpnet/helpers" "${launcher_tmp}/scp.log"

: > "${launcher_tmp}/ssh.log"
: > "${launcher_tmp}/scp.log"

if FAKE_SSH_LOG="${launcher_tmp}/ssh.log" \
  FAKE_SCP_LOG="${launcher_tmp}/scp.log" \
  PATH="${launcher_tmp}/fakebin:${PATH}" \
  REMOTE_HOST="remote6000" \
  REMOTE_PORT="6600" \
  REMOTE_ROOT="/srv/remote6000" \
  REMOTE_ENV="/envs/chrombpnet" \
  REMOTE_PYTHON="/envs/chrombpnet/bin/python" \
  CHROMBPNET_OFFICIAL_ROOT="/official/chrombpnet" \
  DATASETS="GM12878 K562" \
  THREADS="8" \
  NICE_LEVEL="3" \
  RUN_TAG="unit_unsafe_6000" \
  bash "${start_dataset_prep_6000}" > "${launcher_tmp}/unsafe6000.stdout" 2> "${launcher_tmp}/unsafe6000.stderr"; then
  echo "ERROR: start_6000 unexpectedly accepted unsafe dataset input" >&2
  exit 1
fi
if ! grep -q "unsafe value for DATASETS" "${launcher_tmp}/unsafe6000.stderr"; then
  echo "ERROR: start_6000 did not emit the expected unsafe-value error" >&2
  exit 1
fi
assert_log_empty "${launcher_tmp}/ssh.log"
assert_log_empty "${launcher_tmp}/scp.log"

: > "${launcher_tmp}/ssh.log"
: > "${launcher_tmp}/scp.log"

if FAKE_SSH_LOG="${launcher_tmp}/ssh.log" \
  FAKE_SCP_LOG="${launcher_tmp}/scp.log" \
  PATH="${launcher_tmp}/fakebin:${PATH}" \
  REMOTE_HOST="-oProxyCommand=bad" \
  REMOTE_PORT="6600" \
  REMOTE_ROOT="/srv/remote6000" \
  REMOTE_ENV="/envs/chrombpnet" \
  REMOTE_PYTHON="/envs/chrombpnet/bin/python" \
  CHROMBPNET_OFFICIAL_ROOT="/official/chrombpnet" \
  DATASETS="GM12878" \
  THREADS="8" \
  NICE_LEVEL="3" \
  RUN_TAG="unit_unsafe_host_6000" \
  bash "${start_dataset_prep_6000}" > "${launcher_tmp}/unsafe_host_6000.stdout" 2> "${launcher_tmp}/unsafe_host_6000.stderr"; then
  echo "ERROR: start_6000 unexpectedly accepted unsafe REMOTE_HOST input" >&2
  exit 1
fi
if ! grep -q "unsafe value for REMOTE_HOST" "${launcher_tmp}/unsafe_host_6000.stderr"; then
  echo "ERROR: start_6000 did not emit the expected REMOTE_HOST unsafe-value error" >&2
  exit 1
fi
assert_log_empty "${launcher_tmp}/ssh.log"
assert_log_empty "${launcher_tmp}/scp.log"

: > "${launcher_tmp}/ssh.log"
: > "${launcher_tmp}/scp.log"

if FAKE_SSH_LOG="${launcher_tmp}/ssh.log" \
  FAKE_SCP_LOG="${launcher_tmp}/scp.log" \
  PATH="${launcher_tmp}/fakebin:${PATH}" \
  REMOTE_HOST="remote6000" \
  REMOTE_PORT="6600" \
  REMOTE_ROOT="/srv/remote6000" \
  REMOTE_ENV="/envs/chrombpnet" \
  REMOTE_PYTHON="/envs/chrombpnet/bin/python" \
  CHROMBPNET_OFFICIAL_ROOT="/official/chrombpnet" \
  DATASETS="GM12878" \
  THREADS="8" \
  NICE_LEVEL="3" \
  RUN_TAG="../escaped_path" \
  bash "${start_dataset_prep_6000}" > "${launcher_tmp}/unsafe_runtag_6000.stdout" 2> "${launcher_tmp}/unsafe_runtag_6000.stderr"; then
  echo "ERROR: start_6000 unexpectedly accepted path-traversal RUN_TAG input" >&2
  exit 1
fi
if ! grep -q "unsafe value for RUN_TAG" "${launcher_tmp}/unsafe_runtag_6000.stderr"; then
  echo "ERROR: start_6000 did not emit the expected RUN_TAG unsafe-value error" >&2
  exit 1
fi
assert_log_empty "${launcher_tmp}/ssh.log"
assert_log_empty "${launcher_tmp}/scp.log"

: > "${launcher_tmp}/ssh.log"
: > "${launcher_tmp}/scp.log"

if FAKE_SSH_LOG="${launcher_tmp}/ssh.log" \
  FAKE_SCP_LOG="${launcher_tmp}/scp.log" \
  PATH="${launcher_tmp}/fakebin:${PATH}" \
  REMOTE_HOST="remote6000" \
  REMOTE_PORT="6600" \
  REMOTE_ROOT="/srv/remote6000" \
  REMOTE_ENV="/envs/chrombpnet" \
  REMOTE_PYTHON="/envs/chrombpnet/bin/python" \
  CHROMBPNET_OFFICIAL_ROOT="/official/chrombpnet" \
  DATASETS="GM12878" \
  THREADS="8" \
  NICE_LEVEL="3" \
  RUN_TAG="." \
  bash "${start_dataset_prep_6000}" > "${launcher_tmp}/unsafe_runtag_dot_6000.stdout" 2> "${launcher_tmp}/unsafe_runtag_dot_6000.stderr"; then
  echo "ERROR: start_6000 unexpectedly accepted dot RUN_TAG input" >&2
  exit 1
fi
if ! grep -q "unsafe value for RUN_TAG" "${launcher_tmp}/unsafe_runtag_dot_6000.stderr"; then
  echo "ERROR: start_6000 did not reject dot RUN_TAG input as expected" >&2
  exit 1
fi
assert_log_empty "${launcher_tmp}/ssh.log"
assert_log_empty "${launcher_tmp}/scp.log"

: > "${launcher_tmp}/ssh.log"
: > "${launcher_tmp}/scp.log"

if FAKE_SSH_LOG="${launcher_tmp}/ssh.log" \
  FAKE_SCP_LOG="${launcher_tmp}/scp.log" \
  PATH="${launcher_tmp}/fakebin:${PATH}" \
  REMOTE_HOST="remote6002" \
  REMOTE_PORT="6602" \
  REMOTE_KEY="/tmp/fake_6002_key" \
  REMOTE_ROOT="/srv/remote6002" \
  REMOTE_ENV="/envs/transchrombp" \
  REMOTE_PYTHON="/envs/transchrombp/bin/python" \
  SOURCE_HOST="source6000" \
  SOURCE_PORT="6610" \
  SOURCE_OFFICIAL_ROOT="/official/\$unsafe" \
  DATASETS="K562" \
  THREADS="6" \
  NICE_LEVEL="4" \
  RUN_TAG="unit_unsafe_6002" \
  bash "${start_dataset_prep_6002}" > "${launcher_tmp}/unsafe6002.stdout" 2> "${launcher_tmp}/unsafe6002.stderr"; then
  echo "ERROR: start_6002 unexpectedly accepted unsafe source official root input" >&2
  exit 1
fi
if ! grep -q "unsafe value for SOURCE_OFFICIAL_ROOT" "${launcher_tmp}/unsafe6002.stderr"; then
  echo "ERROR: start_6002 did not emit the expected unsafe-value error" >&2
  exit 1
fi
assert_log_empty "${launcher_tmp}/ssh.log"
assert_log_empty "${launcher_tmp}/scp.log"

: > "${launcher_tmp}/ssh.log"
: > "${launcher_tmp}/scp.log"

if FAKE_SSH_LOG="${launcher_tmp}/ssh.log" \
  FAKE_SCP_LOG="${launcher_tmp}/scp.log" \
  PATH="${launcher_tmp}/fakebin:${PATH}" \
  REMOTE_HOST="remote6002" \
  REMOTE_PORT="6602" \
  REMOTE_KEY="/tmp/fake_6002_key" \
  REMOTE_ROOT="/srv/remote6002" \
  REMOTE_ENV="/envs/transchrombp" \
  REMOTE_PYTHON="/envs/transchrombp/bin/python" \
  SOURCE_HOST="-oProxyCommand=bad" \
  SOURCE_PORT="6610" \
  SOURCE_OFFICIAL_ROOT="/official/source_root" \
  DATASETS="K562" \
  THREADS="6" \
  NICE_LEVEL="4" \
  RUN_TAG="unit_unsafe_host_6002" \
  bash "${start_dataset_prep_6002}" > "${launcher_tmp}/unsafe_host_6002.stdout" 2> "${launcher_tmp}/unsafe_host_6002.stderr"; then
  echo "ERROR: start_6002 unexpectedly accepted unsafe SOURCE_HOST input" >&2
  exit 1
fi
if ! grep -q "unsafe value for SOURCE_HOST" "${launcher_tmp}/unsafe_host_6002.stderr"; then
  echo "ERROR: start_6002 did not emit the expected SOURCE_HOST unsafe-value error" >&2
  exit 1
fi
assert_log_empty "${launcher_tmp}/ssh.log"
assert_log_empty "${launcher_tmp}/scp.log"

: > "${launcher_tmp}/ssh.log"
: > "${launcher_tmp}/scp.log"

if FAKE_SSH_LOG="${launcher_tmp}/ssh.log" \
  FAKE_SCP_LOG="${launcher_tmp}/scp.log" \
  PATH="${launcher_tmp}/fakebin:${PATH}" \
  REMOTE_HOST="remote6002" \
  REMOTE_PORT="6602" \
  REMOTE_KEY="/tmp/fake_6002_key" \
  REMOTE_ROOT="/srv/remote6002" \
  REMOTE_ENV="/envs/transchrombp" \
  REMOTE_PYTHON="/envs/transchrombp/bin/python" \
  SOURCE_HOST="source6000" \
  SOURCE_PORT="6610" \
  SOURCE_OFFICIAL_ROOT="/official/source_root" \
  DATASETS="K562" \
  THREADS="6" \
  NICE_LEVEL="4" \
  RUN_TAG="../escaped_path_6002" \
  bash "${start_dataset_prep_6002}" > "${launcher_tmp}/unsafe_runtag_6002.stdout" 2> "${launcher_tmp}/unsafe_runtag_6002.stderr"; then
  echo "ERROR: start_6002 unexpectedly accepted path-traversal RUN_TAG input" >&2
  exit 1
fi
if ! grep -q "unsafe value for RUN_TAG" "${launcher_tmp}/unsafe_runtag_6002.stderr"; then
  echo "ERROR: start_6002 did not emit the expected RUN_TAG unsafe-value error" >&2
  exit 1
fi
assert_log_empty "${launcher_tmp}/ssh.log"
assert_log_empty "${launcher_tmp}/scp.log"

: > "${launcher_tmp}/ssh.log"
: > "${launcher_tmp}/scp.log"

if FAKE_SSH_LOG="${launcher_tmp}/ssh.log" \
  FAKE_SCP_LOG="${launcher_tmp}/scp.log" \
  PATH="${launcher_tmp}/fakebin:${PATH}" \
  REMOTE_HOST="remote6002" \
  REMOTE_PORT="6602" \
  REMOTE_KEY="/tmp/fake_6002_key" \
  REMOTE_ROOT="/srv/remote6002" \
  REMOTE_ENV="/envs/transchrombp" \
  REMOTE_PYTHON="/envs/transchrombp/bin/python" \
  SOURCE_HOST="source6000" \
  SOURCE_PORT="6610" \
  SOURCE_OFFICIAL_ROOT="/official/source_root" \
  DATASETS="K562" \
  THREADS="6" \
  NICE_LEVEL="4" \
  RUN_TAG=".." \
  bash "${start_dataset_prep_6002}" > "${launcher_tmp}/unsafe_runtag_dotdot_6002.stdout" 2> "${launcher_tmp}/unsafe_runtag_dotdot_6002.stderr"; then
  echo "ERROR: start_6002 unexpectedly accepted dotdot RUN_TAG input" >&2
  exit 1
fi
if ! grep -q "unsafe value for RUN_TAG" "${launcher_tmp}/unsafe_runtag_dotdot_6002.stderr"; then
  echo "ERROR: start_6002 did not reject dotdot RUN_TAG input as expected" >&2
  exit 1
fi
assert_log_empty "${launcher_tmp}/ssh.log"
assert_log_empty "${launcher_tmp}/scp.log"

step3_smoke_tmp="${tmpdir}/step3_smoke"
mkdir -p "${step3_smoke_tmp}/fakebin" "${step3_smoke_tmp}/out"
cat > "${step3_smoke_tmp}/fakebin/bedtools" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
subcommand="${1:?missing subcommand}"
shift
case "${subcommand}" in
  slop)
    input=""
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -i)
          input="$2"
          shift 2
          ;;
        -g|-b)
          shift 2
          ;;
        *)
          shift
          ;;
      esac
    done
    cat "${input}"
    ;;
  sort)
    cat
    ;;
  merge)
    cat
    ;;
  *)
    echo "unexpected bedtools subcommand: ${subcommand}" >&2
    exit 1
    ;;
esac
EOF
chmod +x "${step3_smoke_tmp}/fakebin/bedtools"

cat > "${step3_smoke_tmp}/stub_make_gc_matched_negatives.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$*" > "$4/helper_args.txt"
cat > "$4/negatives.bed" <<'BED'
chr1	100	200
BED
EOF
chmod +x "${step3_smoke_tmp}/stub_make_gc_matched_negatives.sh"

cat > "${step3_smoke_tmp}/reference.fa" <<'EOF'
>chr1
ACGT
EOF
cat > "${step3_smoke_tmp}/chrom.sizes" <<'EOF'
chr1	1000
EOF
cat > "${step3_smoke_tmp}/blacklist.bed" <<'EOF'
chr1	10	20
EOF
cat > "${step3_smoke_tmp}/overlap.bed" <<'EOF'
chr1	30	40
EOF
cat > "${step3_smoke_tmp}/genomewide_gc.bed" <<'EOF'
chr1	0	100
EOF
cat > "${step3_smoke_tmp}/fold.json" <<'EOF'
{"fold": 0}
EOF

PATH="${step3_smoke_tmp}/fakebin:${PATH}" \
MAKE_GC_MATCHED_NEGATIVES_SCRIPT="${step3_smoke_tmp}/stub_make_gc_matched_negatives.sh" \
bash "${tutorial_step3}" \
  "${step3_smoke_tmp}/reference.fa" \
  "${step3_smoke_tmp}/chrom.sizes" \
  "${step3_smoke_tmp}/blacklist.bed" \
  "${step3_smoke_tmp}/overlap.bed" \
  2114 \
  "${step3_smoke_tmp}/genomewide_gc.bed" \
  "${step3_smoke_tmp}/out" \
  "${step3_smoke_tmp}/fold.json"

if [[ ! -s "${step3_smoke_tmp}/out/negatives_with_summit.bed" ]]; then
  echo "ERROR: step3 smoke test did not produce negatives_with_summit.bed" >&2
  exit 1
fi
if ! grep -q $'^chr1\t100\t200' "${step3_smoke_tmp}/out/negatives_with_summit.bed"; then
  echo "ERROR: step3 smoke test produced unexpected negatives_with_summit.bed contents" >&2
  exit 1
fi
if ! grep -F -q -- "${step3_smoke_tmp}/overlap.bed ${step3_smoke_tmp}/out/exclude.bed 2114 ${step3_smoke_tmp}/out" "${step3_smoke_tmp}/out/helper_args.txt"; then
  echo "ERROR: step3 smoke test helper was not invoked with the expected arguments" >&2
  exit 1
fi

step3_official_tmp="${tmpdir}/step3_official_smoke"
mkdir -p "${step3_official_tmp}/official_root/chrombpnet/helpers/make_gc_matched_negatives" "${step3_official_tmp}/out"
cp "${step3_smoke_tmp}/fakebin/bedtools" "${step3_official_tmp}/bedtools"
mkdir -p "${step3_official_tmp}/fakebin"
mv "${step3_official_tmp}/bedtools" "${step3_official_tmp}/fakebin/bedtools"
cat > "${step3_official_tmp}/official_root/chrombpnet/helpers/make_gc_matched_negatives/make_gc_matched_negatives.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$*" > "$4/helper_args.txt"
cat > "$4/negatives.bed" <<'BED'
chr1	300	400
BED
EOF
chmod +x "${step3_official_tmp}/official_root/chrombpnet/helpers/make_gc_matched_negatives/make_gc_matched_negatives.sh" "${step3_official_tmp}/fakebin/bedtools"

PATH="${step3_official_tmp}/fakebin:${PATH}" \
CHROMBPNET_OFFICIAL_ROOT="${step3_official_tmp}/official_root" \
bash "${tutorial_step3}" \
  "${step3_smoke_tmp}/reference.fa" \
  "${step3_smoke_tmp}/chrom.sizes" \
  "${step3_smoke_tmp}/blacklist.bed" \
  "${step3_smoke_tmp}/overlap.bed" \
  2114 \
  "${step3_smoke_tmp}/genomewide_gc.bed" \
  "${step3_official_tmp}/out" \
  "${step3_smoke_tmp}/fold.json"

if ! grep -q $'^chr1\t300\t400' "${step3_official_tmp}/out/negatives_with_summit.bed"; then
  echo "ERROR: step3 official-root smoke test produced unexpected negatives_with_summit.bed contents" >&2
  exit 1
fi

full_workflow_tmp="${tmpdir}/full_workflow_smoke"
mkdir -p "${full_workflow_tmp}/case_a" "${full_workflow_tmp}/case_b" "${full_workflow_tmp}/case_c" "${full_workflow_tmp}/fakebin"
make_fake_full_workflow_tools "${full_workflow_tmp}/fakebin"

(
  cd "${full_workflow_tmp}/case_a"
  PATH="${full_workflow_tmp}/fakebin:${PATH}" \
  CHROMBPNET_OFFICIAL_ROOT="/custom/root" \
  FAKE_FULL_WORKFLOW_STEP3_ENV_FILE="${full_workflow_tmp}/case_a/step3_env.txt" \
  bash "${full_workflow_test}" 0 > "${full_workflow_tmp}/case_a/stdout" 2> "${full_workflow_tmp}/case_a/stderr"
)

if [[ "$(<"${full_workflow_tmp}/case_a/step3_env.txt")" != "/custom/root" ]]; then
  echo "ERROR: tests/full_workflow.sh did not pass through the pre-set CHROMBPNET_OFFICIAL_ROOT" >&2
  exit 1
fi

mkdir -p "${full_workflow_tmp}/fallback_root"
(
  cd "${full_workflow_tmp}/case_c"
  env -u CHROMBPNET_OFFICIAL_ROOT \
    PATH="${full_workflow_tmp}/fakebin:${PATH}" \
    CHROMBPNET_OFFICIAL_ROOT_FALLBACK="${full_workflow_tmp}/fallback_root" \
    FAKE_FULL_WORKFLOW_STEP3_ENV_FILE="${full_workflow_tmp}/case_c/step3_env.txt" \
    bash "${full_workflow_test}" 0 > "${full_workflow_tmp}/case_c/stdout" 2> "${full_workflow_tmp}/case_c/stderr"
)

if [[ "$(<"${full_workflow_tmp}/case_c/step3_env.txt")" != "${full_workflow_tmp}/fallback_root" ]]; then
  echo "ERROR: tests/full_workflow.sh did not pass through the existing fallback official root" >&2
  exit 1
fi

if (
  cd "${full_workflow_tmp}/case_b"
  env -u CHROMBPNET_OFFICIAL_ROOT \
    PATH="${full_workflow_tmp}/fakebin:${PATH}" \
    CHROMBPNET_OFFICIAL_ROOT_FALLBACK="${full_workflow_tmp}/missing_default_root" \
    FAKE_FULL_WORKFLOW_STEP3_ENV_FILE="${full_workflow_tmp}/case_b/step3_env.txt" \
    bash "${full_workflow_test}" 0 > "${full_workflow_tmp}/case_b/stdout" 2> "${full_workflow_tmp}/case_b/stderr"
); then
  echo "ERROR: tests/full_workflow.sh unexpectedly succeeded without CHROMBPNET_OFFICIAL_ROOT and without a fallback root" >&2
  exit 1
fi
if ! grep -q "set CHROMBPNET_OFFICIAL_ROOT" "${full_workflow_tmp}/case_b/stderr"; then
  echo "ERROR: tests/full_workflow.sh did not emit the expected missing official root guidance" >&2
  exit 1
fi
if [[ -e "${full_workflow_tmp}/case_b/step3_env.txt" ]]; then
  echo "ERROR: tests/full_workflow.sh reached step 3 despite missing helper root configuration" >&2
  exit 1
fi

run_remote_smoke_tmp="${tmpdir}/run_remote_success_smoke"
mkdir -p "${run_remote_smoke_tmp}/official_mode" "${run_remote_smoke_tmp}/gc_helper_mode"
make_fake_run_remote_prep_fixture "${run_remote_smoke_tmp}/official_mode"
make_fake_official_helper_root "${run_remote_smoke_tmp}/official_mode/official_root"

bash "${run_remote_dataset_prep}" \
  --root "${run_remote_smoke_tmp}/official_mode/root" \
  --env-dir "${run_remote_smoke_tmp}/official_mode/env" \
  --python-bin "${run_remote_smoke_tmp}/official_mode/env/bin/python" \
  --datasets GM12878 \
  --official-root "${run_remote_smoke_tmp}/official_mode/official_root" \
  > "${run_remote_smoke_tmp}/official_mode/stdout" \
  2> "${run_remote_smoke_tmp}/official_mode/stderr"

for path in \
  "${run_remote_smoke_tmp}/official_mode/root/chrombpnet_refs/genomewide_gc_hg38_stride_1000_inputlen_2114.bed" \
  "${run_remote_smoke_tmp}/official_mode/root/chrombpnet_datasets/GM12878/merged.bam" \
  "${run_remote_smoke_tmp}/official_mode/root/chrombpnet_datasets/GM12878/merged.bam.bai" \
  "${run_remote_smoke_tmp}/official_mode/root/chrombpnet_datasets/GM12878/prep_v1/bigwig/merged_unstranded.bw" \
  "${run_remote_smoke_tmp}/official_mode/root/chrombpnet_datasets/GM12878/prep_v1/background/tutorial_folds/nonpeaks_tutorial_folds.bed"; do
  if [[ ! -s "${path}" ]]; then
    echo "ERROR: run_remote official-root smoke missing expected output: ${path}" >&2
    exit 1
  fi
done

make_fake_run_remote_prep_fixture "${run_remote_smoke_tmp}/gc_helper_mode"
make_fake_staged_gc_helper_dir "${run_remote_smoke_tmp}/gc_helper_mode/staged_gc_helpers"

bash "${run_remote_dataset_prep}" \
  --root "${run_remote_smoke_tmp}/gc_helper_mode/root" \
  --env-dir "${run_remote_smoke_tmp}/gc_helper_mode/env" \
  --python-bin "${run_remote_smoke_tmp}/gc_helper_mode/env/bin/python" \
  --datasets GM12878 \
  --gc-helper-dir "${run_remote_smoke_tmp}/gc_helper_mode/staged_gc_helpers" \
  > "${run_remote_smoke_tmp}/gc_helper_mode/stdout" \
  2> "${run_remote_smoke_tmp}/gc_helper_mode/stderr"

for path in \
  "${run_remote_smoke_tmp}/gc_helper_mode/root/chrombpnet_refs/genomewide_gc_hg38_stride_1000_inputlen_2114.bed" \
  "${run_remote_smoke_tmp}/gc_helper_mode/root/chrombpnet_datasets/GM12878/merged.bam" \
  "${run_remote_smoke_tmp}/gc_helper_mode/root/chrombpnet_datasets/GM12878/merged.bam.bai" \
  "${run_remote_smoke_tmp}/gc_helper_mode/root/chrombpnet_datasets/GM12878/prep_v1/bigwig/merged_unstranded.bw" \
  "${run_remote_smoke_tmp}/gc_helper_mode/root/chrombpnet_datasets/GM12878/prep_v1/background/tutorial_folds/nonpeaks_tutorial_folds.bed"; do
  if [[ ! -s "${path}" ]]; then
    echo "ERROR: run_remote gc-helper-dir smoke missing expected output: ${path}" >&2
    exit 1
  fi
done

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
