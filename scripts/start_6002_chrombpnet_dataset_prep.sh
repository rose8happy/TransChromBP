#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

REMOTE_HOST="${REMOTE_HOST:-zhengwei@127.0.0.1}"
REMOTE_PORT="${REMOTE_PORT:-6002}"
REMOTE_KEY="${REMOTE_KEY:-/home/zhengwei/.ssh/codex_6002_ed25519}"
REMOTE_ROOT="${REMOTE_ROOT:-/home/zhengwei/bylw_atac}"
REMOTE_ENV="${REMOTE_ENV:-/home/zhengwei/bylw_atac/.mamba/envs/transchrombp}"
REMOTE_PYTHON="${REMOTE_PYTHON:-$REMOTE_ENV/bin/python}"
SOURCE_HOST="${SOURCE_HOST:-zhoujiazhen@127.0.0.1}"
SOURCE_PORT="${SOURCE_PORT:-6000}"
SOURCE_OFFICIAL_ROOT="${SOURCE_OFFICIAL_ROOT:-${CHROMBPNET_OFFICIAL_ROOT:-/data1/zhoujiazhen/bylw_atac/chrombpnet_official}}"
DATASETS="${DATASETS:-GM12878,K562}"
THREADS="${THREADS:-4}"
NICE_LEVEL="${NICE_LEVEL:-10}"
RUN_TAG="${RUN_TAG:-chrombpnet_dataset_prep_6002_$(date +%Y%m%d_%H%M%S)}"
LOCAL_STAGE_DIR="$(mktemp -d)"

REMOTE_JOB_DIR="$REMOTE_ROOT/.codex_jobs/chrombpnet_dataset_prep/$RUN_TAG"
REMOTE_LOG="$REMOTE_ROOT/logs/${RUN_TAG}.log"

trap 'rm -rf "${LOCAL_STAGE_DIR}"' EXIT

fail() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

contains_unsafe_remote_chars() {
  local value="$1"
  [[ "$value" == *"'"* || "$value" == *'"'* || "$value" == *'$'* || "$value" == *'`'* || "$value" == *'\'* || "$value" == *$' '* || "$value" == *$'\t'* || "$value" == *$'\n'* || "$value" == *$'\r'* ]]
}

validate_host_value() {
  local name="$1"
  local value="$2"
  if [[ -z "$value" || "$value" == -* || ! "$value" =~ ^[A-Za-z0-9_.@-]+$ ]]; then
    fail "unsafe value for ${name}: ${value}"
  fi
}

validate_remote_value() {
  local name="$1"
  local value="$2"
  if contains_unsafe_remote_chars "$value"; then
    fail "unsafe value for ${name}: ${value}"
  fi
}

validate_run_tag() {
  local value="$1"
  if [[ -z "$value" || ! "$value" =~ ^[A-Za-z0-9._-]+$ ]]; then
    fail "unsafe value for RUN_TAG: ${value}"
  fi
}

validate_host_value "REMOTE_HOST" "$REMOTE_HOST"
validate_host_value "SOURCE_HOST" "$SOURCE_HOST"
validate_remote_value "REMOTE_ROOT" "$REMOTE_ROOT"
validate_remote_value "REMOTE_ENV" "$REMOTE_ENV"
validate_remote_value "REMOTE_PYTHON" "$REMOTE_PYTHON"
validate_remote_value "SOURCE_OFFICIAL_ROOT" "$SOURCE_OFFICIAL_ROOT"
validate_remote_value "DATASETS" "$DATASETS"
validate_remote_value "THREADS" "$THREADS"
validate_remote_value "NICE_LEVEL" "$NICE_LEVEL"
validate_run_tag "$RUN_TAG"

ssh_base=(
  ssh
  -i "$REMOTE_KEY"
  -p "$REMOTE_PORT"
  "$REMOTE_HOST"
)

scp_base=(
  scp
  -i "$REMOTE_KEY"
  -P "$REMOTE_PORT"
)

source_scp_base=(
  scp
  -P "$SOURCE_PORT"
)

"${ssh_base[@]}" "mkdir -p '$REMOTE_JOB_DIR' '$REMOTE_ROOT/logs' '$REMOTE_ROOT/chrombpnet_refs'"

"${source_scp_base[@]}" \
  "$SOURCE_HOST:$SOURCE_OFFICIAL_ROOT/chrombpnet/helpers/make_gc_matched_negatives/get_gc_content.py" \
  "$SOURCE_HOST:$SOURCE_OFFICIAL_ROOT/chrombpnet/helpers/make_gc_matched_negatives/get_gc_matched_negatives.py" \
  "$SOURCE_HOST:$SOURCE_OFFICIAL_ROOT/chrombpnet/helpers/make_gc_matched_negatives/get_genomewide_gc_buckets/get_genomewide_gc_bins.py" \
  "$LOCAL_STAGE_DIR/"

"${scp_base[@]}" \
  "$REPO_ROOT/scripts/run_remote_chrombpnet_dataset_prep.sh" \
  "$LOCAL_STAGE_DIR/get_gc_content.py" \
  "$LOCAL_STAGE_DIR/get_gc_matched_negatives.py" \
  "$LOCAL_STAGE_DIR/get_genomewide_gc_bins.py" \
  "$REMOTE_HOST:$REMOTE_JOB_DIR/"

"${ssh_base[@]}" "chmod +x '$REMOTE_JOB_DIR/run_remote_chrombpnet_dataset_prep.sh'"

pid="$("${ssh_base[@]}" \
  "nohup bash '$REMOTE_JOB_DIR/run_remote_chrombpnet_dataset_prep.sh' \
    --root '$REMOTE_ROOT' \
    --env-dir '$REMOTE_ENV' \
    --python-bin '$REMOTE_PYTHON' \
    --gc-helper-dir '$REMOTE_JOB_DIR' \
    --datasets '$DATASETS' \
    --threads '$THREADS' \
    --nice-level '$NICE_LEVEL' \
    > '$REMOTE_LOG' 2>&1 < /dev/null & echo \$!")"

echo "started pid=$pid"
echo "remote_job_dir=$REMOTE_JOB_DIR"
echo "remote_log=$REMOTE_LOG"
echo "check: ssh -i $REMOTE_KEY -p $REMOTE_PORT $REMOTE_HOST 'tail -n 40 $REMOTE_LOG'"
