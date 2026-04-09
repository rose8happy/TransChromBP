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
OFFICIAL_HOST="${OFFICIAL_HOST:-zhoujiazhen@127.0.0.1}"
OFFICIAL_PORT="${OFFICIAL_PORT:-6000}"
OFFICIAL_ROOT="${CHROMBPNET_OFFICIAL_ROOT:-/data1/zhoujiazhen/bylw_atac/chrombpnet_official}"
DATASETS="${DATASETS:-GM12878,K562}"
THREADS="${THREADS:-4}"
NICE_LEVEL="${NICE_LEVEL:-10}"
RUN_TAG="${RUN_TAG:-chrombpnet_dataset_prep_6002_$(date +%Y%m%d_%H%M%S)}"

REMOTE_JOB_DIR="$REMOTE_ROOT/.codex_jobs/chrombpnet_dataset_prep/$RUN_TAG"
REMOTE_LOG="$REMOTE_ROOT/logs/${RUN_TAG}.log"
LOCAL_STAGE="$(mktemp -d)"
trap 'rm -rf "$LOCAL_STAGE"' EXIT

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

"${ssh_base[@]}" "mkdir -p '$REMOTE_JOB_DIR' '$REMOTE_ROOT/logs' '$REMOTE_ROOT/chrombpnet_refs'"

scp -P "$OFFICIAL_PORT" \
  "$OFFICIAL_HOST:$OFFICIAL_ROOT/chrombpnet/helpers/make_gc_matched_negatives/get_gc_content.py" \
  "$OFFICIAL_HOST:$OFFICIAL_ROOT/chrombpnet/helpers/make_gc_matched_negatives/get_gc_matched_negatives.py" \
  "$OFFICIAL_HOST:$OFFICIAL_ROOT/chrombpnet/helpers/make_gc_matched_negatives/get_genomewide_gc_buckets/get_genomewide_gc_bins.py" \
  "$LOCAL_STAGE/"

"${scp_base[@]}" \
  "$REPO_ROOT/scripts/run_remote_chrombpnet_dataset_prep.sh" \
  "$LOCAL_STAGE/get_gc_content.py" \
  "$LOCAL_STAGE/get_gc_matched_negatives.py" \
  "$LOCAL_STAGE/get_genomewide_gc_bins.py" \
  "$REMOTE_HOST:$REMOTE_JOB_DIR/"

"${ssh_base[@]}" "chmod +x '$REMOTE_JOB_DIR/run_remote_chrombpnet_dataset_prep.sh'"

pid="$("${ssh_base[@]}" \
  "nohup bash '$REMOTE_JOB_DIR/run_remote_chrombpnet_dataset_prep.sh' \
    --root '$REMOTE_ROOT' \
    --env-dir '$REMOTE_ENV' \
    --python-bin '$REMOTE_PYTHON' \
    --datasets '$DATASETS' \
    --threads '$THREADS' \
    --nice-level '$NICE_LEVEL' \
    --gc-helper-dir '$REMOTE_JOB_DIR' \
    > '$REMOTE_LOG' 2>&1 < /dev/null & echo \$!")"

echo "started pid=$pid"
echo "remote_job_dir=$REMOTE_JOB_DIR"
echo "remote_log=$REMOTE_LOG"
echo "check: ssh -i $REMOTE_KEY -p $REMOTE_PORT $REMOTE_HOST 'tail -n 40 $REMOTE_LOG'"
