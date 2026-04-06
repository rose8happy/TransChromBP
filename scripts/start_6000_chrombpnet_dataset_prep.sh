#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

REMOTE_HOST="${REMOTE_HOST:-zhoujiazhen@127.0.0.1}"
REMOTE_PORT="${REMOTE_PORT:-6000}"
REMOTE_ROOT="${REMOTE_ROOT:-/data1/zhoujiazhen/bylw_atac}"
REMOTE_ENV="${REMOTE_ENV:-/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet}"
REMOTE_PYTHON="${REMOTE_PYTHON:-$REMOTE_ENV/bin/python}"
OFFICIAL_ROOT="${CHROMBPNET_OFFICIAL_ROOT:-/data1/zhoujiazhen/bylw_atac/chrombpnet_official}"
DATASETS="${DATASETS:-GM12878,K562}"
THREADS="${THREADS:-4}"
NICE_LEVEL="${NICE_LEVEL:-10}"
RUN_TAG="${RUN_TAG:-chrombpnet_dataset_prep_6000_$(date +%Y%m%d_%H%M%S)}"

REMOTE_JOB_DIR="$REMOTE_ROOT/.codex_jobs/chrombpnet_dataset_prep/$RUN_TAG"
REMOTE_LOG="$REMOTE_ROOT/logs/${RUN_TAG}.log"

ssh_base=(
  ssh
  -p "$REMOTE_PORT"
  "$REMOTE_HOST"
)

scp_base=(
  scp
  -P "$REMOTE_PORT"
)

"${ssh_base[@]}" "mkdir -p '$REMOTE_JOB_DIR' '$REMOTE_ROOT/logs' '$REMOTE_ROOT/chrombpnet_refs'"

"${scp_base[@]}" \
  "$REPO_ROOT/scripts/run_remote_chrombpnet_dataset_prep.sh" \
  "$REMOTE_HOST:$REMOTE_JOB_DIR/"

"${ssh_base[@]}" "chmod +x '$REMOTE_JOB_DIR/run_remote_chrombpnet_dataset_prep.sh'"

pid="$("${ssh_base[@]}" \
  "nohup bash '$REMOTE_JOB_DIR/run_remote_chrombpnet_dataset_prep.sh' \
    --root '$REMOTE_ROOT' \
    --env-dir '$REMOTE_ENV' \
    --python-bin '$REMOTE_PYTHON' \
    --official-root '$OFFICIAL_ROOT' \
    --datasets '$DATASETS' \
    --threads '$THREADS' \
    --nice-level '$NICE_LEVEL' \
    > '$REMOTE_LOG' 2>&1 < /dev/null & echo \$!")"

echo "started pid=$pid"
echo "remote_job_dir=$REMOTE_JOB_DIR"
echo "remote_log=$REMOTE_LOG"
echo "check: ssh -p $REMOTE_PORT $REMOTE_HOST 'tail -n 40 $REMOTE_LOG'"
