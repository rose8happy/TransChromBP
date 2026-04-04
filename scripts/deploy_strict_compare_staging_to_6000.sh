#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

REMOTE_USER="zhoujiazhen"
REMOTE_HOST="127.0.0.1"
REMOTE_PORT="6000"
REMOTE_ROOT="/data1/zhoujiazhen/bylw_atac/TransChromBP"
SSH_BASE=(ssh -p "${REMOTE_PORT}" "${REMOTE_USER}@${REMOTE_HOST}")
SCP_BASE=(scp -P "${REMOTE_PORT}")

LOCAL_TRAIN_CFG_DIR="${REPO_ROOT}/vendor/transchrombp/transchrombp/configs/train"
LOCAL_SCRIPT_DIR="${REPO_ROOT}/vendor/transchrombp/scripts"

"${SSH_BASE[@]}" "mkdir -p '${REMOTE_ROOT}/configs/train' '${REMOTE_ROOT}/scripts'"

"${SCP_BASE[@]}" \
  "${LOCAL_TRAIN_CFG_DIR}/train_tutorial_corrected_b_strict_compare_6000.yaml" \
  "${LOCAL_TRAIN_CFG_DIR}/train_gm12878_corrected_b_strict_compare_6000.yaml" \
  "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_ROOT}/configs/train/"

"${SCP_BASE[@]}" \
  "${LOCAL_SCRIPT_DIR}/select_best_epoch.py" \
  "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_ROOT}/scripts/"

echo "Deployed strict-compare staged files to ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_ROOT}"
