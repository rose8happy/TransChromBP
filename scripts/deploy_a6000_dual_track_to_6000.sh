#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

REMOTE_HOST="${REMOTE_HOST:-zhoujiazhen@127.0.0.1}"
REMOTE_PORT="${REMOTE_PORT:-6000}"
REMOTE_ROOT="${REMOTE_ROOT:-/data1/zhoujiazhen/bylw_atac/TransChromBP}"

require_file() {
    local path="$1"
    if [ ! -f "${path}" ]; then
        echo "[error] Missing local file: ${path}" >&2
        exit 1
    fi
}

deploy_files() {
    local destination="$1"
    shift
    scp -P "${REMOTE_PORT}" "$@" "${REMOTE_HOST}:${destination}"
}

LOCAL_MODEL_DIR="${REPO_ROOT}/vendor/transchrombp/transchrombp/models"
LOCAL_CONFIG_DIR="${REPO_ROOT}/vendor/transchrombp/transchrombp/configs/model"
LOCAL_RUNTIME_SCRIPT_DIR="${REPO_ROOT}/vendor/transchrombp/transchrombp/scripts"
LOCAL_ALPHA_DIR="${REPO_ROOT}/scripts/alphagenome_pilot"

MODEL_FILES=(
    "${LOCAL_MODEL_DIR}/profile_decoder.py"
    "${LOCAL_MODEL_DIR}/transchrombp.py"
    "${LOCAL_MODEL_DIR}/__init__.py"
)

CONFIG_FILES=(
    "${LOCAL_CONFIG_DIR}/transchrombp_teacher_v2_center_pool_msdls_v2.yaml"
)

RUNTIME_SCRIPT_FILES=(
    "${LOCAL_RUNTIME_SCRIPT_DIR}/run_msdls_v2_gate.sh"
)

ALPHAGENOME_FILES=(
    "${LOCAL_ALPHA_DIR}/run_alphagenome_pilot.py"
    "${LOCAL_ALPHA_DIR}/merge_locus_totals.py"
    "${LOCAL_ALPHA_DIR}/build_matched_panel_v2.py"
    "${LOCAL_ALPHA_DIR}/regions_k562_tutorial_matched_panel_v2.csv"
)

for path in \
    "${MODEL_FILES[@]}" \
    "${CONFIG_FILES[@]}" \
    "${RUNTIME_SCRIPT_FILES[@]}" \
    "${ALPHAGENOME_FILES[@]}"; do
    require_file "${path}"
done

ssh -p "${REMOTE_PORT}" "${REMOTE_HOST}" \
    "mkdir -p \
        '${REMOTE_ROOT}/src/transchrombp/models' \
        '${REMOTE_ROOT}/configs/model' \
        '${REMOTE_ROOT}/scripts' \
        '${REMOTE_ROOT}/scripts/alphagenome_pilot'"

deploy_files "${REMOTE_ROOT}/src/transchrombp/models/" "${MODEL_FILES[@]}"
deploy_files "${REMOTE_ROOT}/configs/model/" "${CONFIG_FILES[@]}"
deploy_files "${REMOTE_ROOT}/scripts/" "${RUNTIME_SCRIPT_FILES[@]}"
deploy_files "${REMOTE_ROOT}/scripts/alphagenome_pilot/" "${ALPHAGENOME_FILES[@]}"

echo "[ok] deployed dual-track assets to ${REMOTE_HOST}:${REMOTE_ROOT}"
