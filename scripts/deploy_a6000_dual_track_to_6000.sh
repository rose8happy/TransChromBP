#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

REMOTE_HOST="${REMOTE_HOST:-zhoujiazhen@127.0.0.1}"
REMOTE_PORT="${REMOTE_PORT:-6000}"
REMOTE_ROOT="${REMOTE_ROOT:-/data1/zhoujiazhen/bylw_atac/TransChromBP}"
DRY_RUN="${DRY_RUN:-0}"
STAGE_PARENT="${REMOTE_ROOT}/scripts/.deploy_stage"
STAGE_ID="a6000_dual_track_$(date +%Y%m%d_%H%M%S)_$$"
REMOTE_STAGE_DIR="${STAGE_PARENT}/${STAGE_ID}"

REMOTE_SSH=(ssh -p "${REMOTE_PORT}" "${REMOTE_HOST}")
REMOTE_SCP=(scp -P "${REMOTE_PORT}")

LOCAL_MODEL_DIR="${REPO_ROOT}/vendor/transchrombp/transchrombp/models"
LOCAL_CONFIG_DIR="${REPO_ROOT}/vendor/transchrombp/transchrombp/configs/model"
LOCAL_RUNTIME_SCRIPT_DIR="${REPO_ROOT}/vendor/transchrombp/transchrombp/scripts"
LOCAL_ALPHA_DIR="${REPO_ROOT}/scripts/alphagenome_pilot"

LOCAL_FILES=(
    "${LOCAL_MODEL_DIR}/profile_decoder.py"
    "${LOCAL_MODEL_DIR}/transchrombp.py"
    "${LOCAL_MODEL_DIR}/__init__.py"
    "${LOCAL_CONFIG_DIR}/transchrombp_teacher_v2_center_pool_msdls_v2.yaml"
    "${LOCAL_RUNTIME_SCRIPT_DIR}/run_msdls_v2_gate.sh"
    "${LOCAL_ALPHA_DIR}/run_alphagenome_pilot.py"
    "${LOCAL_ALPHA_DIR}/merge_locus_totals.py"
    "${LOCAL_ALPHA_DIR}/build_matched_panel_v2.py"
    "${LOCAL_ALPHA_DIR}/regions_k562_tutorial_matched_panel_v2.csv"
)

STAGE_DESTS=(
    "src/transchrombp/models/profile_decoder.py"
    "src/transchrombp/models/transchrombp.py"
    "src/transchrombp/models/__init__.py"
    "configs/model/transchrombp_teacher_v2_center_pool_msdls_v2.yaml"
    "scripts/run_msdls_v2_gate.sh"
    "scripts/alphagenome_pilot/run_alphagenome_pilot.py"
    "scripts/alphagenome_pilot/merge_locus_totals.py"
    "scripts/alphagenome_pilot/build_matched_panel_v2.py"
    "scripts/alphagenome_pilot/regions_k562_tutorial_matched_panel_v2.csv"
)

cleanup_remote_stage() {
    if [ -n "${REMOTE_STAGE_DIR:-}" ]; then
        "${REMOTE_SSH[@]}" "rm -rf '${REMOTE_STAGE_DIR}'" >/dev/null 2>&1 || true
    fi
}

require_file() {
    local path="$1"
    if [ ! -f "${path}" ]; then
        echo "[error] Missing local file: ${path}" >&2
        exit 1
    fi
}

assert_remote_root() {
    "${REMOTE_SSH[@]}" "
        set -euo pipefail
        test -d '${REMOTE_ROOT}'
        test -d '${REMOTE_ROOT}/.git'
        test -d '${REMOTE_ROOT}/src/transchrombp'
        test -f '${REMOTE_ROOT}/src/transchrombp/models/transchrombp.py'
        test -d '${REMOTE_ROOT}/configs/model'
        test -d '${REMOTE_ROOT}/scripts/alphagenome_pilot'
    "
}

prepare_remote_stage() {
    "${REMOTE_SSH[@]}" "
        set -euo pipefail
        mkdir -p \
            '${REMOTE_STAGE_DIR}/src/transchrombp/models' \
            '${REMOTE_STAGE_DIR}/configs/model' \
            '${REMOTE_STAGE_DIR}/scripts/alphagenome_pilot'
    "
}

upload_to_stage() {
    local index
    local local_path
    local rel_path
    for index in "${!LOCAL_FILES[@]}"; do
        local_path="${LOCAL_FILES[index]}"
        rel_path="${STAGE_DESTS[index]}"
        "${REMOTE_SCP[@]}" "${local_path}" "${REMOTE_HOST}:${REMOTE_STAGE_DIR}/${rel_path}"
    done
}

verify_stage() {
    local quoted_paths=()
    local rel_path
    for rel_path in "${STAGE_DESTS[@]}"; do
        quoted_paths+=("'${REMOTE_STAGE_DIR}/${rel_path}'")
    done

    "${REMOTE_SSH[@]}" "
        set -euo pipefail
        for path in ${quoted_paths[*]}; do
            test -f \"\${path}\"
        done
    "
}

promote_from_stage() {
    local quoted_paths=()
    local rel_path
    local dry_run_flag="${DRY_RUN}"

    for rel_path in "${STAGE_DESTS[@]}"; do
        quoted_paths+=("'${rel_path}'")
    done

    "${REMOTE_SSH[@]}" "
        set -euo pipefail
        dry_run='${dry_run_flag}'
        for rel in ${quoted_paths[*]}; do
            staged='${REMOTE_STAGE_DIR}/'\${rel}
            final='${REMOTE_ROOT}/'\${rel}
            final_dir=\$(dirname \"\${final}\")
            final_tmp=\"\${final}.deploy_tmp\"
            test -f \"\${staged}\"
            if [ \"\${dry_run}\" = '1' ]; then
                continue
            fi
            mkdir -p \"\${final_dir}\"
            cp \"\${staged}\" \"\${final_tmp}\"
            cmp -s \"\${staged}\" \"\${final_tmp}\"
            mv \"\${final_tmp}\" \"\${final}\"
            cmp -s \"\${staged}\" \"\${final}\"
        done
    "
}

main() {
    local path
    local rel_path

    trap cleanup_remote_stage EXIT

    for path in "${LOCAL_FILES[@]}"; do
        require_file "${path}"
    done

    assert_remote_root
    prepare_remote_stage
    upload_to_stage
    verify_stage
    promote_from_stage

    if [ "${DRY_RUN}" = "1" ]; then
        echo "[ok] dry-run staged and verified deploy payload at ${REMOTE_HOST}:${REMOTE_STAGE_DIR}"
        return
    fi

    for rel_path in "${STAGE_DESTS[@]}"; do
        "${REMOTE_SSH[@]}" "
            set -euo pipefail
            test -f '${REMOTE_ROOT}/${rel_path}'
            cmp -s '${REMOTE_STAGE_DIR}/${rel_path}' '${REMOTE_ROOT}/${rel_path}'
        "
    done

    echo "[ok] deployed dual-track assets to ${REMOTE_HOST}:${REMOTE_ROOT}"
}

main "$@"
