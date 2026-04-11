#!/usr/bin/env bash

set -euo pipefail

if [ -z "${BASH_VERSION:-}" ]; then
    echo "[error] Please run this script with bash." >&2
    exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QUEUE_CONFIG="${QUEUE_CONFIG:-${ROOT_DIR}/configs/queues/factor_ladder_unattended_20260411.yaml}"
STATE_DIR="${STATE_DIR:-${ROOT_DIR}/outputs/queue/factor_ladder_unattended_20260411}"
VENV_DIR="${VENV_DIR:-/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b}"

export PATH="${VENV_DIR}/bin:${PATH:-}"
export LD_LIBRARY_PATH="${VENV_DIR}/lib:${LD_LIBRARY_PATH:-}"
if [ -d "${ROOT_DIR}/src/transchrombp" ]; then
    PACKAGE_IMPORT_ROOT="${ROOT_DIR}/src"
elif [ -d "${ROOT_DIR}/transchrombp" ]; then
    PACKAGE_IMPORT_ROOT="${ROOT_DIR}"
else
    echo "[error] Could not locate transchrombp package import root from ${ROOT_DIR}" >&2
    exit 1
fi
export PYTHONPATH="${PACKAGE_IMPORT_ROOT}:${PYTHONPATH:-}"

exec python -m transchrombp.orchestration.factor_ladder_unattended_queue \
    --queue-config "${QUEUE_CONFIG}" \
    --state-dir "${STATE_DIR}" \
    "$@"
