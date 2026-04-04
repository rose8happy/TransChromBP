#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_DIR="${1:-$ROOT_DIR/.venv-report}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
BOOTSTRAP_PYTHON="${BOOTSTRAP_PYTHON:-$ROOT_DIR/.venv/bin/python}"

echo "[info] root: $ROOT_DIR"
echo "[info] env:  $ENV_DIR"
echo "[info] python: $PYTHON_BIN"

if "$PYTHON_BIN" -m venv "$ENV_DIR" 2>/tmp/chrombpnet_report_env.err; then
  echo "[info] created env with stdlib venv"
else
  echo "[warn] stdlib venv failed, falling back to virtualenv"
  cat /tmp/chrombpnet_report_env.err
  if [ -x "$BOOTSTRAP_PYTHON" ]; then
    echo "[info] bootstrap python: $BOOTSTRAP_PYTHON"
    "$BOOTSTRAP_PYTHON" -m ensurepip --upgrade >/dev/null 2>&1 || true
    "$BOOTSTRAP_PYTHON" -m pip install virtualenv
    "$BOOTSTRAP_PYTHON" -m virtualenv "$ENV_DIR"
  else
    "$PYTHON_BIN" -m pip install --user --break-system-packages virtualenv
    "$PYTHON_BIN" -m virtualenv "$ENV_DIR"
  fi
fi

source "$ENV_DIR/bin/activate"

python -m pip install --upgrade pip setuptools wheel
python -m pip install -r "$ROOT_DIR/requirements-report.txt"

if [ "${INSTALL_NOTEBOOK:-0}" = "1" ]; then
  python -m pip install -r "$ROOT_DIR/requirements-report-notebook.txt"
  if python -m ipykernel install --user --name chrombpnet-report --display-name "Python (chrombpnet-report)" >/dev/null 2>&1; then
    echo "[info] registered Jupyter kernel: chrombpnet-report"
  else
    echo "[warn] failed to register Jupyter kernel; env is still usable"
  fi
fi

echo "[done] report env ready"
echo "[hint] activate with: source \"$ENV_DIR/bin/activate\""
echo "[hint] optional notebook extras: INSTALL_NOTEBOOK=1 bash scripts/setup_report_env.sh"
