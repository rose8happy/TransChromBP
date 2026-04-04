#!/usr/bin/env bash
set -euo pipefail

ROOT=/data1/zhoujiazhen/bylw_atac
MODEL_ID=InstaDeepAI/nucleotide-transformer-v2-500m-multi-species
MODEL_ROOT="$ROOT/foundation_models/nucleotide_transformer"
MODEL_DIR="$MODEL_ROOT/nt-v2-500m-multi-species"
LOG_ROOT="$ROOT/foundation_models/logs"
HF_HOME="$ROOT/.cache/huggingface"
GENOS_ENV="$ROOT/.venvs/genos-1.2b"

mkdir -p "$MODEL_ROOT" "$MODEL_DIR" "$LOG_ROOT" "$HF_HOME"

if [[ ! -x "$GENOS_ENV/bin/python" ]]; then
	echo "missing required env: $GENOS_ENV" >&2
	exit 1
fi

source "$GENOS_ENV/bin/activate"

export HF_ENDPOINT="${HF_ENDPOINT:-https://hf-mirror.com}"
export HF_HUB_DOWNLOAD_TIMEOUT="${HF_HUB_DOWNLOAD_TIMEOUT:-30}"
export HF_HUB_DISABLE_TELEMETRY=1
export HF_HOME
export HUGGINGFACE_HUB_CACHE="$HF_HOME/hub"
export TRANSFORMERS_CACHE="$HF_HOME/transformers"

echo "[$(date '+%F %T')] starting nucleotide transformer download"
echo "model_id=$MODEL_ID"
echo "model_dir=$MODEL_DIR"
echo "hf_endpoint=$HF_ENDPOINT"
echo "python=$(python -V 2>&1)"
echo "huggingface_hub=$(python -c 'import huggingface_hub; print(huggingface_hub.__version__)' 2>&1)"

hf download "$MODEL_ID" --local-dir "$MODEL_DIR"

echo "[$(date '+%F %T')] download finished"
du -sh "$MODEL_DIR"
find "$MODEL_DIR" -maxdepth 1 -type f | sort
