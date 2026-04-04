#!/usr/bin/env bash
set -euo pipefail

REMOTE_HOST="zhoujiazhen@127.0.0.1"
REMOTE_PORT="6000"
ROOT="/data1/zhoujiazhen/bylw_atac"
LOCAL_REPO="$(cd "$(dirname "$0")/.." && pwd)"
REMOTE_REPO="$ROOT/chromBPNet"
TRANSCHROMBP_ROOT="$ROOT/TransChromBP"
LOG_ROOT="$ROOT/logs"
NT_ENV="$ROOT/.mamba/envs/nucleotide-transformer-py311"
MODEL_DIR="$ROOT/foundation_models/nucleotide_transformer/nt-v2-500m-multi-species"
RUN_TS="$(date +%Y%m%d_%H%M%S)"
RUN_NAME="nt_v2_repair_smoke_probe_${RUN_TS}"
REMOTE_SCRIPT="/tmp/${RUN_NAME}.sh"
REMOTE_LOG="${LOG_ROOT}/${RUN_NAME}.log"
PROBE_DIR="$TRANSCHROMBP_ROOT/outputs/nt_v2_probe/${RUN_TS}"

scp -P "${REMOTE_PORT}" \
  "${LOCAL_REPO}/scripts/nt_v2_smoke.py" \
  "${LOCAL_REPO}/scripts/nt_v2_probe.py" \
  "${REMOTE_HOST}:${REMOTE_REPO}/scripts/"

ssh -p "${REMOTE_PORT}" "${REMOTE_HOST}" "cat > '${REMOTE_SCRIPT}' <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

ROOT='${ROOT}'
REMOTE_REPO='${REMOTE_REPO}'
TRANSCHROMBP_ROOT='${TRANSCHROMBP_ROOT}'
LOG_ROOT='${LOG_ROOT}'
NT_ENV='${NT_ENV}'
MODEL_DIR='${MODEL_DIR}'
RUN_NAME='${RUN_NAME}'
PROBE_DIR='${PROBE_DIR}'
STAGING_MODEL_DIR="\$ROOT/foundation_models/nucleotide_transformer/nt-v2-500m-multi-species_repair_\${RUN_NAME}"
MODEL_ID='InstaDeepAI/nucleotide-transformer-v2-500m-multi-species'

mkdir -p \"\$LOG_ROOT\" \"\$PROBE_DIR\"

MICROMAMBA_BIN=\"\$(command -v micromamba)\"
if [[ -z \"\$MICROMAMBA_BIN\" ]]; then
  echo 'micromamba not found' >&2
  exit 1
fi
export MAMBA_ROOT_PREFIX=\"\$ROOT/.mamba\"
export HF_ENDPOINT='https://hf-mirror.com'
export HF_HUB_DISABLE_TELEMETRY=1
export HF_HUB_DOWNLOAD_TIMEOUT=30

if [[ ! -x \"\$NT_ENV/bin/python\" ]]; then
  \"\$MICROMAMBA_BIN\" create -y -p \"\$NT_ENV\" python=3.11 pip
fi

\"\$MICROMAMBA_BIN\" run -p \"\$NT_ENV\" pip install --upgrade pip setuptools wheel
\"\$MICROMAMBA_BIN\" run -p \"\$NT_ENV\" pip install --upgrade 'numpy<2.3' scipy scikit-learn pyfaidx pyBigWig safetensors huggingface_hub 'transformers>=4.57.0,<5'
\"\$MICROMAMBA_BIN\" run -p \"\$NT_ENV\" pip install --upgrade --index-url https://download.pytorch.org/whl/cu126 torch

mkdir -p \"\$STAGING_MODEL_DIR\"
if [[ -d \"\$MODEL_DIR\" ]]; then
  rsync -a \"\$MODEL_DIR/\" \"\$STAGING_MODEL_DIR/\"
fi
if ! \"\$MICROMAMBA_BIN\" run -p \"\$NT_ENV\" hf download \"\$MODEL_ID\" --local-dir \"\$STAGING_MODEL_DIR\"; then
  echo '[warn] hf download failed, continuing with staged local copy'
fi

SMOKE_JSON=\"\$PROBE_DIR/smoke_summary.json\"
EXTRACT_DIR=\"\$PROBE_DIR/extract\"
ANALYZE_DIR=\"\$PROBE_DIR/analyze\"

PEAK_CENTER=\"\$(awk 'BEGIN{FS=\"\t\"} \$1 !~ /^#/ {start=\$2+0; end=\$3+0; summit=(NF>9 && \$10 != \"\") ? int(\$10+0) : int((start+end)/2)-start; print start + summit; exit}' \"\$TRANSCHROMBP_ROOT/outputs/preprocessing/tutorial_canonical_v1/step4_filtering/filtered.peaks.bed\")\"
PEAK_CHROM=\"\$(awk 'BEGIN{FS=\"\t\"} \$1 !~ /^#/ {print \$1; exit}' \"\$TRANSCHROMBP_ROOT/outputs/preprocessing/tutorial_canonical_v1/step4_filtering/filtered.peaks.bed\")\"

\"\$MICROMAMBA_BIN\" run -p \"\$NT_ENV\" python \"\$REMOTE_REPO/scripts/nt_v2_smoke.py\" \
  --model-dir \"\$STAGING_MODEL_DIR\" \
  --genome-fasta \"\$ROOT/hg38.fa\" \
  --chrom \"\$PEAK_CHROM\" \
  --center \"\$PEAK_CENTER\" \
  --device cuda \
  --output-json \"\$SMOKE_JSON\"

\"\$MICROMAMBA_BIN\" run -p \"\$NT_ENV\" python \"\$REMOTE_REPO/scripts/nt_v2_probe.py\" extract \
  --model-dir \"\$STAGING_MODEL_DIR\" \
  --genome-fasta \"\$ROOT/hg38.fa\" \
  --bigwig \"\$TRANSCHROMBP_ROOT/outputs/preprocessing/tutorial_canonical_v1/step2_bigwig/merged_unstranded.bw\" \
  --peaks-bed \"\$TRANSCHROMBP_ROOT/outputs/preprocessing/tutorial_canonical_v1/step4_filtering/filtered.peaks.bed\" \
  --nonpeaks-bed \"\$TRANSCHROMBP_ROOT/outputs/preprocessing/tutorial_canonical_v1/step4_filtering/filtered.nonpeaks.bed\" \
  --folds-json \"\$ROOT/chrombpnet_tutorial/data/folds.json\" \
  --split valid \
  --sample-size-per-class 500 \
  --layers auto_quartiles \
  --feature-types global_mean,bins4_mean \
  --batch-size 8 \
  --device cuda \
  --output-dir \"\$EXTRACT_DIR\"

python_bin='${ROOT}/.mamba/envs/transchrombp/bin/python'
\"\$python_bin\" \"\$REMOTE_REPO/scripts/nt_v2_probe.py\" analyze \
  --probe-dir \"\$EXTRACT_DIR\" \
  --baseline-checkpoint \"\$TRANSCHROMBP_ROOT/outputs/checkpoints/v2fix_20260320_cpool_s42/best.pt\" \
  --baseline-model-config \"\$TRANSCHROMBP_ROOT/configs/model/transchrombp_teacher_v2_center_pool.yaml\" \
  --transchrombp-root \"\$TRANSCHROMBP_ROOT\" \
  --batch-size 64 \
  --device cuda \
  --output-dir \"\$ANALYZE_DIR\"

echo \"[done] \$(date '+%F %T') \$RUN_NAME\"
EOF
chmod +x '${REMOTE_SCRIPT}'
nohup bash '${REMOTE_SCRIPT}' > '${REMOTE_LOG}' 2>&1 & echo \$!"

echo "run_name=${RUN_NAME}"
echo "remote_log=${REMOTE_LOG}"
echo "probe_dir=${PROBE_DIR}"
