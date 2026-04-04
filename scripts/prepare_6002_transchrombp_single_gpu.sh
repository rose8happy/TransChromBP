#!/usr/bin/env bash
set -euo pipefail

REMOTE_HOST="zhengwei@127.0.0.1"
REMOTE_PORT="6002"
REMOTE_KEY="/home/zhengwei/.ssh/codex_6002_ed25519"
SOURCE_HOST="zhoujiazhen@127.0.0.1"
SOURCE_PORT="6000"

REMOTE_ROOT="${REMOTE_ROOT:-/home/zhengwei/bylw_atac}"
REMOTE_TRANS_ROOT="$REMOTE_ROOT/TransChromBP"
REMOTE_ENV="$REMOTE_ROOT/.mamba/envs/transchrombp"
SYNC_BIGWIG="${SYNC_BIGWIG:-0}"

ssh_6002=(
  ssh
  -i "$REMOTE_KEY"
  -p "$REMOTE_PORT"
  "$REMOTE_HOST"
)

ssh_6000=(
  ssh
  -p "$SOURCE_PORT"
  "$SOURCE_HOST"
)

echo "[1/4] ensure remote directories"
"${ssh_6002[@]}" "mkdir -p \
  '$REMOTE_TRANS_ROOT/outputs/preprocessing/tutorial_canonical_v1/step2_bigwig' \
  '$REMOTE_TRANS_ROOT/outputs/preprocessing/tutorial_canonical_v1/step4_filtering' \
  '$REMOTE_TRANS_ROOT/configs/data' \
  '$REMOTE_TRANS_ROOT/configs/train' \
  '$REMOTE_TRANS_ROOT/scripts' \
  '$REMOTE_TRANS_ROOT/src' \
  '$REMOTE_ROOT/logs'"

echo "[2/4] sync tutorial preprocessing small files and required source tree from 6000"
small_sync_paths=(
  "TransChromBP/outputs/preprocessing/tutorial_canonical_v1/step4_filtering/filtered.peaks.bed"
  "TransChromBP/outputs/preprocessing/tutorial_canonical_v1/step4_filtering/filtered.nonpeaks.bed"
  "TransChromBP/src/transchrombp"
)

if [[ "$SYNC_BIGWIG" == "1" ]]; then
  small_sync_paths+=(
    "TransChromBP/outputs/preprocessing/tutorial_canonical_v1/step2_bigwig/merged_unstranded.bw"
  )
fi

"${ssh_6000[@]}" "cd /data1/zhoujiazhen/bylw_atac && tar -cf - ${small_sync_paths[*]}" \
| "${ssh_6002[@]}" "cd '$REMOTE_ROOT' && tar -xf -"

echo "[3/4] write 6002-specific configs"
"${ssh_6002[@]}" "python3 - <<'PY'
from pathlib import Path
import yaml

root = Path('/home/zhengwei/bylw_atac/TransChromBP')
data_cfg = root / 'configs/data/data_tutorial_canonical_v1_6002.yaml'
bias_cfg = root / 'configs/train/train_bias_tutorial_teacher_v2_6002_single.yaml'
main_cfg = root / 'configs/train/train_tutorial_teacher_v2_main_6002_single.yaml'

with (root / 'configs/data/data_tutorial_canonical_v1.yaml').open() as f:
    data = yaml.safe_load(f)
with (root / 'configs/train/train_bias_tutorial_teacher_v2.yaml').open() as f:
    bias = yaml.safe_load(f)
with (root / 'configs/train/train_tutorial_teacher_v2_main.yaml').open() as f:
    main = yaml.safe_load(f)

prefix_old = '/data1/zhoujiazhen/bylw_atac'
prefix_new = '/home/zhengwei/bylw_atac'

def rewrite(obj):
    if isinstance(obj, dict):
        return {k: rewrite(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [rewrite(v) for v in obj]
    if isinstance(obj, str):
        return obj.replace(prefix_old, prefix_new)
    return obj

data = rewrite(data)
bias = rewrite(bias)
main = rewrite(main)

for cfg in (bias, main):
    cfg.setdefault('trainer', {})['grad_accum_steps'] = 2
    cfg.setdefault('data', {})['config_path'] = 'configs/data/data_tutorial_canonical_v1_6002.yaml'
    cfg.setdefault('logging', {})['output_dir'] = '/home/zhengwei/bylw_atac/TransChromBP/outputs'

data_cfg.write_text(yaml.safe_dump(data, sort_keys=False), encoding='utf-8')
bias_cfg.write_text(yaml.safe_dump(bias, sort_keys=False), encoding='utf-8')
main_cfg.write_text(yaml.safe_dump(main, sort_keys=False), encoding='utf-8')
PY"

echo "[4/4] write single-GPU launcher"
"${ssh_6002[@]}" "cat > '$REMOTE_TRANS_ROOT/scripts/launch_teacher_v2_6002_single.sh' <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=\"/home/zhengwei/bylw_atac/TransChromBP\"
ENV_DIR=\"/home/zhengwei/bylw_atac/.mamba/envs/transchrombp\"
LOG_DIR=\"/home/zhengwei/bylw_atac/logs\"

PIPELINE_NAME=\"\${PIPELINE_NAME:-teacher_v2_6002_single_\$(date +%Y%m%d_%H%M%S)}\"
MASTER_PORT_BASE=\"\${MASTER_PORT_BASE:-29561}\"
BIAS_BATCH_SIZE_PER_GPU=\"\${BIAS_BATCH_SIZE_PER_GPU:-16}\"
MAIN_BATCH_SIZE_PER_GPU=\"\${MAIN_BATCH_SIZE_PER_GPU:-16}\"
BIAS_NUM_WORKERS=\"\${BIAS_NUM_WORKERS:-2}\"
MAIN_NUM_WORKERS=\"\${MAIN_NUM_WORKERS:-2}\"

mkdir -p \"\$LOG_DIR\"
export PATH=\"\$ENV_DIR/bin:\$PATH\"
export LD_LIBRARY_PATH=\"\$ENV_DIR/lib:\${LD_LIBRARY_PATH:-}\"
export PYTHONPATH=\"\$ROOT_DIR/src:\${PYTHONPATH:-}\"

cd \"\$ROOT_DIR\"

NPROC_PER_NODE=1 \
MASTER_PORT_BASE=\"\$MASTER_PORT_BASE\" \
PIPELINE_NAME=\"\$PIPELINE_NAME\" \
PIPELINE_VARIANTS=\"\${PIPELINE_VARIANTS:-learnable}\" \
BIAS_BATCH_SIZE_PER_GPU=\"\$BIAS_BATCH_SIZE_PER_GPU\" \
MAIN_BATCH_SIZE_PER_GPU=\"\$MAIN_BATCH_SIZE_PER_GPU\" \
BIAS_NUM_WORKERS=\"\$BIAS_NUM_WORKERS\" \
MAIN_NUM_WORKERS=\"\$MAIN_NUM_WORKERS\" \
/usr/bin/bash \"\$ROOT_DIR/scripts/run_tutorial_bias_to_main_pipeline.sh\" \
  \"\$ROOT_DIR/configs/model/transchrombp_teacher_v2.yaml\" \
  \"\$ROOT_DIR/configs/train/train_bias_tutorial_teacher_v2_6002_single.yaml\" \
  \"\$ROOT_DIR/configs/train/train_tutorial_teacher_v2_main_6002_single.yaml\"
EOF
chmod +x '$REMOTE_TRANS_ROOT/scripts/launch_teacher_v2_6002_single.sh'"

echo "[verify] show generated files"
"${ssh_6002[@]}" "sed -n '1,80p' '$REMOTE_TRANS_ROOT/configs/data/data_tutorial_canonical_v1_6002.yaml'"
echo "---"
"${ssh_6002[@]}" "sed -n '1,120p' '$REMOTE_TRANS_ROOT/scripts/launch_teacher_v2_6002_single.sh'"
