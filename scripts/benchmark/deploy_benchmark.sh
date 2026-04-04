#!/usr/bin/env bash
# ==============================================================================
# 将 benchmark 文件部署到 6000 服务器
# 用法: bash scripts/benchmark/deploy_benchmark.sh
# ==============================================================================
set -euo pipefail

LOCAL_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
REMOTE_HOST="zhoujiazhen@127.0.0.1"
REMOTE_PORT=6000
REMOTE_TRANSCHROMBP="/data1/zhoujiazhen/bylw_atac/TransChromBP"

echo "=== Deploying benchmark files to 6000 server ==="
echo "Local root: $LOCAL_ROOT"
echo "Remote: $REMOTE_HOST:$REMOTE_PORT:$REMOTE_TRANSCHROMBP"

# 1. 部署 benchmark 配置
echo ""
echo "[1/4] Deploying benchmark config..."
scp -P "$REMOTE_PORT" \
    "$LOCAL_ROOT/.codex_remote_edit/TransChromBP/configs/train/train_benchmark_base.yaml" \
    "$REMOTE_HOST:$REMOTE_TRANSCHROMBP/configs/train/train_benchmark_base.yaml"

# 2. 部署 train_ddp_bench.py
echo "[2/4] Deploying train_ddp_bench.py..."
scp -P "$REMOTE_PORT" \
    "$LOCAL_ROOT/.codex_remote_edit/TransChromBP/src/transchrombp/training/train_ddp_bench.py" \
    "$REMOTE_HOST:$REMOTE_TRANSCHROMBP/src/transchrombp/training/train_ddp_bench.py"

# 3. 部署运行脚本
echo "[3/4] Deploying benchmark scripts..."
ssh -p "$REMOTE_PORT" "$REMOTE_HOST" "mkdir -p $REMOTE_TRANSCHROMBP/scripts/benchmark"
for script in run_training_benchmark.sh run_profiler.sh parse_benchmark_logs.py; do
    scp -P "$REMOTE_PORT" \
        "$LOCAL_ROOT/scripts/benchmark/$script" \
        "$REMOTE_HOST:$REMOTE_TRANSCHROMBP/scripts/benchmark/$script"
done

# 4. 远端验证
echo "[4/4] Verifying on remote..."
ssh -p "$REMOTE_PORT" "$REMOTE_HOST" bash -c "'
    echo \"--- Deployed files ---\"
    ls -la $REMOTE_TRANSCHROMBP/configs/train/train_benchmark_base.yaml
    ls -la $REMOTE_TRANSCHROMBP/src/transchrombp/training/train_ddp_bench.py
    ls -la $REMOTE_TRANSCHROMBP/scripts/benchmark/

    echo \"\"
    echo \"--- Syntax check train_ddp_bench.py ---\"
    export PYTHONPATH=$REMOTE_TRANSCHROMBP/src:\${PYTHONPATH:-}
    export PATH=/data1/zhoujiazhen/bylw_atac/.mamba/envs/transchrombp/bin:\$PATH
    python3 -m py_compile $REMOTE_TRANSCHROMBP/src/transchrombp/training/train_ddp_bench.py && echo \"OK\" || echo \"FAILED\"

    echo \"\"
    echo \"--- Syntax check benchmark config ---\"
    python3 -c \"import yaml; yaml.safe_load(open(\\\"$REMOTE_TRANSCHROMBP/configs/train/train_benchmark_base.yaml\\\"))\" && echo \"OK\" || echo \"FAILED\"

    echo \"\"
    echo \"--- Bash syntax check ---\"
    bash -n $REMOTE_TRANSCHROMBP/scripts/benchmark/run_training_benchmark.sh && echo \"OK\" || echo \"FAILED\"
'"

echo ""
echo "=== Deployment complete ==="
echo ""
echo "To run benchmark on 6000:"
echo "  ssh -p $REMOTE_PORT $REMOTE_HOST"
echo "  cd $REMOTE_TRANSCHROMBP"
echo "  nohup bash scripts/benchmark/run_training_benchmark.sh > outputs/logs/benchmark/full_run.log 2>&1 &"
echo ""
echo "To check progress:"
echo "  tail -f $REMOTE_TRANSCHROMBP/outputs/logs/benchmark/full_run.log"
