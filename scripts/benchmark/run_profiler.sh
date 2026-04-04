#!/usr/bin/env bash
# ==============================================================================
# TransChromBP PyTorch Profiler 独立运行脚本
# 用途：Profile 20 步训练，找出 CUDA 时间瓶颈
# ==============================================================================
set -euo pipefail

export PATH=/data1/zhoujiazhen/bylw_atac/.mamba/envs/transchrombp/bin:$PATH
export LD_LIBRARY_PATH=/data1/zhoujiazhen/bylw_atac/.mamba/envs/transchrombp/lib:${LD_LIBRARY_PATH:-}
export PYTHONPATH=/data1/zhoujiazhen/bylw_atac/TransChromBP/src:${PYTHONPATH:-}

ROOT_DIR="/data1/zhoujiazhen/bylw_atac/TransChromBP"
MODEL_CONFIG="$ROOT_DIR/configs/model/transchrombp_base.yaml"
TRAIN_CONFIG="$ROOT_DIR/configs/train/train_benchmark_base.yaml"
OUTPUT_DIR="$ROOT_DIR/outputs"
LOG_DIR="$OUTPUT_DIR/logs/benchmark"
PROFILER_DIR="$OUTPUT_DIR/profiler_traces"
MASTER_PORT="${MASTER_PORT:-29511}"
NPROC="${NPROC_PER_NODE:-2}"
BS="${BATCH_SIZE:-16}"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

cd "$ROOT_DIR"
mkdir -p "$LOG_DIR" "$PROFILER_DIR"

echo "=== PyTorch Profiler: ${NPROC}GPU bs=${BS} ==="

# 使用带 profiler 的 bench 模块运行 25 步 (wait=5 + warmup=5 + active=10 + margin=5)
ENABLE_PROFILER=1 \
PROFILER_OUTPUT_DIR="$PROFILER_DIR/trace_${NPROC}gpu_bs${BS}_${TIMESTAMP}" \
torchrun \
    --standalone \
    --nproc_per_node="$NPROC" \
    --master_port="$MASTER_PORT" \
    -m transchrombp.training.train_ddp_bench \
    --model-config "$MODEL_CONFIG" \
    --train-config "$TRAIN_CONFIG" \
    --run-name "profiler_${NPROC}gpu_bs${BS}" \
    --dry-run-steps 25 \
    --batch-size-per-gpu "$BS" \
    2>&1 | tee "$LOG_DIR/profiler_${NPROC}gpu_bs${BS}_${TIMESTAMP}.log"

echo ""
echo "Profiler traces saved to: $PROFILER_DIR/trace_${NPROC}gpu_bs${BS}_${TIMESTAMP}"
echo "View with: tensorboard --logdir=$PROFILER_DIR"
