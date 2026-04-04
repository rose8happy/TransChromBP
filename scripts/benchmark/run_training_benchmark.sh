#!/usr/bin/env bash
# ==============================================================================
# TransChromBP 训练性能诊断主脚本
# 用途：在 6000 服务器上运行，诊断训练性能瓶颈
# 前提：主实验已结束，GPU 空闲
#
# 使用方法：
#   bash scripts/benchmark/run_training_benchmark.sh
#   DRY_RUN_STEPS=500 bash scripts/benchmark/run_training_benchmark.sh
# ==============================================================================
set -euo pipefail

# ── 环境配置 ──────────────────────────────────────────────────────
export PATH=/data1/zhoujiazhen/bylw_atac/.mamba/envs/transchrombp/bin:$PATH
export LD_LIBRARY_PATH=/data1/zhoujiazhen/bylw_atac/.mamba/envs/transchrombp/lib:${LD_LIBRARY_PATH:-}
export PYTHONPATH=/data1/zhoujiazhen/bylw_atac/TransChromBP/src:${PYTHONPATH:-}

ROOT_DIR="/data1/zhoujiazhen/bylw_atac/TransChromBP"
MODEL_CONFIG="$ROOT_DIR/configs/model/transchrombp_base.yaml"
TRAIN_CONFIG="$ROOT_DIR/configs/train/train_benchmark_base.yaml"
OUTPUT_DIR="$ROOT_DIR/outputs"
LOG_DIR="$OUTPUT_DIR/logs/benchmark"
PROFILER_DIR="$OUTPUT_DIR/profiler_traces"
MASTER_PORT="${MASTER_PORT:-29510}"

DRY_RUN_STEPS="${DRY_RUN_STEPS:-300}"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
SUMMARY_FILE="$LOG_DIR/benchmark_summary_${TIMESTAMP}.txt"

cd "$ROOT_DIR"
mkdir -p "$LOG_DIR" "$PROFILER_DIR"

header() {
    echo ""
    echo "========================================================"
    echo " $1"
    echo "========================================================"
}

header "TransChromBP Training Performance Benchmark"
echo " Timestamp:      $TIMESTAMP"
echo " DRY_RUN_STEPS:  $DRY_RUN_STEPS"
echo " ROOT_DIR:       $ROOT_DIR"
echo " TRAIN_CONFIG:   $TRAIN_CONFIG"
echo " Summary:        $SUMMARY_FILE"

# ── 实验 0: 环境信息 ──────────────────────────────────────────────
header "[Exp 0] Environment Info"
ENV_LOG="$LOG_DIR/env_info_${TIMESTAMP}.log"

python3 -c "
import torch
print(f'PyTorch version: {torch.__version__}')
print(f'CUDA version: {torch.version.cuda}')
print(f'cuDNN version: {torch.backends.cudnn.version()}')
print(f'GPU count: {torch.cuda.device_count()}')
for i in range(torch.cuda.device_count()):
    p = torch.cuda.get_device_properties(i)
    print(f'  GPU {i}: {p.name} ({p.total_memory / 1e9:.1f} GB)')
print(f'torch.backends.cuda.matmul.allow_tf32: {torch.backends.cuda.matmul.allow_tf32}')
print(f'torch.backends.cudnn.allow_tf32: {torch.backends.cudnn.allow_tf32}')
try:
    print(f'NCCL version: {torch.cuda.nccl.version()}')
except Exception as e:
    print(f'NCCL version: unavailable ({e})')
" 2>&1 | tee "$ENV_LOG"

nvidia-smi 2>&1 | tee -a "$ENV_LOG"
echo "" | tee -a "$ENV_LOG"
echo "[env] Log saved: $ENV_LOG"

# ── 辅助函数 ──────────────────────────────────────────────────────
# run_bench <nproc> <bs> <run_name> [extra_env_vars] [dry_run_steps_override]
run_bench() {
    local nproc="$1"
    local bs="$2"
    local run_name="$3"
    local extra_env="${4:-}"
    local dry_steps="${5:-$DRY_RUN_STEPS}"
    local log_file="$LOG_DIR/${run_name}_${TIMESTAMP}.log"
    local port=$((MASTER_PORT + RANDOM % 100))

    echo ""
    echo "--- Running: $run_name (nproc=$nproc, bs=$bs, dry_steps=$dry_steps, port=$port) ---"
    echo "    Log: $log_file"

    # GPU 基线状态
    echo "=== GPU status before run ===" > "$log_file"
    nvidia-smi --query-gpu=index,memory.used,memory.total,utilization.gpu,power.draw \
        --format=csv 2>/dev/null >> "$log_file" || true
    echo "" >> "$log_file"

    # 在后台采样 GPU 状态 (每 5 秒)
    local gpu_sample_file="$LOG_DIR/${run_name}_gpu_samples_${TIMESTAMP}.csv"
    echo "timestamp,gpu_idx,mem_used_MiB,mem_total_MiB,gpu_util_pct,power_W" > "$gpu_sample_file"
    (
        while true; do
            nvidia-smi --query-gpu=index,memory.used,memory.total,utilization.gpu,power.draw \
                --format=csv,noheader,nounits 2>/dev/null | head -"$nproc" | \
                while IFS= read -r line; do
                    echo "$(date +%s),$line"
                done >> "$gpu_sample_file"
            sleep 5
        done
    ) &
    local gpu_monitor_pid=$!

    # 启动训练 (使用 bench 模块，带 timing + 可选 profiler)
    local exit_code=0
    env $extra_env \
    torchrun \
        --standalone \
        --nproc_per_node="$nproc" \
        --master_port="$port" \
        -m transchrombp.training.train_ddp_bench \
        --model-config "$MODEL_CONFIG" \
        --train-config "$TRAIN_CONFIG" \
        --run-name "$run_name" \
        --dry-run-steps "$dry_steps" \
        --batch-size-per-gpu "$bs" \
        2>&1 | tee -a "$log_file" || exit_code=$?

    # 停止 GPU 监控
    kill "$gpu_monitor_pid" 2>/dev/null || true
    wait "$gpu_monitor_pid" 2>/dev/null || true

    # GPU 结束状态
    echo "" >> "$log_file"
    echo "=== GPU status after run ===" >> "$log_file"
    nvidia-smi --query-gpu=index,memory.used,memory.total,utilization.gpu,power.draw \
        --format=csv 2>/dev/null >> "$log_file" || true

    if [ "$exit_code" -ne 0 ]; then
        echo "    [WARN] $run_name exited with code $exit_code (possible OOM)"
        echo "$run_name: FAILED (exit_code=$exit_code)" >> "$SUMMARY_FILE"
        return "$exit_code"
    fi

    echo "    Done: $run_name"

    # 提取 TIMING SUMMARY 到汇总文件
    echo "" >> "$SUMMARY_FILE"
    echo "--- $run_name (nproc=$nproc, bs=$bs) ---" >> "$SUMMARY_FILE"
    grep -A 20 "TIMING SUMMARY" "$log_file" >> "$SUMMARY_FILE" 2>/dev/null || echo "  (no timing summary found)" >> "$SUMMARY_FILE"
    grep "GPU MEMORY SUMMARY" -A 10 "$log_file" >> "$SUMMARY_FILE" 2>/dev/null || true

    # 提取 GPU 采样峰值
    if [ -f "$gpu_sample_file" ]; then
        local peak_mem
        peak_mem=$(tail -n +2 "$gpu_sample_file" | cut -d',' -f3 | sort -n | tail -1)
        echo "  peak_gpu_mem_used: ${peak_mem:-N/A} MiB" >> "$SUMMARY_FILE"
    fi
}

# ── 初始化汇总文件 ──────────────────────────────────────────────
echo "TransChromBP Training Performance Benchmark Summary" > "$SUMMARY_FILE"
echo "Timestamp: $TIMESTAMP" >> "$SUMMARY_FILE"
echo "DRY_RUN_STEPS: $DRY_RUN_STEPS" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

# ── 实验 1: 单卡基线 (bs=16) ─────────────────────────────────────
header "[Exp 1] Single GPU Baseline (bs=16)"
run_bench 1 16 "bench_1gpu_bs16"
sleep 3

# ── 实验 2: 双卡基线 (bs=16) ─────────────────────────────────────
header "[Exp 2] Dual GPU Baseline (bs=16)"
run_bench 2 16 "bench_2gpu_bs16"
sleep 3

# ── 实验 3: Batch Size Sweep (双卡) ──────────────────────────────
header "[Exp 3] Batch Size Sweep (2 GPUs)"

for bs in 24 32 48; do
    sleep 3
    run_bench 2 "$bs" "bench_2gpu_bs${bs}" || true
done

# ── 实验 4: PyTorch Profiler (双卡, bs=16) ───────────────────────
header "[Exp 4] PyTorch Profiler (2 GPUs, bs=16, ~20 steps)"

PROF_TRACE_DIR="$PROFILER_DIR/trace_2gpu_bs16_${TIMESTAMP}"
run_bench 2 16 "bench_profiler_2gpu_bs16" \
    "ENABLE_PROFILER=1 PROFILER_OUTPUT_DIR=$PROF_TRACE_DIR" \
    25 || true

# ── 最终汇总 ─────────────────────────────────────────────────────
header "Benchmark Complete"
echo " Timestamp:    $TIMESTAMP"
echo " Logs:         $LOG_DIR"
echo " Summary:      $SUMMARY_FILE"
echo " Profiler:     $PROFILER_DIR"
echo ""
echo " ── Summary File Contents ──"
cat "$SUMMARY_FILE"
echo ""
echo " Next steps:"
echo "   1. Review summary above (or cat $SUMMARY_FILE)"
echo "   2. If profiler ran, view with: tensorboard --logdir=$PROFILER_DIR"
echo "   3. Run parse script: python3 scripts/benchmark/parse_benchmark_logs.py $LOG_DIR/*_${TIMESTAMP}.log"
echo "   4. Write final report based on docs/plan/training_perf_benchmark_plan.md"
