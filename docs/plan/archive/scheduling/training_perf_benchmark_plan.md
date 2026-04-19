# TransChromBP 训练性能诊断方案

> 创建日期: 2026-03-18
> 状态: 待执行（等当前 teacher_v2_20260317 训练结束后运行）

## 一、背景与目标

**当前主实验**: teacher_v2_20260317 (teacher_v2_main_learnable)
**硬件**: 2 × RTX A6000 48GB (6000 服务器)
**现象**:
- 两张卡显存都只用了约 4.2 / 49 GB
- GPU util 长期 97%–99%，功耗接近满载
- 训练一轮约 4–5 小时

**已知配置**:
```yaml
precision: bf16
grad_accum_steps: 1
batch_size_per_gpu: 16
profile_bin_size: 1
input_len: 2114
num_workers: 2
```

**目标**：回答 5 个核心问题:
1. 当前训练慢是不是正常现象
2. 双卡扩展是否有效
3. 时间主要花在 dataloader、forward/backward、optimizer 还是通信
4. batch_size_per_gpu 还有没有必要继续增大
5. input_len=4096 应不应该单独作为新实验线

## 二、诊断实验设计

### 实验 0: 建立 benchmark 基线配置

从 `train_tutorial_teacher_v2_main.yaml` 复制出 benchmark 版本，改动点：
- `max_epochs: 999` (由 `--dry-run-steps` 控制)
- `validate_every_epochs: 999` (禁用 validation)
- `checkpoint_every_epochs: 999` (禁用 checkpoint)
- `early_stop_patience: 0` (禁用 early stop)
- `log_every_steps: 10` (更频繁 logging)
- `run_name: benchmark_xxx`

使用已有的 `--dry-run-steps 300` 机制，每组实验只跑 300 步。

**配置文件**: `configs/train/train_benchmark_base.yaml` (将部署到远端)
**与主实验 diff**: 仅关闭 validation/checkpoint/early stop，其他参数完全一致

### 实验 1: Step Time 分解 (timing instrumentation)

需要在 `train_ddp.py` 训练循环中加入轻量级 timing，统计：
- `data_time`: 取 batch 的时间（DataLoader.__next__）
- `forward_backward_time`: forward + loss + backward
- `optimizer_time`: optimizer.step() + zero_grad + scheduler.step()
- `step_time`: 单步总耗时

**实现方式**: 新建 `train_ddp_bench.py`（复制 train_ddp.py + timing patch），不修改原文件。
见下文 [timing patch 设计](#timing-patch)。

### 实验 2: 1 卡 vs 2 卡对比

| 配置 | NPROC_PER_NODE | batch_size_per_gpu | dry_run_steps |
|------|---------------|-------------------|---------------|
| A: 单卡 | 1 | 16 | 300 |
| B: 双卡 | 2 | 16 | 300 |

指标：
- avg step_time
- samples/s (= batch_size_per_gpu × world_size / step_time)
- speedup = throughput_2gpu / throughput_1gpu
- scaling efficiency = speedup / 2

判断标准：
- scaling efficiency > 0.85: 正常
- 0.7–0.85: 可接受但有优化空间
- < 0.7: 有瓶颈需排查

**启动方式**: `torchrun --standalone --nproc_per_node=N`
**分布式模式**: 确认是 DDP (DistributedDataParallel)，不是 DataParallel

### 实验 3: 运行时环境确认

在 benchmark 脚本开头打印并记录：
```python
import torch
print(f"PyTorch version: {torch.__version__}")
print(f"CUDA version: {torch.version.cuda}")
print(f"cuDNN version: {torch.backends.cudnn.version()}")
print(f"bf16 autocast enabled: True (config)")
print(f"torch.backends.cuda.matmul.allow_tf32: {torch.backends.cuda.matmul.allow_tf32}")
print(f"torch.backends.cudnn.allow_tf32: {torch.backends.cudnn.allow_tf32}")
try:
    print(f"NCCL version: {torch.cuda.nccl.version()}")
except:
    print("NCCL version: unavailable")
```

### 实验 4: Batch Size Sweep

| bs_per_gpu | 其他参数 | dry_run_steps |
|-----------|---------|---------------|
| 16 (基线) | 不变 | 300 |
| 24 | 不变 | 300 |
| 32 | 不变 | 300 |
| 48 | 不变 | 300 |

使用 `--batch-size-per-gpu N` 命令行 override。

收集指标：
- OOM 与否
- avg_step_time
- samples/s
- loss 是否异常抖动
- GPU 显存峰值 (`nvidia-smi`)
- GPU util/功耗

输出汇总表：
```
| bs_per_gpu | OOM | avg_step_time | samples/s | peak_mem | gpu_util | comment |
```

### 实验 5: PyTorch Profiler (20 步)

选双卡 + bs16 配置，用 `torch.profiler` profile 20 步。

```python
from torch.profiler import profile, ProfilerActivity, schedule

with profile(
    activities=[ProfilerActivity.CPU, ProfilerActivity.CUDA],
    schedule=schedule(wait=5, warmup=5, active=10, repeat=1),
    on_trace_ready=torch.profiler.tensorboard_trace_handler("./profiler_logs"),
    record_shapes=True,
    profile_memory=True,
    with_stack=False,
) as prof:
    for step_idx, batch in enumerate(train_loader, start=1):
        # ... training step ...
        prof.step()
        if step_idx >= 20:
            break

print(prof.key_averages().table(sort_by="cuda_time_total", row_limit=20))
```

分析重点：
- CUDA 时间主要花在哪类算子
- all-reduce / NCCL 通信耗时占比
- 是否有大量零碎小 kernel
- backward vs forward 时间比
- DataLoader/CPU 是否成为瓶颈

### 实验 6: input_len=4096 性质判断

**分析依据**（无需实验运行）：

从代码和配置看，`input_len` 的改变影响：
1. **输入序列长度**: `seq` tensor shape 从 `[B, 2114, 4]` → `[B, 4096, 4]`
2. **output_len / supervised_bp**: 目前固定 1000bp，需确认 input_len 变化后是否自动调整
3. **模型 max_len**: `transchrombp_base.yaml` 中 `max_len: 2114`，需改为 4096
4. **Transformer 计算量**: Self-attention 是 O(n²)，2114→4096 约增加 3.75× 计算量
5. **数据窗口**: 更长输入意味着捕获更远的调控信号，改变了任务定义

**预期结论**: input_len=4096 是**新实验线**，理由：
- 改变了模型的感受野和输入上下文
- 可能需要调整 output_len / supervised_bp 以保持一致
- Benchmark 口径改变，结果不可直接与 2114 比较
- 计算量质变（~3.75×），不是简单优化

## 三、Timing Patch 设计 {#timing-patch}

在训练循环中加入最小侵入 timing，方案如下（伪代码）：

```python
import time

# 在 epoch 循环开始前初始化
timing_stats = {"data": [], "fwd_bwd": [], "opt": [], "step": []}

# 在 batch 循环中
t_step_start = time.perf_counter()

# --- data loading ---
t_data_end = time.perf_counter()
timing_stats["data"].append(t_data_end - t_step_start)

# --- forward + backward ---
torch.cuda.synchronize()
t_fwd_start = time.perf_counter()

with no_sync_ctx:
    with torch.autocast(...):
        outputs = model(seq)
        loss, metrics = compute_losses(...)
        loss = loss / grad_accum_steps
    loss.backward()

torch.cuda.synchronize()
t_fwd_end = time.perf_counter()
timing_stats["fwd_bwd"].append(t_fwd_end - t_fwd_start)

# --- optimizer ---
t_opt_start = time.perf_counter()
if should_step:
    optimizer.step()
    optimizer.zero_grad(set_to_none=True)
    scheduler.step()
torch.cuda.synchronize()
t_opt_end = time.perf_counter()
timing_stats["opt"].append(t_opt_end - t_opt_start)

t_step_end = time.perf_counter()
timing_stats["step"].append(t_step_end - t_step_start)

# 每 log_every 步打印 rolling average
if step_idx % log_every == 0:
    print(f"[timing] data={mean(timing_stats['data'][-log_every:]):.4f}s "
          f"fwd_bwd={mean(timing_stats['fwd_bwd'][-log_every:]):.4f}s "
          f"opt={mean(timing_stats['opt'][-log_every:]):.4f}s "
          f"step={mean(timing_stats['step'][-log_every:]):.4f}s "
          f"samples/s={batch_size * world_size / mean(timing_stats['step'][-log_every:]):.1f}")

# dry_run 结束后输出总平均
print(f"[timing-summary] avg_data={mean(...):.4f}s avg_fwd_bwd={mean(...):.4f}s "
      f"avg_opt={mean(...):.4f}s avg_step={mean(...):.4f}s "
      f"avg_samples_per_sec={...:.1f}")
```

**重要**: 必须在 forward/backward 前后使用 `torch.cuda.synchronize()`，否则 CUDA 异步操作会导致 timing 不准。

## 四、执行脚本

所有脚本将放在 `scripts/benchmark/` 目录下，部署到远端 6000 服务器运行。

### 4.1 Benchmark 配置文件

文件: `configs/train/train_benchmark_base.yaml`（远端部署路径）

### 4.2 主诊断脚本

文件: `scripts/benchmark/run_training_benchmark.sh`

功能：
1. 打印环境信息
2. 单卡 300 步 benchmark
3. 双卡 300 步 benchmark
4. batch sweep (bs=16,24,32,48) × 300 步
5. profiler run (双卡 bs16, 20步)
6. 汇总输出

### 4.3 Timing 注入

方案 A（推荐）: 创建 `train_ddp_bench.py`，从 `train_ddp.py` 复制并注入 timing
方案 B（备选）: 通过环境变量 `BENCH_TIMING=1` 在原文件中控制 timing 代码路径

推荐方案 A，避免任何风险。

## 五、预期交付

运行完成后输出结构化报告 `reports/training_perf_benchmark_report.md`，包含：

1. **执行摘要**
2. **实验设置**（配置路径、命令、硬件、软件版本）
3. **基线结果**（含 step time 分解: data / fwd_bwd / opt）
4. **1 卡 vs 2 卡对比**（throughput, speedup, scaling efficiency）
5. **Batch Sweep 结果**（汇总表）
6. **Profiler 结果**（top-5 热点 + 分析）
7. **对 input_len=4096 的判断**
8. **最终建议**（按优先级排序）

必须明确回答：
- [ ] 当前 2×A6000 训练 4–5 小时/epoch，是正常耗时还是实现有问题
- [ ] 当前瓶颈主要在哪
- [ ] 双卡是否真的有效
- [ ] batch_size_per_gpu 是否该从 16 提到 24/32/48
- [ ] input_len=4096 是否必须单列为新实验

## 六、风险与注意事项

1. **不修改正在运行的主实验**：所有 benchmark 使用独立配置和输出目录
2. **GPU 抢占**：benchmark 需等主实验结束后再跑，或确保主实验不在跑
3. **OOM**: batch_size=48 可能 OOM，正常记录即可，不视为失败
4. **Profiler 开销**: profiler 本身有约 10–20% 开销，不要用 profiler 数据衡量绝对性能
5. **seed 固定**: benchmark 使用 `seed: 1234`（与主实验一致），保证可复现
