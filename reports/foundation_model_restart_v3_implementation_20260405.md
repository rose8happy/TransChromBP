# foundation-model restart v3 实施记录（2026-04-05）

## 结论

- 这轮不是只写了计划或脚本占位，而是已经把 foundation-model restart v3 的第一条训练链路真正起跑。
- 当前状态可明确分成三层：
  - 基础设施已落地：训练、评估、数据、缓存、配置、launcher 都已补齐。
  - 动态可用性已验证：`residual_head` 与 `cross_attention` 在 6000 上都通过了真实 forward + `compute_losses`。
  - 首条正式任务已启动，且已根据用户要求升级为“cache 也走双卡”的 restart 版本：`NT v2 cached residual short10` 当前正在 6000 以双卡 shard-cache 方式执行。

## 已完成实现

### 1. 统一 foundation_model 训练接口

- 新增 `foundation_model.enabled/source/mode` 配置入口。
- 支持三种 mode：
  - `cross_attention`
  - `residual_head`
  - `distill_only`
- 模型前向现在可消费两类外部输入：
  - `foundation_tokens`
  - `foundation_summary`

### 2. 新接法

- `FoundationCrossAttentionAdapter`
  - 用于 late token-level cross-attention
  - 输出投影零初始化
  - residual gate 负偏置初始化
- `FoundationResidualHead`
  - 只在 debiased signal 上预测 profile/count residual
  - 输出头零初始化
  - 支持 coarse profile bins -> full profile delta 的上采样

### 3. 数据与 loss plumbing

- dataset 新增：
  - `foundation_cache_dir`
  - `foundation_cache_features`
  - `teacher_cache_dir`
  - `teacher_target_names`
- `train_ddp.py` 新增：
  - cached foundation features 读取
  - teacher targets 读取
  - distillation losses
  - residual-head freeze/unfreeze policy

### 4. 新脚本与配置

- 新脚本：
  - `scripts/build_foundation_cache.py`
  - `scripts/foundation_model_probe.py`
  - `scripts/run_ntv2_residual_short10.sh`
- 新配置：
  - `transchrombp_teacher_v2_center_pool_ntv2_residual.yaml`
  - `transchrombp_teacher_v2_center_pool_ntv2_cross_attention.yaml`
  - `train_tutorial_foundation_short10.yaml`

## 验证证据

### 静态验证

- 本地 `py_compile` 通过。
- 6000 远端 `py_compile` 通过。
- launcher `bash -n` 通过。

### 动态验证

- 6000 环境：
  - `/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b`
- synthetic dry-run 结果：
  - residual path `profile=(2,1000)`、`count=(2,1)`
  - cross-attention path `profile=(2,1000)`、`count=(2,1)`

### 缓存 dry-run

- tutorial canonical:
  - train `426545`
  - valid `23417`
- NT v2 cache 计划：
  - `layer_07__global_mean`
  - `layer_07__bins4_mean`
  - `layer_14__global_mean`
  - `layer_14__bins4_mean`
- 总 fp16 空间预算：
  - 约 `8.136 GiB`

### 双卡 cache smoke

- 在正式替换首版 run 前，已额外做过一次小规模双卡 smoke：
  - `world_size=2`
  - `max_records=64`
  - `split=valid`
  - `layer_07__global_mean`
- 结果：
  - 最终 `valid_layer_07__global_mean.f16.npy` 正常生成
  - `manifest_valid.json` 正常写出，`world_size=2`
  - `records_valid.jsonl` 行数正确为 `64`
  - shard 临时文件已被合并并清理

## 已启动任务

- run name：
  - 首版：`ntv2_residual_short10_s42_20260405`
  - 当前：`ntv2_residual_short10_s42_20260405_dualcache`
- 机器：
  - `6000`
- GPU：
  - 首版启动时：`GPU0 active, GPU1 idle`
  - 当前 dual-cache restart：`GPU0 active, GPU1 active`
- 日志：
  - 首版：`/data1/zhoujiazhen/bylw_atac/logs/ntv2_residual_short10_20260405_6000.log`
  - 当前：`/data1/zhoujiazhen/bylw_atac/logs/ntv2_residual_short10_dualcache_20260405_6000.log`
- 进程：
  - 首版 launcher PID `3798798`（已停止）
  - 当前 launcher PID `3803855`
- 启动后首个确认：
  - 首版日志已写出 launcher header，但因用户明确要求“cache 阶段也别浪费”而在训练开始前被替换
  - 当前 dual-cache restart 日志已写出 `CACHE_GPU_IDS=0,1`、`CACHE_NPROC_PER_NODE=2`
  - 当前已进入 `torchrun --nproc_per_node=2` 的 cache build
  - 两张 A6000 都已被实际占用

## 当前判断

- 旧的 NT v2 sampled probe 已经说明它“有独立信号但未显出互补性”。
- 这次 residual-head short10 是必要的下一步，因为它第一次真正回答：
  - 在不直接扰动主 backbone 的前提下，
  - 让 NT v2 只学习 debiased residual，
  - 是否仍然拿不到正收益。
- 若这条线仍无收益，NT 路线的否证强度会明显提高；若它有收益，就能直接进入第二个 seed，而不需要再回到旧的 concat / FiLM recipe。

## 下一步

- 等待 cache build 与 short10 训练自动完成。
- 收口后优先看 valid `JSD/count_r` 是否优于当前 corrected-B baseline。
- 如果有正迹象，下一步最值得做的是补 `seed1234`；如果无明显增益，则优先切到 `cross_attention` 或 teacher/distill 线，而不是继续给旧 recipe 加算力。
