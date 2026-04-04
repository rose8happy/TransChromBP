# foundation-model restart v3 执行计划（2026-04-05）

## 1. 目标

- 把“接入 foundation model 可能带来正收益”从 probe 级假设推进到真正可训练、可验证的代码链路。
- 这轮不再停留在 `Genos/Caduceus/NT v2` 的 sampled probe，而是先打通：
  - dataset-aligned cache
  - unified `foundation_model` 训练接口
  - `residual_head` / `cross_attention` 两种新接法
  - teacher/distill 所需的 batch / loss / cache plumbing
- 第一条实际执行线固定为：
  - `corrected-B / center-pool` baseline
  - `NT v2 cached residual_head`
  - tutorial canonical short10
  - seed42

## 2. 本轮已实施代码

### 2.1 训练与评估主链路

- 在 `src/transchrombp/models/transchrombp.py` 中新增通用 `foundation_model` 接口：
  - `foundation_tokens`
  - `foundation_summary`
  - `mode=cross_attention | residual_head | distill_only`
- 新增 `src/transchrombp/models/foundation_adapter.py`：
  - `FoundationCrossAttentionAdapter`
  - `FoundationResidualHead`
- 在 `src/transchrombp/training/train_ddp.py` 中补齐：
  - `extract_foundation_feature_kwargs`
  - `extract_teacher_targets`
  - `apply_foundation_freeze_policy`
  - `distill_profile/count/rank` 三类辅助损失
  - `foundation_model` 的 config 校验与 dataloader 透传
- 在 `src/transchrombp/evaluation/evaluate_checkpoint.py` 中同步支持 cached foundation features 与 teacher targets。

### 2.2 数据侧

- `src/transchrombp/data/real_data.py` 现已支持三类外部缓存：
  - `genos_cache_dir`
  - `foundation_cache_dir`
  - `teacher_cache_dir`
- cache/teacher manifest 统一按 `record_sha1 + n_records + split` 校验，防止 dataset 顺序错位。

### 2.3 脚本与配置

- 新增训练对齐缓存构建脚本：
  - `scripts/build_foundation_cache.py`
- 新增统一 probe 脚本：
  - `scripts/foundation_model_probe.py`
- 升级 `scripts/nt_v2_probe.py`，使其同时能输出 `true_profile16`、manifest 和训练侧可复用的 cache 文件。
- 新增配置：
  - `configs/model/transchrombp_teacher_v2_center_pool_ntv2_residual.yaml`
  - `configs/model/transchrombp_teacher_v2_center_pool_ntv2_cross_attention.yaml`
  - `configs/train/train_tutorial_foundation_short10.yaml`
- 新增 launcher：
  - `scripts/run_ntv2_residual_short10.sh`

## 3. 已完成验证

### 3.1 本地静态验证

- `python3 -m py_compile` 已通过：
  - `build_foundation_cache.py`
  - `foundation_model_probe.py`
  - `nt_v2_probe.py`
  - `foundation_adapter.py`
  - `transchrombp.py`
  - `real_data.py`
  - `train_ddp.py`
  - `evaluate_checkpoint.py`
- `bash -n scripts/run_ntv2_residual_short10.sh` 已通过。

### 3.2 6000 远端最小动态验证

- 环境：`/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b`
- 结果：
  - `residual_head` 路径可真实 forward + `compute_losses`
  - `cross_attention` 路径可真实 forward + `compute_losses`
  - synthetic dry-run 输出 shape：
    - `profile=(2, 1000)`
    - `count=(2, 1)`

### 3.3 dataset-aligned cache dry-run

- data config：`configs/data/data_tutorial_canonical_v1.yaml`
- train config：`configs/train/train_tutorial_foundation_short10.yaml`
- backend：`nt_v2`
- layers/features：`7,14 × global_mean,bins4_mean`
- dry-run 结果：
  - train `n_records=426545`
  - valid `n_records=23417`
  - 总 fp16 体积约 `8.136 GiB`
  - 输出目录剩余空间约 `9715.7 GiB`

## 4. 已启动后台任务

- 启动时间：
  - 首版：`2026-04-05 01:01 CST`
  - dual-cache restart：`2026-04-05 01:16 CST`
- 机器/GPU：`6000 / A6000 x2`
- launcher：
  - `bash scripts/run_ntv2_residual_short10.sh`
- run name：
  - 首版：`ntv2_residual_short10_s42_20260405`
  - 当前：`ntv2_residual_short10_s42_20260405_dualcache`
- launcher PID：
  - 首版：`3798798`（已主动停止）
  - 当前：`3803855`
- 日志：
  - 首版：`/data1/zhoujiazhen/bylw_atac/logs/ntv2_residual_short10_20260405_6000.log`
  - 当前：`/data1/zhoujiazhen/bylw_atac/logs/ntv2_residual_short10_dualcache_20260405_6000.log`
- 当前状态：
  - 首版 run 在还停留于单卡 cache 阶段时，被用户要求“cache 也别浪费”后主动终止
  - 当前 dual-cache restart 已进入 `[Step 0] Building dataset-aligned NT v2 cache on CACHE_GPU_IDS=0,1 with 2 workers...`
  - 当前即时观测为 `GPU0 ~3540 MiB / util 94%`、`GPU1 ~3540 MiB / util 93%`
- 当前串行步骤：
  1. 双卡分片生成 NT v2 train/valid cache
  2. 跳过 probe（本次未传 `BASELINE_CHECKPOINT`）
  3. 自动进入双卡 short10 DDP 训练

## 5. 运行后判定规则

- 这条首跑只承担“训练接入是否可行 + 是否值得扩 seed”的 gate，不承担最终路线定论。
- run 结束后先看：
  - `best.pt`
  - `epoch_metrics.jsonl`
  - valid `peak.profile_target_jsd_full_mean`
  - valid `count_pearson_full`
- 若相对当前 baseline 有明确正向迹象：
  - 立即补 `seed1234`
- 若接近零增益但不坏：
  - 评估是否改为 `cross_attention`
- 若明显变坏：
  - NT 路线本阶段按 negative outcome 收口，优先转 teacher/distill 线

## 6. 下一步

- 等 cache build 结束并确认 short10 自动起跑。
- run 收口后优先给出：
  - 指标对比
  - 是否晋级 `seed1234`
  - 是否继续 `cross_attention`
  - 是否转入 `gReLU / AlphaGenome` teacher 线
