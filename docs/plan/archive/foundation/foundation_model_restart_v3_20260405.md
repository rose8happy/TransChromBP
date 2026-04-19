# foundation-model restart v3 执行计划（2026-04-05）

> 更新（`2026-04-06 03:30 CST`）：本计划对应的默认决策链已经全部走完。
> `restart v3` 自身 held-out 触发 `fail-or-unsafe`，之后仅在用户显式授权下补了
> `short10 matched no-foundation control` 与 `NT v2 bins16 center-aligned residual`
> 两步受控重入；两者都未给出 clean gain。本文档现在保留为“restart v3 如何被启动并判负”的
> 历史执行计划，当前默认权威动作以
> `docs/plan/archive/foundation/post_chatgpt_pro_priority_execution_20260405.md`、
> `docs/plan/archive/foundation/a6000_foundation_reentry_candidates_20260405.md` 和 `TRACKING.md` 为准。

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

## 4. 本轮执行与最新状态

- 启动时间：
  - 首版：`2026-04-05 01:01 CST`
  - dual-cache restart：`2026-04-05 01:16 CST`
- 机器/GPU：`6000 / A6000 x2`
- launcher：
  - `bash scripts/run_ntv2_residual_short10.sh`
- run name：
  - 首版：`ntv2_residual_short10_s42_20260405`
  - crash 版：`ntv2_residual_short10_s42_20260405_dualcache`
  - smoke：`ntv2_residual_short10_fixsmoke_20260405`
  - 当前：`ntv2_residual_short10_s42_20260405_dualcache_fixddpval`
- launcher PID：
  - 首版：`3798798`（已主动停止）
  - crash 版：`3803855`（已退出）
  - 当前：`3851275`（已结束）
- 日志：
  - 首版：`/data1/zhoujiazhen/bylw_atac/logs/ntv2_residual_short10_20260405_6000.log`
  - crash 版：`/data1/zhoujiazhen/bylw_atac/logs/ntv2_residual_short10_dualcache_20260405_6000.log`
  - 当前：`/data1/zhoujiazhen/bylw_atac/logs/ntv2_residual_short10_dualcache_fixddpval_20260405_6000.log`
- 当前状态：
  - 首版 run 在还停留于单卡 cache 阶段时，被用户要求“cache 也别浪费”后主动终止
  - dual-cache crash 版已完成 `[Step 0] Building dataset-aligned NT v2 cache on CACHE_GPU_IDS=0,1 with 2 workers...`
  - `manifest_train.json` 与 `manifest_valid.json` 已写出，说明 dataset-aligned cache 构建成功
  - crash 版在 `[Step 2]` 进入 `[foundation] epoch=1: residual head only` 后，于 `2026-04-05 03:07 CST` 因 DDP `Expected to have finished reduction in the prior iteration` 报错退出
  - 随后已完成两处代码热修：`residual_head` 动态 freeze 场景下自动开启 `find_unused_parameters`，以及 validation 显式补传 `model_cfg`
  - 6000 上真实 `2-step` smoke `ntv2_residual_short10_fixsmoke_20260405` 已跑通 train + val + checkpoint + `dry run limit reached: 2 steps`
  - 正式 rerun `ntv2_residual_short10_s42_20260405_dualcache_fixddpval` 于 `2026-04-05 15:12 CST` 启动，并在 `2026-04-05 16:41 CST` 左右正常收口到 `epoch 10`
  - `best.pt`、`epoch_010.pt`、`epoch_metrics.jsonl`、`run_meta.json` 已全部写出；`run_meta.json` 记录 `best_epoch=10`、`stopped_early=false`
  - 最新/最佳 validation 指标为：peak `profile_target_jsd_full_mean=0.37030`、peak `count_pearson_full=0.74198`、overall `profile_target_jsd_full_mean=0.39012`、overall `count_pearson_full=0.74679`
  - `2026-04-05 17:01 CST` 首次 held-out `test-full` 评估尝试失败，报错为 `Foundation cache manifest not found: .../manifest_test.json`；根因是 launcher 当前只构建了 `train/valid` cache，没有为 held-out test 预生成 foundation cache
  - 因此已于 `2026-04-05 17:05 CST` 启动 sidecar：先双卡补 `split=test + nonpeak_ratio=1.0` 的 full test cache，再单卡对 `best.pt` 跑 held-out `test-full`
  - sidecar 信息：
    - launcher：`/tmp/ntv2_testfull_eval_20260405.sh`
    - PID：`3872655`
    - 日志：`/data1/zhoujiazhen/bylw_atac/logs/ntv2_residual_short10_testfull_sidecar_20260405_170520.log`
    - 临时 runtime：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/runtime/ntv2_residual_short10/train_tutorial_foundation_short10_testfull_runtime.yaml`
    - 关键配置差异：仅把 `nonpeak_ratio` 从 `0.1` 改成 `1.0`
    - 机器/GPU：`6000 / GPU0,1`（cache），随后 `GPU0`（eval）
    - 预计结束：`2026-04-05 17:25-17:35 CST`
- 当前串行步骤：
  1. 双卡分片生成 NT v2 train/valid cache（已完成）
  2. 跳过 probe（本次未传 `BASELINE_CHECKPOINT`，已完成）
  3. 自动进入双卡 short10 DDP 训练（已完成并正常收口）
  4. 补 `test-full` 所需的 full test cache，再跑 held-out eval（进行中）

## 5. 运行后判定规则（现已完成）

- 这条首跑首先承担“训练接入是否真正可行”的 gate；当前这道 gate 已通过，说明 cache / launcher / train / validation 这条最小闭环已经真正跑通。
- 当前已知事实是：
  - dual-cache cache 路径可用
  - launcher 自动串联可用
  - `residual_head-only` 的双卡训练路径原先存在 DDP unused-parameter / reduction 问题，现已通过 trainer 护栏修复
  - validation 还暴露出一个 `model_cfg` 未传入 `run_validation()` 的 plumbing bug，也已修复
  - held-out 评估阶段还暴露出第三个更晚才会触发的缺口：若只生成 `train/valid` foundation cache，`evaluate_checkpoint --split test --nonpeak-ratio 1.0` 会在 dataset 初始化阶段因缺 `manifest_test.json` 直接失败
- 这些缺口随后都已补齐：
  - `best.pt`
  - `epoch_metrics.jsonl`
  - best checkpoint 对应的 held-out `test-full` 指标
- 现有 validation 读数已经说明这条 run 至少没有表现出“明显优于 corrected-B 锚点”的表象；但正式 `fail-or-unsafe / near-zero / clean gain` 归类仍只看 held-out，不用 validation 先行替代。
- held-out 结果出来后，统一按 `docs/plan/archive/foundation/post_chatgpt_pro_priority_execution_20260405.md` 的决策树执行：
  - 若触发 `fail-or-unsafe`，foundation 线本阶段停表，不扩 `seed1234`、不切 `cross_attention`
  - 若未触发红线，只允许先开 matched `no-foundation short10 control`
  - 只有当前 run 相对 matched control 给出 `clean gain`，才允许进入“中心对齐、高带宽 residual head A/B”

补充结果（`2026-04-06`）：

- `restart v3` 的 held-out 已触发 `fail-or-unsafe`
- 显式授权下补跑的 matched `short10 no-foundation control` 证明 `short10` 预算本身仍健康
- 随后唯一 genuinely new-hypothesis 的 `bins16 center-aligned residual` 也仍明显落后于 matched control
- 因此当前 stable 结论已收口为：`summary / token-fusion / coarse residual / bins16 center-aligned residual`
  这些已实测 recipe family 均未显示 clean gain，foundation 线默认停表

## 6. 当前后续动作

- 本文件不再承载新的待执行步骤，只保留为历史计划与链路修复记录。
- 默认后续动作已经切回论文主线、claim matrix 和归档整理，而不是继续扩这条 foundation 支线。
- 若未来还要重开 foundation，必须先提出一个不同于当前 residual short10 family 的新 hypothesis，并单独说明它为什么不违反现有停表规则。
