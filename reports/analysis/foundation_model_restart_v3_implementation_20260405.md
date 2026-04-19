# foundation-model restart v3 实施记录（2026-04-05）

> 更新说明（`2026-04-05 18:01 CST`）：本报告中的“fail 后默认停表”属于当时默认执行建议；若用户显式授权受控重入，以
> `docs/plan/archive/foundation/a6000_foundation_reentry_candidates_20260405.md` 和 `TRACKING.md` 最新状态为准。

> 更新说明（`2026-04-06 03:30 CST`）：随后按显式授权执行的唯一受控重入
> `ntv2_bins16_centerres_short10_s42_20260406_dual` 已完成 full held-out。它相对 matched
> `short10_nofoundation_control` 仍明显更差：overall `count_pearson_full=0.7727` vs `0.8457`、
> `profile_target_jsd_full_mean=0.4625` vs `0.4395`；peak `count_pearson_full=0.7516` vs
> `0.8298`、`profile_target_jsd_full_mean=0.3588` vs `0.3193`。因此当前结论已从“旧 coarse-summary
> residual 已 fail-or-unsafe”收紧为“当前 NT v2 residual short10 family（含 bins16 center-aligned
> 重入）均未显示 clean gain，foundation 线默认停表”。

## 结论

- 这轮不是只写了计划或脚本占位，而是已经把 foundation-model restart v3 的第一条训练链路真正起跑。
- 当前状态可明确分成三层：
  - 基础设施已落地：训练、评估、数据、缓存、配置、launcher 都已补齐。
  - 动态可用性已验证：`residual_head` 与 `cross_attention` 在 6000 上都通过了真实 forward + `compute_losses`。
  - 首条正式任务不仅经历了“双卡 cache 完成 -> 双卡训练起跑 -> DDP/validation bug 修复 -> held-out gate 补齐”，其后续的 matched control 与唯一 bins16 受控重入也都已执行完成；当前已经有足够结果把这条 family 收口为默认停表。

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

## 已执行任务与最新状态

- run name：
  - 首版：`ntv2_residual_short10_s42_20260405`
  - crash 版：`ntv2_residual_short10_s42_20260405_dualcache`
  - smoke：`ntv2_residual_short10_fixsmoke_20260405`
  - 当前：`ntv2_residual_short10_s42_20260405_dualcache_fixddpval`
- 机器：
  - `6000`
- GPU：
  - 首版启动时：`GPU0 active, GPU1 idle`
  - dual-cache cache 阶段：`GPU0 active, GPU1 active`
  - 当前：run 已结束，`GPU0/GPU1` 均空闲
- 日志：
  - 首版：`/data1/zhoujiazhen/bylw_atac/logs/ntv2_residual_short10_20260405_6000.log`
  - crash 版：`/data1/zhoujiazhen/bylw_atac/logs/ntv2_residual_short10_dualcache_20260405_6000.log`
  - 当前：`/data1/zhoujiazhen/bylw_atac/logs/ntv2_residual_short10_dualcache_fixddpval_20260405_6000.log`
- 进程：
  - 首版 launcher PID `3798798`（已停止）
  - crash 版 launcher PID `3803855`（已退出）
  - 当前 launcher PID `3851275`（已结束）
- 启动后首个确认：
  - 首版日志已写出 launcher header，但因用户明确要求“cache 阶段也别浪费”而在训练开始前被替换
  - 当前 dual-cache restart 日志已写出 `CACHE_GPU_IDS=0,1`、`CACHE_NPROC_PER_NODE=2`
  - `Step 0` 双卡 cache 已完成，`manifest_train.json` 与 `manifest_valid.json` 均已写出
  - `Step 1` probe 因未设置 `BASELINE_CHECKPOINT` 被明确跳过
  - `Step 2` 双卡训练已实际起跑，并写出 `[foundation] epoch=1: residual head only`
  - dual-cache crash 版随后在 `2026-04-05 03:07 CST` 因 DDP `Expected to have finished reduction in the prior iteration` 错误退出
  - 之后已完成两处代码修复：
    - trainer 在 `foundation residual_head + freeze_backbone_until_epoch>0` 下自动开启 `find_unused_parameters`
    - validation 路径显式补传 `model_cfg`，修复 `run_validation()` 的 `NameError`
  - 6000 上真实 `2-step` smoke 已跑通 train + val + checkpoint + meta 写回，确认两处热修都打中了真实故障路径
  - 正式 rerun 已于 `2026-04-05 15:12 CST` 重启，并在 `2026-04-05 16:41 CST` 左右正常收口到 `epoch 10`
  - `best.pt`、`epoch_010.pt`、`epoch_metrics.jsonl`、`run_meta.json` 都已写出；`run_meta.json` 记录 `best_epoch=10`、`stopped_early=false`
  - 当前 best/最新 validation 为：peak `profile_target_jsd_full_mean=0.37030`、peak `count_pearson_full=0.74198`、overall `profile_target_jsd_full_mean=0.39012`、overall `count_pearson_full=0.74679`
  - `2026-04-05 17:01 CST` 第一次 held-out `test-full` 尝试已真实触发一个新的 pipeline 缺口：`evaluate_checkpoint` 在 dataset 初始化阶段因缺 `.../foundation_cache/ntv2_tutorial_canonical_v1/manifest_test.json` 失败
  - 这说明当前 restart v3 的 launcher 虽然足以支持 train/validation，却还不足以支持 post-run held-out gate，因为它只生成了 `train/valid` foundation cache
  - 因此已于 `2026-04-05 17:05 CST` 启动 sidecar `full test cache -> held-out test-full`：
    - PID：`3872655`（已结束）
    - 日志：`/data1/zhoujiazhen/bylw_atac/logs/ntv2_residual_short10_testfull_sidecar_20260405_170520.log`
    - 临时 runtime：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/runtime/ntv2_residual_short10/train_tutorial_foundation_short10_testfull_runtime.yaml`
    - 唯一配置差异：`nonpeak_ratio=1.0`
    - 资源：`6000 / GPU0,1` 先补 full test cache，之后 `GPU0` 跑 eval
  - `2026-04-05 17:22 CST` 左右 sidecar 完成并写出 held-out JSON：
    - `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/ntv2_residual_short10_s42_20260405_dualcache_fixddpval_best_test_full_20260405_170520.json`
    - overall：`count_r=0.6800`、`jsd=0.4612`、`prof_r=0.6606`、`auroc=0.7950`、`auprc=0.7983`、`f1=0.7387`
    - peak：`count_r=0.7729`、`jsd=0.3560`、`prof_r=0.7896`
    - nonpeak：`count_r=0.1211`、`jsd=0.5664`、`prof_r=0.5316`

## 当前判断

- 旧的 NT v2 sampled probe 已经说明它“有独立信号但未显出互补性”。
- 这次 residual-head short10 仍然是必要的下一步，因为它第一次真正回答：
  - 在不直接扰动主 backbone 的前提下，
  - 让 NT v2 只学习 debiased residual，
  - 是否仍然拿不到正收益。
- 当前这轮新增的稳定结论是：
  - cache / launcher / train / validation 这条最小闭环已经真正跑通
  - 首个 crash 不是 NT v2 指标层面的负结果，而是 DDP 与 validation plumbing bug
  - 两个 bug 修掉后，正式 rerun 已完整收口，说明“可训练”这道 gate 已不再是 blocker
- held-out gate 已经补齐，而且结果不是 `near-zero`，而是直接触发 `fail-or-unsafe`：
  - peak `profile_target_jsd_full_mean=0.3560`，显著高于停表阈值 `0.3197`
  - peak `count_pearson_full=0.7729`，显著低于停表阈值 `0.8400`
- 这说明当前 `FoundationResidualHead + NT v2 cached short10` 方案在 tutorial canonical 上不仅没有 clean gain，连 paper-facing 安全线都没有守住。
- 因此这轮新增的最终结论是：
  - “训练链路可跑通”与“foundation 路线值得继续”是两回事
  - restart v3 已经完成它应回答的问题；在默认规则下它应直接触发停表
  - 随后仅在显式授权下补的 matched control 与 `bins16 center-aligned residual` 也都没有扭转结论，因此当前 family 的 stop-rule 已进一步收紧
  - A6000 后续应回到 paper-aligned 的高价值工作，而不是继续烧在当前 foundation 方案上

## 下一步

- 把 `fail-or-unsafe` 结论写回 `TRACKING.md`、gate 文档和论文主线，不再把这轮 run 记作“待 held-out 判定”。
- foundation 线当前默认停表：不扩 `seed1234`、不切 `cross_attention`、不重开旧 `Genos/Caduceus` recipe，也不再把已完成的 matched control 或 `bins16` 重入当作可重复默认动作。
- 结合本轮工作目录整理，把 6000/6002 的 foundation 产物明确收口为“已判负的归档输出”，为后续 paper-aligned GPU 工作让路。
