# 当前项目审查（计划 + 代码，2026-04-05）

## 审查范围

- 状态与决策文档：`TRACKING.md`、`docs/plan/*.md`、`reports/foundation_model_restart_v3_implementation_20260405.md`
- 训练与启动代码：`vendor/transchrombp/transchrombp/training/train_ddp.py`、`vendor/transchrombp/transchrombp/scripts/run_ntv2_residual_short10.sh`、`vendor/transchrombp/transchrombp/scripts/run_short10_no_foundation_control.sh`
- 语法/静态检查：`python -m py_compile`、`bash -n`

## 处理进展（同日更新）

- `2026-04-05` 已落地：`run_ntv2_residual_short10.sh` 新增 `Step 3/4`，默认补齐 held-out `test` cache 并自动执行 `evaluate_checkpoint`，以输出 `test-full` gate JSON。
- `2026-04-05` 已落地：launcher 进一步补齐 `TRAIN_DRY_RUN_STEPS` 与 `HELDOUT_MAX_REGIONS`，现在可以做 bounded smoke，而不必总是跑完整训练 / 完整 held-out。
- `2026-04-05` 已落地：`build_foundation_cache.py` 现在会把 `dataset_split`、`region_source`、`nonpeak_ratio`、`max_records` 写入 manifest，并按 split 读取 `max_{split}_regions` / `{split}_region_source`，避免 test cache 复用判定过粗。
- `2026-04-05` 已落地：`docs/plan/post_chatgpt_pro_priority_execution_20260405.md` 与 `reports/foundation_model_restart_v3_implementation_20260405.md` 已补充“默认停表 vs 显式授权受控重入”的优先级说明。
- `2026-04-05 20:16 CST` 已在 6000 挂起等待式 smoke wrapper：当前双卡长跑结束后，会自动执行 `TRAIN_DRY_RUN_STEPS=2`、`HELDOUT_MAX_REGIONS=64` 的 bounded gate 验证。
- `2026-04-05 20:55-20:57 CST` 已新增最小回归测试 `tests/test_foundation_cache_alignment.py`：先在 6000 稳定复现 `test` split 下 cache builder / held-out evaluator 的 `record_sha1` 分叉，以及 `ChromBPNetBigWigDataset` 未暴露 `max_records` 导致 manifest 元数据恒为 `0`。
- `2026-04-05 20:57-20:58 CST` 已落地根因修复：`real_data.py` 新增共享的 `resolve_dataset_seed()`，并把 `max_records` 挂到 dataset 实例；`build_foundation_cache.py` 与 `evaluate_checkpoint.py` 改为共用该 helper。修复后 6000 上的回归测试已转绿。
- `2026-04-05 20:59 CST` 已在 6000 用新 run `ntv2_residual_short10_gate_smoke_20260405_fix1` 复跑 bounded smoke：`manifest_test.json` 已正确写出 `max_records=64`，held-out gate 成功产出 JSON，不再报 `record_sha1 mismatch`。
- `2026-04-05 21:02-21:05 CST` 已补跑 `short10_nofoundation_control_s42_20260405_dual` 的 full held-out：输出 `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/short10_nofoundation_control_s42_20260405_dual_best_test_full_20260405_2100.json`，核心结果为 `count_pearson_full=0.8457`、`count_pearson_debiased=0.8419`、`profile_target_jsd_full_mean=0.4395`、`profile_pearson_full_mean=0.7116`；日志 summary 里的分类指标为 `auroc=0.8677`、`auprc=0.8647`、`f1=0.7902`。
- `2026-04-05` 同轮结果收口：已读取 `ntv2_residual_short10_gate_smoke_20260405_fix1` 的 held-out JSON；在 `n=64` 的 bounded gate 上，其 `count_pearson_full=-0.0980`、`profile_pearson_full_mean=-0.0142`、`profile_target_jsd_full_mean=0.6668`、`peak_auroc=0.5543`。与 matched no-foundation control 对照后，可以把当前 `FoundationResidualHead + bins4/coarse-summary residual` family 明确收口为 `unsafe / not-worth-expanding`。
- `2026-04-06 03:25 CST` 已完成用户显式授权的 genuinely new-hypothesis 重入 `ntv2_bins16_centerres_short10_s42_20260406_dual`：full held-out 为 overall `count_pearson_full=0.7727`、`profile_target_jsd_full_mean=0.4625`、`profile_pearson_full_mean=0.6580`，peak `count_pearson_full=0.7516`、`profile_target_jsd_full_mean=0.3588`、`profile_pearson_full_mean=0.7837`，仍明显落后于 matched no-foundation control。
- 当前剩余：不再是补链路，也不是继续替当前 family 找借口，而是把“链路修复已完成 + 旧 bins4 residual 判负 + bins16 center-aligned residual 仍判负”统一写进 gate 文档和论文口径。

## 发现（按严重度）

### [高][已落地修复] `run_ntv2_residual_short10.sh` 曾未覆盖 held-out gate 所需的 `test` foundation cache

- 证据：
  - launcher 仅构建 `train valid`（`--splits train valid`）：`vendor/transchrombp/transchrombp/scripts/run_ntv2_residual_short10.sh:179` 与 `:193`
  - 当前计划文档已记录过一次真实故障：缺 `manifest_test.json` 导致 held-out `test-full` 失败：`docs/plan/foundation_model_restart_v3_20260405.md:126`
- 影响：
  - 重跑同脚本时，训练可完成但 held-out gate 仍可能再次失败；
  - gate 闭环依赖人工 sidecar，降低可复现性与可操作性。
- 建议：
  - 把 `test` 纳入 launcher 的可配置 cache split（至少在需要 `test-full` gate 时自动包含）；
  - 或在 launcher 末尾显式检测 `manifest_test.json` 并给出阻断式提示，避免“训练成功但 gate 缺口”。

### [高][部分缓解] foundation 停表规则存在多文档口径冲突

- 证据：
  - `post_chatgpt_pro_priority_execution` 定义：`fail-or-unsafe` 后应停表，不继续 foundation：`docs/plan/post_chatgpt_pro_priority_execution_20260405.md:77-90`
  - `restart_v3 implementation` 也写明“不跑 matched no-foundation control”：`reports/foundation_model_restart_v3_implementation_20260405.md:173`
  - 但 `TRACKING` 显示已在 fail 后按“新重入规则”启动 matched control，并把上述两份文档同时列为关键路径：`TRACKING.md:18-19`
- 影响：
  - 同一任务可被两套规则解释，后续执行人可能按不同文档做出相反动作。
- 建议：
  - 指定单一权威决策文档（建议直接在 `TRACKING.md` 首行声明）；
  - 在被替代文档顶部加 `Superseded by ... (date)`，并去掉旧规则中的“固定动作”表述，避免再次冲突。

### [高][结果层面已收口] 当前 `bins4/coarse-summary residual` family 不再值得继续扩线

- 证据：
  - matched `short10 no-foundation control` 的 full held-out 仍然健康：`count_pearson_full=0.8457`、`profile_pearson_full_mean=0.7116`、`profile_target_jsd_full_mean=0.4395`
  - 对照之下，`ntv2_residual_short10_gate_smoke_20260405_fix1` 的 bounded gate 在 `n=64` 上仅有 `count_pearson_full=-0.0980`、`profile_pearson_full_mean=-0.0142`、`profile_target_jsd_full_mean=0.6668`、`peak_auroc=0.5543`
- 影响：
  - 不能再把当前失败主要归因于 `short10` 预算本身；
  - 若继续为旧 `bins4` residual 配方分配 A6000，只会扩大负结果而不是增加信息增益。
- 建议：
  - 默认停止当前 `FoundationResidualHead + layer_07__bins4_mean + feature_tokens=4 + profile_bin_count=16` family；
  - 若用户显式授权继续 foundation，只允许进入 `bins16 center-aligned / higher-bandwidth residual` 这种 genuinely new-hypothesis 路线。

### [中][已落地修复] 当前 gate 流程曾缺少脚本级“一键闭环”保障

- 证据：
  - `run_ntv2_residual_short10.sh` 当前流程止于训练结束，没有把 held-out `test-full` gate 串成标准步骤：`vendor/transchrombp/transchrombp/scripts/run_ntv2_residual_short10.sh:219-235`
  - 本轮实际执行依赖额外 sidecar 脚本补齐 held-out。
- 影响：
  - 运行成功与“可判定”并不等价，容易出现“有 checkpoint 但无最终 gate 结论”。
- 建议：
  - 增加 `--run-heldout-gate` 开关，训练结束后自动执行 `test-full` 评估并产出统一 JSON；
  - 至少在日志中强提醒“未跑 held-out gate 不可用于 stop-rule 判定”。

## 本轮已做检查

- `python -m py_compile vendor/transchrombp/transchrombp/training/train_ddp.py`：通过
- `bash -n vendor/transchrombp/transchrombp/scripts/run_ntv2_residual_short10.sh`：通过
- `bash -n vendor/transchrombp/transchrombp/scripts/run_short10_no_foundation_control.sh`：通过
- 6000 / `python -m pytest tests/test_foundation_cache_alignment.py -q`：先红灯稳定复现，再在修复后转绿（`2 passed`）
- 6000 / bounded smoke `ntv2_residual_short10_gate_smoke_20260405_fix1`：通过，成功写出 `manifest_test.json` 与 held-out JSON `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/ntv2_residual_short10_gate_smoke_20260405_fix1_best_test_full_20260405_205905.json`
- 6000 / `short10_nofoundation_control_s42_20260405_dual` full held-out：通过，日志 `/data1/zhoujiazhen/bylw_atac/logs/short10_nofoundation_control_testfull_20260405_6000.log` 已收口

## 总结

这轮代码风险已经从“闭环链路本身不通”收敛到“当前结果足不足以支持继续扩线”。`NT v2` smoke 的根因不是模型效果，而是 `test` split seed 规则在 cache builder 与 held-out evaluator 间不一致；修复后闭环已被真实 run 证实可用。结果层面，matched no-foundation baseline 已证明 `short10` 预算本身仍然健康，而后续的两条 NT v2 residual gate 也都给出了负面 stop/go 信号：旧 `bins4/coarse-summary residual` 在 bounded smoke 下已是 `unsafe / not-worth-expanding`，新的 `bins16 center-aligned residual` full held-out 也仍明显落后于 matched baseline。当前没有继续给 residual short10 family 分配 A6000 的依据。
