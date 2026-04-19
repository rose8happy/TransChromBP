# 外发复用前复审记录（2026-04-07）

## 1. 目的

在继续复用 `reports/external/chatgpt_bundle_project_handoff_20260405/` 之前，再做一轮独立核实，避免把“文档互相转述”误当成证据闭环。

本轮复审覆盖三类证据：

1. 本地归档代码是否真的包含文档所述的关键实现；
2. 6000 远端实际运行仓与实验产物是否支持当前 stop-rule；
3. 本地归档仓与 6000 实际运行仓是否存在会影响复用判断的偏差。

## 2. 结论

### 2.1 总结

当前**没有发现会推翻“复用现有外发包”判断的阻塞问题**。

复用判断成立，原因是：

- 本地归档代码与 6000 实际运行代码在本轮核验的 4 个关键文件上哈希一致；
- 6000 真实运行产物、日志与 `TRACKING.md` / 计划文档中的核心数字一致；
- foundation 相关的关键补丁（`center_token_bins` residual、cache/eval split seed 对齐、DDP `find_unused_parameters` 自动兜底）都能被代码和远端自检独立支撑。

### 2.2 本轮唯一需要显式提醒的风险

外发或后续内部复核时，**不要把本地归档仓的路径/环境假定直接套到 6000 运行仓上**：

- 本地归档快照主路径是 `vendor/transchrombp/...`
- 6000 实际运行仓主路径是 `src/transchrombp/...`
- 6000 上真正可跑 foundation 训练/自检的是 `.venvs/genos-1.2b/bin/python`
- `chrombpnet` 环境适合官方 ChromBPNet / `modisco` / `bedGraphToBigWig` 等工具，但**不含 `torch`，也不适合直接拿来跑 foundation 训练代码核验**

这不是结果层面的 blocker，但如果不写清楚，下一轮审计很容易在“代码路径不对 / 环境不对”上浪费时间。

## 3. 代码核验

### 3.1 本地归档仓

本地归档代码与当前文档叙述一致：

- `FoundationResidualHead` 已支持 `alignment_mode=center_token_bins`，且在 token 对齐路径下要求 `foundation_tokens`、按中心区域池化本地 token，再预测局部 residual bins 和 count residual。
- `evaluate_checkpoint.py` 与 `build_foundation_cache.py` 都已使用同一套 `resolve_dataset_seed(...)` 与 split/max-records 解析逻辑。
- `train_ddp.py` 已在 `foundation residual_head + freeze_backbone_until_epoch>0` 时自动开启 `find_unused_parameters=true`。
- bins16 的模型配置也明确写成：
  - `feature_name: layer_07__bins16_mean`
  - `feature_layout: token`
  - `alignment_mode: center_token_bins`
  - `aligned_token_count: 16`

### 3.2 本地归档仓与 6000 运行仓的一致性

本轮直接比对了 4 个关键文件的 sha1：

- `foundation_adapter.py`
- `evaluate_checkpoint.py`
- `build_foundation_cache.py`
- `train_ddp.py`

结果：**本地归档仓与 6000 运行仓哈希完全一致**。

因此，本地 `vendor/` 中用于外发/审查的关键 foundation 代码快照，和 6000 实际跑出结论的代码不是“两套不同实现”。

## 4. 远端运行环境与最小自检

### 4.1 运行环境事实

本轮在 6000 实测得到：

- `/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin/python`
  - 有 `pyBigWig`、`pyfaidx`
  - **没有** `torch`
  - **没有** `yaml`
- `/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b/bin/python`
  - 有 `torch`
  - 有 `pyBigWig`
  - 有 `pyfaidx`
  - 有 `yaml`

因此 foundation 相关核验以 `.venvs/genos-1.2b` 为准。

### 4.2 最小自检结果

在 6000 实际运行仓 `/data1/zhoujiazhen/bylw_atac/TransChromBP` 中，用 `.venvs/genos-1.2b/bin/python` 做了等价自检：

1. toy dataset 上验证 `build_dataset_for_cache('test', ...)` 与 `build_eval_dataset(..., split='test')` 产出的记录完全一致；
2. 验证 `ChromBPNetBigWigDataset.max_records` 会正确暴露；
3. 验证 `FoundationResidualHead(center_token_bins)` 的输出 shape 为 `(B,1000)` 和 `(B,1)`，且零初始化下初始 residual 为 0；
4. 对关键 Python 文件跑 `py_compile`；
5. 对 `scripts/run_ntv2_residual_short10.sh` 跑 `bash -n`。

结果：**全部通过**。

## 5. 实验产物核验

### 5.1 `bins16 center-aligned residual` 与 baseline 的真实对照

本轮直接读取 6000 上的 held-out JSON 与日志，得到：

- baseline `short10_nofoundation_control`
  - overall `count_pearson_full=0.8457`
  - overall `profile_target_jsd_full_mean=0.4395`
  - peak `count_pearson_full=0.8298`
  - peak `profile_target_jsd_full_mean=0.3193`
- `bins4` residual
  - overall `count_pearson_full=0.6800`
  - overall `profile_target_jsd_full_mean=0.4612`
  - peak `count_pearson_full=0.7729`
  - peak `profile_target_jsd_full_mean=0.3560`
- `bins16 center-aligned residual`
  - overall `count_pearson_full=0.7727`
  - overall `profile_target_jsd_full_mean=0.4625`
  - peak `count_pearson_full=0.7516`
  - peak `profile_target_jsd_full_mean=0.3588`

相对 baseline：

- `bins4`：
  - overall count `-0.1657`
  - overall JSD `+0.0217`
  - peak count `-0.0569`
  - peak JSD `+0.0367`
- `bins16`：
  - overall count `-0.0730`
  - overall JSD `+0.0230`
  - peak count `-0.0782`
  - peak JSD `+0.0395`

这说明：

> `bins16 center-aligned residual` 虽然比旧 `bins4` 路线更像“真正的新假设”，但它并没有把结果拉回 matched baseline 之上，反而继续落后。

### 5.2 日志与 run_meta / manifest 的一致性

本轮还核对了：

- `ntv2_bins16_centerres_short10_20260406_6000.log`
- `short10_nofoundation_control_testfull_20260405_6000.log`
- `outputs/logs/ntv2_bins16_centerres_short10_s42_20260406_dual/run_meta.json`
- `outputs/foundation_cache/ntv2_tutorial_canonical_bins16_v1/manifest_test.json`

核对结果：

- held-out log 尾部 summary 数字与 JSON 一致；
- `run_meta.json` 明确记录了 bins16 模型配置、world size、best epoch、best checkpoint；
- `manifest_test.json` 明确记录了
  - `dataset_split=test`
  - `feature_types=["bins16_mean"]`
  - `features=["layer_07__bins16_mean"]`
  - `n_records=62342`
  - `nonpeak_ratio=1.0`
  - `region_source=both`

因此，当前 `TRACKING.md` 中关于 bins16 gate 判负的高层结论，确实有真实远端产物支撑。

## 6. 对外发复用的影响

当前最稳妥的做法仍然是：

1. 继续复用 `reports/external/chatgpt_bundle_project_handoff_20260405/`；
2. 新增一个 `post-20260406 delta`，只补：
   - `bins16 center-aligned residual` 也判负；
   - foundation 线已从“保留一条窄门”收紧为“当前 measured family 默认停表”；
   - 当前工作重心已经切回 `paper-facing` 收口；
3. 如需在 bundle 里补运行说明，只写一句：
   - “6000 实际运行仓为 `src/transchrombp`，foundation 相关复核应使用 `.venvs/genos-1.2b`，不要机械沿用本地归档仓的 `vendor/` 路径或 `chrombpnet` 环境。”

## 7. 一句话结论

本轮复审支持继续复用现有外发包；当前需要补的是 `post-20260406` 的结果增量和运行仓/环境说明，而不是重做一整套研究路线裁决材料。
