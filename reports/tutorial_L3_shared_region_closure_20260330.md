# Tutorial L3 Shared-Region 收口报告（2026-03-30）

## 结论

tutorial `Layer 3` 的 shared-region system compare 现已收口。结论明确：

- `TransChromBP` 在 held-out test 上同时优于 official `ChromBPNet`
- 优势同时出现在 `profile` 与 `count` 两侧
- 本轮结果可以作为 `Layer 2` 结论的 shared-region robustness check
- 这仍然是 `shared-region system comparison`，不是 `architecture-only` 归因

## 主表：Held-Out Test 主指标

| Arm | best epoch | peak mean JSD | peak median JSD | peak count_r | overall count_r | 产物 |
|---|---:|---:|---:|---:|---:|---|
| Official ChromBPNet controlled L3 | 29 | 0.33853 | 0.33989 | 0.69958 | 0.72215 | `/data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/runs/tutorial_official_controlled_s42_L3/fold_0/seed_42/chrombpnet/external_best_test_clsfix_20260330_2200/chrombpnet.epoch_029_metrics.json` |
| TransChromBP corrected-B controlled L3 | 47 | 0.31319 | 0.31547 | 0.84016 | 0.85468 | `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/tutorial_corrected_b_strict_compare_L3_s42_ep047_test_full_20260330_214519.json` |

按 peak 指标计算：

- `mean JSD` 改善 `0.02534`
- `median JSD` 改善 `0.02442`
- `count_r` 提升 `0.14058`

这说明 shared-region 约束并没有抹掉优势；相反，在相同 peaks / nonpeaks / bigWig 的设定下，`TransChromBP` 仍稳定领先。

## 附表：ATAC 分类指标

本轮按“先放附表、不并入主表”的口径处理。

| Arm | AUROC | AUPRC | best F1 | 备注 |
|---|---:|---:|---:|---|
| Official ChromBPNet controlled L3 | 0.82295 | 0.83148 | 0.75248 | 来自 `classification_metrics.peaks_vs_nonpeaks.logcounts`；该 run 中 count/logcount 数值一致，仅 threshold 不同 |
| TransChromBP corrected-B controlled L3 | 0.87861 | 0.87530 | 0.80991 | 来自 `results.overall.peak_{auroc,auprc,best_f1}_count_full` |

补充说明：

- `TransChromBP` 当前 production JSON 已稳定带出分类指标
- official 侧本轮最初缺失分类指标，根因不是命令参数，而是 6000 远端 `chrombpnet/training/predict.py` 与 `metrics.py` 仍是旧版
- 2026-03-30 22:00 CST 左右已将本地带 `classification_metrics` 的实现同步到 6000，随后对 official `epoch 29` 做了单 checkpoint 最小重评估

## 复核与一致性

为避免把“补字段”误变成“改结果”，本轮对 official `epoch 29` 做了新旧对比：

| 指标 | 旧 test JSON | 新 clsfix JSON | 差值 |
|---|---:|---:|---:|
| peak mean JSD | 0.33853003 | 0.33853003 | 0.00000000 |
| peak median JSD | 0.33989056 | 0.33989056 | 0.00000000 |
| peak count_r | 0.69958020 | 0.69958020 | 0.00000000 |

因此本轮修复只补上了 official 分类指标字段，没有改变原先已经成立的 profile/count 结论。

## 对论文与主线的影响

- `Layer 3` 现在可以写成已完成的 `shared-region robustness check`
- 主文 strict-compare 主表仍建议保留 `profile/count`
- `AUROC/AUPRC/F1` 可放 supplementary / 附表，不必挤进主表
- 下一步不再是继续补 strict-compare 技术链路，而是把这组数字并入论文主表、claim matrix 与结果叙事
