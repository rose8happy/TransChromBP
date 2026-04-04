# 6002 独立真实数据线总结（2026-03-23）

## 目标

用 6002 单卡 RTX 3080 在两套独立真实数据上做最小闭环验证，回答两个问题：

1. `GM12878/K562` 这条与 Genos 无关的真实数据线是否能稳定训练。
2. validation 指标是否能在 held-out test 上站住，而不是只在单机单 split 上好看。

评估口径统一为：

- checkpoint：各 run 的 `best.pt`
- split：`test`
- nonpeak_ratio：`1.0`
- 脚本：`python -m transchrombp.evaluation.evaluate_checkpoint`

## 结果总表

| run | 状态 | best epoch | best val peak JSD | best val peak count_r | last val peak JSD | last val peak count_r | test peak JSD | test peak count_r | 判断 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `gm12878_v2_6002_smoke_retry_20260322_003645` | finished | 3 | 0.4351 | 0.7383 | 0.4371 | 0.7970 | - | - | smoke 证明链路可跑，但只够做起点，不够收口 |
| `gm12878_v2_6002_pilot20_20260322_191746` | finished | 19 | 0.4302 | 0.7998 | 0.4312 | 0.8008 | 0.4226 | 0.8040 | held-out 比 val 略强，`GM12878` 线稳定成立 |
| `k562_v2_6002_smoke_20260322_161602` | finished | 9 | 0.6331 | 0.8396 | 0.6354 | 0.8411 | - | - | smoke 已接近平台，但没有 held-out，不能单独下结论 |
| `k562_v2_6002_single_20260323_060026` | finished | 35 | 0.6281 | 0.8604 | 0.6289 | 0.8587 | 0.6124 | 0.8586 | held-out 与 val 一致，`K562` 线稳定成立 |

## 关键判断

- `GM12878 pilot20` 已经完成这条线真正需要的验证。它不是只在 validation 上小幅改善，而是在 held-out test 上也保持甚至略优于 validation，说明没有出现明显的 split-specific 失真。
- `K562 single` 的 held-out 结果同样没有反转。`test peak JSD=0.6124` 比 best val `0.6281` 更低，`test peak count_r=0.8586` 也基本贴住 best val `0.8604`，说明这条线在第二个真实数据集上也站住了。
- 因为 `GM12878` 和 `K562` 两套 held-out 都已补齐，6002 这条线已经完成了“独立 sanity check”的职责。它现在不是缺更多训练，而是缺收口和归档。
- `GM12878` 的主要收益仍然集中在 count 分支校准；`K562` 则更像是在 smoke 已经较强的前提下，把稳定性和最终 held-out 口径补齐。

## 当前建议

1. 把 6002 独立线视为阶段性完成，准备迁入 `TRACKING_archive.md`。
2. 6002 暂不继续追加新结构或新口径实验，避免和 6000 主线比较口径混乱。
3. 如果后续确实还要占用 6002，最合理的任务不是新结构，而是补一个 `GM12878 pilot20` 的新 seed，用来估计单卡真实数据线的方差。

## 额外验证

`2026-03-23 17:48 CST` 已在 6002 用当前部署的 `transchrombp` 环境对 `GM12878/K562` 执行 `build_eval_dataset(valid/test, max_regions=32)` loader smoke：

- 四个 dataset 均可正常构造
- 样本键为 `seq / profile_counts / region_source`
- 形状为 `seq=(2114, 4)`、`profile_counts=(1000,)`

这意味着 6002 上的 `prep_v1` 产物不只是“文件存在”，而是已被当前 loader 实际消费。

## 证据路径

- `GM12878 smoke`：`/home/zhengwei/bylw_atac/logs/gm12878_v2_6002_smoke_retry_20260322_003645.log`
- `GM12878 pilot20`：`/home/zhengwei/bylw_atac/logs/gm12878_v2_6002_pilot20_20260322_191746.log`
- `GM12878 pilot20 test`：`/home/zhengwei/bylw_atac/logs/gm12878_v2_6002_pilot20_20260322_191746_test_full_20260323_165903.log`
- `GM12878 pilot20 test json`：`/home/zhengwei/bylw_atac/TransChromBP/outputs/metrics/gm12878_v2_6002_pilot20_20260322_191746_best_test_full_20260323_165903.json`
- `K562 smoke`：`/home/zhengwei/bylw_atac/logs/k562_v2_6002_smoke_20260322_161602.log`
- `K562 single`：`/home/zhengwei/bylw_atac/logs/k562_v2_6002_single_20260323_060026.log`
- `K562 single test`：`/home/zhengwei/bylw_atac/logs/k562_v2_6002_single_20260323_060026_test_full_20260323_172622.log`
- `K562 single test json`：`/home/zhengwei/bylw_atac/TransChromBP/outputs/metrics/k562_v2_6002_single_20260323_060026_best_test_full_20260323_172622.json`
