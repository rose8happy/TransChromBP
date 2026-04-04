# 2026-03-23 实验收口分析

## 1. 6000 Online Genos Phase 1

held-out test 评估口径：

- split=`test`
- nonpeak_ratio=`1.0`
- 评估脚本：`python -m transchrombp.evaluation.evaluate_checkpoint`

| run | best epoch | val peak JSD | val peak count_r | test peak JSD | test peak count_r | test peak count_r_debiased | 结论 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `genos_20260321_baseline_s42` | 20 | 0.3337 | 0.7998 | 0.3163 | 0.8410 | 0.8396 | test 比 val 更强，基线稳定 |
| `genos_20260321_gate_s42` | 10 | 0.3360 | 0.8010 | 0.3388 | 0.4689 | 0.4127 | val 看似接近 baseline，但 test count 明显塌陷 |
| `genos_20260322_mean_s42` | 5 | 0.3414 | 0.7955 | 0.3387 | 0.6360 | 0.5612 | 比 `G1` 稍好，但仍明显落后 baseline |

关键判断：

- `G0 baseline` 仍是唯一同时在 val / test 上都稳定的方案。
- `G1 gate` 的 held-out test 明显失真：`count_r` 相比 val 掉了 `0.3321`，说明 validation 已不足以作为 online Genos recipe 的 gate。
- `G2 mean` 也没有兑现 validation 上的接近度；虽然 `JSD` 与 val 接近，但 `count_r` 仍比 baseline 低很多。
- `G1/G2` 在 test 上的 full / debiased count 相关性差距都显著放大，提示它们对 bias 分支或 split-specific 统计存在更强依赖。

当前结论：

- 暂停继续放大 online Genos `G1/G2` 训练。
- 后续 Genos 方向只保留 cached-fusion 主线：`summary cache -> probe -> P0/P1/P2`。

## 2. 6002 单卡真实数据线

当前只有 validation 级分析；held-out test 评估尚未补跑。

| run | status | best epoch | best peak JSD | best peak count_r | 备注 |
| --- | --- | ---: | ---: | ---: | --- |
| `gm12878_v2_6002_smoke_retry_20260322_003645` | finished | 3 | 0.4351 | 0.7383 | epoch 7 early-stop |
| `gm12878_v2_6002_pilot20_20260322_191746` | finished | 19 | 0.4302 | 0.7998 | 相比 smoke，JSD 仅小幅改善，但 count_r 提升明显；曲线在 epoch 19-20 附近基本饱和 |
| `k562_v2_6002_smoke_20260322_161602` | finished | 9 | 0.6331 | 0.8396 | epoch 10 末略有回摆 |
| `k562_v2_6002_single_20260323_060026` | running | 10 so far | 0.6339 | 0.8390 | `epoch 11` validation 回到 `0.6352 / 0.8361`，当前更像接近平台而非持续上升 |

关键判断：

- `GM12878 pilot20` 的增益主要体现在 count 分支校准，而不是 profile JSD 的大幅下降。
- `K562 single` 目前没有显示出明显超越昨日 smoke best 的趋势；至少到 `epoch 11`，它更像是在最优点附近震荡。
- 在 6002 当前长跑结束前，不建议再为这条单卡线追加新的长训练。

## 3. 建议动作

1. 保留 6002 当前 `K562 single`，等更多 validation 点或 early-stop 再决定是否补 test。
2. `GM12878 pilot20` 的 held-out test 评估等 6002 GPU0 空出来后再补。
3. 6000 资源不要再投给 online Genos `G1/G2`；把后续试验集中到 cached-fusion。

