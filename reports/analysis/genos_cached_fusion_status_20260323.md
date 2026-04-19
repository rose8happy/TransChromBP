# Genos Cached-Fusion 状态总结（2026-03-23）

## 目标

整理 `G0/G1/G2` 与 `cached-fusion P0/P2` 的当前状态，明确：

1. 哪些路线已经可以停。
2. 哪些路线还能继续。
3. 为什么后续必须坚持 matched `2-rank DDP`，不能再用单卡 fallback 混入口径。

## 判断原则

- gate 以 held-out test 为准，不以 validation 单独定生死。
- `P0/P2/P1` 必须保持 matched recipe：
  - `nproc_per_node=2`
  - `batch_size_per_gpu=10`
  - `global_batch=20`
- 训练拓扑不能成为实验变量。因此 `P2` 曾经的单卡 fallback 已中止，不纳入任何比较。

## 已完成结果

| run | 类别 | val peak JSD | val peak count_r | test peak JSD | test peak count_r | 结论 |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| `genos_20260321_baseline_s42` | online Genos `G0` baseline | 0.3337 | 0.7998 | 0.3163 | 0.8410 | 现阶段最稳 baseline，val 与 test 一致且 test 更强 |
| `genos_20260321_gate_s42` | online Genos `G1` gate | 0.3360 | 0.8010 | 0.3388 | 0.4689 | 明显泛化失败，validation 误导性很强 |
| `genos_20260322_mean_s42` | online Genos `G2` mean | 0.3414 | 0.7955 | 0.3387 | 0.6360 | 比 `G1` 略好，但仍远低于 baseline |
| `genos_cached_P0_baseline_ddp2_20260323_121650` | cached-fusion `P0` matched baseline | 0.3353 | 0.8033 | 0.3179 | 0.8384 | 路线稳定，已通过“可以继续放大”的 gate，但尚未超过 `G0` |

## 最新状态

| run | 类别 | 当前状态 | 最新已知指标 | 判断 |
| --- | --- | --- | --- | --- |
| `genos_cached_P2_count_only_ddp2_20260323_173318` | cached-fusion `P2` count-only | held-out test finished | val best 为 `epoch 9` 的 `0.3344 / 0.7984`；held-out peak 为 `0.3174 / 0.6037` | profile 基本追回，但 count 分支显著崩塌；不通过继续推进 `P1` 的 gate |

## probe 与路线选择

前置 probe 的结论没有变：

- `genos_global_mean -> count residual` 有小但稳定的解释力：`R^2=0.0167`，`Pearson r=0.1987`
- `genos_only` 对 peak / nonpeak 有中等区分度：`AUC=0.6997`
- 但与 encoded 直接拼接后 `AUC_concat=0.9018`，反而低于 `encoded_only=0.9219`

这意味着：

- Genos summary 里确实有独立信号
- 但“有信号”不等于“直接拼接就一定有收益”
- 所以后续必须靠 `P0 -> P2 -> P1` 的受控实验，而不是再相信 validation 上的直觉

## 关键判断

- online Genos `G1/G2` 这条线已经可以停止。它们最大的问题不是 JSD 一定坏，而是 validation 与 held-out test 的关系不可信。
- cached-fusion `P0` 已经证明这条路线至少是稳定可复现的，而且 held-out test 没有崩。这就是它可以继续推进 `P2` 的理由。
- 但截至 `P0`，cached-fusion 还不能宣称优于当前最好 baseline：`P0 test = 0.3179 / 0.8384`，而 `G0 test = 0.3163 / 0.8410`。差距很小，但方向上仍是 baseline 略优。
- `P2` 的 held-out test 已经给出明确负面结果：虽然 peak `JSD=0.3174` 与 `P0/G0` 相距不大，但 peak `count_r=0.6037`，相比 `P0=0.8384` 和 `G0=0.8410` 出现了不可接受的断崖式下滑。
- 这说明 `count-only` 这条 cached-fusion 接法没有把 probe 里的“count residual 有信号”转化成可用的下游收益，反而破坏了主任务的 count 读出。
- 因此，`P2` 不只是“没有超过 `P0/G0`”，而是已经足够差到不该再用它作为继续烧 `P1` 的理由。

## 关于训练拓扑

`2026-03-23` 当天曾一度误判 `GPU1` 被持续打满，尝试过 `GPU0` 单卡 fallback 的 `P2`。这个尝试已经中止，不纳入结果。

原因很简单：

- 单卡 fallback 会改变训练拓扑
- 训练拓扑变化会改变实验语义
- 与 `P0` 的 matched `2-rank DDP` 比较就不再干净

因此后续规则应固定为：

- `P0/P2/P1` 主线一律按 matched `2-rank DDP`
- 如果资源紧张，宁可等待或并发调度，也不要偷偷切成单卡口径

## 当前建议

1. 不再推进 `P1 global_late_film`，当前 cached-fusion 主线到 `P2` 为止应停止。
2. 把 `P0/P2` 与 `G0/G1/G2` 的结论归档，作为“Genos cached summary 在当前 recipe 下未形成正收益”的阶段性结论。
3. 只有在后续明确更换问题设定或融合位置时，才值得重新开新的 Genos 主线；不要沿着 `P0 -> P2 -> P1` 继续顺推。

## 证据路径

- 总结性前置分析：`reports/analysis/genos_and_6002_run_analysis_20260323.md`
- `G0` 训练日志：`/data1/zhoujiazhen/bylw_atac/logs/genos_20260321_baseline_s42.log`
- `G1` 训练日志：`/data1/zhoujiazhen/bylw_atac/logs/genos_20260321_gate_s42.log`
- `G2` 训练日志：`/data1/zhoujiazhen/bylw_atac/logs/genos_20260322_mean_s42.log`
- `P0` 训练日志：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/logs/genos_cached_P0_baseline_ddp2_20260323_121650.log`
- `P0` held-out test 日志：`/data1/zhoujiazhen/bylw_atac/logs/genos_cached_P0_baseline_ddp2_20260323_121650_test_full_20260323_170352.log`
- `P0` held-out test json：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/genos_cached_P0_baseline_ddp2_20260323_121650_best_test_full_20260323_170352.json`
- `P2` 训练日志：`/data1/zhoujiazhen/bylw_atac/logs/genos_cached_P2_count_only_ddp2_20260323_173318.log`
- `P2` held-out test 日志：`/data1/zhoujiazhen/bylw_atac/logs/genos_cached_P2_count_only_ddp2_20260323_173318_test_full_20260323_224747.log`
- `P2` held-out test json：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/genos_cached_P2_count_only_ddp2_20260323_173318_best_test_full_20260323_224747.json`
