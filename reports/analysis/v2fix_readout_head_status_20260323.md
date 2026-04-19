# 2026-03-23 V2fix 读出头实验收口报告

> 更新时间：2026-03-23 22:35 CST
> 范围：`s42` 单 seed 的 `A/B/F/G` 四条 V2fix 结构实验
> 口径：validation 取训练期 `best_epoch` 的 `peak.profile_target_jsd_full_mean`；held-out test 统一为 `split=test`、`nonpeak_ratio=1.0`

## 1. 当前结果表

参考 baseline 取 `V2-full` 的历史 held-out test：

- `full_s42`: peak `JSD=0.31425`, `count_r=0.83536`
- `full(mean±std)`: peak `JSD=0.31419±0.00012`, `count_r=0.83676±0.00531`

| 方案 | run | best epoch | best val peak JSD | best val peak count_r | held-out test peak JSD | held-out test peak count_r | 当前判断 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| baseline | `full_s42` | 29 | - | - | 0.31425 | 0.83536 | 历史参考锚点 |
| A: freeze fix | `v2fix_20260320_freeze_s42` | 40 | 0.3318 | 0.8042 | 0.3148 | 0.8446 | held-out 已证实正向，但整体弱于 `B/F` |
| B: center pool | `v2fix_20260320_cpool_s42` | 34 | 0.3319 | 0.8079 | 0.3147 | 0.8503 | 当前单 seed 最优，下一轮优先级升高 |
| F: attention pool | `v2fix_20260321_attnpool_s42` | 40 | 0.3318 | 0.7991 | 0.3148 | 0.8492 | 仍是正向候选，但不再是唯一赢家 |
| G: profile refine | `v2fix_20260321_profref_s42` | 15 | 0.4425 | 0.8081 | 0.4260 | 0.8460 | count 不差，但 profile 明显落后，不建议继续 |

注：

- `A/B/F` 的 held-out peak `JSD` 全都卡在 `0.3147-0.3148`，相对 baseline 只差 `+0.00045 ~ +0.00055`，profile 主指标基本同档。
- 真正把三者拉开的是 held-out peak `count_r`：`B (0.8503) > F (0.8492) > A (0.8446) > baseline (0.8354)`。
- `B` 不只在 peak count 上领先；在当前单 seed held-out summary 里，`overall` 与 `nonpeak` 口径也都略优于 `A/F`。

## 2. 关键判断

### 2.1 `B: center pool` 现在是当前单 seed 最强候选

`B_s42` 的 held-out test peak 为 `0.3147 / 0.8503`。

相对 `full_s42` baseline：

- peak `JSD` 只差 `+0.00045`
- peak `count_r` 提升 `+0.01494`

相对 `F_s42`：

- peak `JSD` 略好 `-0.0001`
- peak `count_r` 略好 `+0.0011`

虽然这两个差距都很小，但方向是一致的，而且 `B` 在 `overall/nonpeak` summary 上也没有露出新的副作用。基于目前已经补齐的 held-out 证据，如果只能选一个“下一条新 seed”，`B` 应该排在最前。

### 2.2 `F: attention pool` 仍然成立，但“唯一正向方案”的判断被推翻了

`F_s42` 的 held-out test peak 为 `0.3148 / 0.8492`，依然明显优于 baseline 的 count 分支，而且 profile 没有实质退化。

需要修正的是结论口径：

- `F` 仍是正向候选
- 但它不再是唯一在 held-out 上给出正向信号的读出头改动
- 在当前单 seed 上，`B` 甚至略微压过了 `F`

不过 `F_s1234=attnpool` 已经在 6002 跑起来了，这条不应该中途打断。更合理的做法是先让它正常收口，再决定是否需要给 `B` 做 matched 多 seed。

### 2.3 `A: freeze fix` 通过了 held-out，但优先级应低于 `B/F`

`A_s42` 的 held-out test peak 为 `0.3148 / 0.8446`。

这说明 `A` 并没有像 `G` 那样在真实 held-out 上翻车；它也确实比 baseline 更好：

- peak `JSD` 仍与 baseline 基本同档
- peak `count_r` 提升 `+0.00924`

但和 `B/F` 相比，`A` 的增益幅度更小，当前没有证据支持再给它单独扩 seed。若算力预算紧，`A` 应该先降级到观察位，而不是继续抢第一批资源。

### 2.4 `G: profile refine` 可以停线

`G_s42` 的 held-out test peak 为 `0.4260 / 0.8460`。

虽然 count 相关性不差，但 profile JSD 相比 baseline 偏差过大，说明这条 profile refine 方案没有兑现为真实峰形收益。现在 `A/B` 的 held-out 也已经补齐，因此 `G` 更没有继续扩 seed 的必要。

### 2.5 validation 依然不够用，但 `B` 的 val 提示并非完全噪声

`A/B/F` 的 best val peak `JSD` 几乎完全挤在 `0.3318-0.3319`，单看 validation 主指标仍然分不出胜负。

但这轮 held-out 补齐之后，可以看到：

- `B` 的 val count 本来就略高于 `A/F`
- 最终 held-out 也是 `B` 最强
- `F` 的 val count 最低，但 held-out 仍然很好

这意味着当前判断应当是：

- validation 不能单独做 gate
- 但它也不是纯噪声，尤其 `B` 的轻微领先在 held-out 上保留下来了
- 因此后续如果要压缩实验数，最好直接比 matched held-out，而不是再看一轮 val 再决定

### 2.6 `40 epoch` 对所有头部方案并不等价

best epoch 分布：

- `A = 40`
- `B = 34`
- `F = 40`
- `G = 15`

这说明：

- `A/F` 的确吃到了 `30 -> 40 epoch` 的尾部收益
- `B` 在更早就已经收口，而且最后 held-out 仍然最好
- `G` 的问题不是训练不够久，而是结构方向本身不对

所以当前没有证据支持“所有头部方案都要统一继续补到 50 epoch”；更合理的是按结构分别判断。

## 3. 建议

如果继续占用 `6002` 做 sidecar：

1. 让已经在跑的 `F_s1234` 正常收口，不要中途打断。
2. 下一条新的 matched seed 预算，优先留给 `B`，而不是默认继续押 `F`。
3. 如果算力只够保留两条候选线，就收敛到 `B/F`，不要再给 `A` 新 seed。
4. `G` 可以准备迁出 live 清单，不再追加结构投入。

更简短地说：

- `B` 现在是当前单 seed 领跑者
- `F` 仍值得看完正在运行的多 seed
- `A` 是正向但次优
- `G` 应该停止

## 4. 关键文件

- baseline 参考表：`reports/assets/ablation_tf_20260318/summary_table.csv`
- `A_s42` test log：`/data1/zhoujiazhen/bylw_atac/logs/v2fix_20260320_freeze_s42_test_full_20260323_182611.log`
- `A_s42` test JSON：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/v2fix_20260320_freeze_s42_best_test_full_20260323_182611.json`
- `B_s42` test log：`/data1/zhoujiazhen/bylw_atac/logs/v2fix_20260320_cpool_s42_test_full_20260323_182611.log`
- `B_s42` test JSON：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/v2fix_20260320_cpool_s42_best_test_full_20260323_182611.json`
- `F_s42` test log：`/data1/zhoujiazhen/bylw_atac/logs/v2fix_20260321_attnpool_s42_test_full_20260323_180905.log`
- `F_s42` test JSON：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/v2fix_20260321_attnpool_s42_best_test_full_20260323_180905.json`
- `G_s42` test log：`/home/zhengwei/bylw_atac/logs/v2fix_20260321_profref_s42_test_full_20260323_180431.log`
- `G_s42` test JSON：`/home/zhengwei/bylw_atac/TransChromBP/outputs/metrics/v2fix_20260321_profref_s42_best_test_full_20260323_180431.json`
- `A_s42` epoch metrics：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/logs/v2fix_20260320_freeze_s42/epoch_metrics.jsonl`
- `B_s42` epoch metrics：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/logs/v2fix_20260320_cpool_s42/epoch_metrics.jsonl`
- `F_s42` epoch metrics：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/logs/v2fix_20260321_attnpool_s42/epoch_metrics.jsonl`
- `G_s42` epoch metrics：`/home/zhengwei/bylw_atac/TransChromBP/outputs/logs/v2fix_20260321_profref_s42/epoch_metrics.jsonl`
