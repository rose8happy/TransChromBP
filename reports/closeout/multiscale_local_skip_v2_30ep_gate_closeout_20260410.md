# A6000 Formal Gate Closeout: `teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1`（2026-04-10）

## Verdict

- `fail`
- 这是一次 **A6000 formal gate** 的正式判读，不是新模型晋级，也不是 benchmark 胜利。
- 该 run 只应被记录为“已完成但未通过 formal gate”。

## Run Identity

- run name：`teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1`
- machine / GPU：`6000 / A6000 x2`
- runtime repo：`/data1/zhoujiazhen/bylw_atac/TransChromBP`
- log：`/data1/zhoujiazhen/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1.log`
- output dir：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/logs/teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1`
- completion time：`2026-04-10 10:10:10 CST` 左右
- early stop：`true`
- stop epoch：`28`
- best epoch：`22`

## Evidence

### 1. 运行已完成，且 closeout 所需的结束标记齐全

- 日志与输出目录已经收口，没有活跃训练进程继续挂在这条 run 上。
- `run_meta.json` 给出的关键字段为：
  - `best_metric_name=peak.loss_total`
  - `best_metric_value=16.912663229837005`
  - `stopped_early=true`

### 2. 最优 epoch 与最终 epoch 都没有达到 corrected B comparator

本次 formal gate 只与历史 `corrected B` 比较。仓内已稳定引用的 two-seed mean comparator 为：

- peak `profile_target_jsd_full_mean` ≈ `0.3146`
- peak `count_pearson_full` ≈ `0.8496`

`best epoch=22` 的读数：

- peak `profile_target_jsd_full_mean=0.3337811803218466`
- peak `count_pearson_full=0.7948227405497477`
- peak `profile_full_debiased_jsd=0.00042317341300050147`
- overall `profile_target_jsd_full_mean=0.35636016641071955`
- overall `count_pearson_full=0.8012534669611618`

相对 `corrected B` two-seed mean，best epoch 的差值约为：

- peak JSD `+0.01918`
- peak count `-0.05478`

`final epoch=28` 的读数：

- peak `profile_target_jsd_full_mean=0.3330134000718902`
- peak `count_pearson_full=0.798037612304924`
- peak `profile_full_debiased_jsd=0.0004180776680376493`
- overall `profile_target_jsd_full_mean=0.3556691982690527`
- overall `count_pearson_full=0.8040587947255652`

终点虽然比 best epoch 略有回弹，但仍然明显落后于 corrected B comparator。

### 3. Gate thresholds 视角下的判读很清楚

按本轮设计/spec gate 的口径，只比较 corrected B 时：

- peak JSD 没有改进，反而变差了远超过 `0.003`
- peak count 也没有守住，下降幅度远超过 `0.01`

因此这不是边缘摇摆，也不是“差一点过线”，而是明确不过 gate。

## Decision

- 最终判定：`fail`
- 通过含义：没有通过 A6000 formal gate，不能进入 promotion / follow-up seed / benchmark 扩面。
- 唯一合理的后续口径：把它当作已完成但未通过的 formal gate 存档；如果未来再看这条方向，必须是新的显式 hypothesis，不是当前 run 的延续。

## Closeout Conclusion

`teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1` 已完成并正式判定为 `fail`。它没有达到 `corrected B` comparator，因此不应被写成性能胜利，也不应被自然升级为下一阶段 benchmark。对当前仓库而言，这条线到此收口。
