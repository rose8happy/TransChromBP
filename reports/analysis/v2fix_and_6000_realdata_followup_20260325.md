# 2026-03-25 最新收口实验分析：V2fix `B_s1234` 与 6000 独立数据线

> 更新时间：2026-03-25 00:55 CST
> 范围：`v2fix_20260324_005834_cpool_s1234_gpu0`、`v2fix_20260323_attnpool_s1234`、`k562_v2_6000_single_20260324_184918`、`gm12878_v2_6000_pilot20_20260323_235205`，以及对应的 6002 参考 run
> 直接结论：P0 held-out 评估已补齐；`B_s1234` 在 held-out 上确认优于 `F_s1234`，`K562_6000` 的峰区弱势也已在 held-out 上复现

## 1. 直接回答

### 1.1 是否还需要跑评估

主线所需的评估已经补齐：

1. `B_s1234` 的 held-out test：已完成
2. `F_s1234` 的 held-out test：已完成
3. `K562_6000` 的 held-out test：已完成

因此，当前已经不需要为了主线判断再补新的 validation，也不需要继续开新的 matched seed。

现在只剩次优先级评估：

1. `GM12878_6000` 的 held-out test
2. `GM12878_6002_s1234` 的 held-out test

它们主要用于估计 seed / 机器 / batch 差异，不再是当前 `v2fix` 主线的必要 gate。

### 1.2 现在不值得优先做什么

- 不需要再补新的 validation 评估；当前 validation 曲线信息已经够了。
- 不建议继续给 `B/F` 追加新 seed，再用训练替代已经拿到的 held-out 证据。
- 不建议立刻重跑 `K562_6000`；这次偏弱已在 held-out 上复现，下一步应先决定要不要做 batch / epoch 诊断。

## 2. 当前结果表

| run | 状态 | best epoch | best val peak JSD | best val peak count_r | 当前判断 |
| --- | --- | ---: | ---: | ---: | --- |
| `v2fix_20260324_005834_cpool_s1234_gpu0` | finished + tested | 39 | 0.33185 | 0.79763 | held-out peak `0.3145 / 0.8488`，与 `B_s42` 的 `0.3147 / 0.8503` 同档 |
| `v2fix_20260323_attnpool_s1234` | early-stop + tested | 19 | 0.44338 | 0.81099 | held-out peak `0.4269 / 0.8457`，profile 端仍明显落后，不能作为默认方案 |
| `k562_v2_6000_single_20260324_184918` | early-stop + tested | 10 | 0.63365 | 0.82014 | held-out peak `0.6186 / 0.8287`，峰区弱势不是纯 val 偏差 |
| `k562_v2_6002_single_20260323_060026` | finished + tested | 35 | 0.62809 | 0.86031 | held-out 为 `0.6124 / 0.8586`，说明 K562 这条线本身是稳定的 |
| `gm12878_v2_6000_pilot20_20260323_235205` | finished | 18 | 0.43097 | 0.81172 | 与 6002 s1234 的 JSD 很接近，但 count 更高；需要 held-out 才能判断是否真优 |
| `gm12878_v2_6002_pilot20_s1234_20260324_105529` | finished | 18 | 0.42852 | 0.79372 | validation 不差，但目前缺 held-out，暂时只能作为 side evidence |

## 3. 主线一：`B_s1234` vs `F_s1234`

这组对照现在已经不只是 validation 信号，而是 held-out 已经给出明确结论：在 matched seed=1234、单卡、`grad_accum_steps=2`、教程数据同口径下，`B=center pool` 在真实 held-out 上继续优于 `F=attention pool`。

关键事实：

- `B_s1234`：`best_epoch=39`，`best val peak JSD=0.33185`，held-out peak 为 `0.3145 / 0.8488`。
- `F_s1234`：`best_epoch=19`，`best val peak JSD=0.44338`，held-out peak 为 `0.4269 / 0.8457`。
- 两条 run 的训练日程口径基本一致，`B_s1234` 配置文件里还明确写了“尽量贴近 6002 的 `F_s1234` 语义”。因此当前差异更像是结构信号，而不是随机波动。
- 更关键的是，`B_s1234` 的 held-out 结果几乎贴住 `B_s42` 的 `0.3147 / 0.8503`；而 `F_s1234` 即使到了 held-out，也没有追回 `F_s42` 的 `0.3148 / 0.8492`。

这意味着：

1. 当前已经没有充分理由继续把 `B/F` 视作并列默认候选。
2. `B` 现在不仅是 validation 领跑者，也是 held-out 领跑者。
3. `F` 仍然能给出不错的 count，但 profile 端在 `seed=1234` 上明显失守，因此更适合作为补充对照，而不是默认方案。

因此，`v2fix` 主线现在已经可以收敛到 `B=center pool`。  
如果后续还想继续补，只需要把这个结论并入总报告或论文表格，不需要再为 `B/F` 开新的训练来做 gate。

## 4. 主线二：`K562_6000` 为什么现在也不能简单归因为“机器差”

`K562_6000` 这次的 held-out 也弱于 6002 参考：

- 6000 held-out peak：`0.6186 / 0.8287`
- 6002 held-out peak：`0.6124 / 0.8586`

因此，“这次只是 validation 偏差”的解释已经站不住。  
不过仍然不能直接把它解释成“6000 机器不行”，原因有两个。

第一，训练口径并非完全等价。  
6000 的 `train_k562_v2_6000_single.yaml` 为了留显存余量，把 `batch_size_per_gpu` 从 6002 的 `8` 降到了 `6`。这已经足以改变优化轨迹。

第二，这条曲线本身表现出明显的 profile/count tradeoff。  
`K562_6000` 在 `epoch 10` 得到最优 profile JSD，但后面 count_r 还在继续上升：

- `epoch 10`: `0.63365 / 0.82014`
- `epoch 16`: `0.63487 / 0.84726`
- `epoch 19`: `0.63648 / 0.84974`
- `epoch 20`: `0.63651 / 0.84238`

这说明它不是简单“整体训坏”，而是：

1. profile 分支在中前期最优
2. count 分支在后期还在继续变强
3. 当前 early-stop 规则只盯 `peak.profile_target_jsd_full_mean`，因此最终 `best.pt` 会落在更偏 profile 的 epoch 10

现在 held-out 已经补上，结论是：

1. 6000 这次 `K562` 的 peak 指标确实弱于 6002 参考。
2. 但它的 overall count 并不差，说明差异更像是“峰区优化没有追到 6002 那么好”，而不是整条链路失效。
3. `2026-04-03 20:38-20:40 CST` 又补跑了 `epoch_019.pt` 的 held-out test，结果是 peak `0.6213 / 0.8414`。相对历史 `best.pt` 的 peak `0.6186 / 0.8287`，它把 count_r 再抬高了约 `+0.0127`，但 profile JSD 也恶化了约 `+0.0027`。这证明后期 checkpoint 确实在继续强化 count 分支，但不足以把它改写成一个比 `best.pt` 更均衡的默认选择。

## 5. 主线三：`GM12878` 现在更像 side evidence，而不是主 gate

`GM12878_6000` 和 `GM12878_6002_s1234` 都在 `epoch 18` 触底，validation 上已经比较接近：

- 6000：`0.43097 / 0.81172`
- 6002 s1234：`0.42852 / 0.79372`

这里的 profile JSD 差距很小，count_r 差距稍大，但没有 held-out 就很难分清：

- 是 seed 差异
- 是机器差异
- 还是 `batch / IO / dataloader` 的副作用

因此 `GM12878` 现在的优先级低于 `B/F` 和 `K562_6000`。  
它值得补，但更适合作为“机器 / seed 方差估计”的第二层证据，而不是当前主结论的第一道 gate。

## 6. 建议的后续评估矩阵

| 优先级 | run | 需要程度 | 目的 | 产出后能回答什么 |
| --- | --- | --- | --- | --- |
| P1 | `gm12878_v2_6000_pilot20_20260323_235205` | 推荐 | 补 6000 侧 held-out | 6000 `GM12878` 是否与历史 6002 结果一致 |
| P1 | `gm12878_v2_6002_pilot20_s1234_20260324_105529` | 推荐 | 补 6002 s1234 held-out | `GM12878` 的 seed 方差大概有多大 |
| P2 | `k562_v2_6000_single_20260324_184918/epoch_019.pt` | 已补（2026-04-03） | 诊断 profile/count tradeoff | 结果为 peak `0.6213 / 0.8414`，验证了 tradeoff，但不改写 `best.pt` 的 profile-first 选择 |

## 7. 后续评估

P0 评估已经完成。  
如果后续还要继续，只需要补 `GM12878` 两侧 held-out，或者做 `K562_6000` 的 `epoch_019.pt` 诊断。

### 7.1 通用 held-out test 模板

```bash
CUDA_VISIBLE_DEVICES=0 nohup python -u -m transchrombp.evaluation.evaluate_checkpoint \
  --checkpoint <best.pt> \
  --split test \
  --nonpeak-ratio 1.0 \
  --output <metrics.json> \
  > <log>.log 2>&1 &
```

### 7.2 当前如果还要继续评估，最值得先跑的 checkpoint

```bash
# 6000: GM12878
/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/checkpoints/gm12878_v2_6000_pilot20_20260323_235205/best.pt

# 6002: GM12878 s1234
/home/zhengwei/bylw_atac/TransChromBP/outputs/checkpoints/gm12878_v2_6002_pilot20_s1234_20260324_105529/best.pt

# 6000: K562 epoch 19
/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/checkpoints/k562_v2_6000_single_20260324_184918/epoch_019.pt
```

### 7.3 2026-04-03 更新：`GM12878_6000` held-out 已补跑

`2026-04-03 18:38-18:40 CST` 在 `6000 GPU1` 重新补跑了 `gm12878_v2_6000_pilot20_20260323_235205/best.pt` 的 full held-out：

- log: `/data1/zhoujiazhen/bylw_atac/logs/gm12878_v2_6000_pilot20_best_test_full_20260403_183834.log`
- json: `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/gm12878_v2_6000_pilot20_best_test_full_20260403_183834.json`
- overall: `count_r=0.7746`、`profile_jsd=0.5517`
- peak: `count_r=0.8179`、`profile_jsd=0.4231`
- nonpeak: `count_r=0.0762`、`profile_jsd=0.6802`

这组数字与 6000 上历史 `best_test_full` 记录实质一致，说明 `GM12878_6000` 的 held-out 偏弱不是瞬时波动，而是当前这条训练配置在 6000 上可复现的结果。换句话说，这条 P1 follow-up 已经回答了“6000 侧结果是否偶然”这个问题；如果还要继续追，下一步更值得补的是 `K562_6000 epoch_019` 或 `GM12878_6002 s1234`，而不是再重复同一个 6000 checkpoint 的 best test。

## 8. 证据路径

### 8.1 `B_s1234`

- log: `/data1/zhoujiazhen/bylw_atac/logs/v2fix_20260324_005834_cpool_s1234_gpu0.log`
- epoch metrics: `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/logs/v2fix_20260324_005834_cpool_s1234_gpu0/epoch_metrics.jsonl`
- run meta: `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/logs/v2fix_20260324_005834_cpool_s1234_gpu0/run_meta.json`
- checkpoint: `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/checkpoints/v2fix_20260324_005834_cpool_s1234_gpu0/best.pt`
- held-out log: `/data1/zhoujiazhen/bylw_atac/logs/v2fix_20260324_005834_cpool_s1234_gpu0_test_full_20260325_002234.log`
- held-out json: `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/v2fix_20260324_005834_cpool_s1234_gpu0_best_test_full_20260325_002234.json`

### 8.2 `F_s1234`

- log: `/home/zhengwei/bylw_atac/logs/v2fix_20260323_attnpool_s1234.log`
- epoch metrics: `/home/zhengwei/bylw_atac/TransChromBP/outputs/logs/v2fix_20260323_attnpool_s1234/epoch_metrics.jsonl`
- run meta: `/home/zhengwei/bylw_atac/TransChromBP/outputs/logs/v2fix_20260323_attnpool_s1234/run_meta.json`
- checkpoint: `/home/zhengwei/bylw_atac/TransChromBP/outputs/checkpoints/v2fix_20260323_attnpool_s1234/best.pt`
- held-out log: `/home/zhengwei/bylw_atac/logs/v2fix_20260323_attnpool_s1234_test_full_20260325_002233.log`
- held-out json: `/home/zhengwei/bylw_atac/TransChromBP/outputs/metrics/v2fix_20260323_attnpool_s1234_best_test_full_20260325_002233.json`

### 8.3 `K562_6000`

- log: `/data1/zhoujiazhen/bylw_atac/logs/k562_v2_6000_single_20260324_184918.log`
- epoch metrics: `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/logs/k562_v2_6000_single_20260324_184918/epoch_metrics.jsonl`
- run meta: `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/logs/k562_v2_6000_single_20260324_184918/run_meta.json`
- checkpoint: `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/checkpoints/k562_v2_6000_single_20260324_184918/best.pt`
- held-out log: `/data1/zhoujiazhen/bylw_atac/logs/k562_v2_6000_single_20260324_184918_test_full_20260325_002234.log`
- held-out json: `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/k562_v2_6000_single_20260324_184918_best_test_full_20260325_002234.json`
- epoch_019 held-out log: `/data1/zhoujiazhen/bylw_atac/logs/k562_v2_6000_single_epoch019_test_full_20260403_203835.log`
- epoch_019 held-out json: `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/k562_v2_6000_single_epoch019_test_full_20260403_203835.json`

### 8.4 参考报告

- `reports/analysis/v2fix_readout_head_status_20260323.md`
- `reports/analysis/6002_independent_datasets_summary_20260323.md`
