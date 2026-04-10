# U-Net-lite v1 实验严谨性复核（2026-04-10）

## 一句话结论

截至 `2026-04-10 03:34:42 CST`，`teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409` 这组尝试已经完成 `r4` 确认性复跑。`r1` 与 `r2` 仍然是基础设施 / 配置失败，`r3` 与 `r4` 是两条有效 cheap-screen / confirmation run；其中 `r4` 没有把这条线推上 shortlist，反而把原先的 provisional 判断收得更硬。当前最严谨的结论是：

- `U-Net-lite v1` 已经有两条有效 short10 结果；
- 这两条结果都明显不占优，因此**不具备晋级资格**；
- 这条 family 现在可以收成 `no-go / stop`，但不应写成正向 shortlist。

## 1. 复核问题

要回答的问题不是“`r3` 看起来差不差”，而是：

1. `r1/r2/r3` 里哪些 run 真正有效；
2. `r3` 是否跑在预期配置上；
3. 现有证据是否足以支持“正式判负”，还是只能支持“暂不晋级”。

## 2. 逐条 run 复核

### 2.1 `r1` 不是有效负样本

- 日志：`/home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r1.log`
- 根因：
  - `ChromBPNetBigWigDataset.__init__() got an unexpected keyword argument 'foundation_cache_dir'`
- 结论：
  - 这是 loader / config contract 没对齐导致的启动失败；
  - `r1` 不能算一次有效训练，更不能当作负向重复。

### 2.2 `r2` 不是有效负样本

- 日志：`/home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r2.log`
- 根因：
  - bigWig 路径仍指向 `6000` 的 `/data1/.../merged_unstranded.bw`
  - `pyBigWig.open(...)` 直接报错
- 结论：
  - 这是数据路径没有完全 remap 到 `6002` 导致的启动失败；
  - `r2` 同样不能算有效训练。

### 2.3 `r3` 是当前唯一有效 run

- 日志：`/home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r3.log`
- 关键事实：
  - log header 明确写的是 `Data config: /home/zhengwei/bylw_atac/TransChromBP/configs/data/data_tutorial_canonical_v1_6002.yaml`
  - seed 为 `42`
  - `BATCH_SIZE_PER_GPU=16`
  - `epoch=5` 触发 early-stop 并正常收口
  - checkpoint 目录存在 `best.pt` 和 `epoch_001` 到 `epoch_005`
- best 指标：
  - `epoch=2`
  - peak `profile_target_jsd_full_mean=0.45525`
  - peak `count_pearson_full=0.7711`
- 结论：
  - `r3` 是一条技术上有效、配置口径基本正确的 short10 cheap-screen。

### 2.4 `r4` 是有效确认性复跑，但仍不构成晋级

- 日志：`/home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r4.log`
- 完成时间：`2026-04-10 03:34:42 CST`
- 关键事实：
  - `best_epoch=6`
  - `stopped_early=true`
  - `best peak profile_target_jsd_full_mean=0.44839715448596246`
  - `best peak count_pearson_full=0.7989915325863597`
  - `best overall profile_target_jsd_full_mean=0.4681349953322501`
  - `best overall count_pearson_full=0.8049882518128701`
  - final epoch `9` 的 peak 仍是 `0.45226327185176557 / 0.7932712502144175`
- 结论：
  - `r4` 比 `r3` 更严谨，但没有改变这条 family 的方向；
  - 它确认了 `U-Net-lite v1` 不是“只差一次好 seed 就能晋级”的那类候选。

## 3. 与当前 shortlist 的对照

### 3.1 对 `skipprobe_v1_wide`

- 日志：`/home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_skipprobe_v1_wide_short10_s42_6002_r1.log`
- 同口径关键信息：
  - 同样是 `6002`
  - seed `42`
  - `BATCH_SIZE_PER_GPU=16`
  - data config 也是 `data_tutorial_canonical_v1_6002.yaml`
- 代表性读数：
  - best peak 约 `profile_target_jsd_full_mean=0.44754`
  - 对应 peak `count_pearson_full=0.8024`

相较之下，`U-Net-lite r3` 的 `0.45525 / 0.7711` 在 `JSD` 和 `count_r` 两侧都更差。

`r4` 虽然相较 `r3` 有所改善，但仍然没有追上 `skipprobe_v1_wide`：

- `r4` best peak `0.448397 / 0.798992`
- `skipprobe_v1_wide` best peak `0.447541 / 0.802370`

这说明 `r4` 只能把 `U-Net-lite v1` 从“更差”推到“仍然不够好”，不能把它推成 shortlist。

### 3.2 对 matched `short10_nofoundation_control`

- 现有锚点：
  - peak `profile_target_jsd_full_mean=0.3193`
  - peak `count_pearson_full=0.8298`

相较这个 matched baseline，`U-Net-lite r3` 也明显落后，距离不属于“边缘摇摆”。

## 4. 这是否足以支持“正式判负”

### 4.1 足以支持的结论

现有证据足以支持：

1. `U-Net-lite v1` 当前**不具备晋级资格**
2. `U-Net-lite v1 r4` 只是把原先的 provisional no-go 再确认了一次
3. 不能把这条结果写成“readout family 的正向候选”
4. 这条 family 可以停表，不需要再按当前同配方继续 cheap rerun

### 4.2 还不足以支持的结论

现有证据还**不足以**支持：

1. 把 `U-Net-lite family` 写成“正向 shortlist”
2. 把 `r1/r2/r3/r4` 写成需要继续扩线的信号
3. 用 `r4` 的小幅改进去主张 family-level promotion

核心原因只有一个：

- 即使加上 `r4`，这条 family 也只是“有效但不够好”，不是“值得晋级”。

## 5. 当前最严谨的状态标记

因此，当前最严谨的状态应写成：

- `U-Net-lite v1`: `no-go / stop`
- 含义：`r4` 已经把它从“暂不晋级”推进到“无需继续按当前配方再跑”

## 6. 建议的下一步

当前不建议再为 `U-Net-lite v1` 追加同配方 run。

如果未来要重开，只能是一个显式的新 hypothesis，而不是延续当前 `r3/r4` 配方。
