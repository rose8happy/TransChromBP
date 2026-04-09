# U-Net-lite v1 实验严谨性复核（2026-04-09）

## 一句话结论

截至 `2026-04-09 22:40 CST`，`teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409` 这组尝试里，只有 `r3` 算真正有效的 cheap-screen run；`r1` 与 `r2` 都是基础设施 / 配置失败，不能当作负向重复。因此当前最严谨的结论不是“`U-Net-lite family` 已正式判负”，而是：

- `U-Net-lite v1` 当前只有一个有效 short10 结果；
- 这个有效结果相对现有 shortlist 明显不占优，因此**不具备晋级资格**；
- 但若要把它写成更强的“正式停表 / family 判负”，还缺一次确认性复核。

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

### 3.2 对 matched `short10_nofoundation_control`

- 现有锚点：
  - peak `profile_target_jsd_full_mean=0.3193`
  - peak `count_pearson_full=0.8298`

相较这个 matched baseline，`U-Net-lite r3` 也明显落后，距离不属于“边缘摇摆”。

## 4. 这是否足以支持“正式判负”

### 4.1 足以支持的结论

现有证据足以支持：

1. `U-Net-lite v1` 当前**不具备晋级资格**
2. `U-Net-lite v1 r3` 没有达到原先 cheap-screen gate
3. 不能把这条结果写成“readout family 的正向候选”

### 4.2 还不足以支持的结论

现有证据还**不足以**支持：

1. 把 `U-Net-lite family` 整体写成强意义上的“正式判负”
2. 把 `r1/r2/r3` 写成三次负向重复
3. 用“已经充分复现”来证明这条 family 不值得再看

核心原因只有一个：

- 目前只有 `r3` 是有效 run，`r1/r2` 都是启动失败，不是有效重复。

## 5. 当前最严谨的状态标记

因此，当前最严谨的状态应写成：

- `U-Net-lite v1`: `provisional no-go`
- 含义：单个有效 run 已经显示它不值得直接晋级，但还没到“family 级正式停表”的证据强度

## 6. 建议的下一步

若要把这条线收得更硬，最小额外动作应该是二选一：

1. **推荐：同配方确认性 cheap rerun**
   - 在 `6002` 重新跑一次完全同口径 `U-Net-lite v1`
   - 目的不是“救这条线”，而是确认 `r3` 不是单次偶然波动

2. **备选：只做 best checkpoint 的补充 closeout**
   - 明确承认“当前只有一个有效 run”
   - 把结论限制为“v1 single-run no-go，不晋级，不足以 family-level stop”

在没有执行上述任一补强前，不建议把它写成“正式判负已成立”。
