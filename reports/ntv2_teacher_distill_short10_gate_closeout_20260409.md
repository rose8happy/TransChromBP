# NT v2 Teacher-Distill Short10 Gate Closeout（2026-04-09）

## 一句话结论

`ntv2_teacher_distill_short10_s42_6000_20260409_r2` 已于 `2026-04-09 22:43 CST` early-stop 收口，并且按原计划定义的 short10 gate **明确不过线**。因此：

- `NT v2 teacher-distill tutorial` 这条线当前不升 `full-budget distill`
- 不扩 seed
- 按 stop-rule 直接停表

## 1. Gate 规则

根据 [dual_track_ntv2_distill_unet_lite_20260409.md](/home/zhengwei/.config/superpowers/worktrees/chromBPNet/dual-track-20260409/docs/plan/dual_track_ntv2_distill_unet_lite_20260409.md) 的原始定义，`teacher-distill short10` 相对 matched `short10_nofoundation_control`，必须同时满足：

1. `peak JSD` 改善 `>= 0.002`
2. `peak count_r` 不下降超过 `0.005`
3. profile / count 的优势不能靠单侧明显塌陷换来

## 2. 当前 run 结果

run：

- `ntv2_teacher_distill_short10_s42_6000_20260409_r2`

关键事实：

- `best_epoch=2`
- `stopped_early=true`
- `early_stop_reason=no improvement in peak.profile_target_jsd_full_mean for 3 validations`

best `val:peak`：

- `profile_target_jsd_full_mean=0.345324`
- `count_pearson_full=0.5781`

日志 / 元数据来源：

- `/data1/zhoujiazhen/bylw_atac/TransChromBP/logs/ntv2_teacher_distill_short10_s42_6000_20260409_r2.log`
- `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/logs/ntv2_teacher_distill_short10_s42_6000_20260409_r2/run_meta.json`

## 3. 同口径对照：matched no-foundation short10 control

为避免把 `valid` 与 `held-out` 混比，本次 gate 使用 matched control 的 **同口径 `valid peak`**：

run：

- `short10_nofoundation_control_s42_20260405_dual`

best `val:peak`：

- `profile_target_jsd_full_mean=0.336346`
- `count_pearson_full=0.8006`
- `best_epoch=9`

来源：

- `/data1/zhoujiazhen/bylw_atac/logs/short10_nofoundation_control_s42_20260405_dual_6000.log`
- `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/logs/short10_nofoundation_control_s42_20260405_dual/run_meta.json`

## 4. Gate 判读

逐条对照：

1. `peak JSD` 改善 `>= 0.002`
   - control：`0.336346`
   - distill：`0.345324`
   - 结果：没有改善，反而变差约 `+0.00898`

2. `peak count_r` 不下降超过 `0.005`
   - control：`0.8006`
   - distill：`0.5781`
   - 结果：下降约 `-0.2225`，远超允许范围

3. profile / count 不靠单侧塌陷换来
   - 这条无需进一步判，因为前两条已经同时明显失败

结论：

- 这不是边缘摇摆，也不是“差一点过线”
- 这是**明确不过 gate**

## 5. Stop-Rule 执行结果

按原 stop-rule，本轮执行结果固定为：

1. 不升 `full-budget distill`
2. 不扩 seed
3. `teacher-distill tutorial` 线当前直接停表

## 6. 对当前双机队列的影响

这次 gate 判负的直接影响是：

- `6000` 当前不再有合法的 `teacher-distill tutorial` 后续
- 但这**不等于**默认切回“先写论文”
- 这份 closeout 只裁决 `teacher-distill tutorial` 线本身；当前双机默认动作以 charter / TRACKING 为准
- `6000` 不再因为这份 closeout 而等待 `6002`
- 若要提 `6002 r4`，只能把它当作当时并行存在的现场事实，不能当默认串行下一步

因此当前默认下一步变成：

1. 先按 charter / TRACKING 继续当前双机各自的 live 队列
2. 完成 `U-Net-lite v1` 的正式判读
3. 若后续要重新排 `6002` 的 `r4` 或其他分支，再按当时现场事实单独判断，不要从这份 closeout 里倒推出默认串行等待链
