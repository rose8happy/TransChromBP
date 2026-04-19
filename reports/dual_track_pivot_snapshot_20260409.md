# Dual-Track Pivot Snapshot（2026-04-09）

## 一句话结论

`dual-track-20260409` 保存的是从 `msdec_v1 / skipprobe` 切到 `NT v2 teacher-distill + U-Net-lite readout` 的本地 pivot snapshot。现在这条 branch 不再需要 live mounted worktree；恢复入口改成 `snapshot/dual-track-20260409/20260409` 与 branch `dual-track-20260409`。

## 1. 这条 snapshot 冻结了什么

- snapshot commit: `4ffc6e7`
- snapshot tag: `snapshot/dual-track-20260409/20260409`
- run tags:
  - `run/6000/ntv2_teacher_distill_short10_s42_6000_20260409_r2`
  - `run/6002/teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r1`
  - `run/6002/teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r2`
  - `run/6002/teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r3`
  - `run/6002/teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r4`

这条 branch 的价值不是“最终胜出代码”，而是把当时的 pivot 思路、配置和 launcher 冻结成一个可恢复现场。

## 2. branch 上最关键的材料

- 计划与切换说明：
  - `docs/plan/dual_track_ntv2_distill_unet_lite_20260409.md`
  - `reports/dual_track_experiment_pivot_20260409.md`
- 关键代码入口：
  - `vendor/transchrombp/transchrombp/evaluation/teacher_cache_export.py`
  - `vendor/transchrombp/transchrombp/models/profile_decoder.py`
  - `vendor/transchrombp/transchrombp/scripts/run_ntv2_teacher_distill.sh`
  - `vendor/transchrombp/transchrombp/scripts/run_unet_lite_decoder_probe.sh`
- 关键测试：
  - `tests/test_teacher_cache_export.py`
  - `tests/test_ntv2_teacher_distill_launcher.py`
  - `tests/test_unet_lite_launcher.py`

## 3. 为什么现在可以卸载 mounted worktree

- 这条 branch 当前是 clean。
- 相关 family 的 run / closeout / snapshot 证据链已经在 canonical registry 和报告里闭环。
- 真正还缺的不是实验事实，而是 master-side 的解释入口；本文正是那层入口。

所以现在的正确保留方式是：

- 保留 branch `dual-track-20260409`
- 保留 `snapshot/dual-track-20260409/20260409`
- 卸载 mounted worktree

## 4. 局限与注意事项

- 这不是完整 runtime reconstruction 包。
- 尤其 `U-Net-lite` 相关实现，只能视作当时本地 archival snapshot；它和 6002 runtime 上实际运行过的实现并不完全等价。这个差异已在 `reports/unet_vs_alphagenome_reassessment_20260411.md` 中记录。
- 因此，若未来真的要重开这条线，不能直接把 snapshot 当成“已验证可运行的当前事实”，而应先区分 archival snapshot 与 runtime 事实。

## 5. 恢复方式

如需重新查看当时现场，优先顺序如下：

1. 先读本报告与相关 closeout：
   - `reports/unet_lite_v1_rigor_review_20260409.md`
   - `reports/ntv2_teacher_distill_short10_gate_closeout_20260409.md`
2. 再从 snapshot 恢复代码现场：
   - `git switch dual-track-20260409`
   - 或基于 `snapshot/dual-track-20260409/20260409` 新开 branch/worktree

当前不建议再为“只是可浏览”而长期保留 mounted worktree。
