# ChromBPNet Official Patch Ledger（2026-04-06）

## 当前官方根

- 6000 official repo: `/data1/zhoujiazhen/bylw_atac/chromBPNet`
- 本仓不再携带完整官方源码树；这里只记录仍会影响操作的 patch/bridge 接触点。

## 仍需关心的官方侧接触点

1. `chrombpnet/helpers/make_gc_matched_negatives/make_gc_matched_negatives.sh`
   - 用途：tutorial step3 / dataset-prep 背景区域生成。
   - 当前验证方式：`chrombpnet_official_step3_bridge_smoke_20260408_110411` 真机 smoke。
   - 当前要求：本仓 launcher 必须显式走 `CHROMBPNET_OFFICIAL_ROOT=/data1/zhoujiazhen/bylw_atac/chromBPNet` 或等价 `--official-root`。

2. `chrombpnet/training/predict.py`
   - 用途：paper-aligned strict compare / selector / best-epoch 评估。
   - 当前要求：不得再从本地 retired payload 导入；应通过远端 official root 调用。

3. `chrombpnet/training/metrics.py`
   - 用途：official compare 指标链。
   - 当前要求：任何对官方指标语义的讨论都应明确“来自 6000 official repo”，不要把本仓 `vendor/transchrombp/` 与官方指标实现混成一套。

## 本仓对应桥接入口

- `scripts/run_remote_chrombpnet_dataset_prep.sh`
- `scripts/start_6000_chrombpnet_dataset_prep.sh`
- `scripts/start_6002_chrombpnet_dataset_prep.sh`
- `workflows/tutorial/step3_get_background_regions.sh`
- `tests/test_chrombpnet_official_externalization.sh`

## 解释边界

- 本文不是完整 implementation plan。
- 本文也不是结果 closeout；结果判断见 [reports/chrombpnet_official_externalization_closeout_20260408.md](reports/chrombpnet_official_externalization_closeout_20260408.md)。
- 本文只回答“现在到底还依赖官方仓的哪些触点，以及操作时应该把官方根指向哪里”。
