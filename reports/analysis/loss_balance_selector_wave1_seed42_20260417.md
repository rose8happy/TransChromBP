# Loss Balance Selector Wave 1 Seed42 Summary

## Scope

- Family: `loss_balance_curriculum`
- Run: `teacher_v2_center_pool_lossbal_e0_selector_s42_6000_20260417_r1`
- Machine: `6000` dual A6000
- Objective: quantify the gap between selecting checkpoints by `peak.profile_target_jsd_full_mean` and selecting them by `peak.loss_total` under the current corrected-B loss contract

Remote evidence:

- Log: `/data1/zhoujiazhen/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_lossbal_e0_selector_s42_6000_20260417_r1.log`
- Metrics: `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/logs/teacher_v2_center_pool_lossbal_e0_selector_s42_6000_20260417_r1/epoch_metrics.jsonl`
- Meta: `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/logs/teacher_v2_center_pool_lossbal_e0_selector_s42_6000_20260417_r1/run_meta.json`
- Structured summary: [reports/assets/loss_balance_selector_seed42_summary_20260417.json](assets/loss_balance_selector_seed42_summary_20260417.json)

## Result

The run early-stopped at epoch 29. The JSD-selected checkpoint was epoch 21, while the loss-selected checkpoint was epoch 22.

| Selector | Epoch | `peak.profile_target_jsd_full_mean` | `peak.count_pearson_full` | `peak.profile_full_debiased_jsd` | `peak.count_full_debiased_abs` |
|---|---:|---:|---:|---:|---:|
| best by JSD | 21 | 0.331321 | 0.805554 | 0.00010606 | 0.283627 |
| best by loss | 22 | 0.332464 | 0.804109 | 0.00010614 | 0.211907 |

Derived selector gap:

- `selector_epoch_gap = 1`
- `selector_jsd_gap = 0.001143`

## Interpretation

Seed42 shows a real selector mismatch: choosing by JSD recovers a better main metric than choosing by loss, and the count Pearson is also slightly better. The gain is modest rather than dramatic. Safety on `profile_full_debiased_jsd` is effectively unchanged, while `count_full_debiased_abs` is somewhat worse at the JSD-selected epoch than at the loss-selected epoch.

This means the selector is not a no-op, but one seed is not enough to conclude that selector-only fully explains the gap. The current evidence supports one more seed under the same selector-aligned contract before paying the cost of the 9-point static safe-envelope sweep.

## Next Step

This next-step recommendation was superseded on `2026-04-17`: `teacher_v2_center_pool_lossbal_e0_selector_s1234_6000_20260417_r1` was started but then manually stopped at `2026-04-17 14:55:46 CST`, and the dual-A6000 slot was repurposed to `teacher_v2_center_pool_lossbal_e2_dynamic_count_s42_6000_20260417_r1`.

The current canonical continuation is therefore:

- use this report as the selector-only baseline reference
- use [reports/closeout/loss_balance_dynamic_count_seed42_closeout_20260418.md](loss_balance_dynamic_count_seed42_closeout_20260418.md) as the dynamic-count closeout
- if the family continues, decide whether experiment 1 static safe-envelope is worth paying for
