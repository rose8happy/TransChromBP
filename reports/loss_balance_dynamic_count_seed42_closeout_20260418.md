# Loss Balance Dynamic Count Seed42 Closeout

## Scope

- Family: `loss_balance_curriculum`
- Dynamic run: `teacher_v2_center_pool_lossbal_e2_dynamic_count_s42_6000_20260417_r1`
- Smoke: `teacher_v2_center_pool_lossbal_e2_dynamic_count_smoke_s42_6000_20260417_r1`
- Machine: `6000` dual A6000
- Comparator: `teacher_v2_center_pool_lossbal_e0_selector_s42_6000_20260417_r1`

Remote evidence:

- Log: `/data1/zhoujiazhen/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_lossbal_e2_dynamic_count_s42_6000_20260417_r1.log`
- Metrics: `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/logs/teacher_v2_center_pool_lossbal_e2_dynamic_count_s42_6000_20260417_r1/epoch_metrics.jsonl`
- Meta: `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/logs/teacher_v2_center_pool_lossbal_e2_dynamic_count_s42_6000_20260417_r1/run_meta.json`

## Queue Change

The originally planned selector-only second seed `teacher_v2_center_pool_lossbal_e0_selector_s1234_6000_20260417_r1` was manually stopped at `2026-04-17 14:55:46 CST` before first validation. The dual-A6000 slot was then repurposed to the dynamic-count experiment.

## Dynamic Run Result

The formal dynamic-count run started at `2026-04-17 15:03:06 CST` and completed at `2026-04-17 21:45:24 CST`. It ran the full `30` epochs without early stop. The best checkpoint was the final epoch.

Best validation snapshot:

| Run | Best Epoch | `peak.profile_target_jsd_full_mean` | `peak.count_pearson_full` | `peak.profile_full_debiased_jsd` | `peak.count_full_debiased_abs` |
|---|---:|---:|---:|---:|---:|
| selector-only seed42 | 21 | 0.331321 | 0.805554 | 0.00010606 | 0.283627 |
| dynamic-count seed42 | 30 | 0.331519 | 0.799202 | 0.00009175 | 0.222861 |

Dynamic schedule summary:

- initial `count_weight = 0.10`
- final `count_weight = 0.2318567`
- total non-null updates = `27`
- update reasons:
  - `push_count = 14`
  - `hold = 8`
  - `protect_profile = 5`

The final logged adjustment at epoch 30 was:

- reason: `hold`
- `jsd_delta = -0.0002267`
- `count_delta = +0.006721`
- `count_weight = 0.2318567 -> 0.2318567`

## Interpretation

This run demonstrates that the dynamic schedule was technically active rather than a no-op: `count_weight` more than doubled over training and the epoch payloads recorded the adjustment history.

But under the current contract, it did not produce a clean win over the selector-only seed42 baseline:

- JSD was slightly worse: `0.331519` vs `0.331321`
- count Pearson was also worse: `0.799202` vs `0.805554`

The safety metrics were not catastrophic and in fact looked somewhat cleaner on the final chosen epoch, but that is not enough to justify saying dynamic-count improved the main trade-off.

## Next Step

The current evidence supports this closeout:

- selector mismatch is real but modest
- dynamic-count schedule is implementable and active
- dynamic-count seed42 did not show a clear gain over selector-only seed42

If this family continues, the next justified move is the static safe-envelope sweep from experiment 1, not a larger dynamic-count sweep and not reopening selector-only `s1234` by default.
