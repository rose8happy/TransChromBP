# TransChromBP 项目实验全景图

这份文档尽量按 **family** 而不是按单条 run 做总览。

它的目的不是替代所有 closeout，而是帮你快速回答：

1. 这个项目到底做过哪些主要实验线？
2. 每条线改了什么？
3. 哪些变成当前主线，哪些停在补证、边界或负结果？

## 1. 先给一张总判断

### 已经进入当前默认口径的 family

- backbone ablation：证明 Transformer 有真实收益
- readout design：证明 `center pool` 是默认 readout
- shortcut clean matrix：把叙事从“Transformer 特有 shortcut”改写成“bias-safe framework + full/debiased diagnostics”
- strict compare `L3 shared-region`：当前最强 external evidence

### 只保留为补证或边界的 family

- no-bias：supplementary boundary result
- AlphaGenome matched raw-track：external sidecar/reference

### 已完成但没有晋级的 family

- Genos integration
- Caduceus token fusion
- NT v2 probe / residual gates
- U-Net-lite v1
- multiscale local skip formal gate
- recent loss-balance dynamic-count note

## 2. 主线 family

| 问题 | 改动点 | 代码入口 | 配置入口 | 代表 run | 关键指标 | 结论 | 详细报告 |
|---|---|---|---|---|---|---|---|
| Transformer 到底有没有真实收益 | 只切 `sequence_encoder.enabled`，做 `V2-full` vs `V2-noTF` matched ablation | `../../vendor/transchrombp/transchrombp/models/transchrombp.py`、`../../vendor/transchrombp/transchrombp/models/transformer_encoder.py` | `../../vendor/transchrombp/transchrombp/configs/model/ablations/ablation_no_transformer.yaml`、`../../vendor/transchrombp/transchrombp/configs/train/train_ablation_v2_main_profile_select.yaml` | `V2-full` / `V2-noTF` 三 seed 汇总 | `V2-full=0.31419±0.00012 / 0.83676±0.00531`；`V2-noTF=0.32393±0.00021 / 0.82503±0.00489` | Transformer 的收益已稳住，是主 claim 之一 | [reports/transchrombp_internal_design_and_experiment_history_20260331.md](../../reports/transchrombp_internal_design_and_experiment_history_20260331.md) |
| count 读出头该怎么选 | 比较 `center pool` / `attention pool` / 其他读出头 | `../../vendor/transchrombp/transchrombp/models/transchrombp.py` 中的 `count_pool_mode` | `../../vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool.yaml`、`../../vendor/transchrombp/transchrombp/configs/model/v2fix_attn_pool.yaml` | `B_s42`、`B_s1234` | `B center pool=0.31460±0.00010 / 0.84958±0.00103`；`F_s1234` 退化到 `0.4269 / 0.8457` | `center pool` 是当前默认，不再把 `B/F` 当并列主候选 | [reports/v2fix_readout_head_status_20260323.md](../../reports/v2fix_readout_head_status_20260323.md)、[reports/v2fix_and_6000_realdata_followup_20260325.md](../../reports/v2fix_and_6000_realdata_followup_20260325.md) |
| 早期 shortcut 叙事到底哪里稳、哪里不稳 | 系统扫描 `sg`、`deb2`、TF/noTF、center pool 等组合 | `../../vendor/transchrombp/transchrombp/training/train_ddp.py`、`../../vendor/transchrombp/transchrombp/models/transchrombp.py` | clean-matrix 相关 ablation 配置 | clean matrix A/C/corrected-B | `A=0.31410 / 0.85064`；`C(two-seed)=0.31730±0.00014 / 0.83671±0.00389`；corrected-B 的 `full-debiased gap=0.00268` | 当前没有证据支持“Transformer 特有 shortcut”；更稳的是 `deb2 + full/debiased gap` 叙事 | [reports/profile_shortcut_revalidation_summary_20260326.md](../../reports/profile_shortcut_revalidation_summary_20260326.md)、[reports/paper_claim_evidence_matrix_20260326.md](../../reports/paper_claim_evidence_matrix_20260326.md) |
| 相对 official 的外部比较能不能站稳 | strict compare 一路推进到 `L3 shared-region` | `../../vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py`、`../../scripts/paper_aligned_repro/select_best_epoch.py` | `../../vendor/transchrombp/transchrombp/configs/train/train_tutorial_corrected_b_strict_compare_L3_6000.yaml` | `tutorial_corrected_b_strict_compare_L3_s42` | ours `0.31319 / 0.84016` vs official `0.33853 / 0.69958` | 当前最强 external evidence，不应再只盯旧 tutorial baseline | [reports/tutorial_L3_shared_region_closure_20260330.md](../../reports/tutorial_L3_shared_region_closure_20260330.md) |

## 3. 边界 / 补证 family

| 问题 | 改动点 | 代码入口 | 配置入口 | 代表 run | 关键指标 | 结论 | 详细报告 |
|---|---|---|---|---|---|---|---|
| bias branch 是不是绝对必要 | 直接去掉 bias branch 看 boundary | `../../vendor/transchrombp/transchrombp/models/transchrombp.py` | `../../vendor/transchrombp/transchrombp/configs/model/ablations/ablation_no_bias.yaml` | `paper_no_bias_s42_2gpu` | peak `JSD=0.31496`，peak `count_r=0.84978` | 这条线只能放 supplementary，不改主线 | [reports/paper_writeup_flow_20260330.md](../../reports/paper_writeup_flow_20260330.md) |
| AlphaGenome 外部坐标对照能不能跑通 | 不接入训练主干，只做 matched raw-track sidecar | sidecar pipeline | `regions_k562_tutorial_matched_panel_v2.csv` 等 | `alphagenome_matched_raw_track_slice_v2_20260410` | `16` 个 loci 全部完成；每个位点 `num_tracks_used=1` | technical/reference `pass`，不是模型质量胜利 | [reports/alphagenome_matched_raw_track_slice_v2_closeout_20260410.md](../../reports/alphagenome_matched_raw_track_slice_v2_closeout_20260410.md) |

## 4. foundation-model family

| 问题 | 改动点 | 代码入口 | 配置入口 | 代表 run | 关键指标 | 结论 | 详细报告 |
|---|---|---|---|---|---|---|---|
| Genos 的 frozen summary 能不能补当前 backbone | online gate / mean、cached summary、count-only 注入 | `../../vendor/transchrombp/transchrombp/models/genos_adapter.py`、`../../vendor/transchrombp/transchrombp/models/transchrombp.py` | `v2fix_genos_gate.yaml`、`v2fix_genos_mean.yaml`、`train_genos_cached_short10.yaml` | `G1`、`G2`、`P2`、`G0/P0` | `G1=0.3388 / 0.4689`；`G2=0.3387 / 0.6360`；`P2=0.3174 / 0.6037`；baseline `G0=0.3163 / 0.8410`；probe `genos_only AUC=0.6997`, `concat AUC=0.9018 < encoded_only 0.9219` | 高质量负结果；模型本体没坏，但 integration recipe 不成立 | [reports/genos_no_positive_gain_analysis_20260323.md](../../reports/genos_no_positive_gain_analysis_20260323.md) |
| token-level frozen fusion 会不会更自然 | 在主干前接 `Caduceus-PS` token hidden states | `../../vendor/transchrombp/transchrombp/models/caduceus_adapter.py`、`../../vendor/transchrombp/transchrombp/training/train_ddp.py` | `../../vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_caduceus_ps.yaml` | matched A/B 单 seed | A `0.3132558744 / 0.8388160470`；B `0.3126602104 / 0.8466657043`；JSD 只改善 `0.00060` | `near-null / marginal positive`，不值得自然扩线 | [docs/plan/caduceus_tutorial_ab_online_fusion_20260331.md](../plan/caduceus_tutorial_ab_online_fusion_20260331.md) |
| NT v2 是否比 Genos 更像有效特征源 | 先 probe，再做 residual gate | `../../vendor/transchrombp/transchrombp/models/foundation_adapter.py`、`../../vendor/transchrombp/transchrombp/scripts/build_foundation_cache.py` | `transchrombp_teacher_v2_center_pool_ntv2_residual.yaml` 等 | `nt_v2_repair_smoke_probe_20260404_211117` + residual short10 family | probe: `nt_only AUC=0.742424`，`concat AUC=0.857792 < encoded_only 0.868116`；held-out residual: `0.3560 / 0.7729` 与 `0.3588 / 0.7516`，均落后 control `0.3193 / 0.8298` | 有信号，但仍不互补；measured family 该停 | [reports/nt_v2_probe_20260404.md](../../reports/nt_v2_probe_20260404.md)、[reports/paper_claim_evidence_matrix_20260326.md](../../reports/paper_claim_evidence_matrix_20260326.md) |

## 5. AlphaGenome-like / next-gen structure family

| 问题 | 改动点 | 代码入口 | 配置入口 | 代表 run | 关键指标 | 结论 | 详细报告 |
|---|---|---|---|---|---|---|---|
| U-Net-lite 这种更像 encoder-decoder 的结构会不会直接更强 | 轻量 U-Net 风格 cheap-screen / confirmation | 对应 worktree 的 experimental implementation | `unet_lite_v1` family configs | `teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r4` | best peak `0.448397 / 0.798992`；matched no-foundation control `0.3193 / 0.8298` | `no-go / stop`，当前配方不具备晋级资格 | [reports/unet_lite_v1_rigor_review_20260409.md](../../reports/unet_lite_v1_rigor_review_20260409.md) |
| multiscale local skip decoder 能不能在 A6000 formal gate 过线 | 更显式的 multiscale / local skip readout | 对应 `msdls_v2` experimental modules | `msdls_v2` formal gate configs | `teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1` | best peak `0.333781 / 0.794823`；历史 corrected-B comparator 约 `0.3146 / 0.8496` | 正式 `fail`，不是边缘摇摆 | [reports/multiscale_local_skip_v2_30ep_gate_closeout_20260410.md](../../reports/multiscale_local_skip_v2_30ep_gate_closeout_20260410.md) |

## 6. 近期维护 / 训练语义 family

| 问题 | 改动点 | 代码入口 | 配置入口 | 代表 run | 关键指标 | 结论 | 详细报告 |
|---|---|---|---|---|---|---|---|
| 只改 selector 会不会比按总 loss 选模更好 | selector-only，不改 loss 合同 | `train_ddp.py` 的 `best_metric` 语义；外部 selector 脚本 | selector-aligned train config | `teacher_v2_center_pool_lossbal_e0_selector_s42_6000_20260417_r1` | best-by-JSD `0.331321 / 0.805554`；best-by-loss `0.332464 / 0.804109`；`selector_jsd_gap=0.001143` | selector mismatch 是真的，但幅度不大 | [reports/loss_balance_selector_wave1_seed42_20260417.md](../../reports/loss_balance_selector_wave1_seed42_20260417.md) |
| 动态调 `count_weight` 会不会改善 JSD/count trade-off | warmup 后按 validation 调整 `count_weight` | 当前归档快照未保留稳定代码入口；以 closeout 为准 | dynamic-count local patch config | `teacher_v2_center_pool_lossbal_e2_dynamic_count_s42_6000_20260417_r1` | selector-only `0.331321 / 0.805554`；dynamic `0.331519 / 0.799202`；`count_weight 0.10 -> 0.2318567` | 动态器是活的，但当前 contract 下没有 clear gain | [reports/loss_balance_dynamic_count_seed42_closeout_20260418.md](../../reports/loss_balance_dynamic_count_seed42_closeout_20260418.md) |

## 7. 如果你只想记住“哪些线值得先看”

推荐优先级：

1. backbone ablation
2. readout / `center pool`
3. clean matrix
4. `L3 shared-region`
5. no-bias
6. foundation negative results
7. U-Net / multiscale / loss-balance 这些后续维护和下一代探索

## 8. 一句话结论

当前项目最值钱的不是“做过很多实验”本身，而是：

> 多数 family 已经被问到一个很清楚的结论层级上了：哪些构成今天的默认主线，哪些只是补证，哪些是高质量负结果，哪些只能作为 sidecar/reference 保留。
