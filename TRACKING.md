# 项目追踪清单（工作清单 + 下载资源清单）

本文件用于减少支线过多导致的信息遗漏。每次开始新任务或汇报前，先查看并按最新状态更新本文件。

> 已完成 / 已合并事项已迁至 [TRACKING_archive.md](TRACKING_archive.md)。
> 第一节原则上只保留 `进行中` / `待处理` / `待验证` 条目；阶段性完成或已并入其它主线的事项转入 `TRACKING_archive.md`。
> 当前没有 active 的 `6000` training run。`loss_balance_curriculum` 的最新 terminal run 是 `teacher_v2_center_pool_lossbal_e2_dynamic_count_s42_6000_20260417_r1`（dual A6000，`2026-04-17 15:03:06-21:45:24 CST`）：formal run 跑满 `epoch=30`，best checkpoint 在 `epoch=30`，`peak.profile_target_jsd_full_mean=0.33152`、`peak.count_pearson_full=0.79920`、`peak.profile_full_debiased_jsd=9.18e-05`、`peak.count_full_debiased_abs=0.22286`；动态 `count_weight` 从 `0.10` 增至 `0.23186`，共记录 `27` 次更新（`push_count=14 / hold=8 / protect_profile=5`）。原计划的 selector-only `teacher_v2_center_pool_lossbal_e0_selector_s1234_6000_20260417_r1` 已在 `2026-04-17 14:55:46 CST` 手动停止并让位给动态实验。`alphagenome_factor_ladder` 的 unattended queue 已在 `2026-04-12 08:46 CST` 停在 `export_e2_teacher_cache`，错误为 `model_teacher_cache_export.py: error: unrecognized arguments: --targets profile16 logcount`。
> 实验 family / run 的 canonical 索引统一看 [docs/experiments/registry.md](docs/experiments/registry.md) 和 [docs/experiments/runs.csv](docs/experiments/runs.csv)。

状态约定：
- `已完成`：可直接使用
- `进行中`：任务仍在跑或尚未收尾
- `待处理`：明确要做但还没开始
- `待验证`：文件/目录已在，但需要实际加载或运行验证

## 一、在做事情清单（实时更新）

> 双机实验规则源统一看 [docs/plan/2026-04-09_dual_machine_experiment_charter.md](docs/plan/2026-04-09_dual_machine_experiment_charter.md)。
> 当前仓库 / worktree / 双机运行快照统一看 [reports/repository_status_handoff_20260409.md](reports/repository_status_handoff_20260409.md)；本节只保留结论级 live 状态。
> family/workstream 与 run 的检索入口统一看 [docs/experiments/registry.md](docs/experiments/registry.md) 和 [docs/experiments/runs.csv](docs/experiments/runs.csv)。

| 事项 | 当前状态 | 当前结论 / 进度 | 下一步 | 关键路径 |
|---|---|---|---|---|
| 动态多 loss 协调外发专题包（2026-04-16） | 待验证 | 已新建 `reports/chatgpt_bundle_loss_balance_20260416/`，把本轮外发主问题从“总体技术路线”收紧到“动态多 loss 协调是否值得迁入当前 `profile/count/full-debiased` 框架”。当前包内已整理：任务与 loss contract、loss 协调相关证据链、面向 `GPT Pro` 的主问题清单、关键指标简表，以及完整代码 / 配置快照 `06_transchrombp_current_model.py`、`07_train_ddp_current.py`、`08_train_teacher_v2_main.yaml`、`09_train_teacher_v2_main_profile_select.yaml`。这次外发明确依赖 `GPT Pro` 自行联网检索近年多目标 / 多 loss 协调论文，而不是仓内手写综述；当前最关键的内部锚点是：`debiased_profile_weight=2.0` 比 `stop-gradient` 更像主效应、`peak.loss_total` 与 `peak.profile_target_jsd_full_mean` selector 不一致、以及 `chrombpnet_auto` 风格 count weighting 不能直接迁到 `learnable_scales` 语境。 | 用 `00/10` 的新 prompt 把第一轮文件发给 `GPT Pro`，让其先基于现有代码事实和搜索结果回答“哪类动态 loss 协调方法值得迁入、最小验证闭环是什么”；若对方要求核对实现，再补传 `06/07`。拿到外部建议后，再决定是否立一条显式新 hypothesis 和最小实验计划。 | `reports/chatgpt_bundle_loss_balance_20260416/00_readme_upload_order_and_prompt.md`、`reports/chatgpt_bundle_loss_balance_20260416/03_questions_for_external_models_dynamic_loss.md`、`reports/chatgpt_bundle_loss_balance_20260416/10_ready_to_send_chatgpt_pro_prompt.md` |
| loss balance curriculum（2026-04-17） | 进行中 | selector-only seed42 baseline `teacher_v2_center_pool_lossbal_e0_selector_s42_6000_20260417_r1` 已收口：best-by-JSD 在 `epoch=21`，`peak.profile_target_jsd_full_mean=0.3313`、`peak.count_pearson_full=0.8056`；best-by-loss 在 `epoch=22`，`JSD=0.3325`、`count_r=0.8041`，`selector_jsd_gap=0.00114`。随后原计划的 selector-only `s1234` 在 `2026-04-17 14:55:46 CST` 手动停表，把双卡让给 dynamic-count formal `teacher_v2_center_pool_lossbal_e2_dynamic_count_s42_6000_20260417_r1`。这条动态 run 已于 `2026-04-17 21:45:24 CST` 跑满 `epoch=30` 并收口，best 在 `epoch=30`，`peak.profile_target_jsd_full_mean=0.33152`、`peak.count_pearson_full=0.79920`、`peak.profile_full_debiased_jsd=9.18e-05`、`peak.count_full_debiased_abs=0.22286`；`count_weight` 从 `0.10` 增到 `0.23186`，共发生 `27` 次调整。相对 selector-only seed42，dynamic-count 没拿到更好的 JSD，也没有更好的 count。 | 先把这条 closeout 视为“dynamic count 在当前 contract 下没有给出 clear gain”的证据；若还要继续 loss-balance 线，优先决定是否进入实验 1 的静态 safe-envelope，而不是重开 selector-only `s1234`。除非出现显式新 hypothesis，否则暂不继续扩 dynamic-count 的 sweep。 | `reports/loss_balance_selector_wave1_seed42_20260417.md`、`reports/loss_balance_dynamic_count_seed42_closeout_20260418.md`、`docs/experiments/runs.csv` |
| AlphaGenome-like factor ladder（E1/E2/E3） | 待处理 | `alphagenome_factor_ladder` 的 unattended queue 在 `2026-04-12 08:46 CST` 停在 `export_e2_teacher_cache`。已完成的 terminal training runs 是 `teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1` 与 `teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1`；queue state dir=`/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411/outputs/queue/factor_ladder_unattended_20260411`。halt log=`/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411/logs/export_e2_teacher_cache.log`，报错 `model_teacher_cache_export.py: error: unrecognized arguments: --targets profile16 logcount`，因此 export 没有成功落盘，后续 `E3 distill` 也没有继续推进。 | 先修 exporter 参数漂移：对齐 `model_teacher_cache_export.py` 当前 CLI，解决 `--targets profile16 logcount` 不再被接受的问题；修复并补一条最小 smoke 前，未修复前不要继续推进 E3 distill，也不要重新启动 factor-ladder unattended queue。 | `docs/superpowers/specs/2026-04-11-alphagenome-like-factor-ladder-design.md`、`reports/unet_vs_alphagenome_reassessment_20260411.md`、`docs/experiments/registry.md` |
| U-Net vs AlphaGenome 架构复核（2026-04-11） | 待处理 | 截至 `2026-04-11 17:07:56 CST`，已完成本地 canonical 文档、`dual-track-20260409` 冻结 worktree、`6002` runtime 代码与 `r4` 日志、以及 AlphaGenome 官方论文的交叉复核。当前最稳妥的结论是：本仓现有负向证据**足以否掉**“在 `corrected B` 上仅把 `2114bp -> 1000bp` 的 debiased profile readout 换成 `U-Net-lite / multiscale-local-skip decoder` 就能带来 clean gain”这类 drop-in family；但它**不足以否掉** AlphaGenome 那种 `1 Mb` 长上下文、全模型 `U-Net-style encoder + transformer tower + decoder`、多模态多任务、且带 teacher/distillation 的更大范式。额外发现：`6002` runtime 上确实存在并运行了 `UNetLiteProfileDecoder`，但当前冻结 worktree 的 `vendor/.../profile_decoder.py` / `transchrombp.py` 没把这段实现完整归档，说明事后复盘时必须区分 runtime 事实和 archival snapshot。基于这层复核，当前已新增一份正式设计 spec，把后续高价值实验收成 `AlphaGenome-like factor ladder`。 | 若要因为 AlphaGenome 重新考虑 U-Net，先决定是否立一个**显式新 hypothesis**：至少要把问题改写成“更接近 AlphaGenome 训练前提的全模型/长上下文/teacher-or-distill 路线是否值得做最小验证”，而不是继续追加当前 `unet_lite_v1` / `msdls_v2` recipe 的 cheap rerun；若真要重开，先补 canonical archive 与 runtime 差异记录，并以新 spec 为准进入实现计划。 | `docs/superpowers/specs/2026-04-11-alphagenome-like-factor-ladder-design.md`、`reports/unet_vs_alphagenome_reassessment_20260411.md`、`reports/unet_lite_v1_rigor_review_20260409.md`、`reports/multiscale_local_skip_v2_30ep_gate_closeout_20260410.md`、`reports/alphagenome_matched_raw_track_slice_v2_closeout_20260410.md`、`/home/zhengwei/.config/superpowers/worktrees/chromBPNet/dual-track-20260409/vendor/transchrombp/transchrombp/models/profile_decoder.py`、`/home/zhengwei/.config/superpowers/worktrees/chromBPNet/dual-track-20260409/vendor/transchrombp/transchrombp/models/transchrombp.py` |
| 仓库任务盘点（2026-04-09） | 进行中 | 截至 `2026-04-10 10:10:10 CST` 的 worktree cleanup 确实把旧 `6000/6002` training run 清成了“全部收口”的终态：`master` 仍是 clean 的 canonical 档案仓，冻结支线继续保留；`6000` 的 A6000 formal gate `teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1` 已收口并判为 `fail`；并行 sidecar `alphagenome_matched_raw_track_slice_v2_20260410` 已收口并判为 `pass`；`6002` 的 `U-Net-lite v1` 也已完成 `r4` 确认性复跑并收成 `no-go / stop`。后续新增的 `alphagenome_factor_ladder` 与 `loss_balance_curriculum` 进展，现已分别转入各自 family 条目与 canonical manifest，不再在本行重复维护 active run。 | 继续按 charter 维护其余 live 项；family / run 的当前状态统一以下方各 family 行和 `docs/experiments/registry.md` / `docs/experiments/runs.csv` 为准。未来若再重开其他 family，仍必须先以显式新 hypothesis 重新准入，并同步更新 registry / runs manifest。 | `docs/plan/2026-04-09_dual_machine_experiment_charter.md`、`docs/experiments/registry.md`、`docs/experiments/runs.csv`、`reports/multiscale_local_skip_v2_30ep_gate_closeout_20260410.md`、`reports/alphagenome_matched_raw_track_slice_v2_closeout_20260410.md`、`reports/unet_lite_v1_rigor_review_20260409.md`、`reports/repository_status_handoff_20260409.md` |
| 官方 ChromBPNet 外置化实施 | 进行中 | 截至 `2026-04-08 11:21 CST`，`6000` 的真实 bridge smoke `chrombpnet_official_step3_bridge_smoke_20260408_110411` 已收口并**通过**。fresh evidence：`pid=486268` 已退出，日志 `/data1/zhoujiazhen/bylw_atac/logs/chrombpnet_official_step3_bridge_smoke_20260408_110411.log` 末尾出现 `Completed execution`；scratch 输出 `/data1/zhoujiazhen/bylw_atac/.codex_jobs/chrombpnet_official_step3_bridge_smoke_20260408_110411/output/negatives_with_summit.bed` 已落盘，大小约 `21M`，`wc -l=508429`；子进程命令链明确打到了 `/data1/zhoujiazhen/bylw_atac/chromBPNet/chrombpnet/helpers/make_gc_matched_negatives/make_gc_matched_negatives.sh`。这说明 helper bridge 已被真实执行验证，当前文档里原先写的 `/data1/zhoujiazhen/bylw_atac/chrombpnet_official` 仍不是已验证存在的实机路径，真实可用官方根仍是 `/data1/zhoujiazhen/bylw_atac/chromBPNet`。 | 不再追这条 smoke；下一步改成“统一 canonical root 口径”。二选一：要么把仓库文档里仍残留的 `chrombpnet_official` 统一改成 `/data1/zhoujiazhen/bylw_atac/chromBPNet`，要么先在 6000 上补一个经过验证的 `chrombpnet_official` alias，再做文档统一。完成后这条即可迁入归档。 | `tests/test_chrombpnet_official_externalization.sh`、`docs/plan/chrombpnet_official_externalization_design_20260406.md`、`docs/plan/chrombpnet_official_externalization_implementation_plan_20260406.md`、`workflows/tutorial/step3_get_background_regions.sh` |
| 论文主稿重写与 supporting writeup | 进行中 | 中文主稿已经收口到 `bias-safe framework + full/debiased diagnostics + stable readout`；与 foundation / shortcut / L3 / 分类附表相关的稳定证据现已不再单列 live 条目，统一回收到论文主稿与 supporting writeup。 | 继续把同一口径同步到英文 draft / supplementary，并压实 `L3 shared-region + 分类指标附表` 的最终表述。 | `reports/transchrombp_paper_cn_v1.tex`、`reports/transchrombp_paper_draft_v1.tex`、`reports/paper_claim_evidence_matrix_20260326.md` |
| v2fix 与 6000/6002 真实数据对照 | 待处理 | 当前稳定结论是 `B=center pool` 在多 seed 上稳定、`F` 明显失守；`6000` 的 `K562 single` held-out 弱于 `6002` 参考，但暂不能简单归因为机器差。 | 若继续追因，优先补 `GM12878` 双侧 held-out 或 `K562_6000 epoch_019` test；否则按“峰区弱于 6002”收口。 | `reports/v2fix_and_6000_realdata_followup_20260325.md` |
| 学习文档与内部长报告导航 | 待验证 | `docs/learning/` 的新分层入口已搭好，但还没做一轮真实通读校验。 | 按新导航通读一轮，确认是否仍有断链、重复入口或过时口径。 | `docs/learning/`、`reports/transchrombp_internal_design_and_experiment_history_20260331.md` |
| AI 指令体系 + Git 主档案清理 | 进行中 | 在先前 `master`/worktree 收口的基础上，本轮把“多真源心智”压成显式 contract：`docs/env/repository_governance.md` 已统一定义 local canonical trunk、6000 official lookup root、6000 runtime、6002 runtime mirror 四个角色；`scripts/sync_project.sh` 已改成 `publish-runtime-{6000,6002}`、`pull-results-{6000,6002}`、`status-all` 的显式接口，并让旧 `deploy` / `download_results` 直接报 deprecated；根目录零散单文件也已归位到 `reports/`、`docs/research/`、`docs/learning/assets/`、`scripts/`。随后又从保护分支 `master-wip-20260419-pre-governance-merge` 只筛回真正有价值的内容：factor-ladder / loss-balance 的 canonical manifest、`reports/chatgpt_bundle_loss_balance_20260416/`、`reports/repository_hygiene_cleanup_20260417.md`、bias 小资产，以及 `summarize_loss_balance_selectors.py`；同时修掉了 `.gitignore` 对 `docs/env/` 的误伤，`docs/env/transchrombp_genos_env.md` 已恢复为可版本化文档并正式纳入 Git，保护分支也已删除。 | 如果继续做整理，下一步最值得做的不是再扫一次同步脚本，而是对当前 mounted worktree 做 registry / closeout / tag 对账，决定哪些 family 已满足收口条件可以卸载；若仍想把本地仓物理目录从 `chromBPNet` 改成更直观的名字，也应单开一轮 rename，和绝对路径/链接修复一起做。 | `docs/env/repository_governance.md`、`docs/env/transchrombp_genos_env.md`、`reports/repository_governance_cleanup_20260419.md`、`docs/plan/git_cleanup_and_archive_plan_20260321.md`、`docs/plan/workspace_cleanup_plan_20260405.md` |
| foundation models 下载完整性验证 | 进行中 | `Genos-1.2B` 已完成本地加载与前向验证；其它 foundation model 仍按“最小加载可用性”原则核验。 | 网络恢复后优先补 `Genos` 的 HF 连通性 / Docker 兜底，其余模型按当前主线优先级推进。 | `docs/env/transchrombp_genos_env.md`、`docs/research/assets/nt_v2_download_20260324/README.md` |

## 二、下载资源清单（模型与数据）

说明：
- 这里记录"下载/准备了什么、从哪来、放在哪、用于什么"。
- 复现产物（如训练好的 `chrombpnet.h5`）虽不是"下载"，也记录在清单里，便于后续直接复用。

### A. 参考基因组与教程数据（ChromBPNet）

| 名称 | 服务器路径 | 来源（在哪里下载/获得） | 用途 | 当前状态 |
|---|---|---|---|---|
| hg38 参考基因组 | `/data1/zhoujiazhen/bylw_atac/hg38.fa` | 本地工作目录已有（教程中复用） | ChromBPNet 训练/解释、后续自定义模型数据处理 | 已完成 |
| hg38.fa 索引 | `/data1/zhoujiazhen/bylw_atac/hg38.fa.fai` | 本地工作目录已有 | 支持按染色体随机访问 | 已完成 |
| hg38 chrom sizes | `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes` | 教程 step1 从 ENCODE 下载 | `chrombpnet_makebigwig`、`prep nonpeaks`、`prep splits` 的必需输入 | 已完成 |
| ENCODE blacklist | `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/blacklist.bed.gz` | 教程 step1 从 ENCODE 下载 `ENCFF356LFX` | 生成 `nonpeaks` 时作为排除区 | 已完成 |
| 教程数据目录 | `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial` | `chromBPNet` tutorial workflow | 官方教程复现基线 | 已完成 |
| 教程 folds 配置 | `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/folds.json` | 教程流程生成/准备 | 染色体划分（train/valid/test） | 已完成 |
| 预生成 genomewide GC bins | `/data1/zhoujiazhen/bylw_atac/chrombpnet_refs/genomewide_gc_hg38_stride_1000_inputlen_2114.bed` | 本地预生成 | 频繁制备 nonpeaks 时节省时间 | 待处理 |

### A2. 6002 直连官方缓存（TransChromBP 准备中）

| 名称 | 服务器路径 | 来源 | 用途 | 当前状态 |
|---|---|---|---|---|
| 6002 参考基因组与基础 refs | `/home/zhengwei/bylw_atac/hg38.fa` 等 | ENCODE 公开下载 | 6002 上共用参考输入 | 进行中（hg38.fa 已重建完成） |
| 6002 教程原始 peak 输入 | `/home/zhengwei/bylw_atac/chrombpnet_tutorial/data/overlap.bed.gz` | ENCODE `ENCFF333TAT` | tutorial preprocessing 输入 | 已完成 |
| 6002 示例 BAM/peaks/nonpeaks | `/home/zhengwei/bylw_atac/chrombpnet_tutorial/data/ENCSR868FGK_merged.bam` | `chrombpnet_data` 公共示例 | 快速 smoke test | 已完成（文件在，但 BAM 不完整，不用于 tutorial 正式预处理） |
| 6002 tutorial raw BAM（rep1/2/3 + merged） | `/home/zhengwei/bylw_atac/chrombpnet_tutorial/data/rep{1,2,3}.bam`、`merged.bam` | ENCODE tutorial step1 原始 BAM | 复现 A6000 tutorial preprocessing / 生成 `merged_unstranded.bw` | 已完成 |
| 6002 tutorial bigWig | `/home/zhengwei/bylw_atac/TransChromBP/outputs/preprocessing/tutorial_canonical_v1/step2_bigwig/merged_unstranded.bw` | 6002 本地由 `rep1/2/3 -> merged.bam -> bedGraphToBigWig` 生成 | 单卡 tutorial 训练直接输入 | 已完成 |
| 6002 folds / bias model 公共缓存 | `/home/zhengwei/bylw_atac/chrombpnet_tutorial/data/folds.json`、`bias.h5` 等 | Zenodo 小文件改由 6000 补齐 | 快速复用 split 与预训练 bias model | 已完成 |

### B. ChromBPNet 复现产物（可直接复用）

| 名称 | 服务器路径 | 来源（如何得到） | 用途 | 当前状态 |
|---|---|---|---|---|
| 自定义训练可复用 bias 模型 | `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/bias/models/bias.h5` | 教程 bias 训练（step4） | `chrombpnet pipeline -b` 输入候选 | 已完成 |
| ChromBPNet 主模型 | `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/chrombpnet/chrombpnet.h5` | 教程训练（step6b 双卡） | 教程评估、解释 | 已完成 |
| ChromBPNet 无 bias 模型 | `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/chrombpnet/chrombpnet_nobias.h5` | 教程训练自动导出 | 解释和对比分析 | 已完成 |
| 缩放后 bias 模型 | `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/chrombpnet/models/bias_model_scaled.h5` | 教程 bias 训练与预处理 | ChromBPNet 因子化训练输入 | 已完成 |
| SHAP（counts/profile） | `.../evaluation/corrected.{counts,profile}_scores.h5` | 双卡并行解释脚本 | TF-MoDISco / motif 分析 | 已完成 |
| MoDISco 结果+报告（counts/profile） | `.../evaluation/modisco_{results,reports}_{counts,profile}/` | MoDISco 运行 | motif 聚类与可视化 | 已完成 |

### C. 自定义 ChromBPNet 训练候选数据

| 名称 | 服务器路径 | 来源（在哪里下载） | 用途 | 当前状态 |
|---|---|---|---|---|
| `chrombpnet_datasets/GM12878` | `/data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/GM12878` | ENCODE/BAM + overlap peaks | 后续自定义训练 | 已完成（rep1/2/3.bam + overlap.bed.gz） |
| `chrombpnet_datasets/K562` | `/data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/K562` | ENCODE/BAM + overlap peaks | 后续自定义训练 | 已完成（rep1/2.bam + overlap.bed.gz） |
| `chrombpnet_datasets/HeLa-S3` | `/data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/HeLa-S3` | SRA `prefetch` | 需自行对齐/建 bigwig/出 peaks | 部分完成（.sra/.sralite 在，未见 BAM） |
| `TransChromBP` | `/data1/zhoujiazhen/bylw_atac/TransChromBP` | 手动创建 | 自研模型统一输出目录 | 已完成 |
| `6002/chrombpnet_datasets/GM12878` | `/home/zhengwei/bylw_atac/chrombpnet_datasets/GM12878` | ENCODE 公开下载 | 6002 缓存 | 已完成（raw BAM/peaks 与 `prep_v1` 产物已被 6002 当前 loader smoke 实际消费） |
| `6002/chrombpnet_datasets/K562` | `/home/zhengwei/bylw_atac/chrombpnet_datasets/K562` | ENCODE 公开下载 | 6002 缓存 | 已完成（raw BAM/peaks 与 `prep_v1` 产物已被 6002 当前 loader smoke 实际消费） |

### D. 通用 ATAC 数据（自建数据池）

| 名称 | 服务器路径 | 来源 | 用途 | 当前状态 |
|---|---|---|---|---|
| `ATACseq/peaks` | `/data1/zhoujiazhen/bylw_atac/ATACseq/peaks` | `download_atac_datasets.sh` | 候选 peak 集合 | 已完成（~7 个峰文件） |
| `ATACseq/fastq` | `/data1/zhoujiazhen/bylw_atac/ATACseq/fastq` | SRA (`prefetch`/`fasterq-dump`) | 自主对齐 | 已完成（~18 个 FASTQ 文件） |
| `6002/ATACseq/peaks` | `/home/zhengwei/bylw_atac/ATACseq/peaks` | 同上 | 6002 缓存 | 已完成（2026-03-23 17:49 CST 已补齐 6 个候选 peak 文件） |

### E. 基因组基础模型（foundation_models）

| 名称 | 服务器路径 | 来源 | 用途 | 当前状态 |
|---|---|---|---|---|
| Caduceus-PS | `/data1/zhoujiazhen/bylw_atac/foundation_models/caduceus/` | HuggingFace `kuleshov-group/caduceus-ps_*` + GitHub `kuleshov-group/caduceus` | tutorial matched A/B 的首轮在线 token-level foundation model | 进行中（权重和运行环境已就绪，matched A/B 与 held-out `test-full` 也已完成；最终结论固定为 tutorial canonical 单 seed `near-null / marginal positive`，当前不再直接扩到 GM12878，后续只有在改变注入方式或模型族时才值得重开） |
| SegmentNT | `/data1/zhoujiazhen/bylw_atac/foundation_models/segmentnt/` | InstaDeepAI / HuggingFace | 单碱基分辨率 segmentation 备用路线 | 待处理 |
| gReLU / Borzoi / ATAC model zoo | `/data1/zhoujiazhen/bylw_atac/foundation_models/grelu/` | Genentech gReLU model zoo / HuggingFace | 中期 profile foundation / teacher 路线；含 `borzoi-model` 与 `human-atac-catlas-model` | 待处理 |
| AlphaGenome 代码+权重 | `/data1/zhoujiazhen/bylw_atac/foundation_models/alphagenome/` | GitHub + HuggingFace / API | AlphaGenome 推理/对比 | 已完成（代码与 `all-folds` 目录已在；`2026-04-07` 已通过 `ALPHAGENOME_API_KEY` 认证与 `output_metadata()` 实机 smoke；`2026-04-10` 的 `alphagenome_matched_raw_track_slice_v1_20260410` 与 `alphagenome_matched_raw_track_slice_v2_20260410` 均已完成 technical / external-coordinate closeout，并以 `pass` 收口。当前没有 active AlphaGenome 槽位；若未来重开，必须是显式新 hypothesis，且仍不升级成大 benchmark。） |
| Genos (1.2B/10B/10B-v2) | `/data1/zhoujiazhen/bylw_atac/foundation_models/genos/` | HuggingFace + GitHub `BGI-HangzhouAI/Genos` | 基因组基础模型对比 | 进行中（1.2B 已完成关键文件核对、隔离环境安装、`pip check`、本地加载、GPU `flash_attention_2` 前向与真实 usage 验证；剩余主要是 HF 可达性与 Docker 兜底镜像） |
| Nucleotide Transformer v2-500M multi-species | `/data1/zhoujiazhen/bylw_atac/foundation_models/nucleotide_transformer/nt-v2-500m-multi-species/` | HuggingFace `InstaDeepAI/nucleotide-transformer-v2-500m-multi-species` | 基因组基础模型对比、embedding/表征实验 | 进行中（目录、probe、coarse residual gate 与后续 bins16 center-aligned residual gate 都已真实消费。当前高层结论已从“独立信号较强，但还没正式 integration”收口为：probe 信号没有转成 clean integration gain；`summary / coarse residual / bins16 residual` 这条已测 family 默认停表） |
| HyenaDNA（中等体量） | `/data1/zhoujiazhen/bylw_atac/foundation_models/hyenadna/` | HuggingFace / GitHub `HazyResearch/hyena-dna` | 长上下文备用路线 | 待处理 |
| Gengram-10B | `/data1/zhoujiazhen/bylw_atac/foundation_models/gengram/` | HuggingFace | 基因组基础模型对比 | 待验证 |
| AlphaGenome 专用环境 | `/data1/zhoujiazhen/bylw_atac/.mamba/envs/alphagenome` | 6000 本地创建 | Python>=3.10 环境 | 已完成 |
| Nucleotide Transformer 专用环境 | `/data1/zhoujiazhen/bylw_atac/.mamba/envs/nucleotide-transformer-py311` | 6000 本地创建 | Python 3.11 + HF/Transformers + smoke/probe | 待验证（目录已开始创建，但本轮实际 smoke/probe 已证明现有 `genos-1.2b` 环境足以完成加载与特征提取，因此 `py311` 环境暂时不再是 blocker；仅在后续需要更干净的长期维护环境时再补齐） |

### F. Transformer 学习代码仓库

| 名称 | 服务器路径 | 来源 | 当前状态 |
|---|---|---|---|
| nanochat | `/data1/zhoujiazhen/bylw_atac/nanochat` | GitHub `karpathy/nanochat` | 已完成（仓库已克隆） |
| modded-nanogpt | `/data1/zhoujiazhen/bylw_atac/modded-nanogpt` | GitHub `KellerJordan/modded-nanogpt` | 已完成（浅克隆，依赖未装） |

## 三、已知问题与判断规则（避免重复踩坑）

- 自定义训练不要为 `bigWig` 单独安排"下载任务"；`chrombpnet pipeline` 会从 `bam/fragment/tagAlign` 现场生成。真正值得提前缓存的是 `hg38.fa(+.fai)`、`hg38.chrom.sizes`、`blacklist.bed.gz`、`folds.json`、可复用 `bias.h5`。
- AlphaGenome 适合做外部 black-box 对比，不适合把 API 输出拿来训练学生模型；优先做 `1000s` 级别小规模对照。
- 6000 上 `chrombpnet` 环境是 `Python 3.9.23`，与 AlphaGenome SDK 的 `Python >= 3.10` 不兼容；AlphaGenome 使用独立环境 `/data1/zhoujiazhen/bylw_atac/.mamba/envs/alphagenome`。
- `chrombpnet_deepshap`、`chrombpnet_convert_html_to_pdf` 在当前环境里可能不存在；优先使用 Python 模块入口。
- `bias_qc` / `chrombpnet_qc` 历史失败常见原因是环境未激活（`modisco: not found`）；必须先执行 `CLAUDE.md` 中的环境激活命令。
- `foundation_models` 的 HF 下载目录里的 `.lock` 文件不等于下载失败；是否可用以"最小加载验证"为准。
- Genos 官方 `requirements.txt` 当前把 `torchaudio` 写成了 `>=0.22.1,<0.23.0`；按原样安装会直接报 `No matching distribution found`，当前临时修正为 `torchaudio==2.7.1` 继续装，后续如官方更新应回归上游文件为准。
- Genos 安装 `flash-attn` 时，直接按 `pip install -r requirements.txt` 走 build isolation 会在构建前置阶段报 `ModuleNotFoundError: No module named 'torch'`；更稳的处理路径应是先装好 `torch`，再对 `flash-attn` 单独使用 `--no-build-isolation`。

## 四、2026-03-17: 开发环境与版本管理配置

### 完成事项
1. Git 仓库初始化：本地与远程服务器 (6000) 均已关联 GitHub 仓库 `rose8happy/TransChromBP`；`.gitignore` 排除 `*.h5, *.pt`；统一 Git 用户身份 `yangmeisuan <345687960@qq.com>`。
2. 自动化同步脚本：`scripts/sync_project.sh` 支持 `deploy` / `download_results`。
3. 创建了 `DEVELOPMENT.md` 详细说明开发工作流。
