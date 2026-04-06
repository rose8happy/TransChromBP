# 项目追踪清单（工作清单 + 下载资源清单）

本文件用于减少支线过多导致的信息遗漏。每次开始新任务或汇报前，先查看并按最新状态更新本文件。

> 已完成 / 已合并事项已迁至 [TRACKING_archive.md](TRACKING_archive.md)。
> 第一节原则上只保留 `进行中` / `待处理` / `待验证` 条目；阶段性完成或已并入其它主线的事项转入 `TRACKING_archive.md`。

状态约定：
- `已完成`：可直接使用
- `进行中`：任务仍在跑或尚未收尾
- `待处理`：明确要做但还没开始
- `待验证`：文件/目录已在，但需要实际加载或运行验证

## 一、在做事情清单（实时更新）

| 事项 | 当前状态 | 当前结论 / 进度 | 下一步 | 关键路径 |
|---|---|---|---|---|
| 仓库任务盘点（2026-04-06） | 进行中 | 本轮盘点后的 P0 已在同一轮内执行：foundation 停表规则已统一写回 `post_chatgpt_pro_priority_execution`、`foundation_model_restart_v3` 计划/实施记录、claim matrix 与双语论文主稿；英文稿显式 `\todo{}` 也已清零，中文稿作者占位已移除。`2026-04-06 14:33 CST` 再次实时复核 6000 后，两张 A6000 仍为空闲；同日下午 `workspace cleanup` 也已继续推进一轮小清理：删掉了本机 `Zone.Identifier` / `.checked` / 一次性 gate sidecar，并确认 `tmp_remote_edit/transchrombp` 已全部能在 `vendor/` 找到对应副本、6002 继续维持 runtime mirror 角色。随后又形成并写回了新的结构性结论：根仓不应再长期保留官方 `chrombpnet` clone，而应切到“本仓只保留 `TransChromBP` 主线，官方 ChromBPNet 外置到 6000 参考仓”的两阶段迁移设计；用户已确认该设计，且实施计划也已补写完成。当前更适合继续推进的已不再是 foundation 扩线，而是 `workspace cleanup`、双语主稿打磨与少量工程尾项。 | 默认转入 `工作区整理 -> 论文润色/参考文献扩充 -> 小型工程尾项`；cleanup 主线已新增一个上位设计门：先按 `chrombpnet_official_externalization_implementation_plan_20260406` 切断本仓对本地 `chrombpnet/` 的活依赖，再决定何时真正删除本地官方源码树。若用户坚持使用 A6000，当前也只建议做 paper-facing 的短时收口/评估，而不是给 residual short10 family 再开新训练。 | `reports/project_plan_code_review_20260405.md`、`docs/plan/project_roadmap_20260330.md`、`docs/plan/workspace_cleanup_plan_20260405.md`、`docs/plan/chrombpnet_official_externalization_design_20260406.md`、`docs/plan/chrombpnet_official_externalization_implementation_plan_20260406.md` |
| 官方 ChromBPNet 外置化实施 | 待验证 | design and implementation plan executed; active local-official dependencies cut; patch ledger added; local chrombpnet payload and packaging entrypoints removed. | run one real 6000 official compare / dataset-prep smoke to validate the bridged chain; if that passes, archive the row. | `docs/plan/chrombpnet_official_externalization_design_20260406.md`、`docs/plan/chrombpnet_official_externalization_implementation_plan_20260406.md`、`reports/chrombpnet_official_patch_ledger_20260406.md` |
| `NT v2 bins16 center-aligned residual gate`（6000 / A6000×2） | 进行中 | 用户显式授权的唯一 genuinely new-hypothesis 重入已经完成：`ntv2_bins16_centerres_short10_s42_20260406_dual` 于 `2026-04-06 03:25 CST` 左右写出 full held-out JSON `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/ntv2_bins16_centerres_short10_s42_20260406_dual_best_test_full_20260406_032529.json`。结果相对 matched `short10_nofoundation_control` 明显更差：overall `count_pearson_full=0.7727` vs `0.8457`、`profile_target_jsd_full_mean=0.4625` vs `0.4395`、`profile_pearson_full_mean=0.6580` vs `0.7116`；peak `count_pearson_full=0.7516` vs `0.8298`、`profile_target_jsd_full_mean=0.3588` vs `0.3193`、`profile_pearson_full_mean=0.7837` vs `0.8641`。这一结论在本轮已同步写回 gate 文档、claim matrix 与双语论文主稿，因此当前高层结论已稳定为：即便升级到“中心对齐的 `bins16` token residual”，NT v2 residual short10 family 仍未给出 clean gain；6000 两张 A6000 现已空闲。 | 不再把这条 family 当作待执行后续链；后续若仍想重开 foundation，必须先提出不同于当前 `summary / token-fusion / coarse residual / bins16 residual` family 的新 hypothesis。 | `docs/plan/nt_v2_bins16_center_aligned_gate_20260405.md`、`reports/foundation_model_restart_v3_implementation_20260405.md`、`reports/paper_claim_evidence_matrix_20260326.md`、`reports/transchrombp_paper_cn_v1.tex` |
| `NT v2 gate smoke after idle`（6000 / A6000×2） | 进行中 | `record_sha1 mismatch` 已在本轮定位并修复：根因是 `build_foundation_cache.py` 对非 `train` split 用 `seed+10000`，而 `evaluate_checkpoint.py` 的 held-out dataset 用 `seed+20000`，导致 `test + max_test_regions` 子采样记录不一致；同时 `ChromBPNetBigWigDataset` 未暴露 `max_records`，使 manifest 一直把该字段写成 `0`。修复后已在 6000 用新增回归测试通过验证，并于 `2026-04-05 20:59 CST` 用新 run `ntv2_residual_short10_gate_smoke_20260405_fix1` 成功重跑 bounded smoke：`manifest_test.json` 现已写出 `max_records=64`，held-out gate 也成功产出 JSON，不再报 sha1 mismatch。进一步读取结果后确认：这次 bounded smoke 证明的是 gate 链路已恢复，而不是模型信号转正；`n=64` 下 `count_pearson_full=-0.0980`、`profile_pearson_full_mean=-0.0142`、`profile_target_jsd_full_mean=0.6668`、`peak_auroc=0.5543`，足以把当前 `bins4 / coarse residual` family 判为 `unsafe / not-worth-expanding`。 | 把“链路已修 + 当前 family 判负”写回 gate 文档与报告；除非用户显式授权 genuinely new-hypothesis 重入，否则不再为当前配方开新 run。 | `reports/project_plan_code_review_20260405.md`、`docs/plan/a6000_foundation_reentry_candidates_20260405.md`、`vendor/transchrombp/transchrombp/data/real_data.py`、`vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py`、`vendor/transchrombp/transchrombp/scripts/build_foundation_cache.py` |
| 当前项目审查（计划 + 代码，2026-04-05） | 进行中 | 当前这轮审查已经从“发现 launcher 闭环缺口”推进到“完成 root cause 修复并在真实 6000 smoke 上复验”，随后又完成了 result-level stop/go 收口与文档写回。新增证据链是：`tests/test_foundation_cache_alignment.py` 先稳定复现 `test` split cache builder / evaluator 的 `record_sha1` 分叉与 `max_records` 元数据缺失，随后修复 `resolve_dataset_seed` + `self.max_records` 后回归测试转绿，真实 smoke 完整写出 held-out JSON；再与 `short10_nofoundation_control` 的 full held-out 对照后，当前风险已从“链路不通”切换为“当前 foundation recipe 本身判负”，且该结论已同步进报告与论文口径。 | 若继续推进，这条审查线下一步不再是补链路，而是把相关代码重构和 archive/cleanup 继续收口；结果层面的 stop/go 判断已稳定。 | `reports/project_plan_code_review_20260405.md`、`tests/test_foundation_cache_alignment.py`、`vendor/transchrombp/transchrombp/data/real_data.py`、`vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py`、`vendor/transchrombp/transchrombp/scripts/build_foundation_cache.py` |
| A6000 利用率 + foundation-model 停表规则 | 进行中 | 现在资源和结论都已同步收口：`2026-04-06 14:33 CST` 实时 `nvidia-smi` 复核显示 6000 两张 A6000 均为 `0% util / 13 MiB`、仅剩驱动占用，且没有 `zhoujiazhen` 名下训练进程；与此同时，`short10_nofoundation_control` (`0.8457 / 0.4395 / 0.7116`) 与两条 NT v2 residual gate 的对照已经足够形成稳定 stop-rule。旧 `bins4/coarse residual` bounded gate 为 `-0.0980 / 0.6668 / -0.0142`；新的 `bins16 center-aligned residual` full held-out 也只有 `0.7727 / 0.4625 / 0.6580`（peak `0.7516 / 0.3588 / 0.7837`），仍明显落后于 matched baseline。当前合理结论已从“只停旧 family”收紧为“当前 NT v2 residual short10 family 默认停表”。 | 默认不再为当前 foundation family 分配 A6000；若一定要用卡，优先做 paper-facing 的短时收口/评估或其它已确认高价值任务。后续若还想重开 foundation，必须先提出不属于当前 residual short10 family 的全新假设，并在文档里明确为何不违反现有停表规则。 | `docs/plan/post_chatgpt_pro_priority_execution_20260405.md`、`docs/plan/a6000_foundation_reentry_candidates_20260405.md`、`reports/foundation_model_restart_v3_implementation_20260405.md` |
| `short10 matched no-foundation control`（6000 / A6000×2） | 进行中 | full held-out 已完成并写出 `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/short10_nofoundation_control_s42_20260405_dual_best_test_full_20260405_2100.json`；日志 `/data1/zhoujiazhen/bylw_atac/logs/short10_nofoundation_control_testfull_20260405_6000.log` 已收口。核心结果：`n=62342`，`count_pearson_full=0.8457`，`count_pearson_debiased=0.8419`，`profile_target_jsd_full_mean=0.4395`，`profile_pearson_full_mean=0.7116`，分类指标 `auroc=0.8677`、`auprc=0.8647`、`f1=0.7902`。这条 run 已经回答了关键问题：`short10` 预算本身仍能维持健康 baseline，因此它现在是固定校准底座，而不是待继续验证的实验。 | 把这条 matched baseline 固定写进 foundation gate 文档；除非后续真开 `bins16` 新假设，否则不再重跑。 | `docs/plan/a6000_foundation_reentry_candidates_20260405.md`、`vendor/transchrombp/transchrombp/scripts/run_short10_no_foundation_control.sh` |
| ChatGPT Pro 外部分析回流与行动化 | 进行中 | 外部意见已明确支持把论文主线收口到 `bias-safe framework + diagnostics`，并把 foundation 线降为“有停表规则的 side quest”；还特别提醒当前 `FoundationResidualHead` 更像在测“粗 summary residual”，不能过度外推成所有 token/base foundation 路线都被否掉。 | 把这份判断写回论文主线、claim matrix 和 foundation gate 文档；后续所有 foundation 决策都先对照该报告检查是否又回到 dramatic shortcut 或无停表扩线。 | `reports/chatgpt_pro_external_review_20260405.md`、`docs/plan/post_chatgpt_pro_priority_execution_20260405.md`、`reports/chatgpt_bundle_project_handoff_20260405/07_foundation_adapter_current.py` |
| 论文主稿重写与 supporting writeup | 进行中 | 中文主稿已把摘要、贡献点、foundation 讨论与结论进一步收口到 `bias-safe framework + full/debiased diagnostics + stable readout`；foundation 线也已明确降到 appendix / future work 的受控 side quest，且不再只写成 `Genos` 单线负结果。 | 继续把同一口径同步到英文 draft / supplementary，并把 `L3 shared-region + 分类指标附表` 的最终表述压实。 | `reports/transchrombp_paper_cn_v1.tex`、`reports/transchrombp_paper_draft_v1.tex`、`reports/paper_claim_evidence_matrix_20260326.md` |
| foundation model gate 总结（tutorial canonical） | 进行中 | 当前 gate 已固定：`Genos cached P0/P1` 为负结果，`Caduceus-PS` 单 seed A/B 为 `near-null / marginal positive`，`NT v2 restart v3` 的 held-out `test-full` 触发 `fail-or-unsafe`，而新增的 matched `short10 no-foundation control` 又证明短预算本身仍可维持健康 baseline。当前新增结论不是“foundation 仍可继续试”，而是“当前 `FoundationResidualHead + NT v2 cached short10 + bins4/coarse-summary residual` 方案已经完成必要验证并判负；若再进 foundation，只能作为显式授权下的新假设重入，而不是旧 family 的自然延长”。 | 把 foundation 线统一收口成 appendix / side quest 负结果；默认不再为 tutorial canonical 自动扩 matched control、第二个 seed 或 `cross_attention`，只有新假设 `bins16` 路线在显式授权下保留。 | `reports/nt_v2_probe_20260404.md`、`reports/genos_cached_p1_restart_20260403.md`、`reports/chatgpt_pro_external_review_20260405.md`、`docs/plan/a6000_foundation_reentry_candidates_20260405.md` |
| ChromBPNet 官方锚点 + shared-region L3 收口 | 进行中 | strict-compare 的主结论已稳定：tutorial `L3 shared-region` 上，TransChromBP 在 held-out test 的 profile、count 和分类判别三侧都优于 official；剩余工作是把 `L3` 数字和方法口径写回论文主表、claim matrix 和附表。 | 只做写回与表格同步，不再为 tutorial `L3` 额外开新训练。 | `reports/tutorial_L3_shared_region_closure_20260330.md`、`docs/plan/project_roadmap_20260330.md` |
| ATAC 分类指标扩展（AUROC/AUPRC/F1） | 进行中 | peak-vs-nonpeak 分类指标已改为依附外部官方 `predict/metrics` 路径契约与相关报告/脚本工件，不再引用已删除的本地 `chrombpnet/training/predict.py`；当前定位仍是 supplementary 指标，不进入主表。 | 把分类指标补进 strict-compare 附表和论文 supplementary，并继续沿外部官方 predict/metrics contract 校验。 | `docs/plan/chrombpnet_official_externalization_implementation_plan_20260406.md`、`reports/chrombpnet_official_patch_ledger_20260406.md`、`tests/test_chrombpnet_official_externalization.sh` |
| clean matrix / profile shortcut 叙事收口 | 进行中 | clean matrix 双 seed 主证据已齐，结论已固定为“Transformer 特有 shortcut 不成立，风险来自 supervision-driven 的 bias-safe 失稳；conv-only unsafe 风险最高但幅度有 seed 波动”。 | 只做论文与总报告写回，不再额外扩 shortcut 实验。 | `reports/profile_shortcut_revalidation_summary_20260326.md`、`reports/paper_claim_evidence_matrix_20260326.md` |
| v2fix 与 6000/6002 真实数据对照 | 待处理 | 当前可以收敛的结论是：`B=center pool` 在多 seed 上稳定，`F` 明显失守；6000 的 `K562 single` held-out 弱于 6002 参考，但暂不能简单归因为机器差。 | 若继续追因，优先补 `GM12878` 双侧 held-out 或 `K562_6000 epoch_019` test；否则按“峰区弱于 6002”收口。 | `reports/v2fix_and_6000_realdata_followup_20260325.md` |
| 学习文档与内部长报告导航 | 待验证 | `docs/learning/` 已统一到“基础材料补概念，当前模型导读解释主线，内部长报告解释项目如何演化到今天”的分层入口；旧教程入口的重要性已被降级。 | 按新导航实际通读一轮，确认是否还存在断链、重复入口或过时口径。 | `docs/learning/`、`reports/transchrombp_internal_design_and_experiment_history_20260331.md` |
| `chrombpnet-remote` 自定义 skill | 待验证 | 已确认 `~/.codex/skills/chrombpnet-remote/SKILL.md` 存在且当前 Codex 会话已识别到该 skill；其最小内容覆盖 6000/6002 主机映射、6000 环境激活、GPU 余量检查、长任务后台启动、`TRACKING.md` 回写、远端原子写回与写后校验规则。当前剩余缺口不是“是否装上”，而是“是否经得起一次真实远端任务链路前向试用”。 | 在下一次真实 6000/6002 任务里前向试用一轮；若仍频繁手写相同上传/替换流程，再决定是否补 `scripts/` 辅助工具。 | `~/.codex/skills/chrombpnet-remote/SKILL.md`、`AGENTS.md` |
| AI 指令体系 + Git 主档案清理 | 进行中 | `AGENTS.md` 作为主规范的方案已落地；本轮除了收口 foundation 结果，也已在隔离 worktree `~/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-structure` 连续完成三阶段结构整理：先把 `train_ddp.py` / `evaluate_checkpoint.py` 对 cached foundation feature 的解释与校验抽成共享 `foundation cache contract` helper，再把 `build_foundation_cache.py` 与 `run_ntv2_residual_short10.sh` 的 cache request/preflight 接到同一合同层，随后又把 `build_foundation_cache.py` 与 held-out evaluator 的 `split/max_records/region_source/seed` 规则改成共用 helper。与此并行的 workspace cleanup 新进展是：`tmp_remote_edit/transchrombp` 已确认完全镜像到 `vendor/`，但 `tmp_remote_edit/TransChromBP/scripts` 仍有 11 个历史 helper script 未迁移，因此 staging 还不能直接压成 README-only；更上层的 repo-role 设计也已定稿：本仓将不再长期保留官方 `chrombpnet` clone，而是改成“`TransChromBP` 主仓 + 文档档案”，官方 ChromBPNet 外置到 6000。该结构性迁移的实施计划也已补写完成。最小 gate 已通过：`unittest` 8 项、`py_compile` 与 `bash -n` 全绿。 | 工程主线继续给 launcher 增加 manifest 级 preflight，提前阻断“模型配置已切换但旧 cache 仍被误复用”的情况；cleanup 主线则分成两层：一层先处理那 11 个历史 helper script，另一层按 `chrombpnet_official_externalization_implementation_plan_20260406` 切断本仓对本地 `chrombpnet/` 的活依赖，再决定何时彻底删除本地官方源码树。 | `docs/plan/git_cleanup_and_archive_plan_20260321.md`、`docs/plan/workspace_cleanup_plan_20260405.md`、`docs/plan/chrombpnet_official_externalization_design_20260406.md`、`docs/plan/chrombpnet_official_externalization_implementation_plan_20260406.md`、`~/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-structure/docs/plan/foundation_cache_contract_refactor_20260406.md` |
| foundation models 下载完整性验证 | 进行中 | `Genos-1.2B` 已完成本地加载和前向验证；其他 foundation model 仍按“最小加载可用性”原则核验，不以目录是否存在为准。 | 网络恢复后优先补 `Genos` 的 HF 连通性 / Docker 兜底，其余模型按当前主线优先级推进。 | `docs/env/transchrombp_genos_env.md`、`docs/research/assets/nt_v2_download_20260324/README.md` |

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
| AlphaGenome 代码+权重 | `/data1/zhoujiazhen/bylw_atac/foundation_models/alphagenome/` | GitHub + HuggingFace / API | AlphaGenome 推理/对比 | 代码与 `all-folds` 目录已在；API access / 实际加载待验证 |
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
