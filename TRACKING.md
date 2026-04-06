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
| foundation-model restart v3：NT v2 cached residual short10（6000） | 进行中 | 当前正式 run 为 `ntv2_residual_short10_s42_20260405_dualcache`，已从“单卡 cache + 双卡训练”升级为“cache 也双卡”；核心目标是验证 `NT v2 cached residual_head` 在 tutorial short10 上是否能在不牺牲强 baseline 的前提下带来稳定增益。 | 持续跟日志确认 dual-cache 阶段完成、short10 自动起跑且保持 `GPU0,1` 双卡；run 结束后先读 `best.pt + epoch_metrics.jsonl` 再决定是否扩 `seed1234` 或换注入方式。 | `docs/plan/foundation_model_restart_v3_20260405.md`、`reports/foundation_model_restart_v3_implementation_20260405.md`、`/data1/zhoujiazhen/bylw_atac/logs/ntv2_residual_short10_dualcache_20260405_6000.log` |
| foundation model gate 总结（tutorial canonical） | 进行中 | 当前 gate 已基本稳定：`NT v2` 技术可行且独立信号强于 Genos，但与 baseline 互补性不足；`Caduceus-PS` 单 seed A/B 仅 `near-null / marginal positive`；`Genos cached P0/P1` 为负结果。除当前 restart v3 外，不建议再把 tutorial canonical 训练窗口继续投给这三条旧接入线。 | 将 `NT v2`、`Caduceus`、`Genos cached` 的阶段结论统一并入 foundation-model restart bundle / 外部咨询材料。 | `reports/nt_v2_probe_20260404.md`、`docs/plan/caduceus_tutorial_ab_online_fusion_20260331.md`、`reports/genos_cached_p1_restart_20260403.md` |
| 官方 ChromBPNet 外置化实施 | 待验证 | design 与 implementation plan 已执行完成；本仓已切断 active local-official 依赖、补齐 patch ledger、删除本地 `chrombpnet/` 与官方 packaging 入口。 | 在 6000 用一次真实 official compare / dataset-prep smoke 验证桥接链路；若通过，再把该条移入 `TRACKING_archive.md`。 | `docs/plan/chrombpnet_official_externalization_design_20260406.md`、`docs/plan/chrombpnet_official_externalization_implementation_plan_20260406.md`、`reports/chrombpnet_official_patch_ledger_20260406.md` |
| ChromBPNet 官方锚点 + shared-region L3 收口 | 进行中 | strict-compare 的主结论已稳定：tutorial `L3 shared-region` 上，TransChromBP 在 held-out test 的 profile、count 和分类判别三侧都优于 official；剩余工作是把 `L3` 数字和方法口径写回论文主表、claim matrix 和附表。 | 只做写回与表格同步，不再为 tutorial `L3` 额外开新训练。 | `reports/tutorial_L3_shared_region_closure_20260330.md`、`docs/plan/project_roadmap_20260330.md` |
| ATAC 分类指标扩展（AUROC/AUPRC/F1） | 进行中 | peak-vs-nonpeak 分类指标已改为依附外部官方 `predict/metrics` 路径契约、selector 和汇总脚本，不再引用已删除的本地 `chrombpnet/training/predict.py`；当前定位仍是 supplementary 指标，不进入主表。 | 把分类指标补进 strict-compare 附表和论文 supplementary，并继续沿外部官方 `predict/metrics` 契约校验。 | `reports/chrombpnet_official_patch_ledger_20260406.md`、`scripts/paper_aligned_repro/select_best_epoch.py`、`reports/tutorial_L3_shared_region_closure_20260330.md` |
| 论文主稿重写与 supporting writeup | 进行中 | 主线已切到“中文主稿优先 + 双轨主表（历史 tutorial baseline + L3 shared-region）+ 分类指标附表 + clean matrix / no-bias 补证”；当前剩余主要是 overfull、作者单位/致谢/仓库 URL 等写作收口。 | 继续收紧中文稿排版并补齐作者信息；除非 reviewer 明确追问，否则不建议再新开 no-bias 或 tutorial L3 实验。 | `reports/transchrombp_paper_cn_v1.tex`、`reports/transchrombp_paper_draft_v1.tex`、`reports/paper_writeup_flow_20260330.md` |
| clean matrix / profile shortcut 叙事收口 | 进行中 | clean matrix 双 seed 主证据已齐，结论已固定为“Transformer 特有 shortcut 不成立，风险来自 supervision-driven 的 bias-safe 失稳；conv-only unsafe 风险最高但幅度有 seed 波动”。 | 只做论文与总报告写回，不再额外扩 shortcut 实验。 | `reports/profile_shortcut_revalidation_summary_20260326.md`、`reports/paper_claim_evidence_matrix_20260326.md` |
| v2fix 与 6000/6002 真实数据对照 | 待处理 | 当前可以收敛的结论是：`B=center pool` 在多 seed 上稳定，`F` 明显失守；6000 的 `K562 single` held-out 弱于 6002 参考，但暂不能简单归因为机器差。 | 若继续追因，优先补 `GM12878` 双侧 held-out 或 `K562_6000 epoch_019` test；否则按“峰区弱于 6002”收口。 | `reports/v2fix_and_6000_realdata_followup_20260325.md` |
| 学习文档与内部长报告导航 | 待验证 | `docs/learning/` 已统一到“基础材料补概念，当前模型导读解释主线，内部长报告解释项目如何演化到今天”的分层入口；旧教程入口的重要性已被降级。 | 按新导航实际通读一轮，确认是否还存在断链、重复入口或过时口径。 | `docs/learning/`、`reports/transchrombp_internal_design_and_experiment_history_20260331.md` |
| AI 指令体系 + Git 主档案清理 | 进行中 | `AGENTS.md` 作为主规范的方案已落地；当前这轮正在继续做“`vendor/transchrombp/` 作为正式 snapshot、`tmp_remote_edit/` 只保留 staging、报告源码优先、live/archive/report 分层收口”的仓库清理。 | 完成分批提交、归档和 ignore 收紧，确保本地仓库重新成为唯一主档案。 | `reports/ai_instruction_consolidation_20260328.md`、`docs/plan/git_cleanup_and_archive_plan_20260321.md`、`DEVELOPMENT.md` |
| 网页版 ChatGPT 咨询材料整理 | 进行中 | paper-facing 咨询包与 foundation restart 咨询包均已在本地固定；当前缺的是新的外部判断，而不是再堆更多背景文档。 | 等用户上传 `foundation restart` bundle 并回传外部意见，再决定是否正式重启新主线。 | `reports/chatgpt_bundle_foundation_restart_20260331/04_external_analysis_brief.md`、`reports/chatgpt_bundle_20260326/03_current_bottleneck_and_questions.md` |
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
| Nucleotide Transformer v2-500M multi-species | `/data1/zhoujiazhen/bylw_atac/foundation_models/nucleotide_transformer/nt-v2-500m-multi-species/` | HuggingFace `InstaDeepAI/nucleotide-transformer-v2-500m-multi-species` | 基因组基础模型对比、embedding/表征实验 | 进行中（当前目录已被实际加载与 tutorial `valid` probe 消费。高层结论是：NT v2 的独立信号强于 Genos，但直接补当前 corrected-B baseline 的互补性仍不足，本轮不建议直接进入 tutorial A/B 训练接入） |
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
