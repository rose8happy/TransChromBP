# 已完成事项归档

> 从 `TRACKING.md` 第一节"在做事情清单"中迁出的 `已完成` / `已合并` 条目。
> 首次归档时间：2026-03-18
> 最近一次追加归档：2026-04-03

---

## 2026-04-03 第六轮归档（6000 空窗短评估补齐）

| 事项 | 关键结论 | 关键路径 |
|---|---|---|
| `GM12878_6000 pilot20` held-out `test-full` sidecar | `2026-04-03 18:38-18:40 CST` 已完成；overall=`0.5517/0.7746`、peak=`0.4231/0.8179`、nonpeak=`0.6802/0.0762`。结果与历史 6000 记录实质一致，说明 `GM12878_6000` 的 held-out 偏弱不是瞬时波动，而是可复现现象 | `reports/v2fix_and_6000_realdata_followup_20260325.md` |
| `K562_6000 epoch_019` held-out `test-full` sidecar | `2026-04-03 20:38-20:40 CST` 已完成；peak 从历史 `best.pt` 的 `0.6186/0.8287` 变为 `0.6213/0.8414`，count 提升约 `+0.0127`，但 profile JSD 变差约 `+0.0027`。结论是后期 checkpoint 延续了 profile/count tradeoff，但不足以推翻当前 `best.pt` 的 profile-first 选择 | `reports/v2fix_and_6000_realdata_followup_20260325.md` |

---

## 2026-03-28 第五轮归档（clean matrix 补证收口）

| 事项 | 关键结论 | 关键路径 |
|---|---|---|
| `notf_sg1_deb0_s1234` test 补评（6002/3080） | `split=test` 已补齐；peak=`0.43731/0.43841`、gap=`0.01679`、count_r=`0.82869/0.82669`。与 `s42=0.02682` 合并后，conv-only unsafe 双 seed 形成 `0.02682 / 0.01679`，说明它仍是 clean matrix 最高风险格，但幅度存在 seed 波动，不支持“稳定复现强 shortcut”的更强表述 | `reports/profile_shortcut_revalidation_summary_20260326.md` |

## 2026-03-25 第四轮归档（foundation model 缓存验证）

| 事项 | 关键结论 | 关键路径 |
|---|---|---|
| Nucleotide Transformer v2-500M 本地缓存 | 6000 直连 `huggingface.co` 不通，改用“本机解析直链 + 远端下载大文件、小文件本机抓取后写回”；`model.safetensors` 已落盘，并在现有 `genos-1.2b` 环境下跑通本地 `transformers` CPU smoke test | `docs/research/nucleotide-transformer_download_turol.md`, `/data1/zhoujiazhen/bylw_atac/foundation_models/nucleotide_transformer/nt-v2-500m-multi-species/` |

## 2026-03-24 第三轮归档（实验收口与基础设施验证）

| 事项 | 关键结论 | 关键路径 |
|---|---|---|
| TransChromBP validation 指标对齐改造（JSD/count_r/MAE） | 新口径训练器已在 6000 dryrun 跑通，train/val/best-metric/checkpoint 链路可用 | `vendor/transchrombp/transchrombp/training/train_ddp.py` |
| 6002 上历史数据集缓存补齐 | GM12878/K562 已落盘并通过 loader smoke 验证；HeLa-S3/ATAC fastq 暂不推进 | `scripts/start_6002_dataset_cache_downloads.sh` |
| 6000/6002 ChromBPNet 数据集 CPU 预处理（GM12878/K562） | 两边后台作业均已结束，产物已生成；6002 侧已补 loader smoke 验证 | `scripts/run_remote_chrombpnet_dataset_prep.sh` |
| 6000 A6000 vs 6002 RTX3080 训练速度对比 | 历史分析已收口：V2-full 相比 V2-noTF 在 held-out peak 上稳定更优 | `reports/assets/ablation_tf_20260318/summary_table.csv` |
| 40 epoch 尾部收益复盘 | A/F best epoch=40，B=34，G=15；30→40 对 A/F 有收益但不支持统一补到 50 | `reports/v2fix_readout_head_status_20260323.md` |
| TransChromBP tutorial 真实数据基线重训（bs16 + peak-best） | 已被 bias_pretrain → main 主线替代 | — |
| GM12878 自定义 ChromBPNet 复现（全流程） | 该数据集仅作 TransChromBP 输入，不再单独扩 ChromBPNet 复现 | — |
| Genos 融合实验 Phase 0/P2 | P0=0.3179/0.8384 未超 G0；P2 count-only 不通过（count 断崖下滑）；不推进 P1 | `reports/genos_cached_fusion_status_20260323.md`、`reports/genos_no_positive_gain_analysis_20260323.md` |
| Genos OCR 官方基准复现 | layer 12 roc_auc=0.7535 接近官方 0.7569，仅作 sanity check | `docs/research/genos复现ocr.md` |
| nanochat / modded-nanogpt 学习环境 | 仓库已克隆，低优先级暂缓 | — |

## 2026-03-21 第二轮归档（项目管理与运行支撑）

| 事项 | 关键结论 | 关键路径 |
|---|---|---|
| 6000 A6000x2 本地项目工作区整理 | 本地专用工作区 `/home/zhengwei/project/python/server6000_a6000x2` 已按远端一级目录语义整理完成，可作为后续 6000 项目的统一落点 | `/home/zhengwei/project/python/server6000_a6000x2/README.md` |
| 6002 上 TransChromBP 官方输入直连缓存 | tutorial 核心输入、folds 与 bias 模型缓存均已就位，后续新实验可直接复用 | `scripts/start_6002_transchrombp_official_downloads.sh`, `/home/zhengwei/bylw_atac/logs/tutorial_raw_prep_20260318_174729.log` |
| 6002 承担 TransChromBP 实验的单卡准备 | 6002 单卡训练链路、依赖工具和 smoke test 已跑通，3080 可承接短训练、小消融与 smoke run | `scripts/prepare_6002_transchrombp_single_gpu.sh`, `/home/zhengwei/bylw_atac/logs/teacher_v2_6002_bias_smoke_retry_20260318_211241.log` |
| 6002 Bias-Safe 机制消融（2×2 析因：pool\_factor × stop\_gradient） | 4 组实验与分析已收口；`stop_gradient=true` 是关键，`pool_factor` 无显著主效应 | `reports/transchrombp_bias_safe_ablation_20260320.tex` |
| TransChromBP 训练性能诊断（2×A6000 benchmark） | clean run `1GPU bs16 -> 2GPU bs16` 吞吐 `191.0 -> 361.0 samples/s`，scaling efficiency `94.5%`，训练明确 compute-bound | `reports/transchrombp_training_perf_benchmark_20260320.tex`, `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/logs/benchmark/benchmark_summary_20260320_004701.txt` |
| 本地报告/画图独立环境（`.venv-report`） | 独立报告环境已建立并完成核心依赖验证，避免分析依赖污染主训练环境 | `scripts/setup_report_env.sh`, `DEVELOPMENT.md` |
| Ablation 训练轮数复试（30→40 epoch） | 已并入 `V2 代码改进消融实验（修订版）`；第一批实验统一按 40 epoch 执行，不再单列跟踪 | `docs/plan/v2_code_improvement_ablations.md` |

## ChromBPNet 官方复现线

| 事项 | 关键结论 | 关键路径 |
|---|---|---|
| ChromBPNet 教程复现训练（step6b 双卡） | 双卡训练在第 12 轮早停正常结束 | `chrombpnet_tutorial/logs/step6b_train_only_gpu_mgpu.log` |
| ChromBPNet 教程后处理（step7，解释图/TF-MoDISco） | SHAP 打分+MoDISco 完成，motifs.html/pdf 已生成 | `chrombpnet_tutorial/logs/chrombpnet_qc_{counts,profile}_20260202_112230.log` |
| 官方数据线复现验收 | `counts pearsonr=0.7274`, `profile median_jsd=0.3360` | `chrombpnet_tutorial/outputs/chrombpnet/evaluation/chrombpnet_metrics.json` |
| ChromBPNet 论文严格对齐复现（5-fold × 1-seed） | 均值 `counts pearson peaks=0.7007`, `profile median jsd peaks=0.3428` | `chrombpnet_paper_repro/runs/K562_ATAC/summary_fast_seed_1234/` |
| Bias 模型解释图（step5）补跑 | 完整产物已存在 | `chrombpnet_tutorial/outputs/bias/evaluation/` |
| 官方教程线 Markdown 汇报文档（含图） | 可直接浏览 | `chrombpnet_tutorial/TUTORIAL_REPRO_REPORT_20260304.md` |
| 官方教程线"指标解释版"报告（含论文对照） | PDF 已导出 | `chrombpnet_tutorial/TUTORIAL_REPRO_REPORT_EXPLAINED_20260304.pdf` |

## ChromBPNet 前置准备

| 事项 | 关键结论 | 关键路径 |
|---|---|---|
| ChromBPNet 项目目录复制到 Windows | 整仓同步完成，`28656` 文件无差异 | `/mnt/d/project/python/chromBPNet` |
| ChromBPNet 自定义训练前置资源梳理 | 真正需提前准备的是 `genome fasta/chrom.sizes/peaks/bias.h5` | README.md, workflows/ |
| AlphaGenome API 对比实验可行性核实 | 适合做 zero-shot 对比，不能用于训练学生模型；需 Python>=3.10 独立环境 | alphagenomedocs.com |
| AlphaGenome pilot 4 位点对比（6000） | AlphaGenome 严重低估绝对计数（~4.6x），但相对排序正确 | `chromBPNet/outputs/alphagenome_pilot/tutorial_selected_loci_20260315/` |
| AlphaGenome 论文/代码阅读与启发整理 | 吸收三点：长程 Transformer 需验证、distillation 作为方法学变量、scorer 和消融需同步设计 | — |

## TransChromBP 基础设施

| 事项 | 关键结论 | 关键路径 |
|---|---|---|
| 长度×分辨率实验框架 | sweep 脚本+汇总脚本已新增，最小 dry-run 通过 | `TransChromBP/scripts/run_length_resolution_sweep.sh` |
| 模型与实验设计文档化 | 统一设计文档已完成 | `TransChromBP/docs/model_experiment_design.md` |
| 消融实验框架（Transformer × 偏置分解 2x2） | 4 个消融配置+一键运行脚本，dry-run 通过 | `TransChromBP/scripts/run_ablations.sh` |
| 训练策略落地（AdamW 分组 + Warmup Cosine） | 参数分组+cosine 调度+best checkpoint 选择，21 passed | `TransChromBP/src/transchrombp/training/train_ddp.py` |
| 训练语义修正（PID 句柄+epoch nonpeak 重采样+early stop） | 四项修正写回远端并验证，22 passed | `TransChromBP/src/transchrombp/` |
| bias 使用重构 | pretrained_path 修复，fusion 改为 add/logsumexp，19 passed | `TransChromBP/src/transchrombp/models/transchrombp.py` |
| batch size / data loader 吞吐基准 | `bs=16` 最优（364 samples/s），`num_workers>=2` 即可 | `TransChromBP/outputs/logs/batch_size_benchmark/summary.csv` |
| canonical 数据语义固化 | canonical_v1 配置已落地 | `TransChromBP/configs/train/train_tutorial_canonical_v1.yaml` |
| canonical 数据预处理重跑 | step2+step3+hyperparams 完成，manifest.json 已写出 | `TransChromBP/outputs/preprocessing/tutorial_canonical_v1/` |
| preprocessing→teacher 自动衔接 hook | watcher 成功检测并拉起 teacher 流水线 | `TransChromBP/scripts/watch_preprocess_and_launch_teacher.sh` |

## TransChromBP 训练实验

| 事项 | 关键结论 | 关键路径 |
|---|---|---|
| tutorial 真实数据后端接入 | `ChromBPNetBigWigDataset` 完成，13 passed | `TransChromBP/src/transchrombp/data/real_data.py` |
| tutorial 真实数据基线训练（双卡 20 epoch） | `train_loss=5.80`, `val loss_total=5.83` | `TransChromBP/outputs/checkpoints/tutorial_real_baseline_20260313_1740/` |
| tutorial bias→main 正式流水线 | bias 10 epoch + main 20 epoch 完成 | `TransChromBP/outputs/checkpoints/tutorial_bias2main_20260315_023539_main_learnable/` |
| tutorial 独立 test 评估与基线对照 | bias→main 明显优于 baseline：peak count_r 0.8207 vs 0.8048 | `TransChromBP/outputs/metrics/tutorial_test_comparison_20260315.csv` |
| tutorial 深化对照分析（bias-only + locus 图） | bias 在 nonpeak 有信号(count_r=0.4633)，peak 几乎失效(-0.0089) | `reports/transchrombp_tutorial_consolidated_cn.tex` |
| tutorial hybrid_data 两阶段对齐实验 | 流水线完成但定位为"语义绑包+optimization confound 负对照" | `logs/tutorial_bias2main_hybrid_20260315_1032.pipeline.log` |
| tutorial hybrid_data 公平 held-out 评估 | hybrid 落后原主线（count_r 0.8000 vs 0.8347），含 count_weight 混杂 | `TransChromBP/outputs/metrics/tutorial_hybrid_*_test_full.json` |
| canonical teacher 自动首跑 | 高价值失败样本：count_weight×learnable_scales 触发 bias shortcut | `logs/transchrombp_tutorial_teacher_canonical_v1_auto_postfix_20260316_165421.pipeline.log` |
| count_weight × learnable_scales 诊断与 fixcw 对照 | fixcw (cw=0.1) 彻底解决坍塌；peak count_r=0.8205 (+12.8% vs ChromBPNet)；发现 profile bias shortcut | `TransChromBP/outputs/metrics/fixcw_main_learnable_best_test.json` |
| hybrid data semantics 审计 | 定位为"语义绑包+optimization confound 对照"，不再扩训 | `TransChromBP/configs/train/train_tutorial_hybrid_data.yaml` |

## TransChromBP Teacher v2 系列

| 事项 | 关键结论 | 关键路径 |
|---|---|---|
| Teacher v2 设计与实现（修复 profile bias shortcut） | v2 修复成功：debiased JSD 从 0.5633 降至 0.3163；peak count_r 0.8466 | `TransChromBP/outputs/metrics/teacher_v2_main_learnable_best_test.json` |
| Teacher v2 深化分析（full held-out + checkpoint + AlphaGenome shape） | epoch_024 在 test 上更优；4 位点 profile JSD 与 AlphaGenome 持平 | `TransChromBP/outputs/metrics/teacher_v2_main_learnable_epoch024_test_full.json` |
| Teacher v2 late-phase checkpoint sweep（epoch 18–24） | epoch 24 综合最优：peak profile JSD=0.3141, peak count_r=0.8465 | `TransChromBP/outputs/metrics/teacher_v2_main_learnable_epoch0{18..24}_test.json` |
| checkpoint 选择机制优化分析 | best_metric 应从 peak.loss_total 切换到更贴近主表的指标 | `TransChromBP/src/transchrombp/training/train_ddp.py` |

## 报告与文档

| 事项 | 关键结论 | 关键路径 |
|---|---|---|
| Teacher v2 深层机制分析报告 | 10 页 LaTeX：通道非对称性、为何同时改善 count、设计变更分解 | `reports/transchrombp_teacher_v2_deep_analysis_20260318.tex` |
| Teacher v2 代码逐段 LaTeX 讲解 | 面向小白的 LaTeX 讲解稿 | `reports/transchrombp_teacher_v2_code_walkthrough_20260317.tex` |
| fixcw 跨模型综合分析报告 | 9 页 LaTeX：含 AlphaGenome、发现 profile bias shortcut | `reports/fixcw_crossmodel_analysis_20260317.tex` |
| 阶段报告与论文主线重构 | 7 页 LaTeX 源稿，研究主线切换为 Transformer+distillation+bias 关系 | `reports/transchrombp_tutorial_consolidated_cn.tex` |
| 报告冻结快照（远端保护副本） | 快照已保存，硬链接保护 | `TransChromBP/outputs/snapshots/tutorial_report_20260315_frozen/` |
| 数据语义审计文档与冻结归档 | 9 个 checkpoint 快照+18 份评估 JSON | `TransChromBP/outputs/snapshots/data_semantics_audit_20260315_frozen/` |
| 推荐数据处理流程定稿 | count_weight 不再作为 canonical 默认项 | `reports/transchrombp_tutorial_consolidated_cn.tex` |
| 论文 LaTeX 骨架（英文+中文） | 围绕"Transformer×distillation×bias factorization" | `reports/transchrombp_paper_outline.tex`, `reports/transchrombp_paper_outline_cn.tex` |
| 中文项目思路稿 | 收束为长程建模+蒸馏+偏差分解 | `reports/transchrombp_project_note_cn.tex` |
| 本地 WSL LaTeX 环境 | pdflatex/xelatex/latexmk 可用 | — |

## 其他

| 事项 | 关键结论 | 关键路径 |
|---|---|---|
| 6000 远端 CLI 工具补装 | rg/fd/jq/bat/fzf 已装到 ~/.local/bin | `/home/zhoujiazhen/.local/bin/` |
