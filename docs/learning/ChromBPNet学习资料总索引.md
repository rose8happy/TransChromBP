# ChromBPNet / TransChromBP 学习资料总索引

> 最后更新：2026-03-31

这份文件是全量导航版。  
如果你只想按问题找文档，请优先看：

- [学习材料索引.md](学习材料索引.md)

如果你想知道整个目录现在怎么分层使用，请先看：

- [学习资料说明.md](学习资料说明.md)

## 1. 当前目录的完整分层

### A. 内部复盘层

| 文件 | 定位 | 当前用途 |
|---|---|---|
| `reports/transchrombp_internal_design_and_experiment_history_20260331.md` | 长篇内部时间线复盘 | 阶段性回看设计、实验、代码与论文口径演化 |

### B. 当前模型层

| 文件 | 定位 | 当前用途 |
|---|---|---|
| `TransChromBP当前模型学习总览.md` | 当前主线 5 分钟全景图 | 快速跟上当前模型与结论 |
| `TransChromBP当前模型详解.md` | 当前模型、训练、评估、配置、代码入口详解 | 看代码、配实验、想论文 |
| `学习材料索引.md` | 按问题和目标的快速入口 | 边查边读 |

### C. 基础材料层

| 文件 | 定位 | 当前用途 |
|---|---|---|
| `学习资料说明.md` | 目录说明和使用方式 | 理解整个目录该怎么用 |
| `小白学习指南.md` | 从零开始的最短路线 | 新人或自己久别重看时重建认知 |
| `BPNet模型小白教程.md` | profile/count、膨胀卷积、任务定义基础 | 补 BPNet 概念 |
| `BPNet与ChromBPNet模型对比.md` | bias factorization 基础直觉 | 理解为什么要保留 bias branch |

### D. 历史 / 概念补充层

| 文件 | 定位 | 当前用途 |
|---|---|---|
| `ChromBPNet完整项目深度解析.md` | 旧版 ChromBPNet 项目结构与流程 | 历史背景 |
| `ChromBPNet数据处理完全指南.md` | peak/nonpeak、jitter、revcomp 等概念 | 数据语义手册 |

## 2. 最推荐的阅读路径

### 路线 A：阶段性回看项目全貌

1. `reports/transchrombp_internal_design_and_experiment_history_20260331.md`
2. `TransChromBP当前模型学习总览.md`
3. `TransChromBP当前模型详解.md`
4. `reports/paper_claim_evidence_matrix_20260326.md`

### 路线 B：边看代码边想论文

1. `reports/transchrombp_internal_design_and_experiment_history_20260331.md`
2. `TransChromBP当前模型详解.md`
3. `vendor/transchrombp/transchrombp/models/transchrombp.py`
4. `vendor/transchrombp/transchrombp/training/train_ddp.py`
5. `reports/paper_writeup_flow_20260330.md`
6. `reports/transchrombp_paper_cn_v1.tex`

### 路线 C：从零补基础

1. `学习资料说明.md`
2. `小白学习指南.md`
3. `BPNet模型小白教程.md`
4. `BPNet与ChromBPNet模型对比.md`
5. `TransChromBP当前模型学习总览.md`
6. `TransChromBP当前模型详解.md`

### 路线 D：只想快速跟上当前主线

1. `TransChromBP当前模型学习总览.md`
2. `TransChromBP当前模型详解.md`
3. `reports/tutorial_L3_shared_region_closure_20260330.md`
4. `reports/paper_claim_evidence_matrix_20260326.md`

## 3. 当前最常见的 10 个问题

| 问题 | 优先看哪里 |
|---|---|
| 当前默认模型到底是什么？ | `TransChromBP当前模型学习总览.md` |
| 这个项目是怎么一步步走到今天的？ | `reports/transchrombp_internal_design_and_experiment_history_20260331.md` |
| 为什么不再写“Transformer 特有 shortcut”？ | `reports/transchrombp_internal_design_and_experiment_history_20260331.md` |
| `full` / `debiased` 在代码里怎么落地？ | `TransChromBP当前模型详解.md` |
| 为什么 `center pool` 成了默认？ | `TransChromBP当前模型详解.md` |
| strict compare 为什么最后看 `L3`？ | `reports/tutorial_L3_shared_region_closure_20260330.md` |
| no-bias 为什么不进主表？ | `reports/paper_writeup_flow_20260330.md` |
| 当前哪些 claim 还能写？ | `reports/paper_claim_evidence_matrix_20260326.md` |
| 主数字到底从哪里来？ | `reports/assets/paper_metric_source_table_20260326.csv` |
| 我该先看哪段代码？ | `TransChromBP当前模型详解.md` |

## 4. 当前最值得看的代码与配置

### 模型

| 文件 | 作用 |
|---|---|
| `vendor/transchrombp/transchrombp/models/transchrombp.py` | 主模型定义 |
| `vendor/transchrombp/transchrombp/models/transformer_encoder.py` | Transformer 编码器 |
| `vendor/transchrombp/transchrombp/models/bias_branch.py` | bias branch 实现 |

### 训练与评估

| 文件 | 作用 |
|---|---|
| `vendor/transchrombp/transchrombp/training/train_ddp.py` | loss、validation、best metric、训练语义 |
| `vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py` | held-out test 与 classification 指标 |
| `vendor/transchrombp/scripts/select_best_epoch.py` | 自研侧 external selector |
| `scripts/paper_aligned_repro/select_best_epoch.py` | official 侧 external selector |

### 当前最有代表性的配置

| 文件 | 作用 |
|---|---|
| `vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool.yaml` | 当前默认 paper-facing 配置 |
| `vendor/transchrombp/transchrombp/configs/model/ablations/ablation_no_bias.yaml` | no-bias 边界补证 |
| `vendor/transchrombp/transchrombp/configs/train/train_ablation_v2_main_profile_select.yaml` | paper-facing 训练/选模语义 |
| `vendor/transchrombp/transchrombp/configs/train/train_tutorial_corrected_b_strict_compare_L3_6000.yaml` | `L3 shared-region` 主配置 |

## 5. 当前默认口径的最短摘要

今天最稳的目录级共识是：

1. Transformer 的收益成立。
2. `full/debiased gap` 是必要诊断。
3. `debiased_profile_weight=2.0` 比单独 sg 更像主效应。
4. `center pool` 是当前默认 readout。
5. `L3 shared-region` 是当前最强 external evidence。
6. no-bias 只是 supplementary boundary result。

## 6. 如果你要带新人，推荐固定发哪 4 份

1. `学习资料说明.md`
2. `小白学习指南.md`
3. `TransChromBP当前模型学习总览.md`
4. `TransChromBP当前模型详解.md`

如果对方不是要入门，而是要快速理解项目为什么写成今天这样，再补一份：

5. `reports/transchrombp_internal_design_and_experiment_history_20260331.md`

## 7. 和旧版目录最大的不同

以前这个目录更像“ChromBPNet 教学资料 + 若干当前补丁”。  
现在它已经明确改成：

- **内部复盘**：解释项目怎么演化
- **当前模型导读**：解释当前主线是什么
- **基础材料**：解释概念和历史来源

所以今天最不推荐的做法就是：

- 从旧的 `ChromBPNet完整项目深度解析.md` 直接开始，以为那就是当前实现

更好的做法是：

- 想回看全局：先看内部长报告
- 想看当前实现：先看总览和详解
- 想补概念：再回到基础材料
