# ChromBPNet / TransChromBP 学习资料总索引

> 最后更新：2026-04-18

这份文件是全量导航版。

如果你只想按问题快速跳，请先看：

- [学习材料索引.md](学习材料索引.md)

如果你想知道整个目录应该怎么分层使用，请先看：

- [学习资料说明.md](学习资料说明.md)

## 1. 当前目录的完整分层

### A. 入口与导航

| 文件 | 定位 | 当前用途 |
|---|---|---|
| `学习资料说明.md` | 整个目录的分层说明 | 先建立阅读方法 |
| `学习材料索引.md` | 按问题 / 目标跳转 | 快速找答案 |
| `小白学习指南.md` | 从零开始的最短路线 | 新人或久别回看用 |

### B. 补课与当前模型层

| 文件 | 定位 | 当前用途 |
|---|---|---|
| `深度学习补课手册.md` | 深度学习基础概念速成 | 先补 `QKV/JSD/loss/seed` |
| `TransChromBP当前模型学习总览.md` | 当前主线 5 分钟全景图 | 快速跟上当前模型 |
| `TransChromBP当前模型详解.md` | 代码、训练、评估、配置导读 | 看代码和配实验 |

### C. 专题层

| 文件 | 定位 | 当前用途 |
|---|---|---|
| `TransChromBP项目实验全景图.md` | family 级实验总表 | 看项目做过什么、结论怎样 |
| `无大模型基线到当前主线.md` | 当前最佳非 foundation 主线导读 | 只看效果好的主线模型 |
| `基座模型接入路线图.md` | foundation-model 接入专题 | 看 Genos / NT v2 / Caduceus / AlphaGenome |

### D. 基础与历史层

| 文件 | 定位 | 当前用途 |
|---|---|---|
| `BPNet模型小白教程.md` | BPNet 任务和结构基础 | 补 profile/count、膨胀卷积 |
| `BPNet与ChromBPNet模型对比.md` | bias factorization 直觉 | 理解为什么会有 bias branch |
| `ChromBPNet数据处理完全指南.md` | 数据处理和采样语义 | 查 peak/nonpeak、jitter、revcomp |
| `ChromBPNet完整项目深度解析.md` | 旧版 ChromBPNet 工程导读 | 历史背景，不是当前实现 |

### E. 仓库外长报告层

| 文件 | 定位 | 当前用途 |
|---|---|---|
| `../../reports/transchrombp_internal_design_and_experiment_history_20260331.md` | 长篇内部复盘 | 从设计到论文口径的连续时间线 |
| `../../reports/paper_claim_evidence_matrix_20260326.md` | claim / evidence 对照 | 判断什么能写、什么不能写 |
| `../../reports/assets/paper_metric_source_table_20260326.csv` | 主数字真源 | 查论文级核心数字 |

## 2. 最推荐的四条阅读路径

### 路线 A：我不懂深度学习，但想读懂项目

1. `深度学习补课手册.md`
2. `小白学习指南.md`
3. `无大模型基线到当前主线.md`
4. `TransChromBP当前模型学习总览.md`
5. `TransChromBP当前模型详解.md`

### 路线 B：我想看实验全景和版本差异

1. `TransChromBP项目实验全景图.md`
2. `../../reports/transchrombp_internal_design_and_experiment_history_20260331.md`
3. `../../reports/paper_claim_evidence_matrix_20260326.md`

### 路线 C：我想只看最好的几份模型

1. `无大模型基线到当前主线.md`
2. `TransChromBP当前模型学习总览.md`
3. `../../reports/assets/paper_metric_source_table_20260326.csv`

### 路线 D：我想看大模型是怎么接进来的

1. `基座模型接入路线图.md`
2. `TransChromBP项目实验全景图.md`
3. `../../reports/paper_claim_evidence_matrix_20260326.md`

## 3. 当前最常见的 12 个问题

| 问题 | 优先看哪里 |
|---|---|
| 当前默认模型是什么 | `无大模型基线到当前主线.md` |
| 为什么 `center pool` 成为默认 | `无大模型基线到当前主线.md` |
| `full/debiased` 是什么 | `深度学习补课手册.md` + `TransChromBP当前模型详解.md` |
| Transformer 的 `QKV` 怎么理解 | `深度学习补课手册.md` |
| `KL/JSD/NLL/MSE` 在项目里分别出现在哪里 | `深度学习补课手册.md` |
| 做过哪些主要实验 | `TransChromBP项目实验全景图.md` |
| 哪些实验是真正的主线证据 | `TransChromBP项目实验全景图.md` + `../../reports/paper_claim_evidence_matrix_20260326.md` |
| foundation 线为什么没成功升主线 | `基座模型接入路线图.md` |
| Genos / NT v2 / Caduceus / AlphaGenome 分别怎么接 | `基座模型接入路线图.md` |
| no-bias 为什么只放 supplementary | `TransChromBP项目实验全景图.md` |
| 最近 loss balance 做了什么 | `TransChromBP项目实验全景图.md` |
| 真正该先看哪段代码 | `TransChromBP当前模型详解.md` |

## 4. 当前最值得看的代码与配置

### 模型

| 文件 | 作用 |
|---|---|
| `../../vendor/transchrombp/transchrombp/models/transchrombp.py` | 当前主模型定义 |
| `../../vendor/transchrombp/transchrombp/models/transformer_encoder.py` | Transformer 编码器 |
| `../../vendor/transchrombp/transchrombp/models/bias_branch.py` | bias branch 实现 |
| `../../vendor/transchrombp/transchrombp/models/genos_adapter.py` | Genos 支线入口 |
| `../../vendor/transchrombp/transchrombp/models/caduceus_adapter.py` | Caduceus 支线入口 |
| `../../vendor/transchrombp/transchrombp/models/foundation_adapter.py` | foundation 通用挂载入口 |

### 训练与评估

| 文件 | 作用 |
|---|---|
| `../../vendor/transchrombp/transchrombp/training/train_ddp.py` | loss、validation、best metric、训练语义 |
| `../../vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py` | held-out test 与 classification 指标 |
| `../../vendor/transchrombp/scripts/select_best_epoch.py` | 自研侧 external selector |
| `../../scripts/paper_aligned_repro/select_best_epoch.py` | official 侧 selector |

### 当前最有代表性的配置

| 文件 | 作用 |
|---|---|
| `../../vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool.yaml` | 当前默认主线 |
| `../../vendor/transchrombp/transchrombp/configs/model/ablations/ablation_no_bias.yaml` | no-bias 边界补证 |
| `../../vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_caduceus_ps.yaml` | Caduceus token fusion |
| `../../vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_ntv2_residual.yaml` | NT v2 residual 路线 |
| `../../vendor/transchrombp/transchrombp/configs/train/train_tutorial_corrected_b_strict_compare_L3_6000.yaml` | `L3 shared-region` 外部比较主配置 |

## 5. 一句话总结

今天的 `docs/learning/` 最好理解为：

> 用“补课手册 + 三份专题 + 当前模型导读 + 长报告”共同组成的学习闭环，而不是靠某一份旧教程把所有问题全讲完。
