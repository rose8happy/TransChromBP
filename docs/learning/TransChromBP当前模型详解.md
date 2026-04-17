# TransChromBP 当前模型详解

这份文档不再只是“模型结构说明书”，而是服务两个目标：

1. 你想阶段性回到代码细节时，知道该看哪些文件、哪些符号。
2. 你想阶段性思考论文时，知道这些代码入口分别对应哪些实验结论和哪些当前口径。

如果你还没看过总览或内部时间线，建议先看：

- [TransChromBP当前模型学习总览.md](TransChromBP当前模型学习总览.md)
- [reports/transchrombp_internal_design_and_experiment_history_20260331.md](../../reports/transchrombp_internal_design_and_experiment_history_20260331.md)

如果你当前的目标更偏向某个专题，建议先分流：

- 想先补 `QKV / KL / JSD / loss`：看 [深度学习补课手册.md](深度学习补课手册.md)
- 想只看当前最好的非 foundation 主线：看 [无大模型基线到当前主线.md](无大模型基线到当前主线.md)
- 想看所有实验 family：看 [TransChromBP项目实验全景图.md](TransChromBP项目实验全景图.md)
- 想看 Genos / NT v2 / Caduceus / AlphaGenome：看 [基座模型接入路线图.md](基座模型接入路线图.md)

## 1. 当前模型的最短定义

今天的 `TransChromBP` 可以压缩成下面这句话：

> 一个以 `ConvStem + LocalDilatedTower + optional Transformer` 为 signal backbone、以 ChromBPNet-style `bias_branch` 为辅助分支、同时输出 `full` 和 `debiased` 两套 profile/logcount，并以 `center pool` 作为当前默认 count readout 的统一框架。

核心代码入口：

- [vendor/transchrombp/transchrombp/models/transchrombp.py](../../vendor/transchrombp/transchrombp/models/transchrombp.py)

## 2. 模型结构应该怎么读

当前主模型实现里最重要的是 6 个层次。

### 2.1 `ConvStem`

代码入口：

- `ConvStem`：`vendor/transchrombp/transchrombp/models/transchrombp.py`

职责：

- 从 `[B, L, 4]` 的 one-hot DNA 输入提取初级 motif 特征
- 作为整个 signal backbone 的起始层

为什么它重要：

- 当前主线从来都不是“直接把 one-hot DNA 扔给 Transformer”
- 卷积起始层保留了 BPNet / ChromBPNet 在局部 motif 建模上的 inductive bias

### 2.2 `LocalDilatedTower`

代码入口：

- `LocalDilatedTower`：`vendor/transchrombp/transchrombp/models/transchrombp.py`

职责：

- 在 Transformer 前做 BPNet 风格的局部上下文建模
- 保留 dilated residual block 的结构归纳偏置

为什么它重要：

- 当前 full model 不是“纯 Transformer”
- `noTF` 对照组只关闭 `sequence_encoder.enabled`，但保留 `ConvStem + LocalDilatedTower`
- 这让 matched ablation 更公平，也让 backbone 结论更可信

### 2.3 `SequenceTransformerEncoder`

代码入口：

- [transformer_encoder.py](../../vendor/transchrombp/transchrombp/models/transformer_encoder.py)
- 在 [transchrombp.py](../../vendor/transchrombp/transchrombp/models/transchrombp.py) 里由 `sequence_encoder.enabled` 控制是否启用

职责：

- 在 token 级别建模更长距离依赖

当前常见配置：

- `d_model=256`
- `n_heads=8`
- `n_layers=6`
- `ff_mult=4`
- `use_rope=true`
- `use_sdpa=true`

与实验结论的关系：

- 这是 `V2-full` vs `V2-noTF` backbone ablation 的核心开关
- 当前可以保住“Transformer 有真实收益”，就是因为这部分被 matched ablation 稳住了

### 2.4 signal heads

代码入口：

- [transchrombp.py](../../vendor/transchrombp/transchrombp/models/transchrombp.py)

当前 signal 分支显式输出两类量：

1. `profile_signal`
2. `count_signal`

其中：

- profile 头负责位置级别输出
- count 头负责 region-level 标量输出

这两条路径后来分别对应了不同的实验敏感性：

- profile 端更容易暴露 full/debiased gap
- count 端更容易在 readout 或外部 summary 注入时被直接改坏

### 2.5 `count_pool_mode`

代码入口：

- `heads.count_pool_mode`：`vendor/transchrombp/transchrombp/models/transchrombp.py`

当前支持三种聚合方式：

| 模式 | 做法 | 当前地位 |
|---|---|---|
| `full` | 全长度池化 | 历史基线或对照 |
| `center` | 只对中心窗口池化 | 当前默认 |
| `attention` | 学位置权重再聚合 | 研究性分支，不再是默认 |

为什么 `center` 这么重要：

- `B_s42` 和 `B_s1234` 在 held-out 上都稳定
- `F=attention pool` 在第二个 seed 上 profile 明显失守
- 因此今天读代码时看到 `count_pool_mode=center`，应直接把它理解为当前默认 paper-facing readout，而不是一个普通备选项

对应默认配置：

- [transchrombp_teacher_v2_center_pool.yaml](../../vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool.yaml)

### 2.6 `bias_branch`

代码入口：

- [bias_branch.py](../../vendor/transchrombp/transchrombp/models/bias_branch.py)
- [transchrombp.py](../../vendor/transchrombp/transchrombp/models/transchrombp.py)

职责：

- 输出 `profile_bias`
- 输出 `count_bias`

关键配置项：

- `bias_branch.enabled`
- `bias_branch.hidden_channels`
- `bias_branch.n_dil_layers`
- `bias_branch.profile_pool_factor`
- `bias_branch.pretrained_path`
- `bias_branch.freeze_bias_core`

为什么它还重要：

- 项目虽然做了 no-bias 补证，但 no-bias 只形成了 supplementary boundary result
- 当前主线仍然是 bias-factorized 框架，不是完全去掉 bias branch 的框架

## 3. `full` / `debiased` 到底是什么

这是整个项目今天最不应该跳过的地方。

当前模型一次前向会显式返回两套结果：

| 输出 | 含义 |
|---|---|
| `profile_logits_full` / `logcount_full` | signal + bias 融合后的完整输出 |
| `profile_logits_debiased` / `logcount_debiased` | 只看 signal branch 的去偏输出 |

为什么这件事今天这么关键：

- 项目的主诊断不再是某个单独 scale 轨迹
- 而是 `full`、`debiased` 以及它们之间的 gap

因此：

- 只报 full，不报 debiased，不够
- 只看最终 test 指标，不看 gap，也不够

这也是 clean matrix 为什么能改写论文叙事的原因。

## 4. 融合逻辑与 `stop-gradient`

### 4.1 profile 融合

代码入口：

- [transchrombp.py](../../vendor/transchrombp/transchrombp/models/transchrombp.py)

逻辑上可以理解成：

```text
profile_full = profile_signal + profile_scale * profile_bias
```

如果打开：

- `fusion.profile_bias_stop_gradient=true`

那么 bias profile 在参与 full profile 融合前会先 `detach()`。

要怎么理解它：

- bias 仍然参与输出
- 但 full profile loss 不会顺着这条路径回流到 bias branch

当前项目对它的最终评价是：

- 它是重要的辅助性 bias-isolation 设计
- 但 clean matrix 之后，不能再把它写成“唯一关键修复开关”

### 4.2 count 融合

当前 count 端常用的是 log-space 融合语义，目的是让 signal count 和 bias count 的组合更符合 count 的量纲直觉。

这条路径在实验里的意义是：

- count 分支的训练动力学和 profile 分支并不完全一样
- 某些外部 summary 注入（例如 Genos `P2`）会直接在这里把 count 路径改坏

## 5. 当前训练语义应该怎么读

核心代码入口：

- [train_ddp.py](../../vendor/transchrombp/transchrombp/training/train_ddp.py)

最值得直接看的函数有：

| 函数/逻辑 | 为什么重要 |
|---|---|
| `compute_losses()` | 定义 full/debiased 各项 loss 怎么组合 |
| `run_validation()` | 决定 validation 指标怎样聚合、怎样形成 full/debiased gap |
| `resolve_best_metric_name()` | best checkpoint 选模逻辑 |
| `apply_training_mode_defaults()` | 不同训练模式下的默认语义 |
| `validate_semantics_profile()` | 对训练/数据语义做一致性检查 |

### 5.1 当前最关键的 loss 组合

逻辑上可以写成：

```text
total =
  profile_weight * profile_loss(full)
  + count_weight * count_loss(full)
  + debiased_profile_weight * profile_loss(debiased)
  + debiased_count_weight * count_loss(debiased)
```

当前主线推荐值：

| 项 | 当前主线值 |
|---|---:|
| `profile_weight` | `1.0` |
| `count_weight` | `0.1` |
| `debiased_profile_weight` | `2.0` |
| `debiased_count_weight` | `0.0` |

为什么 `debiased_profile_weight=2.0` 今天必须重点记：

- clean matrix 之后，当前最稳的结论是它比单独 `sg` 更像主效应
- 如果你以后开新实验，最不该随手动掉的就是它

### 5.2 当前 best metric 为什么不是随便选的

项目后来越来越强调 paper-facing peak JSD，是因为：

- 只盯 `loss_total` 不够贴近最后要汇报的主指标
- 某些 run 会出现 profile/count tradeoff
- strict compare 和 no-bias 的补证都逐渐走向 external selector + held-out test 的口径

因此今天应该把：

- `trainer.best_metric`
- external selector metric

一起看，而不是只看训练日志里的单一总 loss。

## 6. 当前评估链路应该怎么读

### 6.1 held-out test 与统一评估

代码入口：

- [evaluate_checkpoint.py](../../vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py)

这份文件今天非常关键，因为它已经不仅输出：

- profile 指标
- count 指标

还会输出：

- peak-vs-nonpeak 分类指标

当前最值得盯的点：

| 逻辑 | 位置 |
|---|---|
| 评估数据集如何构建 | `build_eval_dataset()` |
| 输出 JSON 路径如何决定 | `make_output_path()` |
| 分类指标如何计算 | `peak_auroc_*` / `peak_best_f1_*` 一组逻辑 |

这就是为什么：

- `L3` 的 classification supplementary 能闭环
- no-bias 的 held-out test 也能直接产出 AUROC/AUPRC/F1

### 6.2 自研侧 external selector

代码入口：

- [vendor/transchrombp/scripts/select_best_epoch.py](../../vendor/transchrombp/scripts/select_best_epoch.py)

今天读它时最值得关注的不是脚本语法，而是：

1. 它按什么 metric 选 best epoch
2. 它是否同时带出 profile/count/classification 字段
3. 它如何分摊到两张 GPU 做多 checkpoint selector

这份脚本后来对 no-bias 和 strict compare 都很关键。

### 6.3 official 侧 selector

代码入口：

- [scripts/paper_aligned_repro/select_best_epoch.py](../../scripts/paper_aligned_repro/select_best_epoch.py)

它的作用是：

- 对 official ChromBPNet 侧也用统一外部规则选 best epoch
- 后来还补进了 `classification_metrics` 字段的抽取

这就是为什么 `L3 shared-region` 现在可以真正做官方侧和自研侧对齐。

## 7. 实验线与代码 / 配置应该怎么对上

如果你现在想从实验名回到代码，可以先按“改的是哪一层”来找，而不是先背 run 名。

| 实验类别 | 先看哪 1-2 个文件 | 代表配置 | 你只需要先看懂什么 |
|---|---|---|---|
| backbone ablation | `models/transchrombp.py`、`models/transformer_encoder.py` | `ablations/ablation_no_transformer.yaml` | `sequence_encoder.enabled` 是不是打开 |
| readout design | `models/transchrombp.py` | `transchrombp_teacher_v2_center_pool.yaml`、`v2fix_attn_pool.yaml` | `count_pool_mode` 怎么影响 count 头 |
| bias / clean matrix | `models/transchrombp.py`、`training/train_ddp.py` | `transchrombp_teacher_v2_nosg.yaml` 等 clean-matrix configs | `profile_bias_stop_gradient` 和 `debiased_profile_weight` 是怎样一起生效的 |
| loss / selector | `training/train_ddp.py`、`scripts/select_best_epoch.py` | `train_ablation_v2_main_profile_select.yaml` | `best_metric`、`compute_losses()`、`run_validation()` |
| strict compare / held-out test | `evaluation/evaluate_checkpoint.py`、`scripts/paper_aligned_repro/select_best_epoch.py` | `train_tutorial_corrected_b_strict_compare_L3_6000.yaml` | external selector 和 held-out 输出字段 |
| no-bias boundary | `models/transchrombp.py` | `ablations/ablation_no_bias.yaml` | 去掉 bias branch 后还有哪些输出保留 |
| Genos 接入 | `models/genos_adapter.py`、`models/transchrombp.py` | `v2fix_genos_gate.yaml`、`train_genos_cached_short10.yaml` | summary / FiLM / count 注入挂在哪里 |
| Caduceus 接入 | `models/caduceus_adapter.py`、`training/train_ddp.py` | `transchrombp_teacher_v2_center_pool_caduceus_ps.yaml` | token-level 在线融合发生在 Transformer 前面 |
| NT v2 / foundation adapter | `models/foundation_adapter.py`、`scripts/build_foundation_cache.py` | `transchrombp_teacher_v2_center_pool_ntv2_residual.yaml` | cache / residual / adapter 路线怎样挂进主干 |

如果你想按结论而不是按代码看，直接读：

- [TransChromBP项目实验全景图.md](TransChromBP项目实验全景图.md)
- [基座模型接入路线图.md](基座模型接入路线图.md)

## 8. 当前哪些结论已经从“探索”变成“默认”

### 已经是当前默认理解

1. Transformer 是有效 backbone。
2. `full/debiased gap` 是必要诊断。
3. `debiased_profile_weight=2.0` 是当前最该保住的训练语义之一。
4. `center pool` 是当前默认 readout。
5. `L3 shared-region` 是当前 external comparison 主证据。

### 已经不该再这样理解

1. “Transformer 特有 shortcut”
2. “stop-gradient 单独修复了问题”
3. “conv-only 不会 shortcut”
4. “no-bias 证明 bias branch 绝对必要”

## 9. 如果现在要读代码，推荐顺序是什么

### 第一层：先搞清楚模型骨架

1. `vendor/transchrombp/transchrombp/models/transchrombp.py`
2. `vendor/transchrombp/transchrombp/models/transformer_encoder.py`
3. `vendor/transchrombp/transchrombp/models/bias_branch.py`

先回答：

- 当前模型的 forward 到底怎么走
- full/debiased 在 forward 里是怎么分开的
- count pool 的三种模式分别怎么实现

### 第二层：再搞清楚训练语义

1. `vendor/transchrombp/transchrombp/training/train_ddp.py`
2. `vendor/transchrombp/transchrombp/configs/train/train_ablation_v2_main_profile_select.yaml`

重点看：

- `compute_losses()`
- `run_validation()`
- `resolve_best_metric_name()`

### 第三层：最后看评估和 strict compare

1. `vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py`
2. `vendor/transchrombp/scripts/select_best_epoch.py`
3. `scripts/paper_aligned_repro/select_best_epoch.py`

重点看：

- held-out test 的真实输出字段
- selector 用的 metric 路径
- official 与自研是怎么被对齐到同一套外部规则下的

### 第四层：如果你要看 foundation 线

1. `vendor/transchrombp/transchrombp/models/genos_adapter.py`
2. `vendor/transchrombp/transchrombp/models/caduceus_adapter.py`
3. `vendor/transchrombp/transchrombp/models/foundation_adapter.py`

重点看：

- 外部特征是在 backbone 前、backbone 内还是 head 附近被接进来
- 为什么这些模块是“后来挂上去的分支”，不是主干最初就天然需要的一部分
- 哪些接法已经被实验线判成 `no-go` 或 `near-null`

## 10. 如果现在要边看代码边想论文，应该怎么连起来

| 你在看什么代码 | 对应要想的论文问题 |
|---|---|
| `sequence_encoder.enabled` 与 `SequenceTransformerEncoder` | Transformer 的真实收益到底来自哪里 |
| `profile_bias_stop_gradient` 与 bias fusion | 为什么不能再写成“sg 单独修复” |
| `compute_losses()` 里的 `debiased_profile_weight` | 为什么 clean matrix 最终更支持 `deb2` 是主效应 |
| `count_pool_mode=center` | 为什么 `center pool` 会进入默认模型而不是当作支线 |
| `evaluate_checkpoint.py` 的 classification 输出 | 为什么 `AUROC/AUPRC/F1` 现在能稳定进 supplementary |
| self selector / official selector | 为什么 `L3 shared-region` 成了当前 external evidence 主线 |

## 11. 一句话总结

今天最值得掌握的，不是“模型里有 Transformer”这么简单，而是下面这句更完整的话：

> 当前的 TransChromBP 是一套把 BPNet 局部归纳偏置、ChromBPNet bias factorization、可开关 Transformer、full/debiased 双口径训练评估、以及 center-pool count readout 组合在一起的统一框架；而今天对这套框架的理解，已经建立在 clean matrix、`L3 shared-region` 和 no-bias 边界补证之后的新口径上。
