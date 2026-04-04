# TransChromBP 当前模型学习总览

## 1. 先给结论

如果你现在只想知道“我们到底在用什么模型、怎么讲这个项目”，可以先记下面这句话：

> 当前主线已经不是“ChromBPNet 教学仓 + 一点 Transformer 改动”，而是一个以 `teacher v2 + center pool` 为默认 paper-facing 配置的 `TransChromBP` 框架：它保留了 BPNet 风格局部卷积归纳偏置、ChromBPNet-style bias factorization、可开关 Transformer backbone，以及 full/debiased 双口径训练与评估。

更重要的是，今天的项目叙事已经不再围绕“Transformer 特有 shortcut”组织，而是围绕下面三件事组织：

1. Transformer 确实带来真实收益。
2. `full/debiased gap` 是必要诊断。
3. 当前默认模型是 `corrected B = center pool + sg=true + deb2`，而 `L3 shared-region` 保住了外部比较证据。

如果你想知道“这个结论是怎么一步步变成今天这样的”，请直接接着看：

- [reports/transchrombp_internal_design_and_experiment_history_20260331.md](../../reports/transchrombp_internal_design_and_experiment_history_20260331.md)

## 2. 当前默认模型到底是什么

### 2.1 一句话结构图

```text
one-hot DNA
  -> ConvStem
  -> LocalDilatedTower
  -> Transformer encoder
  -> signal profile/count heads
  -> ChromBPNet-style bias branch
  -> full / debiased 双输出
```

当前默认 paper-facing 配置的关键词是：

| 维度 | 当前默认 |
|---|---|
| backbone | `ConvStem + LocalDilatedTower + Transformer` |
| bias 设计 | 显式 `bias_branch`，输出 `profile_bias` / `count_bias` |
| 输出口径 | 同时保留 `full` 与 `debiased` |
| 训练语义 | `debiased_profile_weight=2.0` |
| 默认 readout | `count_pool_mode=center` |

核心模型配置入口：

- [vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool.yaml](../../vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool.yaml)

### 2.2 当前主线里最不该忽略的 5 个开关

| 配置项 | 当前默认值 | 这在项目里意味着什么 |
|---|---:|---|
| `sequence_encoder.enabled` | `true` | 当前是 Transformer 版本；关掉就是 noTF 对照 |
| `fusion.profile_bias_stop_gradient` | `true` | profile 融合时对 bias 分支 stop-gradient |
| `loss.debiased_profile_weight` | `2.0` | 当前最关键的 bias-safe 训练语义 |
| `heads.count_pool_mode` | `center` | 当前最终默认 readout |
| `trainer.best_metric` | `peak.profile_target_jsd_full_mean` 或等价 paper-facing peak JSD 口径 | 选模不再只盯总 loss |

## 3. 这个模型是怎么一步步长成今天这样的

这部分只给短版时间线。要看完整版本，请去长报告。

### 3.1 最初目标：把 Transformer 安全接进 bias factorization

起点不是“做个更复杂的模型”，而是：

- 想保留 ChromBPNet 的 bias factorization
- 同时让 backbone 获得更强的长距离建模能力

因此当前模型从一开始就不是纯 Transformer，而是保留卷积与 local tower。

### 3.2 早期 V1 -> V2 现象很强，但归因不干净

项目早期同时改了：

- `stop-gradient`
- `debiased_profile_weight=2.0`
- v2 的若干 bias/profile 语义

于是曾经形成“Transformer 特有 shortcut 已被发现并修复”的早期叙事。  
后来 clean matrix 证明这条强叙事不稳。

### 3.3 backbone 主线先保住了“Transformer 确实有用”

matched ablation 的 3 seed 结果现在已经比较稳：

- `V2-full`: peak profile `0.31419±0.00012`，peak count `0.83676±0.00531`
- `V2-noTF`: peak profile `0.32393±0.00021`，peak count `0.82503±0.00489`

这一步保住了“Transformer 有真实收益”，也让 noTF 变成后续 clean matrix 的关键对照组。

### 3.4 readout 线最后收敛到 `center pool`

`B/F/G` 对照和后续 matched seed=1234 的 held-out 结果说明：

- `B=center pool` 跨 seed 稳定
- `F=attention pool` 在第二个 seed 上 profile 失守

因此今天已经不该把 `B/F` 当并列默认候选。

### 3.5 clean matrix 改写了整个故事

clean matrix 的真正作用不是“再补一组实验”，而是强迫项目从：

- “shortcut 机制发现论文”

改写为：

- “bias-safe Transformer framework + dual-metric diagnostics”

当前最稳的 clean-matrix 结论是：

- `deb2` 比单独 `sg` 更关键
- `full/debiased gap` 是必要诊断
- 当前没有证据支持 “Transformer 特有 shortcut”

### 3.6 strict compare 直到 `L3 shared-region` 才真正收口

`L2` 方向上已经支持 `TransChromBP` 更强，但还有 candidate region 混杂。  
`L3 shared-region` 统一了 peaks / nonpeaks / bigWig 后，external evidence 才真正站稳。

held-out test 上：

- Official L3: `0.33853 / 0.33989 / 0.69958`
- TransChromBP L3: `0.31319 / 0.31547 / 0.84016`

这就是为什么当前主文外部比较不该再只盯旧 tutorial baseline。

### 3.7 no-bias 没有翻案，只是边界补证

单 seed no-bias 最终结果是：

- peak `mean/median JSD=0.31496/0.31715`
- peak `count_r=0.84978`

它和当前 corrected-B / center-pool 基本同档。  
因此它只能写成 supplementary 的 component-necessity boundary result，而不是新的主证据。

## 4. 当前最应该记住的 6 个结论

### 4.1 Transformer 的收益是真实的

现在可以稳地写：

> Transformer 在当前 factorized ATAC 建模里带来真实收益，而不是只靠 bias leakage 拿到表面改善。

### 4.2 `full/debiased gap` 是方法的一部分，不是附加分析

今天如果有人只给 full 指标，不给 debiased 或 gap，就不能算把这个模型讲完整。

### 4.3 当前没有证据支持“Transformer 特有 shortcut”

这条旧叙事已经失效。  
现在更稳的说法是：**bias-safe 的 supervision 设计，比 backbone 是否是 Transformer 更关键。**

### 4.4 `center pool` 是当前默认 readout

这不是暂时偏好，而是目前最稳的默认。

### 4.5 `L3 shared-region` 是当前最强 external evidence

如果要看对外比较，现在应该优先看 `L3`，而不是只看旧 tutorial baseline。

### 4.6 no-bias 只是 supplementary boundary result

它说明“单 seed 下没有显示清晰额外收益差异”，但不改 clean matrix 和 `L3` 这两条主证据的角色。

## 5. 当前最该看的代码与配置

### 模型代码

1. `vendor/transchrombp/transchrombp/models/transchrombp.py`
2. `vendor/transchrombp/transchrombp/models/transformer_encoder.py`
3. `vendor/transchrombp/transchrombp/models/bias_branch.py`

### 训练与评估

1. `vendor/transchrombp/transchrombp/training/train_ddp.py`
2. `vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py`
3. `vendor/transchrombp/scripts/select_best_epoch.py`
4. `scripts/paper_aligned_repro/select_best_epoch.py`

### 当前最有代表性的配置

1. `vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool.yaml`
2. `vendor/transchrombp/transchrombp/configs/model/ablations/ablation_no_bias.yaml`
3. `vendor/transchrombp/transchrombp/configs/train/train_ablation_v2_main_profile_select.yaml`
4. `vendor/transchrombp/transchrombp/configs/train/train_tutorial_corrected_b_strict_compare_L3_6000.yaml`

## 6. 推荐阅读顺序

### 路线 A：只想尽快理解今天的主线

1. 本文件
2. [TransChromBP当前模型详解.md](TransChromBP当前模型详解.md)
3. [reports/paper_claim_evidence_matrix_20260326.md](../../reports/paper_claim_evidence_matrix_20260326.md)

### 路线 B：想把代码和论文一起重新串起来

1. [reports/transchrombp_internal_design_and_experiment_history_20260331.md](../../reports/transchrombp_internal_design_and_experiment_history_20260331.md)
2. [TransChromBP当前模型详解.md](TransChromBP当前模型详解.md)
3. [reports/paper_writeup_flow_20260330.md](../../reports/paper_writeup_flow_20260330.md)
4. [reports/transchrombp_paper_cn_v1.tex](../../reports/transchrombp_paper_cn_v1.tex)

### 路线 C：想从基础一路学到当前模型

1. [小白学习指南.md](小白学习指南.md)
2. [BPNet模型小白教程.md](BPNet模型小白教程.md)
3. [BPNet与ChromBPNet模型对比.md](BPNet与ChromBPNet模型对比.md)
4. 本文件
5. [TransChromBP当前模型详解.md](TransChromBP当前模型详解.md)

## 7. 一张记忆卡片

如果你需要用最短的话向别人解释当前状态，可以直接说：

> 现在的 TransChromBP 不是“一个加了 Transformer 的 ChromBPNet”，而是一套经过多轮 backbone、readout、clean-matrix、strict-compare 与 supplementary 补证之后收口出来的 bias-safe 框架；当前默认模型是 `teacher v2 + center pool`，当前最强 external evidence 是 `L3 shared-region`，而 no-bias 只保留为边界结果。
