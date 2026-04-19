# 任务定义与当前 Loss Contract

## 1. 任务是什么

这是一个 base-resolution ATAC-seq 预测任务，不是普通分类任务。

- 输入：`2114 bp` DNA 序列，one-hot 编码
- 输出 1：中心 `1000 bp` 的 profile logits
- 输出 2：该区域总 read count 的 `log1p(total_count)` 预测

当前主指标：

- profile：`JSD`，越低越好
- count：`Pearson r`，越高越好

这很重要，因为：

- 只在 OCR / peak classification 上有帮助的特征，不一定对当前主任务有帮助；
- 训练时最小化的 loss，不必然和最终报告的 `JSD` / `Pearson r` 完全单调。

---

## 2. 当前模型返回什么

当前默认主模型不是只返回一套输出，而是显式返回两套：

- `profile_logits_full` / `logcount_full`
- `profile_logits_debiased` / `logcount_debiased`

可以把它理解为：

- `full`：signal 分支和 bias 分支融合后的完整输出
- `debiased`：只看 signal 分支的去偏输出

这意味着：

- `full` 和 `debiased` 不是日志里临时派生的诊断量，而是架构级语义；
- 当前 loss 设计并不是普通的“两头多任务”，而是包含一对结构上强相关、但语义不同的输出目标。

---

## 3. 当前默认模型的最短定义

当前默认最好模型可压缩成：

> `conv stem + local dilated tower + transformer encoder + ChromBPNet-style bias branch + full/debiased outputs + center pool + debiased_profile_weight=2.0 + profile_bias_stop_gradient=true`

与本轮 loss 协调最相关的结构语义是：

1. `full` 和 `debiased` 同时存在。
2. profile 端可以选择是否对 bias 融合路径 stop-gradient。
3. count 端默认不是简单加法，而是 `logsumexp` 融合。
4. count head 当前默认 readout 是 `center pool`。

---

## 4. 当前默认 Loss Contract

当前默认主训练的 loss 计算，在概念上是：

```text
L_total =
  profile_weight * L_profile(full)
  + count_weight * L_count(full)
  + debiased_profile_weight * L_profile(debiased)
  + debiased_count_weight * L_count(debiased)
  + optional distill terms
```

其中当前默认主线最关键的项是：

```yaml
loss:
  profile_weight: 1.0
  count_weight: 0.1
  debiased_profile_weight: 2.0
  debiased_count_weight: 0.0
```

对应的 loss 类型是：

- `L_profile(*)`：multinomial NLL
- `L_count(*)`：`MSE(logcount_pred, log1p(total_count_target))`

这组配置不是“任意初始化”，而是当前主线里反复出现并被保住的训练语义。

---

## 5. 当前默认 Selector Contract

这里有一个对本轮非常关键的事实：

- 默认 `teacher_v2_main` 训练配置里，`best_metric` 仍是 `peak.loss_total`
- 但后来又补了一个只改 selector、不改 loss 的配置，把 `best_metric` 切到 `peak.profile_target_jsd_full_mean`

也就是说，当前已经暴露出：

1. 训练优化目标是一个混合 loss；
2. 但最终更关心的 paper-facing 指标是 `peak profile JSD`；
3. 两者并不总是单调一致。

因此，外部模型不能只把问题理解成“多个 loss 的静态权重是否合理”，还必须回答：

> 训练时的目标组合、验证时的 selector 指标、最终报告的主指标，是否应该一起协调？

---

## 6. 这次外部模型最需要理解的 4 个事实

### 事实 1：`full` / `debiased` 不是普通的两个独立任务

它们共享主干，又通过 bias-branch 融合和去偏语义相互关联。

### 事实 2：`debiased_profile_weight=2.0` 不是随便调出来的

当前证据显示，它更像是保持训练语义健康的重要项，而不只是一个小超参数。

### 事实 3：count 端在结构上更脆弱

它不只受到 `count_weight` 影响，还受到：

- `count_fusion=logsumexp`
- `center pool`
- `track_total_count_target`
- 外部分支 / summary 注入

这些因素共同作用。

### 事实 4：selector 本身已经是问题的一部分

如果 `peak.loss_total` 和 `peak.profile_target_jsd_full_mean` 指向不同 checkpoint，那么“怎么平衡多个 loss”不能只停留在 forward/backward 的标量加权层面。

---

## 7. 当前默认配置摘要

和本轮最相关的默认配置如下。

模型侧：

```yaml
sequence_encoder:
  enabled: true
  d_model: 256
  n_heads: 8
  n_layers: 6

local_tower:
  n_dil_layers: 8

bias_branch:
  enabled: true
  profile_pool_factor: 32

fusion:
  count_fusion: logsumexp
  learnable_scales: true
  profile_bias_stop_gradient: true

heads:
  profile_output_len: 1000
  count_pool_mode: center
```

训练侧：

```yaml
loss:
  profile_weight: 1.0
  count_weight: 0.1
  debiased_profile_weight: 2.0
  debiased_count_weight: 0.0
```

selector 的两个版本：

```yaml
# 旧 / 默认
trainer:
  best_metric: peak.loss_total

# 只改 selector 的版本
trainer:
  best_metric: peak.profile_target_jsd_full_mean
```

---

## 8. 这份文件希望外部模型先回答什么

读完这份文件后，请先形成下面这个判断框架：

1. 我们现在到底在协调哪些目标？
2. 它们是并列目标、主从目标，还是结构耦合目标？
3. 这个问题更像：
   - 静态 loss weight 选得不好
   - 动态多目标协调缺失
   - selector 和最终目标错位
   - 还是三者同时存在

如果不能先回答这三点，后面的“该不该引入动态 loss 平衡方法”就会变成泛泛而谈。
