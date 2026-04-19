# 代码地图：外部模型应重点看什么

## 1. 先读哪些文件

如果你只想最快建立本轮问题图谱，按下面顺序读：

1. `01_task_and_current_loss_contract.md`
2. `02_loss_coordination_evidence_and_failure_modes.md`
3. `08_train_teacher_v2_main.yaml`
4. `09_train_teacher_v2_main_profile_select.yaml`
5. `07_train_ddp_current.py`
6. `06_transchrombp_current_model.py`

其中：

- `08/09` 让你先看配置事实；
- `07` 让你看 loss 和 selector 真正在代码里怎么实现；
- `06` 让你确认 `full/debiased` 的结构语义不是文档口头描述。

---

## 2. `08` 和 `09` 你应重点看什么

### `08_train_teacher_v2_main.yaml`

重点确认：

- `loss.profile_weight=1.0`
- `loss.count_weight=0.1`
- `loss.debiased_profile_weight=2.0`
- `loss.debiased_count_weight=0.0`
- `trainer.best_metric=peak.loss_total`

这代表当前“默认主训练”的事实口径。

### `09_train_teacher_v2_main_profile_select.yaml`

重点确认：

- loss 配置和 `08` 基本一致
- 但 `trainer.best_metric` 改成了 `peak.profile_target_jsd_full_mean`

这代表：

> 我们已经在不改 loss 的前提下，显式承认 selector 和最终关注指标之间存在错位。

---

## 3. `07_train_ddp_current.py` 里最值得搜的函数 / 关键词

请优先搜索这些名字：

### `compute_losses`

这是本轮最核心的函数。

你应该确认：

1. 当前总 loss 是如何由多项相加得到的；
2. `debiased_profile_weight` / `debiased_count_weight` 是否真在代码里参与总 loss；
3. distill 类 loss 是否只是可选项而非当前主线默认项；
4. count target 是不是 `log1p(total_count)`。

### `finalize_selection_metrics`

这是 selector 相关的关键辅助函数。

你应该确认：

1. `profile_target_jsd_full_mean` 等 selector 指标是否真在 validation metrics 中被显式计算；
2. 它们是否只是日志显示项，还是可以被 `best_metric` 真正引用。

### `resolve_best_metric_name`

这是 selector 路由的直接入口。

你应该确认：

1. 如果不显式配置，默认会落到什么指标；
2. 显式写 `peak.profile_target_jsd_full_mean` 时，训练器是否真会按它选 best。

### 其他建议检索的关键词

- `best_metric`
- `loss_total`
- `profile_target_jsd_full_mean`
- `debiased_profile_weight`
- `debiased_count_weight`
- `distill_profile_weight`
- `distill_count_weight`
- `distill_rank_weight`

---

## 4. `06_transchrombp_current_model.py` 里最值得搜的函数 / 关键词

请优先搜索这些名字：

### `forward`

重点确认：

1. `profile_logits_full` / `logcount_full` 如何形成；
2. `profile_logits_debiased` / `logcount_debiased` 是否就是 signal-only 输出；
3. `profile_bias_stop_gradient` 打开时，profile 融合路径是否会对 bias 分支做 `detach()`；
4. count 端是否默认使用 `logsumexp` 融合。

### `profile_bias_stop_gradient`

这是当前 `full/debiased` 语义和 loss 讨论的关键结构开关。

你应该确认：

- 它确实只影响 profile 融合梯度路径，而不是把 bias 分支整体从前向里删掉。

### `count_fusion`

这是理解 count loss 动力学的关键。

你应该确认：

- count 端不是普通双头回归，而是先有 signal / bias，再按融合规则合成 `logcount_full`。

---

## 5. 你读代码时应特别注意的 4 个问题

### 问题 1：`full` 和 `debiased` 在代码里到底是什么关系

如果它们不是标准多任务 learning 里两套独立 head，那么许多常规动态 loss 平衡方法可能并不直接适用。

### 问题 2：loss 的“主从关系”是否已经暗含在当前代码里

例如：

- `debiased_profile_weight=2.0`
- `debiased_count_weight=0.0`

这说明当前系统已经在用人工方式表达“哪些目标更重要”。

### 问题 3：selector 是否已成为训练语义的一部分

如果训练过程中 `best checkpoint` 由另一个指标决定，那么协调方案不能只看 backward loss。

### 问题 4：结构耦合是否会让某些方法产生伪改进

例如某类方法可能只是更快优化了容易下降的 bias-related 路径，而不是改善了真正想要的 debiased signal quality。

---

## 6. 读完代码后，你应能回答什么

理想情况下，读完 `06-09` 之后，你应该能回答：

1. 我们当前真实在优化哪些 loss 项；
2. 它们和最终报告指标之间有哪些错位；
3. `full/debiased` 是否适合被当成普通独立任务；
4. 动态多 loss 协调方法应优先作用于哪几项，而不应盲目包住所有项。
