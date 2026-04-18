# ChatGPT Pro 动态 Loss 协调专题包使用说明

## 这套目录是什么

这是一个给外部大模型做“动态多 loss 协调 / 多目标训练设计”咨询的专题外发包。

它不要求对方先理解我们全部历史路线，也不假定对方能访问本地仓库。它的作用是把**与 loss 协调直接相关**的任务定义、关键证据、关键代码和当前训练配置打包给 `GPT Pro`，让对方结合联网搜索能力去检索近年的多目标 / 多 loss 协调论文，再判断：

> 动态 loss 平衡方法，是否值得迁入我们当前 `profile/count/full-debiased` 框架？

这套包默认不是让外部模型替我们重写整个训练器，也不是让它继续讨论旧的 `Transformer/bias` 误会，而是让它围绕当前 loss contract 和训练动力学做方法判断。

---

## 推荐上传顺序

### 第一轮：先传 00-05 + 08 + 09

1. `00_readme_upload_order_and_prompt.md`
2. `01_task_and_current_loss_contract.md`
3. `02_loss_coordination_evidence_and_failure_modes.md`
4. `03_questions_for_external_models_dynamic_loss.md`
5. `04_code_map_and_what_to_read.md`
6. `05_key_metrics_and_runs.csv`
7. `08_train_teacher_v2_main.yaml`
8. `09_train_teacher_v2_main_profile_select.yaml`

这一轮的目标，是先让对方建立：

- 我们的任务到底是什么；
- 当前 loss 怎么定义；
- 已经暴露过哪些 loss / selector 问题；
- 它应该自己去搜哪类论文，并按什么问题来回答。

### 第二轮：只有在对方明确要核对实现时，再补 06-07

9. `06_transchrombp_current_model.py`
10. `07_train_ddp_current.py`

这两个文件都比较长，但它们足以让对方直接核实：

- `full/debiased` 输出如何形成；
- loss 是如何在代码里相加的；
- selector 指标是如何决定 best checkpoint 的。

---

## 每个文件的职责

- `01`：任务定义、当前输出语义、当前 loss contract、默认训练配置。
- `02`：和 loss 协调直接相关的证据链与失败模式。
- `03`：本轮真正希望外部模型回答的动态 loss 协调问题。
- `04`：读代码路线图，告诉对方看哪些函数、配置键和值得注意的点。
- `05`：可直接抓数的简表，避免对方只读大段叙述。
- `08`：当前默认 `teacher_v2` 主训练配置，保留 `peak.loss_total` selector。
- `09`：JSD-aligned selector 版本配置，保留同一 loss，但把 `best_metric` 切到 `peak.profile_target_jsd_full_mean`。
- `06`：当前主模型源码快照。
- `07`：当前训练主循环和 loss 计算源码快照。
- `10`：可直接复制发送给 `GPT Pro` 的最终提问文本。

---

## 这次外部模型需要做什么

请默认使用 `GPT Pro` 的搜索能力，而不是只根据这组文件里的文字做静态判断。

我们希望对方：

1. 自己检索近年的多目标 / 多 loss 协调方法；
2. 用这些方法回看我们现在的 loss contract；
3. 判断哪些方法值得迁入，哪些大概率会失效；
4. 给出最小但高信息量的验证闭环。

我们不需要对方先复述通用科普，也不希望它把篇幅浪费在“多个 loss 一般怎么加权”这种入门层面。

---

## 建议直接配套使用的提示词

```text
请把这组文件当作一个“研究项目证据包 + 动态多 loss 协调咨询包”来读，而不是只做单点问答。

这次我希望你优先扮演：

- 多目标训练 / loss 协调顾问
- 训练动力学评审者
- 实验设计顾问

请默认使用你的联网搜索能力，主动检索近年与 multiple losses / multi-objective optimization / dynamic loss balancing 相关的代表性方法和论文，而不是只根据我给的文件做常识判断。

这轮默认优先级不是“最小改动、最短周期”，而是：

在当前代码事实和证据约束下，判断哪些动态 loss 协调方法值得迁入我们的框架，哪些不值得。

请先完成下面 7 件事，再给结论：

1. 用你自己的话重建这个项目当前任务、输出语义、默认 loss contract 和 selector 口径。
2. 结合近年的多目标 / 多 loss 协调研究，给出 3 类你认为最值得考虑的方法，并按优先级排序。
3. 判断我们当前是否应该动态协调这些目标：full profile、full count、debiased profile、debiased count，以及可能的 selector 指标。
4. 说明为什么当前框架里 `debiased_profile_weight=2.0` 会表现得像主效应，以及这对动态 loss 平衡设计意味着什么。
5. 解释 `peak.loss_total` 与 `peak.profile_target_jsd_full_mean` selector 不一致的问题，判断这是否只是选模问题，还是说明训练目标本身和最终评估目标错位。
6. 指出哪些动态 loss 协调方法在我们这里最可能失效，尤其要考虑：`full/debiased` 不是独立任务、bias branch 存在、count 端是 logcount MSE、count_fusion=logsumexp、以及主指标是 JSD 而不是 NLL。
7. 给出一个“高信息量、低浪费”的最小验证方案，只允许 1-3 个动作，并明确每一步应该看哪些指标才算成功或失败。

请不要把主要注意力放在旧的 “Transformer 和 bias 机制连接可能有 bug” 叙事上。那不是本轮主问题。

如果你认为还需要核对实现，请再读 `06_transchrombp_current_model.py` 和 `07_train_ddp_current.py`，重点检查：

- `full/debiased` 输出是否真的是主从关系还是并列目标
- loss 是如何汇总的
- selector 是如何决定 best checkpoint 的
```

---

## 你读这些文件时应默认知道的边界

- 这不是普通分类任务，而是 base-resolution ATAC profile/count prediction。
- 当前主模型不是单输出模型，而是显式同时返回 `full` 和 `debiased` 两套输出。
- 当前最关键的问题不是“再找一个新 backbone”，而是当前 loss contract 和最终要优化的指标是否真正对齐。
- `GPT Pro` 应该主动去搜论文，把我们给出的代码事实和最近方法连接起来，而不是等待我把论文综述抄进包里。

---

## 一句话目标

如果你只能回答一个问题，请回答：

> 在我们当前 `profile/count + full/debiased + bias-branch fusion + selector mismatch` 的训练框架下，哪一类动态多 loss 协调方法最值得尝试，最小验证闭环是什么？

---

## 懒人用法

如果你不想从本文件手工摘 prompt，直接打开：

- `10_ready_to_send_chatgpt_pro_prompt.md`

复制其中的“可直接发送的提问文本”即可。
