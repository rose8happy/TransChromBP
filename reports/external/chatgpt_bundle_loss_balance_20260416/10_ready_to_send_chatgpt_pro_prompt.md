# Ready-To-Send ChatGPT Pro Prompt

## 建议上传顺序

第一轮先上传：

1. `00_readme_upload_order_and_prompt.md`
2. `01_task_and_current_loss_contract.md`
3. `02_loss_coordination_evidence_and_failure_modes.md`
4. `03_questions_for_external_models_dynamic_loss.md`
5. `04_code_map_and_what_to_read.md`
6. `05_key_metrics_and_runs.csv`
7. `08_train_teacher_v2_main.yaml`
8. `09_train_teacher_v2_main_profile_select.yaml`

只有当你明确要求核对实现时，我再补：

9. `06_transchrombp_current_model.py`
10. `07_train_ddp_current.py`

---

## 可直接发送的提问文本

```text
请把这组文件当作一个“研究项目证据包 + 动态多 loss 协调咨询包”来读，而不是只做单点问答。

这次我希望你优先扮演：

- 多目标训练 / loss 协调顾问
- 训练动力学评审者
- 最小验证实验设计顾问

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

## 使用备注

- 如果对方开始泛泛讲多任务学习常识，而不回答我们当前框架，请直接追问：
  - “请回到当前代码事实：在 `full/debiased + bias-branch fusion + selector mismatch` 的框架里，哪些动态 loss 协调方法值得迁入？”
- 如果对方只给抽象方法名，不给失败模式和最小实验，可以继续追问：
  - “请把你的推荐压成 1-3 个最小验证动作，并明确每一步看哪些指标、什么结果算停表。”
- 如果对方忽略 selector 问题，可以继续追问：
  - “请显式回答：`peak.loss_total` 和 `peak.profile_target_jsd_full_mean` 不一致时，你建议如何同时设计训练目标和 best-checkpoint selector？”
