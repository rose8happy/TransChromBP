# Ready-To-Send ChatGPT Pro Prompt

## 建议上传顺序

第一轮先上传：

1. `00_readme_upload_order_and_prompt.md`
2. `01_project_state_claims_and_open_questions.md`
3. `02_architecture_code_map_and_best_model.md`
4. `03_experiment_history_and_evidence.md`
5. `04_foundation_attempts_and_current_restart.md`
6. `05_key_metrics_and_runs.csv`
7. `08_post_20260406_delta.md`
8. `09_technical_questions_for_external_models.md`

只有当你明确要求真实代码时，我再补：

9. `06_transchrombp_current_model.py`
10. `07_foundation_adapter_current.py`

---

## 可直接发送的提问文本

```text
请把这组文件当作一个“研究项目证据包 + 技术路线咨询包”来读，而不是只做单点问答。

这次我希望你优先扮演：

- 技术路线顾问
- 下一代架构/实验路线评审者
- 对比设计顾问

而不是先做 reviewer 式论文润色建议。

这轮默认优先级不是“最小改动、最短周期”，而是：

在当前证据约束下，优先判断最有潜力的下一代方向，即使工程代价较大。

请先完成下面 7 件事，再给结论：

1. 用你自己的话重建这个项目当前已经稳定成立的事实、硬边界和 stop-rule。
2. 给出 3 条你认为最有潜力的下一代提升路线，并按优先级排序。
3. 判断 AlphaGenome / U-Net 风格的多尺度 encoder-decoder / readout 设计，是否可能比我们当前架构更适合 base-resolution ATAC profile/count 任务；如果值得借鉴，最值得先借哪一层。
4. 结合 ATAC 预测染色质可及性这一领域，说明常用哪些指标作为考量；在我们当前 `profile JSD + count Pearson + full/debiased` 之外，最值得补哪些指标，哪些适合主文，哪些适合 supplementary。
5. 结合 Genos、NT、Caduceus 等预训练基因组大模型的一般官方做法，分析它们在下游 ATAC/可及性任务里通常怎么接；为什么我们当前接入看起来没有效果；接下来若继续做，最值得尝试的“真正不同”的接法是什么。
6. 说明我们应该怎样和 AlphaGenome 做对比：应做哪些类型的对比、哪些是公平主对比、哪些是不公平但仍有解释价值的附加对比。
7. 给出一个“高信息量、低浪费”的下一步方案，只允许 1-3 个动作。

请不要把主要注意力放在“Transformer 和 bias 机制的连接可能有 bug / 有问题”这个旧怀疑上。那在当前材料里只应被视为已降级的历史误会，而不是当前主问题。

你可以在最后附带简短指出：
- 这些技术路线建议对论文、答辩和评审防守各意味着什么；
但这不是本轮主任务。

如果你认为还需要看真实代码，再告诉我应该重点看 `06` 还是 `07`，以及你想验证什么。
```

---

## 使用备注

- 如果对方一上来就开始讨论旧的 `Transformer/bias` 误会，而不回答技术路线本身，可以直接追问：
  - “请先回到技术路线主问题：下一代最有潜力的架构、指标和 AlphaGenome 对比该怎么做？”
- 如果对方给的建议太保守、只偏向小修小补，可以追问：
  - “请再给一版不受当前实现约束、但仍受当前证据和 stop-rule 约束的下一代路线建议。”
