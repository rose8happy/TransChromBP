# ChatGPT Pro 外发总档案使用说明

## 这套目录是什么

这是一个给外部大模型做“技术路线咨询”的压缩外发包，不是假定对方能访问本地仓库。

当前使用方式已更新到：`2026-04-07`。

这套包保留了 `2026-04-05` 主包的大部分证据底座，但本轮主目标已经切换成：

1. 请外部模型判断当前最有潜力的下一代技术路线；
2. 请外部模型帮我们评估 `AlphaGenome / U-Net` 风格、指标补充、预训练基因组大模型下游接法与对比策略；
3. reviewer / 论文收口问题只作为副轨约束，不再是这轮的主任务。

基础主包仍然形成于 `2026-04-05`，但当前外发时应**强制附带**
`08_post_20260406_delta.md` 和 `09_technical_questions_for_external_models.md`。

---

## 推荐上传顺序

### 第一轮：先传 00-05 + 08 + 09

1. `00_readme_upload_order_and_prompt.md`
2. `01_project_state_claims_and_open_questions.md`
3. `02_architecture_code_map_and_best_model.md`
4. `03_experiment_history_and_evidence.md`
5. `04_foundation_attempts_and_current_restart.md`
6. `05_key_metrics_and_runs.csv`
7. `08_post_20260406_delta.md`
8. `09_technical_questions_for_external_models.md`

这一轮已经足够让对方做高质量技术路线咨询。

### 第二轮：只有在对方明确要求真实代码时，再补 06-07

9. `06_transchrombp_current_model.py`
10. `07_foundation_adapter_current.py`

---

## 每个文件的职责

- `01`：项目当前状态、有效结论、当前高价值开放问题。
- `02`：任务定义、当前最佳模型、关键训练语义、代码地图。
- `03`：实验历史和证据链，重点说明哪些结论被保留、哪些被推翻。
- `04`：foundation model 支线，包括 Genos、Caduceus、NT v2 和当前 measured family 的边界。
- `05`：高信号运行和指标表，便于对方直接抓数字。
- `08`：`2026-04-06` 之后的关键增量，包括 `bins16 center-aligned residual` 判负、stop-rule 收紧，以及 6000 实际运行仓/环境说明。
- `09`：本轮真正希望外部大模型重点回答的技术路线问题总览。
- `10`：可直接复制发送给 ChatGPT Pro 的最终提问文本。
- `06`：当前主模型源码快照。
- `07`：当前 foundation adapter 源码快照。

---

## 建议直接配套使用的提示词

```text
请把这组文件当作一个“研究项目证据包 + 技术路线咨询包”来读，而不是只做单点问答。

这次我希望你优先扮演：

- 技术路线顾问
- 下一代架构/实验路线评审者
- 对比设计顾问

而不是先做 reviewer 式论文润色建议。

这轮默认优先级不是“最小改动、最短周期”，而是：

> 在当前证据约束下，优先判断最有潜力的下一代方向，即使工程代价较大。

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

## 你读这些文件时应默认知道的边界

- 这不是普通 OCR/peak classification 任务，而是 base-resolution ATAC profile/count prediction。
- 当前内部 baseline 很强，不应把它当成容易超过的弱起点。
- 本地归档仓的代码主镜像主要在 `vendor/transchrombp/transchrombp/`，但 6000 实际运行仓使用 `src/transchrombp/` 布局；`08` 已补充这条 runtime note。
- 当前 measured foundation family 已经形成 stop-rule；外部建议应优先讨论“下一代不同方向”，而不是默认沿旧 family 微调。
- 这套外发包已经尽量把“历史实验很多”压缩成“当前真正还重要的事实”。

---

## 一句话目标

如果你只能回答一个问题，请回答：

> 在当前证据和 stop-rule 下，这个项目最有潜力的下一代技术路线是什么，我们又该怎样用指标和 AlphaGenome 对比把它验证清楚？

---

## 懒人用法

如果你不想从本文件手工摘 prompt，直接打开：

- `10_ready_to_send_chatgpt_pro_prompt.md`

复制其中的“可直接发送的提问文本”即可。
