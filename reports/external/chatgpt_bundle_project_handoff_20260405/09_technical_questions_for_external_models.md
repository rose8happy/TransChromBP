# 技术路线咨询问题总览

## 这份文件是什么

这份文件是本轮外发的主问题清单。

它的目的不是重复总结项目历史，而是明确告诉外部大模型：

- 这轮真正想回答的是什么；
- 哪些问题现在已经不值得再花主要注意力；
- 什么样的回答对我们最有价值。

本轮默认优先级不是“最小改动、最短周期”，而是：

> 在当前证据和 stop-rule 约束下，优先判断最有潜力的下一代方向，即使工程代价较大。

---

## 1. 你应优先扮演什么角色

请优先把自己当作：

- 技术路线顾问
- 下一代架构评审者
- 指标与对比设计顾问

请不要把自己首先当作：

- 论文语句润色器
- 单纯 reviewer
- 只会给“最小改动建议”的保守工程师

论文、答辩和评审防守当然重要，但在这轮里，它们是副轨，不是主轨。

---

## 2. 你应该重点回答的 5 组问题

### A. 我们的模型现在还能怎样提升

请回答：

1. 在当前证据下，最有潜力的 3 条下一代提升路线是什么？
2. 这些路线里，哪些最可能真正提高：
   - base-resolution profile
   - count 预测
   - 而不是只提高 classification / OCR 感知
3. 如果你要排序，应如何排：
   - 最值得做
   - 值得观察
   - 当前不值得做

我们不希望只听到“继续在当前 measured family 上微调”这类低信息量建议。

### B. 如果借鉴 AlphaGenome / U-Net 风格，是否可能提升

请回答：

1. AlphaGenome 一类的多尺度 encoder-decoder / U-Net 风格，为什么理论上可能更适合 ATAC profile/count 任务？
2. 与我们当前架构相比，它最可能补强的是：
   - 多尺度特征组织
   - decoder / upsampling 方式
   - readout / head
   - 训练目标
   - 还是别的地方
3. 如果我们要借鉴，最值得先借的是哪一层，而不是整套照搬？
4. 你觉得这条路最大的潜在收益和最大风险分别是什么？

### C. 这个领域常用哪些指标，我们该补哪些

我们当前主看：

- `profile JSD`
- `count Pearson`
- `full/debiased diagnostics`

请回答：

1. 在 ATAC 预测染色质可及性这个领域，常用哪些指标作为考量？
2. 在我们当前任务里，哪些指标最值得补？
3. 哪些指标适合主文，哪些适合 supplementary？
4. 哪些指标虽然常见，但对我们当前任务信息量不高，不值得优先补？

可以讨论但不限于：

- AUROC
- AUPRC
- F1
- Spearman
- calibration
- peak/nonpeak 分层指标
- count 误差类指标
- distribution-aware loss / likelihood 类指标

### D. 预训练基因组大模型官方下游怎么做；为什么我们现在没效果；接下来该怎么做

这点对论文、答辩和评审都很重要。

请回答：

1. 像 Genos、Nucleotide Transformer、Caduceus 这类预训练基因组大模型，在官方论文或常见实践里，通常怎么做下游 ATAC/可及性任务？
2. 它们官方如果有正结果，通常依赖的是什么：
   - probe / linear head
   - frozen feature + shallow head
   - full finetune
   - decoder / adapter
   - distillation / teacher
   - 还是别的路线
3. 为什么我们当前接入看起来没有效果？
4. 更可能的问题是：
   - 任务粒度不匹配
   - readout 结构不匹配
   - profile/count 目标不匹配
   - foundation feature 分辨率不匹配
   - count head 太脆弱
   - 还是别的原因
5. 如果继续做，这条线最值得尝试的“真正不同”的接法是什么？

请特别注意：

我们不希望你只回答“预训练大模型可能没用”；我们更希望你区分：

- 哪些东西已经被当前结果否证；
- 哪些只是当前接法被否证；
- 哪些路线仍值得作为 next-generation hypothesis。

### E. 我们应该怎样和 AlphaGenome 做对比

请回答：

1. 对当前任务而言，什么样的 AlphaGenome 对比是公平主对比？
2. 什么样的对比虽然不完全公平，但仍有解释价值，适合放 supplementary？
3. 对比应至少覆盖哪些维度：
   - profile/count 主指标
   - peak/nonpeak 分层
   - zero-shot / frozen-feature / teacher-like 使用方式
   - black-box 系统对比
4. 对论文和答辩而言，最小但有说服力的 AlphaGenome 对比闭环是什么？

---

## 3. 什么样的回答最有价值

最有价值的回答应该：

1. 先明确当前已经稳定的证据和 stop-rule；
2. 再给出下一代方向判断，而不是把注意力耗在已收口的问题上；
3. 能区分：
   - 当前 measured family 被否证了什么
   - 没被否证什么
4. 对每条建议都说明：
   - 为什么可能有效
   - 为什么可能失败
   - 需要怎样的最小验证
5. 最后能压成 `1-3` 个高信息量动作。

---

## 4. 什么不值得你花主要注意力

请不要把主要注意力放在下面这些地方：

1. 继续把“Transformer 和 bias 机制的连接可能有 bug / 有问题”当作当前主矛盾。
2. 沿当前 measured residual short10 family 再做低信息量微调。
3. 只从 reviewer 语气或论文包装角度给建议，而不回答技术路线本身。

这里第 1 条尤其重要：

> 这条怀疑曾是历史上的一个中期误会，但在当前材料里只应视为已降级的背景说明，而不是当前主问题。

---

## 5. 如果你只能回答一个问题

请优先回答：

> 在当前证据和 stop-rule 下，这个项目最有潜力的下一代技术路线是什么，我们又该怎样用指标和 AlphaGenome 对比把它验证清楚？
