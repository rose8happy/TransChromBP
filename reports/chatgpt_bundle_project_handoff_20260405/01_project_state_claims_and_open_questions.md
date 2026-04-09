# 项目当前状态、有效结论与开放问题

## 1. 项目一句话

我们在做的是：

> 把 Transformer 和更一般的 foundation-model 思路接到 ChromBPNet 风格的 bias-factorized、base-resolution ATAC profile/count 预测框架里，并判断哪些改动是真增益，哪些只是看起来热闹但不成立。

本轮请把这个项目理解成：

- 一个已经有较强 baseline 和较清楚 stop-rule 的系统；
- 当前更需要回答“下一代提升路线是什么”，而不是继续纠缠旧的误会叙事。

---

## 2. 任务定义

- 输入：`2114 bp` DNA 序列，one-hot。
- 输出 1：中心 `1000 bp` 的 base-resolution ATAC profile。
- 输出 2：该区域总 read count/logcount。
- 当前主指标：
  - profile：`JSD`，越低越好
  - count：`Pearson r`，越高越好

这很重要，因为它意味着：

- 这不是普通序列分类任务。
- 只在 OCR/peak-vs-nonpeak 上有用的特征，不一定对当前主任务有用。
- 任何新路线如果只提升分类感知、却不提升 profile/count 主任务，价值都要打折。

---

## 3. 截至 2026-04-07 最稳定的结论

### 可以强判断

1. Transformer 在当前 bias-factorized、base-resolution ATAC 建模里带来真实收益。
2. `full/debiased gap` 是必要的 bias-reliance 诊断口径，单看 full 指标不够。
3. 当前 paper-facing 默认模型 `corrected B` 是稳定且相对 bias-safe 的。
4. `center pool` 是当前最稳的默认 readout。
5. 在 tutorial `L3 shared-region` 的 matched system compare 上，当前系统优于 official ChromBPNet。
6. 模型在独立 GM12878 和独立 K562 上不是完全失效，说明它不只会做 tutorial benchmark。

### 只能弱判断

1. 在当前 clean matrix 下，`debiased profile supervision` 比 `stop-gradient` 更关键。
2. Genos 当前失败，可能主要是任务粒度不匹配，而不只是工程没跑通。
3. NT v2 比 Genos 更像“有独立信号”的候选，但当前已经测过的 measured family 仍未显示 clean gain。

### 当前不该再坚持

1. “我们发现了 Transformer 特有 shortcut。”
2. “纯卷积不会出现 shortcut。”
3. “stop-gradient 单独修复了问题。”
4. “我们已经证明了强而稳定的 shortcut 机制。”
5. “继续沿当前 residual short10 measured family 自然扩线，就大概率能找到正结果。”

---

## 4. 当前项目状态

### 主线 A：当前默认最好模型已经比较稳定

当前最像主线基座的仍然是：

> bias-safe Transformer framework + full/debiased diagnostics + stable readout design + matched external compare。

也就是说，当前我们不是在找“有没有一个明显坏掉的旧 bug 需要继续解释”，而是在问：

> 在这个已经比较稳定的基座上，下一代提升路线应该朝哪里走？

### 主线 B：foundation model 当前 measured family 已默认停表

我们已经做过：

- Genos：较完整接入，但总体负结果
- Caduceus-PS：single-seed `near-null / marginal positive`
- NT v2：probe 有独立信号，但
  - `bins4/coarse residual` 判负
  - `bins16 center-aligned residual` full held-out 也仍落后于 matched no-foundation baseline

因此当前开放问题**不是**：

- “还要不要继续沿这条 measured family 微调？”

而是：

- “如果评审、答辩或后续研究必须面对 foundation-model 潮流，我们下一代真正不同的方向应该是什么？”

---

## 5. 当前最佳模型的简短定义

当前默认最好模型可简写为：

`corrected B = transformer + local tower + ChromBPNet-style bias branch + full/debiased outputs + center pool + debiased_profile_weight=2.0 + profile_bias_stop_gradient=true`

这不是“单纯把 backbone 换成 Transformer”。

---

## 6. 当前最需要外部模型回答的问题

### A. 模型还能怎样提升

1. 在当前证据和 stop-rule 下，最有潜力的 3 条下一代提升路线是什么？
2. 哪些路线最可能真正提高 base-resolution profile/count，而不只是提高分类感知？
3. 哪些方向值得付出中等甚至较大重构代价？

### B. AlphaGenome / U-Net 风格是否可能带来提升

1. AlphaGenome 一类多尺度 encoder-decoder / U-Net 式架构，是否理论上更适合当前任务？
2. 如果值得借鉴，最值得先借的是：
   - 多尺度特征组织
   - decoder / readout
   - head 设计
   - 训练目标
   - 数据/标签组织
3. 对我们当前架构而言，最小但信息量高的 AlphaGenome-inspired 改造会是什么？

### C. ATAC 领域指标应该怎样补

1. 在当前 `profile JSD + count Pearson + full/debiased` 之外，还应补哪些 ATAC 常用指标？
2. 哪些指标适合主文，哪些适合 supplementary？
3. 哪些指标虽然常见，但对我们当前任务信息量不高，不值得优先补？

### D. 预训练基因组大模型为什么当前接入没显出效果

1. Genos / Caduceus / NT 一类模型在官方论文或常见做法里，通常怎样做下游 ATAC/可及性任务？
2. 为什么它们在我们这里看起来没有形成 clean gain？
3. 更可能的问题在于：
   - 任务粒度不匹配
   - readout 结构不匹配
   - 训练语义 / loss 设计
   - 特征分辨率与 profile/count 任务错位
   - 还是别的原因
4. 如果继续做，下一步最值得尝试的“真正不同”的接法是什么？

### E. 我们应该怎样与 AlphaGenome 做对比

1. 应做哪些主对比、附加对比和 black-box 对比？
2. 什么叫对当前任务“公平且有说服力”的 AlphaGenome 对比？
3. 对论文和答辩来说，最小但足够有力的 AlphaGenome 对比闭环是什么？

---

## 7. 当前最不值得外部模型继续花注意力的地方

1. 继续把旧的 `Transformer/bias` 连接问题当作当前主矛盾。
2. 把“有独立信号”直接误判成“有可利用的互补净增益”。
3. 把当前 measured family 的 stop-rule 忽略掉，继续沿旧 family 微调。
4. 用不适合当前主任务的指标或对比，回答不了评审/答辩真正会问的问题。

---

## 8. 如果你只能记住一件事

这个项目当前最值得外部判断的，不再是：

> “旧的 shortcut/bias 误会到底成不成立？”

而是：

> 在已有证据和 stop-rule 下，下一代最有潜力的提升路线是什么，以及我们该怎样用指标和 AlphaGenome 对比把它验证清楚。
