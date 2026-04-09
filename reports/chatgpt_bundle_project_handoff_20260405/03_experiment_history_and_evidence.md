# 实验历史与证据链

## 1. 先给总判断

如果把整个项目压成一句话，当前最稳的表述是：

> 这不是一篇“证明了 Transformer 特有 shortcut”的项目，而是一个把 Transformer 安全接入 ChromBPNet-style bias factorization，并用 full/debiased 双口径做诊断和收口的项目。

这个结论不是一开始就有，而是被实验历史逼出来的。

对当前这轮外部技术路线咨询来说，这份文件里关于 `shortcut / bias-safe` 的内容应主要被读作：

- 历史证据链
- 当前主线与 stop-rule 的边界来源

而不是当前外部模型最该花主要注意力的开放问题。

---

## 2. 早期阶段：先把训练语义和工程地基打稳

在真正比较架构之前，我们先修过一些会直接污染结论的训练语义问题，例如：

- DataLoader worker 里的句柄管理
- nonpeak 每 epoch 重采样
- zero-count profile 监督逻辑
- early stopping / selector 语义

这一步本身不是论文故事，但它决定后续哪些实验值得信。

---

## 3. V1 -> V2：强现象出现，也埋下误判

早期我们同时改了多件事：

- `profile_bias_stop_gradient=true`
- `debiased_profile_weight=2.0`
- v2 的若干 bias/profile 语义

随后看到 debiased JSD 大幅改善，于是早期自然形成了一个很强的叙事：

> “我们发现并修复了 Transformer 特有的 Profile Shortcut。”

后来回看，这一步最大的问题是：

- 现象很强
- 但变量没有拆开
- 所以归因并不干净

这条早期强叙事后来被 clean matrix 明确降级。

---

## 4. backbone ablation：先证明 Transformer 本身有真实收益

在 shortcut 叙事之外，一个更基础的问题必须先回答：

> Transformer 到底有没有带来真实收益？

matched ablation 给出的结论是明确的：

- `V2-full` 三 seed：peak `JSD=0.31419±0.00012`，`count_r=0.83676±0.00531`
- `V2-noTF` 三 seed：peak `JSD=0.32393±0.00021`，`count_r=0.82503±0.00489`

因此：

- Transformer 有真实收益，这条成立
- 后面所有关于 bias-safe 和 diagnostics 的故事，才有继续讨论的价值

---

## 5. clean matrix：推翻“Transformer 特有 shortcut”强叙事

这是整个项目最重要的复核阶段。

我们把 backbone 和 bias-safe 语义拆开看，结果大意如下：

- `A = TF + sg=false + deb2` 很干净
- `C = TF + sg=true + deb0` 出现轻度 gap，并且双 seed 稳定
- `noTF + sg=false + deb2` 也较干净
- `noTF + sg=true + deb0` 反而是全矩阵最高风险格
- 最终默认模型 `corrected B` gap 很小

最重要的推论不是“某个神奇机制被证明了”，而是：

1. 单看 full 指标不够，必须看 `full/debiased gap`
2. 当前矩阵下，`deb2` 比 `stop-gradient` 更像主效应
3. “Transformer 特有 shortcut”不成立
4. “conv-only 不会 shortcut”也不成立

这一步直接重写了整个项目的 paper story。

---

## 6. readout 设计：center pool 留下，attention pool 退场

我们后来又做了 readout 对照，核心候选可以粗看成：

- `B = center pool`
- `F = attention pool`
- `G = profile refine`

关键现象：

- `B` 在两个 seed 上都稳定
- `F` 在单 seed 看着还行，但 `s1234` 的 profile 明显失守
- `G` 直接失败

因此当前结论很干净：

> `center pool` 是默认 readout；`attention pool` 不再是并列主候选。

---

## 7. 外部系统对比：L3 shared-region 让主证据更完整

我们后续把对比收紧到 tutorial `L3 shared-region` 的 matched system compare。

关键结果：

- official controlled L3：peak `mean_jsd=0.33853`，`count_r=0.69958`
- ours corrected-B controlled L3：peak `mean_jsd=0.31319`，`count_r=0.84016`

这条证据的重要性在于：

- 它不是 architecture-only 归因
- 但它是当前最强的 external/system-level compare

因此如果要写论文主证据，这条比早期单点 dramatic 现象更稳。

---

## 8. 独立数据泛化：不是只会做 tutorial

在独立数据上，我们至少能说模型没有彻底失效：

- GM12878：peak `JSD=0.42265`，`count_r=0.80396`
- K562：peak `JSD=0.61235`，`count_r=0.85857`

这些结果不该被包装成“跨数据集大获全胜”，但足以支持：

- 模型不是只会 tutorial benchmark
- 真实数据线仍有一定稳定性

---

## 9. 当前真正留下来的主结论

### 被保留的

1. Transformer backbone 有真实收益
2. `full/debiased gap` 有诊断价值
3. `corrected B` 是当前稳的默认模型
4. `center pool` 是稳的 readout
5. matched external compare 成立

### 被推翻或降级的

1. Transformer 特有 shortcut
2. conv-only 不会 shortcut
3. stop-gradient 单独修复
4. 强 shortcut 机制已经被证明

---

## 10. 这份历史对当前决策真正意味着什么

如果你要给我们下一步建议，最重要的是基于下面这个现实：

- 论文主线已经不再需要用夸张机制叙事来支撑
- foundation model 支线也不能再靠“有一点独立信号”就继续烧算力

因此任何下一步都应该回答：

> 它是在增强已经成立的主线，还是只是在重复制造新的历史噪声？
