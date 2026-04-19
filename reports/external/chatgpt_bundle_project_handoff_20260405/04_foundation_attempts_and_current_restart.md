# foundation model 支线：做过什么、说明了什么、当前在试什么

## 1. 为什么会有这条支线

在核心 Transformer + bias-safe 主线基本站稳之后，我们开始问：

> 预训练基因组大模型，能否给当前 base-resolution ATAC profile/count 任务带来额外收益？

注意这里的门槛很高，因为当前 baseline 已经不弱。

对当前这轮外部技术路线咨询来说，这份文件的主要用途是：

- 说明当前 measured foundation family 已经触到哪些边界；
- 说明哪些“自然扩线”现在已不值得继续；
- 为下一代不同方向的判断提供约束。

---

## 2. 这条线目前已经暴露出的硬约束

1. 当前 baseline 很强，不能用“比一个弱 baseline 好一点”来宣称成功。
2. 当前任务需要接近位置分辨率的有效特征，序列级 summary 往往不够。
3. count head 很脆弱，很多看似无害的外部特征注入方式会先把 count 搞坏。
4. OCR/classification 上有信号，不等于会给当前 profile/count 主任务带来净增益。

---

## 3. Genos：相对完整地做过，但总体是高质量负结果

### 做过什么

- 本地 OCR sanity 基本对齐官方，说明模型本体和提取链路正常
- 做过 online 融合：`G1/G2`
- 做过 cached 融合：`P0/P2`
- 做过 probe：residual ridge、genos_only / encoded_only / concat AUC

### 最重要的结果

- `genos_only` 有一定独立信号，但不强
- `concat(encoded, genos)` 反而比 `encoded_only` 更差
- online/cached integration 都没有形成可靠净收益
- 尤其是把粗 summary 直接打进 count 相关路径时，count 很容易崩

### 最稳的解释

不是“Genos 没跑通”，而是：

- 它确实有信号
- 但当前最自然的信号形态更接近粗 summary
- 这和当前 base-resolution ATAC 主任务不够对齐

---

## 4. Caduceus-PS：不是负爆炸，但也不够强

我们后来做过一次 tutorial canonical 的 matched A/B：

- A：当前 corrected-B baseline
- B：在现有主干前接 `Caduceus-PS` token hidden states

最终 held-out 结论是：

- B 没有触发安全性灾难
- count 还有小幅提升
- 但 peak `mean_jsd` 只比 A 改善大约 `0.00060`
- 明显低于原先设定的推进门槛 `0.002`

因此这条线当前最诚实的表述是：

> near-null / marginal positive

这不足以支持直接扩到新数据集或新大规模主线。

---

## 5. NT v2：比 Genos 更像候选，但仍未证明互补性

### 目前技术上已经知道的事实

- NT v2 在本地可以稳定加载和前向
- 对真实 `2114 bp` window，当前 token 长度稳定是 `355`，不是 base-wise `2114`

这意味着：

- 它更自然像 token/summary 特征源
- 不像可以无痛一比一对齐到 base-resolution 主干里

### probe 最关键结果

- 最好的 `nt_only_auc = 0.742424`
  - 高于 Genos 参考，说明它确实有独立 accessibility 信号
- 但最好的 `concat(encoded, nt)` 仍比 `encoded_only` 差 `0.010324`
- 更关键的是，对 baseline residual 的回归全部失败

因此当前最稳结论是：

> NT v2 有信号，但还不够互补。

它比 Genos 更值得研究，但仍不足以直接进入旧式 tutorial A/B integration。

---

## 6. 当前 restart v3 在试什么

我们没有继续沿用旧的 concat / coarse summary recipe，而是新建了 restart v3：

- 新增 `FoundationCrossAttentionAdapter`
- 新增 `FoundationResidualHead`
- 统一 foundation-model cached feature / teacher target / distill plumbing

当前正式在跑的第一条线是：

- `NT v2 cached residual short10`
- run 名：`ntv2_residual_short10_s42_20260405_dualcache`

这条线真正想回答的问题是：

> 如果不再让 foundation feature 粗暴扰动主 backbone，而只让它学习 debiased residual，它是否仍然拿不到正收益？

这是比“再做一次 concat”信息量更高的问题。

---

## 7. 这条支线现在最需要外部判断什么

1. foundation model 这条线现在是否仍值得继续？
2. 如果值得，最该继续的是：
   - residual/distill 这类保守接法
   - cross-attention 这类 token-level 接法
   - 直接让 foundation model 做主干
   - 还是别的路线
3. Genos / Caduceus / NT v2 目前各自真正否证了什么，又没有否证什么？
4. 在只允许 1-3 个高信息量动作的前提下，下一步最值钱的验证是什么？

---

## 8. 当前我方自己的倾向

### 较倾向接受的

- Genos 旧主线应该降级
- Caduceus 当前 recipe 不值得直接扩线
- NT v2 比 Genos 更像值得继续观察的候选
- 新验证必须优先考虑“互补性”和“是否破坏 count head”，而不是只看独立信号

### 仍然不确定的

- 是否应该继续“外部分支接入”范式本身
- 是否该转向更像 residual/distillation/teacher 的设计
- 是否有比 Genos/Caduceus/NT 更合适、且更贴近当前任务的 foundation 路线

---

## 9. 一句话总结

foundation model 这条线目前最像的状态不是“已经失败”，而是：

> 旧接法大多已经给出负面或近零结果，只有 NT v2 还保留一点值得用更聪明接法继续验证的空间；但继续之前，必须先证明它能补强强 baseline，而不是只证明它自己有点信号。
