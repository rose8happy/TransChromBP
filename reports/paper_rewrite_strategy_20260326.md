# TransChromBP 论文重写策略与补实验优先级（2026-03-26）

## 1. 先给结论

当前这篇论文不应该再写成“我们发现并修复了一个明确、强烈、且 Transformer 特有的 Profile Shortcut”。  
更稳、也更能被当前证据支撑的写法是：

1. `TransChromBP` 是一个将 Transformer 与 bias factorization 结合的碱基分辨率 ATAC-seq 预测模型。
2. 论文的方法学亮点不是“证明了一个强 shortcut 机制”，而是“提出了 bias-safe 的训练/诊断框架，并用受控复核证明最终模型没有明显 bias reliance”。
3. `full/debiased` 双口径指标是本文比单纯换 backbone 更有方法学价值的部分。
4. `debiased profile supervision` 比 `stop-gradient` 更像抑制 bias reliance 的主效应；`stop-gradient` 只能写成辅助性稳定设计，不能再写成唯一关键开关。
5. 最终默认模型 `corrected B = center pool + sg=true + deb2` 在 held-out/test 上没有明显 shortcut，可以继续作为 paper-facing 主模型。

一句话概括新的 paper story：

> 我们提出了一个 bias-safe 的 Transformer ATAC profile/count 建模框架，并用 full/debiased 诊断证明：当显式 debiased supervision 存在时，Transformer 的收益来自真实表征提升而非 bias shortcut；最终 center-pool 读出头进一步改善了 count 预测。

## 2. 当前草稿的核心问题

### 2.1 机制断言过强，已经和 clean matrix 结果不匹配

当前 [transchrombp_paper_cn_v1.tex](/home/zhengwei/project/python/chromBPNet/reports/transchrombp_paper_cn_v1.tex) 里仍有以下强断言：

- “发现了此前未被报道的失效模式 Profile Shortcut”
- “纯卷积架构中不会出现”
- “这是 Transformer 特有现象”
- “stop-gradient 单独修复了该问题”

但 [profile_shortcut_revalidation_summary_20260326.md](/home/zhengwei/project/python/chromBPNet/reports/profile_shortcut_revalidation_summary_20260326.md) 的四条 test 主线不支持这些写法：

- `A = TF + sg=false + deb2`：`profile_full_debiased_jsd = 0.00147`
- `C = TF + sg=true + deb0`：`0.00972 / 0.00997`（双 seed）
- `noTF + sg=false + deb2`：`0.00487`
- `corrected B = TF + center pool + sg=true + deb2`：`0.00268`

最稳的解释是：

- 我们没有复现出“明显且稳健”的强 shortcut。
- 去掉 `deb2` 会增加轻度 bias reliance。
- 当前没有证据支持 “Transformer 特有”。

### 2.2 当前草稿混用了三类证据，层级不清

现在草稿里混着三类东西：

1. 早期探索性现象  
   例如 `v1 -> v2` 的大幅 gap、scale 轨迹、旧版 `2x2` 析因。
2. 后期受控复核  
   例如 A/C/noTF/corrected B 四条 test。
3. 最终模型选择  
   例如 `B=center pool` 相对其他 readout 的优势。

问题在于：第 1 类证据的因果不干净，第 2 类才是当前可以写进论文主结论的证据，第 3 类则回答“最终交付什么模型”。  
如果不把这三层分开，读者会误以为：

- 早期探索现象已经被严格证明；
- 最终主模型就是为“shortcut 修复”而设计；
- 所有后续改动都围绕一个已经证实的机制展开。

这在当前证据下都不够稳。

### 2.3 主结果表、机制表、最终模型表还没有分工

当前最容易写乱的地方是把下面三件事混成一个表：

- `TransChromBP` 相对 `ChromBPNet` 的整体性能提升
- `bias reliance` 的 clean revalidation
- `center pool` 为什么是最终默认 readout

这三件事必须拆成三张表，否则每张表的“回答对象”不清晰。

## 3. 推荐的论文主线

### 3.1 论文定位

建议把论文从“强机制发现论文”改成“模型 + 诊断框架论文”。

更具体地说：

- 不是“我们发现了 Transformer 特有的 shortcut”
- 而是“我们把 Transformer 安全地接入了 bias-factorized base-resolution ATAC modeling，并提出了用于检测 bias reliance 的 full/debiased 诊断口径”

### 3.2 最终应保留的 4 个主 claim

#### Claim 1: Transformer 可以在 bias-factorized ATAC 建模中带来真实收益

证据来源：

- `V2-full` 对 `V2-noTF` 的 matched ablation
- `V2-full` 相对 `ChromBPNet` 的主结果

推荐写法：

> 在 debiased supervision 存在且 full/debiased gap 保持较低的条件下，Transformer backbone 在 tutorial benchmark 上带来稳定的 profile/count 收益。

#### Claim 2: full/debiased 双口径是必要的诊断，而不是锦上添花

证据来源：

- A/C/noTF/corrected B 四条 clean test

推荐写法：

> 单看 full 指标不足以判断模型是否对 bias branch 产生依赖；full/debiased gap 提供了一个低成本、可复用的 bias reliance 诊断口径。

#### Claim 3: debiased profile supervision 是更关键的稳定因素

证据来源：

- A vs C 的 clean 对照

推荐写法：

> 在当前受控矩阵中，去除 debiased profile supervision 会明显抬高 full/debiased gap，而单独关闭 stop-gradient 并不会导致同量级退化。

注意这里应使用 “在当前受控矩阵中观察到” 或 “suggests”，不要写成无条件普适定律。

#### Claim 4: 最终默认模型 corrected B 稳定且可交付

证据来源：

- `B_s42` 和 `B_s1234` held-out/test
- corrected `B` 的 clean shortcut test

目前可直接使用的数字：

- `B_s42`: `0.3147 / 0.8503`
- `B_s1234`: `0.3145 / 0.8488`
- 两 seed 均值约为 `0.3146 ± 0.00014 / 0.8496 ± 0.0011`
- corrected `B` clean test：`profile_full_debiased_jsd = 0.00268`

推荐写法：

> 在当前 readout 设计中，center-aligned count pooling 能在不引入明显 bias reliance 的前提下进一步提升 count 预测，因此被选作最终默认配置。

### 3.3 对外比较证据必须分层表述

当前与 `ChromBPNet` 的对比不应再被统称为“严格控制变量比较”。更稳的分层是：

1. `Layer 1: official/native-pipeline compare`
   - 回答：在社区熟悉的 tutorial 口径下，我们是否能站到外部锚点旁边。
2. `Layer 2: matched-budget system compare`
   - 回答：在共享原始输入、相同硬件、相同训练预算与相同 external best-rule 下，完整系统谁更强。
   - 这一层允许 native preprocessing 与优化器 recipe 不同，因此只能支撑 system-level claim。
3. `Layer 3: shared-region compare`
   - 回答：在统一最终 `peaks/nonpeaks` 后，优势是否仍然成立，从而更接近架构/recipe 归因。

当前 tutorial 主线已经确认两侧 candidate region 集合不同，因此现有结果最多只能归入 `Layer 2`。  
如果正文继续使用这批结果，表注必须明确写成 `shared raw inputs + native preprocessing`，不能再写成 “same frozen dataset” 或 “严格到 region 集完全一致”。

## 4. 论文结构建议

## 4.1 标题

当前标题还能用，但建议副标题或摘要第一句更强调 `bias-safe`：

- `TransChromBP：面向碱基分辨率 ATAC-seq 的 Bias-Safe Transformer 模型`
- 或保留原标题，在摘要中补一句 “a bias-safe Transformer framework”

## 4.2 摘要

摘要应改成三段逻辑：

1. 任务与问题  
   碱基分辨率 ATAC profile/count 预测需要既建模远距离依赖，又控制 Tn5 bias。
2. 方法与核心发现  
   我们提出 TransChromBP，并引入 `full/debiased` 双口径诊断；受控复核显示，显式 debiased supervision 对抑制 bias reliance 更关键。
3. 结果  
   TransChromBP 相对 ChromBPNet 提升主指标；最终 `corrected B` 在 test 上没有明显 shortcut。

摘要里建议删掉：

- “此前未被报道的失效模式”
- “纯卷积不会出现”
- “stop-gradient 单独修复”
- 过度强调 scale 轨迹和跨通路效应

## 4.3 引言

引言不要再把“发现强 shortcut”作为问题提出。  
更好的问题设置是：

1. Transformer 可能带来更强的远程建模能力。
2. 但 bias factorization 让评估不能只看 full output。
3. 因此需要一套 bias-safe 的训练与诊断框架。

引言最后的贡献点建议改成：

1. 提出 `TransChromBP`，将 Transformer 与 bias factorization 结合到 base-resolution ATAC 建模。
2. 提出 `full/debiased` 双口径诊断，用于监测 bias reliance。
3. 受控复核显示 `debiased profile supervision` 是更关键的稳定因素。
4. 提出并验证 `center pool` 读出头，作为最终默认模型。
5. 给出 Genos 融合的负结果分析。

## 4.4 方法

方法节保留：

- 模型结构
- bias branch / signal branch / fusion
- `full` 与 `debiased` 两种评估口径

方法节弱化：

- “stop-gradient 为核心修复”

建议改成：

- `stop-gradient` 是可选的 bias-isolation 设计
- `debiased supervision` 是训练目标中的关键组成部分

## 4.5 实验节建议拆成 5 个子节

### A. 主结果：与 ChromBPNet 的系统级对比

回答问题：

- 这个模型整体上是不是值得做？

如果使用当前这条主线，标题和表注都应明确这是 `matched-budget system compare`，而不是 `shared-region compare`。

建议主表只放：

- ChromBPNet
- TransChromBP 主 scaffold
- 最终 `corrected B`（如果想突出最终交付模型）

如果表太拥挤，就把 `corrected B` 放在 readout 节。

### B. Bias reliance 复核

把现有 “Profile Shortcut 分析” 改名为：

- `Bias reliance 诊断与受控复核`
- 或 `Full/debiased gap 的诊断价值`

只放四条 clean test：

- A
- C
- noTF
- corrected B

主结论写成：

- 没有复现强 shortcut
- `deb2` 更关键
- 当前证据不足以支持 Transformer 特有

### C. Backbone 消融

保留 `V2-full vs V2-noTF`。  
这里回答的是：

- Transformer 是否真的带来有效增益？

要把这节和上面的 bias reliance 复核明确分开：

- 这节讲 “有没有收益”
- 上节讲 “收益是不是建立在 bias reliance 上”

### D. Readout 设计

这里只回答：

- 为什么最终模型用 `center pool`

建议把 `B` 提到主位，把 `F/G` 压成配角。  
如果版面紧张：

- 正文只保留 `baseline vs B vs F`
- `A/G` 放附录或补充材料

### E. Cross-dataset / foundation model

这两节目前都不是本文主冲突。

建议：

- cross-dataset 保留为“外部有效性”证据，但篇幅压缩
- Genos 负结果放正文后半段或附录，不要让它抢走主线

## 5. 哪些现有内容应该删、降级或移附录

### 5.1 建议直接删掉或完全改写

- 任何“纯卷积不会出现 shortcut”的句子
- 任何“Transformer 特有”的句子
- 任何“stop-gradient 单独修复问题”的句子
- 任何把早期 `v1 -> v2` 现象写成已证机制的段落

### 5.2 建议从主文移到附录

- scale 参数轨迹
- 旧版 `2x2` 析因表
- 早期探索阶段的巨大 `full/debiased` gap 叙事

这些内容不是完全不能保留，而是只能作为：

- 开发过程中的启发性观察
- 或历史背景

不能再让它们承担主文结论。

### 5.3 当前 draft 中的一个具体风险

`-BiasBranch` 消融目前还是 `TODO`，但它并不是当前 paper story 的关键路径。  
如果短期内不准备真跑，就不要在主文里把它写成“即将证明的关键问题”；可以直接删掉该小节，或者在附录/未来工作里一笔带过。

## 6. 接下来必须做的分析

这些分析不需要新训练，应该先于任何补实验完成。

### 6.1 建立一份“论文最终数字总表”

把下面几组结果统一抄到一个 `csv` 或 `md` 表里，作为 paper 唯一数字源：

- ChromBPNet 主结果
- `V2-full` / `V2-noTF`
- A/C/noTF/corrected B clean test
- `B_s42` / `B_s1234`
- GM12878 / K562 跨数据集结果

目的：

- 避免正文、表格、摘要分别引用不同版本数字
- 明确哪些是 `val`，哪些是 `test`
- 明确哪些是单 seed，哪些是多 seed

### 6.2 做一份“claim-evidence matrix”

建议逐条列：

- 要写的 claim
- 对应证据文件
- 证据等级：探索 / 受控 / 最终 test
- 允许的措辞强度：`shows` / `suggests` / `is consistent with`

这一步能系统性减少 overstating。

### 6.3 统一数据集叙事

当前草稿容易把：

- tutorial K562 benchmark
- 独立 K562
- GM12878

混写成一个连续故事。  
实际应拆成：

- `tutorial K562`：主 benchmark
- `GM12878 + independent K562`：外部泛化验证

否则读者会疑惑为什么同样叫 K562，数字却差很多。

### 6.4 决定主模型和诊断 scaffold

正文必须明确：

- `A / C / noTF` 是机制诊断 scaffold
- `corrected B` 是最终 paper-facing 默认模型

否则读者会问：为什么 main table 里一个模型，shortcut table 里又是另一个模型。

## 7. 接下来最值得补的实验

## 7.1 P0：不是必须，论文已经可以开始重写

如果接受“弱化 shortcut 叙事”，现在已经足够开写。  
也就是说，当前没有任何一个新实验是“没有它就不能改稿”的。

## 7.2 已完成：`noTF + sg=true + deb0` 第二个 seed 也已落盘

conv-only unsafe 这条 matched control 现在已有双 seed：

- s42: `0.02682`
- s1234: `0.01679`

两者都高于 TF 的 `C = 0.00972 / 0.00997`，说明 architecture-specific 叙事应彻底放弃。  
同时，第二个 seed 没有越过 `0.02` 门槛，也说明更准确的口径不是“稳定复现强 shortcut”，而是：

- conv-only unsafe 仍是 clean matrix 的最高风险格；
- 但风险幅度存在 seed 级波动；
- `full/debiased gap` 的主导因素仍更接近 supervision 设计是否安全。

## 7.3 已完成：`tf_sg1_deb0` 第二个 seed 证实 C 线轻度 gap 稳定

`tf_sg1_deb0` 第二个 seed 的 test 现已补齐：

- s42: `0.00972`
- s1234: `0.00997`

这说明 C 线的轻度 gap 不是单 seed 偏高。  
因此，“`deb2` 比 `sg` 更关键”现在可以继续保留，但仍应写成当前受控矩阵下的稳健经验结论，而不是普适定律。

## 7.4 P2：如果想让最终主模型更漂亮，补 `B` 第三个 seed

这不是科学结论上的硬需求，而是 paper presentation 的需求。

目前 `B` 有两个 held-out/test seed：

- `0.3147 / 0.8503`
- `0.3145 / 0.8488`

已经足够支持 “B 是当前默认 readout”。  
但如果你想让主表也用 `mean ± sd`，并与前文“3 seeds”口径完全统一，那么再补一个 `B_s2024` 是最自然的做法。

## 7.5 当前不推荐补的内容

- `corrected F`
- 全 clean matrix 多 seed 扩全
- 跨机器搬 checkpoint 做评估
- 仅为完成 draft 中的 `-BiasBranch TODO` 而额外开线

这些要么不在当前主路径上，要么 paper 收益明显低于成本。

## 8. 一个可执行的写作顺序

建议按下面顺序推进，而不是直接从头到尾改 `.tex`：

1. 先锁定 claim-evidence matrix
2. 再重写摘要、引言末段、贡献点、讨论、结论
3. 接着重构实验节顺序和三张核心表
4. 最后才决定是否还要补 tutorial `Layer 2` 的 selector / checkpoint follow-up、`Layer 3`，或 `B` 第三个 seed

这样做的好处是：

- 先修“说什么”
- 再修“怎么摆”
- 最后才决定“还缺什么”

## 9. 当前最稳的决策

如果现在就要做项目推进，我建议：

1. 立刻按这份策略重写论文主线，不再等待更多实验。
2. 如果还要补实验，优先考虑 tutorial `Layer 2` 的 selector / checkpoint follow-up 或 `Layer 3`；`B` 第三个 seed 次之。
3. 把 paper 用到的所有数字持续维护在唯一总表里，不再允许“报告数字”和“正文数字”分叉。

换句话说，当前项目不缺“还能不能再跑点东西”，而是更缺“把已经有的证据按正确强度写出来”。
