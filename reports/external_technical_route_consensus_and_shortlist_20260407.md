# 外部技术路线双报告汇总与可执行 Shortlist（2026-04-07）

## 1. 这份文件是什么

这是一份“外部技术路线双报告汇总”，不是新的实验结果文件。

它的作用有四个：

1. 把 Gemini 与 ChatGPT Pro 两份路线咨询报告压缩成一套统一判断，避免后续反复翻原文。
2. 明确两份意见的**共识项**、**分歧项**和**当前最值得执行的 shortlist**。
3. 把“可以直接拿来做决策的结论”和“仍只是外部建议、尚未被仓库证据验证的假设”分开。
4. 为下一步决定是否进入 `multiscale decoder probe`、补指标、或做 AlphaGenome 小样本对照提供单一参照。

若本报告与后续新实验事实冲突，以更新的实验事实为准；若只是与旧叙事冲突，则以本报告、当前 stop-rule 与最新证据链为准。

---

## 2. 汇总范围

本报告只整合两份已经入档的外部路线咨询：

1. [gemini_deep_think_technical_route_review_20260407.md](/home/zhengwei/project/python/chromBPNet/reports/gemini_deep_think_technical_route_review_20260407.md)
2. [chatgpt_pro_technical_route_review_20260407.md](/home/zhengwei/project/python/chromBPNet/reports/chatgpt_pro_technical_route_review_20260407.md)

这两份报告都已经按仓库语境做过一次整理，因此本报告不再重复原文，只做二次归纳与行动化。

---

## 3. 一句话总判断

把两份外部意见压成一句话：

> 当前最值得押注的下一代方向，不是继续在已判负的 foundation adapter family 上打补丁，而是先用 `corrected B` 作为稳定 bias-safe 基座，优先测试“multiscale dense decoder / coarse-to-fine profile readout”是否能突破当前 ceiling；如果未来仍要利用 foundation model，它更适合扮演 `teacher / initialization / 非侵入式辅助`，而不是继续做 summary/token bolt-on。

换句话说：

- **短期最优先**：改 profile readout，而不是再改 adapter。
- **中期更像下一代项目**：sequence-to-track post-training / distillation / 更长上下文的多尺度模型。
- **当前默认不再扩线**：summary fusion、token fusion、coarse residual、`bins16 center-aligned residual` 这一整个 measured adapter family。

---

## 4. 两份外部意见的稳定共识

### 4.1 对当前事实和 stop-rule 的共识

两份报告在以下事实上几乎完全一致：

1. `corrected B` 是当前稳定强基座。
2. `full/debiased gap` 是必须保留的结构化诊断口径。
3. 旧的 `Transformer-specific shortcut` 强叙事已经失效，不应再占主问题。
4. 你们当前任务的难点是 `1-bp dense prediction + 脆弱 count 标度`，不是简单序列分类。
5. 当前 measured foundation adapter family 默认停表，至少包括：
   - summary fusion
   - token fusion
   - coarse residual
   - `bins16 center-aligned residual`

这意味着一个关键结论已经可以写死：

> 后续如果还想继续 foundation 方向，必须提出**不属于当前 adapter family 的新范式**；“再试一个小变体”已经不再是合理默认动作。

### 4.2 对下一代第一优先级的共识

两份报告都把第一优先级压在：

> **multiscale / U-Net-like / dense decoder / skip-connection readout**

而不是：

> 再来一轮 foundation adapter 微调

它们虽然措辞不同，但本质都在强调同一件事：

- 当前最可疑的 ceiling 不在 backbone 本体，而在 **从 bottleneck representation 回到 1-bp profile 的输出形成机制**。
- 当前 profile 支路更像强 backbone 上的直读，而不是显式 coarse-to-fine reconstruction。
- 因此最值得先测的，是把 profile 支路升级成真正的 decoder/readout，而不是继续测试哪种 token 接法也许能带来净增益。

### 4.3 对 foundation model 未来角色的共识

两份报告都没有建议彻底放弃 foundation model，但都明确要求改变它的角色：

1. 不再优先做前向 token/summary 注入；
2. 更适合做 `teacher / distillation / initialization / 非侵入式条件器`；
3. 若有更激进的新接口，也应保证主干继续掌握 1-bp 空间主动权，而不是让粗粒度 token 接管 profile 轨道。

因此，“继续做 foundation”这件事现在的正确含义已经变成：

> 不再问“怎么把 token 塞进去”，而是问“怎么让 foundation 在不破坏主路径的前提下提供辅助信息”。

### 4.4 对 AlphaGenome 对比策略的共识

两份报告都不支持上来就做大而全 benchmark。

共识是：

1. AlphaGenome 更适合做**窄而准**的小样本 external coordinate。
2. 更适合做 matched slice，而不是百万级 exhaustive benchmark。
3. 对比时要区分：
   - 架构内公平主对比；
   - 输出层面的 system-level 外部对比；
   - 不公平但有解释价值的 supplementary 对比。

---

## 5. 两份外部意见的主要分歧

这些分歧不是相互否定，而是帮助我们决定“第二优先级怎么排、先补哪些指标、先做哪类 AlphaGenome slice”。

### 5.1 指标优先级不同

**Gemini 更优先：**

1. `Count Spearman`
2. peak/nonpeak `AUPRC`

它更偏向答辩、生信 reviewer、防守排序稳健性和类别不平衡问题。

**ChatGPT Pro 更优先：**

1. `profile Pearson`
2. `count RMSE / MAE(logcount)`

它更偏向把 dense track 评价补齐到“shape + count calibration/error”。

### 5.2 第二路线不同

**Gemini 更看重：**

- 分布感知 / 生成式 profile+count 建模

**ChatGPT Pro 更看重：**

- sequence-to-track post-training / distillation / AlphaGenome-lite 路线

这里的差异本质上是：

- Gemini 的第二路线更像“loss / output distribution 升级”；
- ChatGPT Pro 的第二路线更像“project-scale 范式升级”。

### 5.3 AlphaGenome 起手对比不同

**Gemini 更看重：**

- bias-only stress test
- 幻觉峰压力测试

**ChatGPT Pro 更看重：**

- matched raw-track slice
- 架构内公平主对比之后，再做外部系统坐标校准

这两者可以并存，但优先级不同：

- Gemini 更强调“机制防守的杀手锏图”；
- ChatGPT Pro 更强调“先把任务层面的 matched 对比站稳”。

---

## 6. 当前最不该做的事

基于两份外部意见和当前 stop-rule，现在最不该做的是：

1. 再开新的 `adapter family` 小变体，希望从 `summary/token/residual` 路线上再挤一点增益。
2. 把“继续做 foundation”误解成“再试一个 cross-attn 或更高带宽 residual”。
3. 一上来就做大规模 AlphaGenome benchmark。
4. 继续把论文和内部思考重心放在旧的 `shortcut / bias 连接误会` 上。
5. 在没有统一指标包的情况下，靠单一 `JSD + count Pearson` 继续猜下一代架构是否有效。

这些动作的问题不只是“可能没用”，而是：

> 它们继续沿着已经被证据明显降权的方向消耗预算和注意力。

---

## 7. 当前最值得执行的统一 Shortlist

如果要把两份外部意见压成一个当前最合理的三步 shortlist，我建议这样排：

### 7.1 动作 1：只做一个 `corrected B + multiscale decoder probe`

这是两份报告的**最强交集**，也应当是当前第一优先级。

约束应写死：

1. 保持 `corrected B` 的 bias branch 不变；
2. 保持 `full/debiased` 输出语义不变；
3. 保持 count 头独立，不把 count 逻辑并进 decoder；
4. 只改 profile 支路，从直读式 readout 升级成 2-3 级 coarse-to-fine decoder；
5. 接入来自 conv/local tower 的 skip features；
6. 显式约束 “先生成 debiased signal，再做 bias fusion”，避免 skip 把 bias 泄漏回 profile 主路径。

这一步的价值最高，因为不论结果正负，解释都很干净：

- 若显著提升：说明当前 ceiling 很可能真在 readout/multiscale。
- 若无提升：则 AlphaGenome/U-Net 借鉴的短期价值会被明显压低，后续应更谨慎转向更重路线。

### 7.2 动作 2：补一个统一的最小指标包

这里不建议二选一，而是直接吸收两份报告的共同高信息量部分。

**主表优先指标：**

1. `profile JSD`
2. `profile Pearson`
3. `count Pearson`
4. `count Spearman`
5. `count RMSE / MAE(logcount)`
6. `full/debiased gap`
7. `peak / nonpeak` 分层

**supplementary 优先指标：**

1. 多分辨率 profile correlation
2. normalized JSD
3. peak/nonpeak `AUPRC`
4. calibration plot / deviance / NLL

这样做的原因是：

- Gemini 关注的“稳健排序 / 类别不平衡”不能丢；
- ChatGPT Pro 关注的“shape 对齐 / count calibration”也不能丢；
- 二者合起来，才足以判断新架构到底提升了什么。

### 7.3 动作 3：做一个窄的 AlphaGenome matched slice，并保留 stress test 作为第二层

这一步也不需要二选一，而是分两层执行：

**第一层：matched raw-track slice**

- 先抽一小批与当前 held-out 对齐的 loci；
- 做 raw track-level compare；
- 不追求大样本，只求快速拿到外部坐标。

**第二层：bias-only stress test**

- 如果第一层值得继续，再补 Gemini 强调的偏置压力测试；
- 把它定位为 supplementary / 机制防守图，而不是最先上的主对比。

这能同时满足两份报告的诉求：

- 先站稳任务层面的 matched 外部坐标；
- 再寻找最有杀伤力的 bias-aware 防守对比。

---

## 8. 如果动作 1 失败，第二层该怎么排

如果 `multiscale decoder probe` 没有给出 clean gain，下一层优先级我建议这样排：

1. **优先考虑 ChatGPT Pro 那条更重的路线：**
   - sequence-to-track post-training / distillation / 更长上下文的多尺度模型

2. **把 Gemini 的生成式 / 分布感知路线保留为探索分支：**
   - 它信息量高，但当前离现有仓库和主线更远，且不如 decoder probe 那样有清晰的最小验证切口

3. **继续 foundation 方向时，只允许新范式：**
   - teacher
   - initialization
   - 非侵入式 conditioner
   - 或其它明确不属于当前 adapter family 的方案

换句话说，动作 1 若失败，不代表项目没路，而代表：

> “轻量 AlphaGenome/U-Net 借鉴”这一层价值下降，应该把注意力转到更上层的 supervised sequence-to-track 项目，而不是退回 adapter 老路。

---

## 9. 对论文、答辩与评审意味着什么

这份统一汇总对论文和答辩的最大价值，不是提供更多口号，而是把故事再向前推了一步：

> 我们已经建立了一个强的 bias-safe 基座，并且有纪律地关掉了无效的 foundation adapter family；因此下一代赌注不再是“再接一个大模型试试看”，而是“把模型真正升级为 multiscale sequence-to-track”。

这会显著改善面对评审和答辩时的解释力：

1. 当被问“为什么不用最新基因组大模型”时，可以回答：不是没试，而是当前 adapter family 已被 stop-rule 判负。
2. 当被问“为什么现在还不直接跟 AlphaGenome 拼大 benchmark”时，可以回答：当前先做窄而准的 matched slice，更符合任务层面比较，也更节省。
3. 当被问“下一步真正准备做什么”时，可以直接回答：先测 multiscale decoder probe，再用统一指标包判断提升来源。

---

## 10. 当前行动化结论

把这份双报告汇总压成一句行动化结论：

> 现在最值得先做的，不是继续 foundation adapter 微调，而是围绕 `corrected B` 做一个最小、干净、bias-safe 的 `multiscale dense decoder probe`，同时补齐能区分 `shape / count calibration / ranking / bias reliance` 的统一指标包；AlphaGenome 对比则先做窄的 matched slice，再决定是否补 bias-only stress test。

如果只允许马上启动一个动作，就选：

> `corrected B + multiscale decoder probe`

如果允许并行准备但不立刻大改代码，则最值得同步准备的是：

1. 指标补算脚本与主表口径；
2. AlphaGenome matched loci shortlist；
3. 后续可能需要的 bias-only stress test 序列构造方案。
