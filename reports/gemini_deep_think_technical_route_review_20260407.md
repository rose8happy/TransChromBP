# Gemini Deep Think 外部技术路线评审回流（2026-04-07）

## 1. 这份文件是什么

这是一份“外部技术路线建议回流报告”，不是新的实验结果文件。

它的作用有四个：

1. 固定本轮 Gemini deep think 对当前项目的技术路线判断，避免后续只剩聊天摘要。
2. 把用户回传的长文本整理成当前仓库可复用的结构化结论。
3. 区分“与现有证据已经对齐的判断”和“仍需验证的新假设”，避免把外部建议误写成已证实事实。
4. 为后续等待中的第二份外部分析（ChatGPT Pro）预留对照基线，便于后续做交集/分歧整理。

若本报告与后续新实验事实冲突，以更新的实验事实为准；若只是与旧叙事冲突，则以本报告、当前 stop-rule 和最新证据链为准。

---

## 2. 来源与整理说明

- 来源：用户于 `2026-04-07` 回传的 Gemini deep think 长文本分析。
- 上下文：Gemini 已被要求以“技术路线顾问和下一代架构评审者”身份，基于当前 `00/01/08/09` 外发包和最新证据回答“下一代最值得尝试什么”。
- 整理原则：
  - 保留原分析的主判断与排序，不逐句复写。
  - 将其中重复出现的 AlphaGenome / U-Net / stress test 段落做去重合并。
  - 将“建议”“推测”“已成立事实”明确拆开，避免后续论文或 `TRACKING.md` 误用。

---

## 3. 归一化后的核心判断

### 3.1 当前稳定事实

Gemini 的前提判断与我们当前仓库证据高度一致：

1. `corrected B` 代表的 `Conv + Local Tower + Transformer + 显式 Bias 分支 + Center Pool` 已经是一个竞争力很强、bias-safe 且具有真实泛化能力的稳定基座。
2. `full/debiased gap` 不是附带分析，而是必须保留的安全诊断口径。
3. 当前任务的难点不是“有没有长程语义”，而是“如何在 1-bp dense prediction 和脆弱 count 标度下不破坏高频空间对齐与数值稳定性”。
4. `NT v2 bins16 center-aligned residual` 判负后，当前已测 `Genos / Caduceus / NT v2 residual short10` family 已足以形成 measured-family stop-rule；继续在这一 family 内做自然扩线不再合理。

### 3.2 Gemini 给出的总判断

Gemini 的总路线判断可以压缩成一句话：

> 当前最值得做的不是继续在已判负的 foundation feature 前向融合 family 里打补丁，而是寻找“能同时保住单碱基高频定位和全局语义”的下一代架构范式。

它实际上把项目从“继续找大模型接法”转成了“寻找更适合 1-bp ATAC dense prediction 的下一代结构”。

---

## 4. Gemini 推荐的 3 条下一代路线

### 4.1 优先级一：多尺度 U-Net / Encoder-Decoder 重构

这是 Gemini 给出的**最高优先级路线**。

其核心逻辑是：

- 你们当前 `flat transformer` 主干虽然有长程感受野，但深层混合容易把高频空间位置“涂抹”掉。
- ATAC profile 预测更像 1D 语义分割，而不是普通序列分类。
- 多尺度 encoder-decoder 可以把“低分辨率全局调控语义”和“高分辨率切割形状细节”解耦。
- 跳跃连接（skip-connection）允许浅层高频特征在最终 readout 前被重新接回，从而减轻 profile 细节丢失。

Gemini 特别强调：**最值得先借的不是完整 AlphaGenome，而是带跳跃连接的 decoder 思想。**

它给出的最小探针建议是：

- 从当前模型里截取 `local_tower` 的高分辨率输出；
- 将该特征越过 transformer，和 transformer 最终输出在通道维做拼接；
- 经过一个简单的 `1x1` 映射后再送入 `profile_signal_head`；
- 先用最小代价验证“补回高频细节”能否击穿当前 JSD 天花板。

### 4.2 优先级二：分布感知 / 形状感知的生成式 profile & count 建模

这是 Gemini 给出的**高潜力但中高代价路线**。

其核心判断是：

- 当前 count 的脆弱性部分来自用简单回归目标去承受极重长尾分布；
- 当前 multinomial NLL 对 profile 的形态先验过弱，不足以显式建模 footprint、峰形和连续波峰结构；
- 若把 count 改成更适合稀疏长尾的分布参数预测（如负二项/ZINB 一类思路），或把 profile 改成显式生成式建模，可能更自然地学习真实信号分布。

这条路线目前**没有仓库内直接证据支撑其必然有效**，但它确实代表一种和现有 readout/loss 族完全不同的思路。

### 4.3 优先级三：把 Foundation Model 变成“非侵入式条件器”或“只留在 loss 里的教师”

Gemini 并没有建议继续沿当前 token residual / concat / coarse summary residual 线扩展，而是建议把 foundation model 的角色降级为：

1. **非侵入式条件器**：例如压成低维环境向量，只通过 FiLM/调制去影响主干，而不污染 1-bp 主路径。
2. **隐式教师**：只用于蒸馏或辅助 loss，不进入主干前向路径。
3. **只读外部知识库**：更激进的版本是 base-to-token cross-attention，让高分辨率主干按需去读冻结 token，而不是让 token 残差强行回写到 1-bp 轨道。

Gemini 在这里的底层逻辑和当前 stop-rule 是一致的：

> 既然 foundation token 一旦进入主路径就容易破坏空间对齐和 count 标度，那么下一步就不该再测试“如何更暴力地塞进去”，而应测试“如何在不接管主路径的前提下提供辅助信息”。

---

## 5. AlphaGenome / U-Net 借鉴建议

### 5.1 值不值得借鉴

Gemini 的结论是：**高度值得借鉴。**

但它主张借的是**物理直觉和结构思想**，不是默认把 AlphaGenome 当成直接可复现白盒基线。

### 5.2 最值得先借什么

Gemini 明确把优先级放在：

1. 多尺度特征组织方式；
2. 带 skip-connection 的 decoder；
3. 用浅层高频特征纠正深层语义特征的 readout 路径。

它不建议先上最重的全量大改，而建议先做一个最小 skip-connection 探针，看当前瓶颈是否真来自“高频细节在 transformer 后被抹平”。

### 5.3 与当前仓库现状的关系

这部分**不是已成立事实**，而是目前最值得验证的新架构假设之一。

当前仓库里已经成立的是：

- transformer 带来真实收益；
- center pool 更稳；
- measured foundation family 判负。

但“当前 JSD 天花板主要由缺少 skip-connection decoder 导致”仍是待验证命题。

---

## 6. 指标建议：该补什么，不该补什么

Gemini 对指标的建议很明确。

### 6.1 最值得补、且适合主文防守的

1. **Count Spearman**
   - 原因：ATAC count 长尾严重，Pearson 容易被离群值绑架；
   - 价值：更能衡量排序能力，适合应对生信 reviewer 对稳健性的追问。

2. **Peak vs. Non-peak AUPRC**
   - 原因：全基因组正负样本极不平衡；
   - 价值：能直接回答“回归模型到底能不能找准开放区”。

### 6.2 更适合 supplementary 的

1. **分层 JSD（例如按 peak/non-peak 或 count 强度分层）**
2. **log-count 上的 MSE / MAE**

Gemini 的逻辑是：这些指标有助于证明模型没有只靠刷平坦背景区来拉低均值，也有助于展示 count calibration。

### 6.3 明确不推荐重点依赖的

1. **AUROC**
   - 原因：在高度不平衡数据上容易虚高，不够有区分度。

2. **严格依赖外部 peak caller 阈值的 overlap 指标**
   - 原因：容易把模型评价绑到外部后处理阈值上。

这部分建议与我们已有 supplementary 方向是兼容的，但是否进入主文仍需结合后续实验和第二份外部分析再做收口。

---

## 7. 对预训练基因组大模型的诊断：为什么官方看起来能做下游，我们这里却没起效

Gemini 的解释链比较完整，可以拆成三层：

### 7.1 官方常做的下游任务并不等价于当前任务

Gemini 认为 Genos、NT 一类模型在官方论文中更常见的是：

- 序列级分类；
- 变异效应预测；
- 粗粒度表征提取；

而不是直接做 `1000bp` 范围内的 1-bp dense profile/count 生成。

### 7.2 当前失败更像“空间分辨率错位 + count 标度脆弱”

Gemini 给出的主因是：

1. foundation token 的天然分辨率过粗；
2. 将粗粒度 token 强行映射回 1-bp 高频网格，会破坏主干已有的精细定位；
3. count head 的对数尺度路径对高方差外来特征非常敏感，因此很容易失稳。

这与当前 measured-family stop-rule 在经验上是高度一致的。

### 7.3 如果未来仍想继续做，这次应当换什么范式

Gemini 更看好：

- 不进入主前向的蒸馏；
- 非侵入式条件调制；
- 或者 base-to-token cross-attention 这类“主干保留空间主动权，只按需读取外部 token”的设计。

需要强调的是：这部分是**未来新 hypothesis 候选**，不是当前已被证明有效的做法。

---

## 8. AlphaGenome 对比建议

Gemini 对 AlphaGenome 的态度是“两层对比，分别服务不同问题”。

### 8.1 主对比思路：黑盒 zero-shot vs. task-specific 专家模型

在无法白盒重训或官方代码/权重不可完全对齐的前提下，Gemini 建议优先做：

- 从当前测试集抽取代表性序列；
- 通过 AlphaGenome API 或开源工具跑小规模 zero-shot 预测；
- 截取对应轨道，与当前 `corrected B` 的 profile/count 指标做非对称黑盒对比。

这条对比的目的不是证明“我们全方面强于巨无霸”，而是回答：

> 在精确 1-bp ATAC 预测这个任务上，task-specific、bias-aware 的白盒模型是否能在关键指标上形成可解释优势。

### 8.2 杀手锏对比：Tn5 偏置压力测试 / 幻觉峰 stress test

Gemini 特别强调这条对比的解释价值：

- 人为构造强 Tn5 切割偏好、但不含真实 TF motif 的假序列；
- 将其输入 AlphaGenome；
- 观察其是否在这些区域产生“幻觉峰”；
- 再用我们的 `debiased` 口径输出做对照。

如果这条压力测试成立，它的意义不只是“多一个附加实验”，而是直接证明：

> 大规模预训练和大模型容量并不能自动替代显式 bias factorization 的结构先验。

### 8.3 需要明确的限制

这部分目前仍是**外部建议**：

- 我们尚未跑通 AlphaGenome 的当前任务 API/工具链对照；
- 也尚未证明 AlphaGenome 一定会在压力测试里出现幻觉峰；
- 因此这些结论只能作为优先实验建议，不能先写成论文事实。

---

## 9. Gemini 建议的 3 个优先动作

Gemini 最后把建议压缩成 3 个动作：

1. **立刻补算 Count Spearman + AUPRC**
   - 直接离线读取已有双 seed / matched baseline 的 full test JSON。
   - 目的：先把论文和答辩最容易被追问的基础指标补齐。

2. **做小规模 AlphaGenome black-box 探测**
   - 不跑大而全，只抽少量关键测试序列和偏置压力测试序列。
   - 目的：尽快摸到行业黑盒天花板和潜在幻觉风险。

3. **做最小 skip-connection 架构探针**
   - 用当前 `local_tower` 高频特征补回 profile head。
   - 目的：低成本验证“多尺度/decoder 方向”是否真能打破当前 profile ceiling。

这三个动作都还**没有在仓库内执行**；当前它们应被视为“Gemini 给出的高信息量下一步 shortlist”。

---

## 10. 与当前仓库证据的对齐情况

### 10.1 已与现有证据对齐的部分

1. 当前 `corrected B` 是强基座。
2. 当前 measured foundation family 应停表。
3. foundation token 的粗粒度表征与 1-bp dense prediction 存在结构性张力。
4. 论文与答辩口径应从“接入 bug / shortcut 玄学”升到“任务物理约束 + 安全架构设计”。

### 10.2 值得验证但尚未成立的部分

1. skip-connection / decoder probe 是否真能显著改善 JSD；
2. Count Spearman / AUPRC 是否会进一步强化当前主线；
3. AlphaGenome zero-shot 在当前测试集上是否真会弱于或暴露不同风险；
4. bias-only stress test 是否能形成真正有杀伤力的对比图；
5. 非侵入式 conditioner / distillation / base-to-token cross-attention 是否比当前 residual family 更值得继续。

### 10.3 当前不应误写成事实的部分

1. “U-Net 一定能提升”
2. “AlphaGenome 一定会幻觉出偏置峰”
3. “生成式 count/profile 建模一定优于当前 NLL + count 回归”
4. “base-to-token cross-attention 已被证明适合本任务”

---

## 11. 对论文、答辩和评审防守意味着什么

Gemini 的建议对论文叙事最大的影响是：

> 项目的主价值不该继续表述成“我们发现并修复了某个接入 bug”，而应升维成“我们识别了 1-bp dense genomic prediction 对大模型粗粒度特征的结构性约束，并据此指向更安全的下一代多尺度架构设计”。

这对答辩和评审的意义在于：

- 当评委质疑“为什么不继续跟最新预训练基因组大模型”时，我们可以用已有 stop-rule 说明：不是不会接，而是当前接法与任务物理分辨率冲突。
- 当 reviewer 追问“为什么不用 AlphaGenome”时，我们可以把问题从“谁更大”改成“谁在该任务上更安全、更符合 bias-aware 机制”。
- 当 reviewer 追问“为什么指标不够全面”时，Gemini 给出了清晰的补指标优先级。

---

## 12. 当前行动化结论

把这份 Gemini 回流压缩成一句行动化结论：

> 当前最值得优先验证的新方向是“多尺度 / skip-connection decoder 探针 + 更完整的 ATAC 指标防守 + 小规模 AlphaGenome 黑盒探测”，而不是继续在已判负的 foundation token 前向融合 family 上扩线。

下一步不建议直接把这份报告单独转写成论文结论。更合理的做法是：

1. 先等待第二份外部分析（ChatGPT Pro）回流；
2. 再把两份外部建议整理成“共识项 / 分歧项 / 可执行 shortlist”；
3. 然后才决定真正要落的实验与论文补充动作。
