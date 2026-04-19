# Loss 协调相关证据与失败模式

## 1. 为什么这次要单独讨论 loss 协调

我们这次不是因为看了一篇文章，就想随手给训练器加一个“动态权重”模块。

真正的原因是：项目里已经出现了多条和 loss / selector 直接相关的证据，说明当前问题并不只是“模型结构还不够强”，而可能也包括：

- 多个目标之间的优化关系不理想；
- 训练时优化的混合 loss 与最终关心的指标不完全一致；
- 某些来自别的框架的 loss 经验值，迁到当前结构后会改变训练动力学。

---

## 2. 证据 1：`debiased_profile_weight=2.0` 比 `stop-gradient` 更像主效应

来自 clean matrix 的关键结果：

| 条件 | `profile_full_debiased_jsd` | 解释 |
|---|---:|---|
| `TF + sg=false + deb2` | `0.00147` | 即使关掉 `sg`，只要保留 `deb2`，依然很干净 |
| `TF + sg=true + deb0` | `0.00972 / 0.00997` | 保留 `sg` 但去掉 `deb2`，gap 明显抬升 |
| `corrected B = TF + center pool + sg=true + deb2` | `0.00268` | 最终默认模型也保持很小 gap |

这说明当前最稳的解释不是：

- `stop-gradient` 单独修复了问题

而是：

- 显式的 debiased profile supervision 更像主效应；
- `stop-gradient` 更像辅助稳态项。

对 loss 协调问题的意义是：

> 当前框架里，loss 项的存在与否本身已经会显著改变训练语义，而不只是轻微地平衡梯度尺度。

---

## 3. 证据 2：`loss_total` 与 `JSD` selector 不一致

当前已经确认：

- `peak.loss_total` 是混合指标
- 它包含 profile NLL、count MSE、debiased profile NLL 等项
- 但最终真正更关心的 paper-facing 指标是 `peak.profile_target_jsd_full_mean`

已有审计结论明确指出：

1. NLL 和 JSD 不完全单调；
2. 对 count 改进实验，`count_weight=0.1` 可能让收益在 `loss_total` 里被掩盖；
3. 对 profile 改进实验，训练后期 profile NLL 和 JSD 的发散会更明显。

因此，这里不是单纯的“选模配置小问题”，而是：

> 当前混合 loss、验证 selector、最终报告主指标之间，已经出现了对齐问题。

如果外部模型只给出“把各 loss 动态调到同量级”这一类建议，而不处理 selector / final metric 对齐问题，就很可能不够。

---

## 4. 证据 3：`count_weight_strategy=chrombpnet_auto` 不能直接迁移

我们已经记录过一个关键训练动力学结论：

- 在原始 ChromBPNet 的刚性 bias 分解语境里，高 `count_weight` 有一套自己的经验值；
- 但在当前 `learnable_scales + additive/logsumexp fusion` 的柔性融合语境下，同样的经验值可能给优化器一条更容易的捷径：
  放大 bias 相关分支，而不是认真优化 debiased signal 分支。

这件事的直接含义是：

1. 当前 loss 权重并不只是“量纲归一化”问题；
2. 同样的权重策略，在不同结构里会改写优化路径；
3. 外部模型如果推荐动态 weight 方法，必须解释它在 `learnable_scales + bias branch` 语境下为何不会把问题变得更糟。

---

## 5. 证据 4：zero-count profile loss 的旧语义会引入额外监督压力

我们已经改过一条高置信训练语义：

- 旧版本里，全零窗口会被构造成近似均匀的 profile target；
- 这会给背景窗口额外施加 profile 监督压力；
- 修正后，全零窗口在 profile loss 里记为 0 贡献。

这条证据的重要性在于：

> “多个 loss 怎么平衡”并不只取决于标量权重，还取决于每一项 loss 自己的 target 语义是否合理。

如果外部模型只谈动态权重，而不检查各 loss 的定义语义，就会漏掉这一层。

---

## 6. 证据 5：`track_total_count_target` 与 auto count-weight 曾经耦合

我们还留过一条专门的训练语义提醒：

- `track_total_count_target`
- `count_weight_strategy=chrombpnet_auto`

这两个变量曾经在同一条 hybrid 路径里一起动，导致很难判断退步来自哪里。

这说明：

1. count 端目标本身不是纯净单变量；
2. 任何动态 loss 协调方案，如果直接包住现有全部 count 相关项，可能会把不同问题缠在一起；
3. 更稳的做法可能是先明确哪些项是核心主目标，哪些项只是辅助训练约束。

---

## 7. 这次外部模型要特别注意的失败模式

### 失败模式 A：把 `full` 和 `debiased` 当成两个普通并列任务

它们共享主干，又通过 bias 融合与去偏语义相互耦合，不是标准多任务学习里的两个独立 head。

### 失败模式 B：只平衡 loss 数值量级，不看最终主指标

即使各 loss 都被调到同量级，也不代表 `JSD` 会更好，尤其当 NLL 与 JSD 不完全单调时。

### 失败模式 C：把来自别的框架的加权经验直接迁过来

ChromBPNet 风格的经验值，在 `learnable_scales` 语境下不必然安全。

### 失败模式 D：只谈训练 loss，不谈 selector

如果最好的 checkpoint 是靠另一个指标定义的，那么“loss 平衡”方案本身就不完整。

### 失败模式 E：动态协调把 bias shortcut 放大

如果某类方法倾向于优先优化更容易下降的项，它可能进一步鼓励模型走 bias-related 的捷径，而不是强化 debiased signal branch。

---

## 8. 这份文件想让外部模型回答什么

读完这份文件后，我们希望外部模型先能回答：

1. 在我们这个框架里，动态 loss 协调最需要解决的主矛盾到底是什么？
2. 这个问题更像“目标冲突”、 “梯度尺度失衡”、 “结构耦合导致的伪冲突”，还是 “selector 与最终指标错位”？
3. 哪一类近期方法最有希望解决这个问题，哪一类方法看起来很合理但其实大概率无效？

如果它不能先回答这三个问题，再多论文名词也没有真正价值。
