# Tutorial L2 Locked Follow-up 复核（Codex，2026-03-29）

## 结论

- 这轮 locked follow-up 是有价值的。它明显削弱了“epoch 44 的 count_r 崩塌代表训练真实退化”这条判断。
- 但 Claude 当前报告里的两句结论仍然写得过满：
  - “崩塌原因很可能就是 valid 只有 1 条染色体”
  - “profile-only selector 当前足够，复合 selector 暂不需要”

更准确的表述应该是：

- 对 **本次 tutorial L2 单 seed run**，valid 选中的 `epoch 44` 在 test 上表现健康，说明这一轮没有出现“选中坏 checkpoint 后测试也崩”的情况。
- 但这更像是“当前 run 的 held-out 表现证明 selector 这次没有踩空”，**还不足以证明 selector 健壮性问题已经解决**；问题只是从“训练后期真崩”转成了“valid/test 对 count 指标的关系不稳定”。

## 主要问题

### F1（P1）：把原因直接归结为“valid 只有 1 条染色体”证据不足

Claude 报告把 epoch 44 的 valid/test 差异，直接解释为：

- `fold_0` 的 valid split 仅 1 条染色体
- 样本量过小
- 所以 valid count_r 估计不稳定

但仓库里的 strict-compare README 明确写了 `fold_0` 是：

- `test chr1`
- `valid chr2`

也就是说，**test 也是单染色体 split**。因此“valid 只有 1 条染色体”本身并不能解释为什么 valid count_r 会掉到 `0.084`，而 test 却能到 `0.838`。要支撑这个因果解释，至少还需要下面之一：

1. valid/test 的实际 `n_examples` 对照
2. chr1/chr2 上 peak/nonpeak 分布差异
3. 多 fold 重复，证明这是 fold-specific 现象而不是偶然

所以目前最多只能写成：

- “这更像 split-specific metric instability / sample-composition effect”

还不能写成：

- “已经确认只是因为 valid 只有 1 条染色体”

### F2（P1）：follow-up 足以降低 selector 改造优先级，但不足以证明“profile-only selector 已经足够”

locked follow-up 的 test 结果说明：

- `epoch 35` 和 `epoch 44` 在 test 上都很好
- `epoch 44` 虽然在 valid count_r 很差，但 test 并没有坏

这能支持的最强结论是：

- **本次 run 中**，profile-only selector 选中的 `epoch 44` 在 held-out test 上没有出问题

但它**不能**自动推出：

- profile-only selector 作为一般规则已经足够
- 复合 selector 没有再实现的价值

因为 selector 的输入仍然是 `valid`，而这次 follow-up 反而提示：**同一个 checkpoint 的 valid/test count 指标关系可能很不稳定**。这说明问题没有消失，只是从“epoch 44 是坏点”变成了“当前 valid count_r 对真实 held-out 表现的代表性不足”。

因此更稳妥的动作应该是：

- 复合 selector 不再是当前 blocker，可以降级
- 但不要把它正式判死；如果后面要补 selector 保护，仍然应以“默认关闭的可选 safeguard”方式保留

## 可以接受的结论

下面这些结论，我认为现在可以接受：

1. tutorial `L2` 这条单 seed matched-budget system compare，TransChromBP 的 held-out test 指标优于 official。
2. `epoch 44` 不应再被描述为“明显坏 checkpoint”。
3. 复合 selector 不再是进入 `L3` 之前的硬 gate。
4. 下一步优先级可以转到 `L3 shared-region compare` 或 `GM12878 L2`。

## 不建议直接写满的结论

下面这些话现在还不建议直接写进主结论：

1. “count_r 崩塌已经确认只是 valid 太小导致的噪声”
2. “profile-only selector 当前已经足够”
3. “复合 selector 暂不需要，因此相关分析可以关闭”

## 建议口径

建议把 locked follow-up 的结论改写为：

> locked follow-up 显示，本次 tutorial `L2` run 中，valid 选中的 `epoch 44` 在 held-out test 上同样健康，说明先前观察到的 valid count_r 异常并未转化为 test 退化。当前更准确的问题表述不再是“训练后期真实崩塌”，而是“valid 与 test 对 count 指标的一致性有限，selector 健壮性仍需在更多 fold / 数据集上继续观察”。

## 下一步建议

1. 可以直接推进 `L3 shared-region compare`。
2. selector 改造不再是当前 blocker；若后续实现，建议保持“默认关闭 + 显式启用”。
3. 若后面还要追这条线，优先补的是：
   - 多 fold / 多 dataset 的 selector 稳定性
   - valid/test 的 `n_examples` 与 region 组成差异

