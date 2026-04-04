# Strict Compare V2 方法论复核（Codex，2026-03-29）

## 结论

- 工程上可做：当前两个 selector 都是“逐 epoch 评估 + 单指标选 best”，增加一个可选的 `--count-r-floor` 开关不复杂。
- 方法学上不建议直接按 Claude 当前文档执行：`count_r floor=0.5` 可以作为 tutorial `Layer 2` 的 robustness follow-up，但当前还不适合升格为新的默认 protocol，更不适合直接平移到 `GM12878` 和多 seed 主线。

## 主要问题

### F1（P1）：`count_r floor=0.5` 是后验阈值，不应被写成默认主 protocol

`docs/plan/strict_compare_v2_updated_methodology_20260329.md` 里已经明确把 floor 的依据写成：

- official best epoch 35 的 `count_r=0.65805`
- 当前这条 TransChromBP run 在 epoch `1-28` 的 `count_r=0.73-0.81`
- 因而推荐 `floor=0.5`

这说明 `0.5` 不是事先冻结的规则，而是看过本轮 tutorial `L2` 结果之后，根据同一批数据反推出来的阈值。它本质上仍然是一个新超参，不能写成“没有引入新超参”或“只是硬安全门槛”。

如果直接把它写进主 protocol，会让外部比较从“统一外部 selector”变成“针对本次异常 run 追加的一条救火规则”。这会削弱主表说服力，尤其是在目前还只有一个 seed、一套 tutorial 数据的情况下。

### F2（P1）：决策顺序过快，跳过了当前最值钱的 locked follow-up

现有执行稿已经把更稳妥的下一步写得比较清楚：先把 `epoch 35` / `epoch 44` 这组 checkpoint 作为 follow-up 写清楚，再决定是补独立 test、调整 selector，还是进入 `L3`。论文 claim matrix 也把这件事列成当前最值得补的实验。

Claude 这版文档把“复合 selector”直接上升为 `Phase 1`，并据此做 `Go/No-go` 到 `L3` 的决策。这一步太快了，因为它把一个单 run 的诊断现象，直接升级成了全流程方法论更新。

更稳妥的顺序应该是：

1. 锁定当前已观察到的 checkpoint（至少 `epoch 35`、`epoch 44`、official `epoch 35`）。
2. 把这些固定 checkpoint 的 `valid/test` 指标补成一张 follow-up 表。
3. 在“不用 test 重新选 best”的前提下，判断 count 崩塌是否只是 `valid` 偶然现象，还是 held-out 上也成立。
4. 只有在这一步成立后，才讨论是否把复合 selector 升级成正式 protocol。

### F3（P2）：协议定义还不够完整，当前更适合“可选 robustness 开关”而不是“默认规则”

如果真的要加 `--count-r-floor`，当前文档还缺几条必须预先锁定的行为定义：

- `count_r` 为 `NaN` 时是否直接视为不通过 floor
- metric 为 `NaN` 但 count 正常时是否参与候选
- 多个 epoch 同时满足 floor 且 metric 并列时如何 tie-break
- 全部 epoch 低于 floor 时，`fallback` 是否仍允许写入主表，还是只作为 warning case 单列
- `GM12878` / 多 seed 是否沿用同一个绝对阈值，还是先做 tutorial 灵敏度分析再决定

从代码现状看，这两个 selector 目前只会输出 `metric/mode/split/n_checkpoints/best`；如果要把 floor gate 作为可复用分析口径，至少还需要把 `count_r_floor`、`count_r_metric`、`excluded_epochs`、`fallback_used` 这类字段稳定写进产物，并同步更新后续汇总表的读取逻辑。否则不同人后续只看 `best_epoch.json`，仍然无法直接判断这次 best 是原始 selector 还是 floor-filtered selector。

## 建议执行顺序

1. 先不要把 composite selector 升格为新的默认 protocol。
2. 先补 tutorial `L2` 的 locked follow-up 表：
   - official `epoch 35`
   - TransChromBP `epoch 35`
   - TransChromBP `epoch 44`
   - `valid/test` 都汇报，但 test 只用于固定 checkpoint 的描述，不参与重新选 best
3. 如果仍想实现 selector 改动，可以先做成“默认关闭、显式传参启用”的可选开关。
4. tutorial 上至少补一轮 floor 灵敏度表（例如 `0.4 / 0.5 / 0.6`），再决定是否推广到 `GM12878` 和多 seed。
5. `L3 shared-region compare` 仍然是提升主结果说服力的硬价值项；复合 selector 不能替代 `L3` 的作用。

## 对实施的直接建议

- 可以写代码，但先不要把它当成“确认后的正式方案”。
- 如果现在就要动 selector，我建议只接受下面这个范围：
  - 增加可选 `--count-r-floor`
  - 默认值保持“关闭”
  - `best_epoch.json` 明确写出 floor 是否启用、排除了哪些 epoch、是否发生 fallback
  - 文档口径明确标成 `robustness / sensitivity analysis`

在这个边界内，代码改动是低风险的；超出这个边界，把它直接写成 `V2` 默认方法论，我不建议现在推进。
