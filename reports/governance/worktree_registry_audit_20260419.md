# Worktree / Registry Audit（2026-04-19）

## 范围

本轮只处理 `master` 已经完成治理收口之后，仍然挂载在本地的实验 worktree 与 `docs/experiments/registry.md` / annotated tag 之间的不一致。

审计前的挂载集合：

- `master`
- `a6000-dual-track-20260410`
- `autonomy/20260406-chrombpnet-externalization`
- `autonomy/20260406-structure`
- `dual-track-20260409`
- `loss-balance-20260417`
- `loss-balance-bundle-20260416`
- `multiscale-decoder-probe-20260407`

## 发现

### 1. 已 clean 且已有 snapshot / closeout 证据链的挂载

以下 worktree 目前都 clean，且已有 registry row 与 annotated tag，对账上没有缺口：

- `a6000-dual-track-20260410`
- `dual-track-20260409`
- `multiscale-decoder-probe-20260407`
- `autonomy/20260406-structure`
- `autonomy/20260406-chrombpnet-externalization`

这些分支本轮不强行卸载；它们已经满足“可追溯”，后续是否继续收挂，取决于是否还需要保留本地可点开的代码现场。

### 2. `loss-balance-bundle-20260416` 是冗余挂载

这条 worktree 不在 registry 中，也没有自己的 snapshot / closeout tag。进一步检查后发现：

- 分支相对 `master` 没有独有提交；
- `reports/external/chatgpt_bundle_loss_balance_20260416/` 的 bundle 文件已经全部在 `master`；
- 本地脏状态只剩一条早期 `TRACKING.md` 插入。

因此本轮直接移除该 worktree，并删除分支 `loss-balance-bundle-20260416`。

### 3. `loss-balance-20260417` 不是“可直接删”的脏挂载

这条 worktree 在审计前虽然没有独有提交，但有一批未提交文件，里面包含：

- `dynamic_count.py`
- `train_ddp.py` 的 dynamic count 接入改动
- 动态/selector 配置
- loss balance 计划、策略报告与测试

也就是说，它不是单纯的“已经合并到 `master` 的旧现场”，而是保存了 **当时实际用于同步到 remote runtime 的本地 WIP 代码**。因此本轮先把它冻结为提交：

- `c9066b6` `wip: snapshot loss balance worktree`

随后补齐 annotated tag：

- `run/6000/teacher_v2_center_pool_lossbal_e0_selector_s42_6000_20260417_r1`
- `run/6000/teacher_v2_center_pool_lossbal_e2_dynamic_count_s42_6000_20260417_r1`
- `snapshot/loss-balance-20260417/20260419`
- `closeout/loss_balance_curriculum/20260418`

在快照和 tag 都落地后，再卸载该 worktree；分支 `loss-balance-20260417` 保留，作为可继续检视该 snapshot commit 的入口。

## 审计后状态

当前挂载集合变为：

- `master`
- `a6000-dual-track-20260410`
- `autonomy/20260406-chrombpnet-externalization`
- `autonomy/20260406-structure`
- `dual-track-20260409`
- `multiscale-decoder-probe-20260407`

也就是说，本轮把两条 `loss-balance*` 挂载中的：

- 一条冗余 bundle 挂载直接移除
- 一条 dirty implementation 挂载先 snapshot、后移除

`loss_balance_curriculum` 的可追溯入口不再依赖 live mounted worktree，而是改由 `master` 上的 canonical report + run/closeout/snapshot tags 共同承担。

## 结论

- 这轮没有盲目删 dirty worktree；真正需要保留的 dynamic-count 本地实现，已经转成显式 snapshot commit + tag。
- 当前剩余挂载都能从 registry row 和现有 annotated tag 反查职责。
- 如果下一轮还要继续收挂，优先讨论的是那 5 条 clean frozen / infra snapshot 是否还需要保留本地可浏览路径，而不是再处理 `loss-balance`。

## 2026-04-19 Follow-up

同日后续的 archive closeout 又继续收了一轮：

- `dual-track-20260409`
- `autonomy/20260406-structure`
- `autonomy/20260406-chrombpnet-externalization`

这三条 branch 现已分别补上 master-side snapshot / closeout 说明，并通过 branch + snapshot/closeout tag 保留恢复入口，不再依赖 live mounted worktree。当前本地挂载只剩 `master`。
