# ChromBPNet Official Externalization Closeout（2026-04-08）

## 一句话结论

`chrombpnet_official_step3_bridge_smoke_20260408_110411` 已在 6000 真机收口并通过；当前唯一经过验证的官方 ChromBPNet 根路径是 `/data1/zhoujiazhen/bylw_atac/chromBPNet`，不是历史文档里残留的 `chrombpnet_official` alias。

## 1. 本轮收口了什么

- 本仓已经退役本地官方 `chrombpnet/` payload 与 packaging 入口，官方源码查阅/复现转到 6000 外部仓。
- helper bridge 的真实 smoke 已在 6000 跑通，run tag 固定在 `run/6000/chrombpnet_official_step3_bridge_smoke_20260408_110411`。
- 这条 infra branch 的冻结快照保留在 `snapshot/autonomy-20260406-chrombpnet-externalization/20260406`。
- 对外操作口径现在应统一成：
  - official root: `/data1/zhoujiazhen/bylw_atac/chromBPNet`
  - canonical archive repo: 本地 `master`
  - runtime repo: `/data1/zhoujiazhen/bylw_atac/TransChromBP`

## 2. 验证证据

- smoke run：`chrombpnet_official_step3_bridge_smoke_20260408_110411`
- 日志：`/data1/zhoujiazhen/bylw_atac/logs/chrombpnet_official_step3_bridge_smoke_20260408_110411.log`
- scratch 输出：`/data1/zhoujiazhen/bylw_atac/.codex_jobs/chrombpnet_official_step3_bridge_smoke_20260408_110411/output/negatives_with_summit.bed`
- 关键事实：
  - 进程已正常退出
  - 日志末尾出现 `Completed execution`
  - `negatives_with_summit.bed` 已落盘，`wc -l=508429`
  - 真正执行到的 helper 是 `/data1/zhoujiazhen/bylw_atac/chromBPNet/chrombpnet/helpers/make_gc_matched_negatives/make_gc_matched_negatives.sh`

## 3. Canonical 档案边界

- implementation plan 继续保留“怎么做”的细节。
- [reports/chrombpnet_official_patch_ledger_20260406.md](reports/chrombpnet_official_patch_ledger_20260406.md) 负责记录仍依赖的官方侧 patch/bridge 点。
- 本文负责“是否已经收口、真实官方根是什么、当前还能不能继续把 implementation plan 当结果报告”的问题。

## 4. 当前状态与下一步

- `chrombpnet_externalization` 仍是一条保留中的 infra branch，因为它还有独有提交和冻结 worktree。
- 但这条线的核心结论已经确定：后续再清理，只剩历史文档里的旧 alias 口径，不再是重新验证 bridge 是否能跑通。
- 若下一轮继续收口，优先做两件事：
  1. 清掉少数历史文档里残留的 `chrombpnet_official` 旧写法。
  2. 评估 `autonomy/20260406-chrombpnet-externalization` 是否已满足卸载 mounted worktree 的条件。
