# Project Broad Audit Cleanup（2026-04-19）

## 本轮目标

对项目做一次“大范围探查 + 可直接落地的整理”，重点不是继续大迁移，而是先把最容易误导维护者的 live 文档、历史入口和 worktree 状态收成一致。

## 本轮实际收口

1. 统一 live 计划命名空间
   - 仍被 `TRACKING.md` 直接引用的 factor-ladder 设计/实现计划已迁回 `docs/plan/`。
   - `docs/superpowers/` 补了 `README`，明确定义为历史 agent 规划档案，不再充当第二套 live planning root。

2. 补 canonical 档案
   - 新增 [chrombpnet_official_externalization_closeout_20260408.md](chrombpnet_official_externalization_closeout_20260408.md)。
   - 新增 [chrombpnet_official_patch_ledger_20260406.md](chrombpnet_official_patch_ledger_20260406.md)。
   - 新增 [foundation_cache_contract_snapshot_20260406.md](foundation_cache_contract_snapshot_20260406.md)。
   - `docs/experiments/registry.md` / `docs/experiments/runs.csv` 已改成引用这些更合适的 canonical 报告路径。

3. 下线明显过时入口
   - `docs/plan/archive/compare/strict_chrombpnet_official_comparison_execution_20260327.md` 现在有显眼的“历史执行稿”提示。
   - 过时的 `scripts/deploy_strict_compare_staging_to_6000.sh` 已退役。
   - `TRACKING.md` 的 2026-03-17 历史条目已改成当前 `publish-runtime-* / pull-results-*` 契约。

4. 收挂已闭环的 clean snapshot
   - `a6000-dual-track-20260410`
   - `multiscale-decoder-probe-20260407`
   - `dual-track-20260409`
   - `autonomy/20260406-structure`
   - `autonomy/20260406-chrombpnet-externalization`

这些 mounted worktree 已卸载；registry 改为通过对应 `branch + snapshot/closeout tag` 恢复，而不是继续依赖本地 live path。

## 当前仍保留的挂载

- 只剩 `master`

## 下一步建议

下一轮最值得做的是目录 rename，把 `chromBPNet` 物理改名和绝对路径修复一起做；如果还要继续减噪，再看历史报告里的绝对路径是否值得统一批量改写。
