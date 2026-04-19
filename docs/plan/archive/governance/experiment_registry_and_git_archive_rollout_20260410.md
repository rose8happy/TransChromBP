# 实验台账与 Git 归档收口方案（2026-04-10）

## 一句话目标

把当前“事实分散在 `master`、worktree、远端运行目录、`TRACKING.md`、`reports/`”的状态，收敛成 `master` 单入口的 canonical 档案体系，避免实验信息、远端路径与结论继续散失。

## 本轮落地范围

- 将 `a6000-dual-track-20260410` 的终态 closeout 文档择净回流到 `master`
- 新增 `docs/experiments/registry.md` 作为 family / workstream 索引
- 新增 `docs/experiments/runs.csv` 作为 run 级 manifest
- 在 `AGENTS.md` 写死 registry、annotated tag 和 worktree 收起前的最小归档规则

## 角色分工

- `TRACKING.md`：只保留 live 状态
- `TRACKING_archive.md`：阶段性完成事项归档
- `docs/experiments/registry.md`：family / workstream 的 canonical 索引
- `docs/experiments/runs.csv`：run 级元信息、远端路径与 verdict manifest
- `reports/`：closeout、gate 判读与可复用证据链

## Git 规则

- `master` 是唯一 canonical 档案入口
- `run/<machine>/<run_id>`、`closeout/<family_id>/<yyyymmdd>`、`snapshot/<branch_slug>/<yyyymmdd>` 一律使用 annotated tag
- 当 6000/6002 运行仓没有可靠 git 历史时，`runs.csv` 的 `commit_sha` 记录本地 canonical archival commit，而不是伪造远端真实 commit
- worktree 只有在 `registry + runs.csv + report + tag` 四项齐全后，才允许被收起

## 当前批次的直接结果

- `A6000 formal gate`、`AlphaGenome v2 sidecar`、`U-Net-lite v1` 的 terminal state 已回流 `master`
- `dual-track-20260409`、`multiscale-decoder-probe-20260407`、`autonomy/*` 这些 worktree 现在都能从 registry 反查其职责与状态
- 后续任何新 run 的最小文档同步要求，变成：`TRACKING.md + runs.csv + report（若收口）`
