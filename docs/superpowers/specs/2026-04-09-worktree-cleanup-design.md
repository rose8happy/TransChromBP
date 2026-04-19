# Worktree Cleanup Design

## Goal

整理 `/home/zhengwei/project/python/TransChromBP` 当前混杂在 `master` 与多个 worktree 里的改动，把它们重新归位到各自负责的分支，并为尚未 ready 的实验线留下本地 `wip` 提交，避免主档案仓继续承担实验 WIP 的混合状态。

## Scope

- 清点 `master`、`dual-track-20260409`、`autonomy/20260406-structure`、`autonomy/20260406-chrombpnet-externalization`、`multiscale-decoder-probe-20260407`
- 把 `master` 中属于实验线或外置化线的代码/脚本改动迁回对应 worktree
- 把主档案仓应保留的 docs / reports / tracking 改动留在 `master`
- 为各 worktree 创建本地提交，未 ready 内容允许用 `wip: ...`

## Non-Goals

- 不 push
- 不强行把 WIP 重构成可合并 feature
- 不重写实验结论，只做归位与可追溯整理

## Target State

- `master` 只保留主规范、长期文档、paper / report / tracking / archive 内容
- `dual-track-20260409` 只保留 `teacher-distill + U-Net-lite` 相关改动
- `autonomy/20260406-structure` 只保留 foundation cache contract 重构
- `autonomy/20260406-chrombpnet-externalization` 只保留 official ChromBPNet 外置化相关改动
- 每条线至少有一个本地提交可追溯当前状态

## Execution Notes

- 优先按文件归属而不是按“最近修改位置”判断归位
- 对于 `master` 与目标 worktree 同路径但内容不一致的文件，以当前 `master` 内容为准迁入对应线，避免丢失新状态
- `TRACKING.md` 以 `master` 为 canonical；worktree 内若仍保留 branch-local `TRACKING.md`，视作该分支的 WIP 快照而非主状态源
