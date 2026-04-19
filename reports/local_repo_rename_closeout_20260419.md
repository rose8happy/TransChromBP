# Local Repository Rename Closeout（2026-04-19）

## 结论

本轮已经完成 pre-rename 文本收口：仓内 repository-local absolute links 已改到新的本地 canonical root，历史与治理文档中的旧路径只保留在受控的 rename narrative / rollback 语境里。`/home/zhengwei/project/python/TransChromBP` 现在是当前对外主入口，但它仍然只是一个临时 forward symlink，指向旧物理目录 `/home/zhengwei/project/python/chromBPNet`，用于在真实物理 rename 前避免 dead links。

真实的物理 rename 还没有执行。下一步应该把目录本体从 `chromBPNet` 改成 `TransChromBP`，然后移除当前临时 forward bridge，并创建最终的兼容 symlink `/home/zhengwei/project/python/chromBPNet -> /home/zhengwei/project/python/TransChromBP`。

## 旧根与新根

- old local canonical root: `/home/zhengwei/project/python/chromBPNet`
- new local canonical root: `/home/zhengwei/project/python/TransChromBP`
- current temporary bridge: `/home/zhengwei/project/python/TransChromBP -> chromBPNet`

## Compatibility Symlink Plan

当前阶段的桥接策略分两层：

1. pre-rename temporary bridge
   - 目的：让已经改写到 `TransChromBP` 的 repository-local absolute links 在物理目录尚未移动前仍然可解析。
   - 现状：`/home/zhengwei/project/python/TransChromBP` 已作为 forward symlink 指向旧目录。
   - 约束：这不是最终兼容方案，必须在真实物理 rename 任务中清理。

2. post-rename compatibility symlink
   - 目的：让旧调用方继续通过 `/home/zhengwei/project/python/chromBPNet` 访问新物理根。
   - 目标：真实 rename 完成后创建 `/home/zhengwei/project/python/chromBPNet -> /home/zhengwei/project/python/TransChromBP`。
   - 责任边界：该兼容 symlink 只服务过渡期，不再是 canonical entrypoint。

## 物理 Rename 前验证

在执行任何物理 rename 之前，已按要求完成治理 / 文档验证：

| 命令 | 结果 |
|---|---|
| `bash -n scripts/sync_project.sh` | PASS |
| `/home/zhengwei/project/python/chromBPNet/.venv/bin/python -m pytest -q tests/test_sync_project_contract.py tests/test_repository_governance_docs.py tests/test_factor_ladder_docs.py` | PASS, `23 passed in 0.19s` |
| `rg -n '/home/zhengwei/\\.config/superpowers/worktrees/TransChromBP' AGENTS.md README.md TRACKING.md docs reports scripts tests workflows` | PASS, no matches after removing the single stale fake-path reference |

## Residual Old-Path Exceptions

仍保留旧路径的地方只有两类：

1. rename narrative / audit trail
   - `docs/superpowers/specs/2026-04-19-local-repo-rename-design.md`
   - `docs/superpowers/plans/2026-04-19-local-repo-rename.md`

   这些文档必须继续描述旧根，因为它们本身就是 rename 设计、验证与回溯记录，不能把历史主语抹掉。

2. bridge / rollback semantics
   - 当前临时 forward symlink 仍然指向旧物理目录，目的是防止 pre-rename 阶段的 dead links。
   - 真实物理 rename 完成后，旧路径会改成最终的 compatibility symlink，而不是继续作为 canonical root。

## 未更改的远端根

本轮只处理本地 canonical root 的文本与审计状态，没有改动远端官方 / runtime roots：

- 6000 official lookup root: `/data1/zhoujiazhen/bylw_atac/chromBPNet`
- 6000 runtime workspace: `/data1/zhoujiazhen/bylw_atac/TransChromBP`
- 6002 runtime mirror: `/home/zhengwei/bylw_atac/TransChromBP`

这些远端根保持原状，是刻意的，不是遗漏。

## 审计结论

`reports/local_repo_rename_closeout_20260419.md` 是 pre-rename 状态的 canonical audit trail。等真实物理 rename 完成后，应该更新 TRACKING、清理临时 forward bridge，并把旧路径切成最终兼容 symlink，再补一轮相同的治理验证。
