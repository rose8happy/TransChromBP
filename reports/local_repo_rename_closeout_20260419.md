# Local Repository Rename Closeout（2026-04-19）

## 结论

本地 canonical root 的文本切换、历史绝对路径批量改写和真实物理 rename 已在同一轮完成。当前物理根是 `/home/zhengwei/project/python/TransChromBP`；旧路径 `/home/zhengwei/project/python/chromBPNet` 已降级为兼容 symlink / rollback entrypoint，不再作为 canonical root。pre-rename 与 post-rename 两阶段验证证据都保留在本报告中。

这次收口只影响本地仓入口和仓内文档。6000 官方 lookup root、6000 runtime workspace 和 6002 runtime mirror 都保持原状。

## 最终路径拓扑

- local canonical physical root: `/home/zhengwei/project/python/TransChromBP`
- local compatibility symlink: `/home/zhengwei/project/python/chromBPNet -> TransChromBP`
- 6000 official lookup root: `/data1/zhoujiazhen/bylw_atac/chromBPNet`
- 6000 runtime workspace: `/data1/zhoujiazhen/bylw_atac/TransChromBP`
- 6002 runtime mirror: `/home/zhengwei/bylw_atac/TransChromBP`

## 实际执行的切换步骤

1. 先完成 live/canonical 文档与历史报告里的 repository-local absolute links 批量改写，把本地主入口统一到 `TransChromBP`。
2. 在真实物理 rename 前，用临时 forward bridge `/home/zhengwei/project/python/TransChromBP -> chromBPNet` 兜住已改写的新链接，避免 dead links。
3. 合并 `local-repo-rename-20260419` 到 `master`，并移除对应实现 worktree / 分支。
4. 清理临时 forward bridge，把物理目录从 `chromBPNet` 改名为 `TransChromBP`。
5. 建立最终兼容 symlink `/home/zhengwei/project/python/chromBPNet -> TransChromBP`。
6. 更新 `TRACKING.md`、`TRACKING_archive.md` 与本报告，并从新根目录复跑治理验证。

## Pre-Rename Verification

真实物理 rename 之前，已先把文本口径与治理约束收口到可安全切换的状态：

| 命令 | 结果 |
|---|---|
| `bash -n scripts/sync_project.sh` | PASS |
| `/home/zhengwei/project/python/chromBPNet/.venv/bin/python -m pytest -q tests/test_sync_project_contract.py tests/test_repository_governance_docs.py tests/test_factor_ladder_docs.py` | PASS, `23 passed in 0.19s` |
| `rg -n '/home/zhengwei/\\.config/superpowers/worktrees/TransChromBP' AGENTS.md README.md TRACKING.md docs reports scripts tests workflows` | PASS, no matches after removing the stale fake-path reference |

## Post-Rename Verification

| 命令 | 结果 |
|---|---|
| `cd /home/zhengwei/project/python/TransChromBP && pwd` | `/home/zhengwei/project/python/TransChromBP` |
| `bash -n scripts/sync_project.sh` | PASS |
| `/home/zhengwei/project/python/TransChromBP/.venv/bin/python -m pytest -q tests/test_sync_project_contract.py tests/test_repository_governance_docs.py tests/test_factor_ladder_docs.py` | PASS, `23 passed in 0.22s` |
| `git rev-parse --show-toplevel` | `/home/zhengwei/project/python/TransChromBP` |
| `ls -ld /home/zhengwei/project/python/chromBPNet /home/zhengwei/project/python/TransChromBP` | PASS, `TransChromBP` 是物理目录，`chromBPNet -> TransChromBP` 是兼容 symlink |

## Final Residual Old-Path List

对 `AGENTS.md README.md TRACKING.md docs reports scripts tests workflows` 复跑
`rg -n '/home/zhengwei/project/python/chromBPNet' ...` 后，最终 residual hits 归为以下几类；每一类都属于计划允许的受控例外：

1. live compatibility / rollback wording
   - `AGENTS.md`
   - `README.md`
   - `TRACKING.md`

   这些命中不是把旧路径当成 canonical root，而是在显式声明“旧路径只剩 compatibility symlink / rollback path”。

2. test contract
   - `tests/test_repository_governance_docs.py`

   这里必须直接保留旧路径字符串，才能验证 live 文档没有把它重新写回 canonical root。

3. rename narrative / audit trail
   - `docs/superpowers/specs/2026-04-19-local-repo-rename-design.md`
   - `docs/superpowers/plans/2026-04-19-local-repo-rename.md`

   这些文档本身就是 rename 设计与执行记录，必须保留旧根作为历史主语。

4. closeout / evidence artifacts
   - `reports/local_repo_rename_closeout_20260419.md`
   - `reports/assets/local_repo_rename_20260419/local_repo_path_residuals_after_rewrite.txt`

   前者需要记录最终 symlink 拓扑与前后验证命令；后者是 pre-rename 批量改写后的冻结 residual 清单，本来就应保留旧路径作为审计快照。

## Closeout Judgment

本地仓 rename 已完成 closeout。后续若再处理路径相关噪音，目标应该只剩批量清理更多历史报告里的旧本地绝对路径，或在兼容期结束后决定何时移除 `chromBPNet -> TransChromBP` symlink；这两件事都不再阻塞当前 canonical repo 使用。
