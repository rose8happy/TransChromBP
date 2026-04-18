# Repository Governance

本文件定义当前仓库的 repository topology、single source of truth 和同步方向。若与 `AGENTS.md` 冲突，以 `AGENTS.md` 为准。

## Canonical Roles

| Root | Role | Can become source of truth? | Notes |
|---|---|---:|---|
| `/home/zhengwei/project/python/chromBPNet` | 本地 canonical trunk | `yes` | 唯一有完整 git / branch / worktree 历史的 TransChromBP 主仓。当前物理目录名暂不改，后续若要 rename，单开任务做。 |
| `/data1/zhoujiazhen/bylw_atac/chromBPNet` | 6000 official lookup / reproduction root | `yes`（仅官方源码） | 只用于官方 ChromBPNet 查阅、复现与对照；不是 TransChromBP runtime。 |
| `/data1/zhoujiazhen/bylw_atac/TransChromBP` | 6000 runtime workspace | `no` | 允许临时热修和运行期文件，但不能直接成为事实真源。 |
| `/home/zhengwei/bylw_atac/TransChromBP` | 6002 runtime mirror | `no` | 无 `.git`，只作为运行/部署镜像。 |

## Sync Model

- 采用“单一真源 + 单向发布”。
- 所有正式代码、脚本、文档、实验摘要，先回收到本地 canonical trunk，再发布到 runtime。
- 远端出现有效临时改动时，默认流程是：
  1. 从 6000/6002 回收改动到本地仓。
  2. 在本地整理、测试、提交。
  3. 通过显式目标命令重新发布到 runtime。
- 禁止把 6000 `TransChromBP` 或 6002 `TransChromBP` 描述成“另一份主仓”或“三边同步节点”。

## Command Contract

标准入口统一通过 `scripts/sync_project.sh`：

- `publish-runtime-6000`
- `publish-runtime-6002`
- `pull-results-6000`
- `pull-results-6002`
- `status-all`

说明：

- `publish-*` 只发布 repo-tracked 代码、脚本和文档；不会覆盖 `logs/`、`outputs/`、数据和大文件。
- `pull-results-*` 只拉回 runtime 下的 `logs/` 和 `reports/`，不会整体同步 `outputs/`。
- 所有发布/回收命令都支持 `--dry-run`。
- 旧命令 `deploy` / `download_results` 视为 deprecated compatibility path，必须报错并提示显式目标命令。

## Directory Hygiene

- `reports/` 只保留可复用 `md/tex/assets` 源文件。
- `references/` 只保留轻量级索引；大文件本体默认进 `references/local-only/`。
- `tmp_remote_edit/` 只做 staging，不做长期归档。
- 根目录只保留稳定入口文件；专题报告、研究笔记、学习材料和脚本都应落回各自命名空间。
