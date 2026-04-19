# Local Repository Rename Design（2026-04-19）

## Summary

本轮目标是把本地 canonical 仓从 `/home/zhengwei/project/python/chromBPNet` 改名为 `/home/zhengwei/project/python/TransChromBP`，并同步收口仓内仍依赖旧本地绝对路径的文档、测试与说明。为降低切换风险，本轮同时保留一个兼容 symlink：`/home/zhengwei/project/python/chromBPNet -> /home/zhengwei/project/python/TransChromBP`。

这次 rename 只处理“本地 canonical 仓”相关路径，不改变 6000 official root `/data1/zhoujiazhen/bylw_atac/chromBPNet` 的角色或命名。

## Goals

- 让本地 canonical 仓的物理目录名与当前项目身份 `TransChromBP` 对齐。
- 把仓内所有仍指向旧本地 canonical 根的 live/canonical 绝对路径改到新路径。
- 批量修正历史报告、历史 spec、测试和 README 中对本地仓文件的绝对链接，避免改名后变成失效链接。
- 保留旧路径 symlink 作为过渡层，避免仓外残留脚本、命令和当前工具会话立刻断开。

## Non-Goals

- 不重命名 6000 official ChromBPNet 根路径 `/data1/zhoujiazhen/bylw_atac/chromBPNet`。
- 不把已卸载 worktree 的历史现场路径机械改成新目录名。
- 不在本轮删除兼容 symlink。
- 不借这次 rename 顺手改实验 family、branch、tag 或远端 runtime 拓扑。

## Path Classes And Rewrite Rules

### 1. Local canonical root

凡是表达“当前本地主仓文件在哪里”的绝对路径，统一从：

- `/home/zhengwei/project/python/chromBPNet`

改为：

- `/home/zhengwei/project/python/TransChromBP`

适用范围：

- `AGENTS.md`
- `README.md`
- `TRACKING.md`
- `docs/env/repository_governance.md`
- `docs/experiments/registry.md`
- `docs/experiments/runs.csv`
- `tests/` 中针对本地 canonical 根的断言
- `reports/`、`docs/plan/`、`docs/superpowers/`、脚本说明文件中指向本地主仓文件的绝对链接
- LaTeX/Markdown 中用于点击跳转本地主仓文件的绝对路径

### 2. Historical worktree paths

以下模式的路径默认不做“同构替换”：

- `/home/zhengwei/.config/superpowers/worktrees/chromBPNet/...`

原因：

- 这些路径描述的是旧 worktree 的历史现场；
- 多数对应 worktree 已经卸载；
- 直接改成 `.../TransChromBP/...` 会制造并不存在的新历史路径。

处理规则：

- live/canonical 文档里，如果这些路径只是恢复入口，应改写为 branch/tag/closeout 说明；
- 如果它们承担“当时现场证据”作用，则保留旧路径，并在文字里明确是 historical path；
- 如果只是普通跳转链接且已无实际目标，应改为当前 canonical 文件链接或删除。

### 3. Remote official/runtime paths

以下路径不属于本地 rename 范围：

- `/data1/zhoujiazhen/bylw_atac/chromBPNet`
- `/data1/zhoujiazhen/bylw_atac/TransChromBP`
- `/home/zhengwei/bylw_atac/TransChromBP`

它们分别对应 6000 official root、6000 runtime、6002 runtime，角色由现有治理文档维持不变。

## Execution Plan

### Phase 1. Canonical documentation first

先把以下文档改到新口径：

- `AGENTS.md`
- `TRACKING.md`
- `docs/env/repository_governance.md`
- 与仓名/本地根路径相关的 README、治理报告、handoff 文档

这样即使后续物理 rename 还没执行，仓内也已经有清晰的 canonical 说明。

### Phase 2. Bulk rewrite repository-local absolute links

批量扫描并修正仓内引用：

- Markdown 链接
- LaTeX `\href{...}{...}`
- shell/README 示例命令
- 测试中的字符串断言

执行时优先区分三类文本：

1. 当前 canonical 说明：直接改到新路径。
2. 纯导航性历史文档：通常也改到新路径，保证点击可用。
3. 历史现场证据：保留旧路径，并增加语义说明，避免被误读成当前 live 路径。

### Phase 3. Validation before physical rename

在实际改目录前完成：

- `rg` 检查仓内是否还残留旧本地 canonical 根路径；
- `rg` 检查是否出现把历史 worktree 路径误改成 `.../worktrees/TransChromBP/...` 的伪历史；
- 运行文档治理相关测试；
- 运行关键脚本的语法检查。

### Phase 4. Physical rename with compatibility symlink

在父目录执行：

1. 把 `/home/zhengwei/project/python/chromBPNet` 改名为 `/home/zhengwei/project/python/TransChromBP`。
2. 创建兼容 symlink：
   `/home/zhengwei/project/python/chromBPNet -> /home/zhengwei/project/python/TransChromBP`
3. 进入新路径后检查：
   - `git rev-parse --show-toplevel`
   - `git status`
   - 关键测试与脚本

### Phase 5. Post-rename closeout

- 更新治理报告，说明新 canonical 物理路径已生效。
- 保留 symlink 作为过渡层，不在同一轮移除。
- 后续单开一轮判断何时删除 symlink，并在那一轮重新扫描仓外依赖。

## Validation

至少执行以下验证：

- `rg -n '/home/zhengwei/project/python/chromBPNet' AGENTS.md README.md TRACKING.md docs reports scripts tests workflows`
  期望结果：命中应只剩本 design/spec、自觉保留的 historical path 说明，或明确不属于本地 rename 范围的上下文；不应再有把旧路径当“当前 canonical 本地仓入口”的普通引用。
- `rg -n '/home/zhengwei/\\.config/superpowers/worktrees/TransChromBP' AGENTS.md README.md TRACKING.md docs reports scripts tests workflows`
  期望结果：不应出现把已卸载历史 worktree 伪造性改写成新目录名的路径。
- `bash -n scripts/sync_project.sh`
- `.venv/bin/python -m pytest -q tests/test_sync_project_contract.py tests/test_repository_governance_docs.py tests/test_factor_ladder_docs.py`

物理 rename 之后，再在新路径下重跑同样的语法/测试验证，并补：

- `pwd`
- `git rev-parse --show-toplevel`
- `test -L /home/zhengwei/project/python/chromBPNet`

## Risks And Mitigations

### Risk: stale external references still use old path

Mitigation:

- 本轮保留兼容 symlink；
- 不在同一轮移除旧路径入口；
- 把“删除 symlink”单独作为后续小任务。

### Risk: accidental rewrite of historical evidence

Mitigation:

- 对旧 worktree 路径不做机械全局替换；
- 用“current canonical link”与“historical path evidence”两种口径区分处理。

### Risk: current shell or tooling loses repo root during rename

Mitigation:

- 先完成仓内文本治理和验证，再做物理 rename；
- rename 后立即以新路径重新验证 git/top-level/test；
- 保留 symlink 兜底当前仍持有旧 cwd 的会话。

## Rollback

若 rename 后发现不可接受的问题：

1. 删除兼容 symlink；
2. 把目录从 `/home/zhengwei/project/python/TransChromBP` 改回 `/home/zhengwei/project/python/chromBPNet`；
3. 回到旧路径重新验证 git 和测试；
4. 保留本轮文档提交，单独修订 spec 与实现计划。

## Success Criteria

- 本地 canonical 仓物理目录名变成 `TransChromBP`。
- 旧路径 `chromBPNet` 作为 symlink 可访问。
- 仓内 live/canonical 绝对路径默认指向新本地根。
- 历史 worktree 路径没有被伪造性替换成不存在的新路径。
- 治理/文档相关测试通过。
- `master` 在 rename 完成后保持干净并可继续工作。
