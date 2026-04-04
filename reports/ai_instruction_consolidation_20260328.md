# AI 指令体系最小整合报告（2026-03-28）

## 背景

仓库内同时存在 `AGENTS.md`、`CLAUDE.md`、`GEMINI.md` 三个 AI 指令入口。此前三者都承载了大量仓库级规则，导致以下问题：

- 同一规则在多处重复，后续修改时容易漂移
- 某些路径和目录描述与当前仓库现实不完全一致
- 代理特有提醒和仓库通用规则混在一起，维护成本偏高

同轮审查还参考了 `repo_ai_instruction_template_exec_v7.md` 与 `repo_ai_instruction_template_user_guide_exec_v7.md`。结论是：模板更适合作为“规范收敛的检查表”，不适合作为当前仓库的目录重构蓝图。

## 当前仓库采用的收敛策略

本轮采用“最小整合”：

- `AGENTS.md` 作为唯一权威源，负责稳定、跨代理通用的仓库规则
- `CLAUDE.md` 和 `GEMINI.md` 收敛为薄包装，只保留代理特有提醒
- 保留现有文档分层，不引入模板默认的 `docs/tasks/active`、`docs/tasks/archive`、`docs/decisions` 等目录

继续沿用的仓库记忆分工：

- `AGENTS.md`：稳定规则、目录图、命令和文档协议
- `TRACKING.md`：live 状态与下一步
- `TRACKING_archive.md`：阶段性完成事项归档
- `reports/`：可复用的完整分析报告
- `docs/plan/`：中大型计划与执行方案

## 已落地修改

### 1. 主规范收敛到 `AGENTS.md`

新增并明确了以下规则：

- `AGENTS.md` 是 canonical repository instruction file
- `CLAUDE.md` 和 `GEMINI.md` 只做代理适配层，不再重复项目全量规则
- 仓库结构、命令或文档协议变化时，优先更新 `AGENTS.md`
- live 状态、完整报告、计划文档分别放在 `TRACKING.md`、`reports/`、`docs/plan/`

同时修正了结构描述中的一个偏差：当前 `docs/` 真实存在的是 `plan/`、`research/`、`env/`、`learning/`，并无常驻 `docs/report/` 目录，因此不再把 `docs/report/` 写成既定结构。

### 2. `CLAUDE.md` / `GEMINI.md` 改为薄包装

两份包装文件现在只保留：

- “先看 `AGENTS.md`，再看 `TRACKING.md`” 的使用顺序
- 本地编辑 / 远端训练的工作边界提醒
- 远端服务器快速入口
- 中文回复要求

这保证代理仍能拿到自己的最小必要上下文，但不会再和主规范发生大面积重复。

## 为什么不照搬模板默认目录

本仓库已经形成可工作的结构：

- 计划主要放 `docs/plan/`
- 长文分析主要放根目录 `reports/`
- 运行态状态集中在 `TRACKING.md`

若直接迁移到模板默认的 `docs/plans/`、`docs/reports/`、`docs/tasks/active/` 等结构，会产生不必要的文档搬迁与链接修复成本。本轮因此只吸收模板的“职责边界”和“主规范优先”思想，不采纳其默认目录蓝图。

## 后续维护规则

- 以后新增仓库级规则时，先改 `AGENTS.md`
- 只有当某个代理存在额外限制或操作差异时，才改对应包装文件
- 如果包装文件再次长成“第二份仓库总规则”，应再次收敛
- 如后续实验线继续增多，再评估是否需要额外引入 task ID 或更细的任务目录
