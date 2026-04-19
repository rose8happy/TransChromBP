# Dual-Machine Experiment Charter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 把双机实验优先级正式收口成一份权威 charter，并同步 `TRACKING.md` 与高风险 handoff/旧计划文档，避免后续默认切回“先写论文”。

**Architecture:** 先复核 6000/6002 的 live 事实，再新建 charter 作为单一规则源；随后把 `TRACKING.md` 顶部 live 区和最容易误导优先级判断的文档改成“规则看 charter、状态看 live/handoff”的分工。

**Tech Stack:** Markdown docs, `ssh`, `nvidia-smi`, `ps`, `tail`, `grep`, `git diff --check`, `rg`

---

### Task 1: Recheck Live Experiment Facts

**Files:**
- Modify: `TRACKING.md`
- Modify: `reports/handoff/repository_status_handoff_20260409.md`

- [ ] **Step 1: Recheck `6000` active run status**

Run:
```bash
ssh zhoujiazhen@127.0.0.1 -p 6000 \
  'nvidia-smi --query-gpu=index,name,utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits && \
   echo --- && \
   ps -eo pid,etimes,cmd | grep transchrombp.training.train_ddp | grep -v grep && \
   echo --- && \
   tail -n 30 /data1/zhoujiazhen/bylw_atac/TransChromBP/logs/ntv2_teacher_distill_short10_s42_6000_20260409_r2.log'
```
Expected: 确认 `6000` 仍在跑 `ntv2_teacher_distill_short10_s42_6000_20260409_r2`，并拿到最新 epoch / step。

- [ ] **Step 2: Recheck `6002` current idle state and `U-Net-lite` evidence**

Run:
```bash
ssh -i /home/zhengwei/.ssh/codex_6002_ed25519 zhengwei@127.0.0.1 -p 6002 \
  'nvidia-smi --query-gpu=index,name,utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits && \
   echo --- && \
   ps -eo pid,etimes,cmd | grep transchrombp.training.train_ddp | grep -v grep || true && \
   echo --- && \
   ls -1t /home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r*.log 2>/dev/null && \
   echo --- && \
   grep -nE "\\[best\\]|Early stopping|completed|U-Net-lite decoder probe completed" \
     /home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r3.log | tail -n 20'
```
Expected: 确认 `6002` 当前空闲，且 `U-Net-lite` 已经实际跑到 `r3`，不是“尚未起跑”。

- [ ] **Step 3: Capture implementation timestamp**

Run:
```bash
date '+%Y-%m-%d %H:%M:%S %Z'
```
Expected: 记录本轮实施的文档时间锚点。

### Task 2: Add The Charter And Sync Live Tracking

**Files:**
- Create: `docs/plan/dual_machine_experiment_charter_20260409.md`
- Modify: `TRACKING.md`

- [ ] **Step 1: Write the charter as the single rule source**

Write:
```markdown
# 双机实验执行宪章（2026-04-09）

## 角色
- 规则源：本文件
- live 状态：`TRACKING.md`
- 仓库/分支快照：`reports/handoff/repository_status_handoff_20260409.md`
```
Expected: charter 明确权威顺序、白名单/黑名单、`6000/6002` 队列、stop-rule、插队规则。

- [ ] **Step 2: Update `TRACKING.md` top note and live rows**

Write:
```markdown
> 实验规则源统一看 `docs/plan/dual_machine_experiment_charter_20260409.md`；
> 仓库 / worktree / 双机运行快照看 `reports/handoff/repository_status_handoff_20260409.md`。
```
Expected: `TRACKING.md` 不再把“论文收口”写成双机默认下一步，并把 `6002` 行改成 `U-Net-lite r1/r2/r3` 已完成待判读。

### Task 3: De-risk High-Visibility Legacy Docs

**Files:**
- Modify: `reports/handoff/repository_status_handoff_20260409.md`
- Modify: `docs/plan/archive/foundation/post_chatgpt_pro_priority_execution_20260405.md`
- Modify: `docs/plan/archive/governance/project_roadmap_20260330.md`

- [ ] **Step 1: Refresh the handoff summary to match current facts**

Run:
```bash
sed -n '1,220p' reports/handoff/repository_status_handoff_20260409.md
```
Expected: 更新 `6002` 不再写成“U-Net-lite 尚未起跑”，并在建议动作里引用 charter。

- [ ] **Step 2: Add supersession note to the old priority plan**

Write:
```markdown
> 更新（2026-04-09）：本文档只保留 `2026-04-05/06` 当时的 foundation 停表背景；
> 当前双机实验优先级与“是否继续实验而非先写论文”的权威规则，改以
> `docs/plan/dual_machine_experiment_charter_20260409.md` 为准。
```
Expected: 旧计划不再被误读成当前机器队列优先级来源。

- [ ] **Step 3: Add a similar note to the long-term paper roadmap**

Write:
```markdown
> 更新（2026-04-09）：本文档是论文/项目总路线图，不负责当前双机实验调度；
> 只要 charter 白名单队列未收口，就不能用本文把默认动作切回“先写论文”。
```
Expected: `project_roadmap` 保留价值，但不再越权覆盖当前实验排程。

### Task 4: Verify And Commit

**Files:**
- Modify: `docs/plan/dual_machine_experiment_charter_20260409.md`
- Modify: `TRACKING.md`
- Modify: `reports/handoff/repository_status_handoff_20260409.md`
- Modify: `docs/plan/archive/foundation/post_chatgpt_pro_priority_execution_20260405.md`
- Modify: `docs/plan/archive/governance/project_roadmap_20260330.md`
- Create: `docs/superpowers/plans/2026-04-09-dual-machine-experiment-charter-implementation.md`

- [ ] **Step 1: Run doc-format verification**

Run:
```bash
git diff --check -- \
  docs/plan/dual_machine_experiment_charter_20260409.md \
  TRACKING.md \
  reports/handoff/repository_status_handoff_20260409.md \
  docs/plan/archive/foundation/post_chatgpt_pro_priority_execution_20260405.md \
  docs/plan/archive/governance/project_roadmap_20260330.md \
  docs/superpowers/plans/2026-04-09-dual-machine-experiment-charter-implementation.md
```
Expected: 无空格或 patch 格式问题。

- [ ] **Step 2: Confirm the new authority links are visible**

Run:
```bash
rg -n "experiment charter|双机实验执行宪章|规则源|先写论文" \
  docs/plan/dual_machine_experiment_charter_20260409.md \
  TRACKING.md \
  reports/handoff/repository_status_handoff_20260409.md \
  docs/plan/archive/foundation/post_chatgpt_pro_priority_execution_20260405.md \
  docs/plan/archive/governance/project_roadmap_20260330.md
```
Expected: 新 charter 被正确引用，旧文档含 supersession 说明。

- [ ] **Step 3: Commit the documentation implementation**

Run:
```bash
git add \
  docs/plan/dual_machine_experiment_charter_20260409.md \
  TRACKING.md \
  reports/handoff/repository_status_handoff_20260409.md \
  docs/plan/archive/foundation/post_chatgpt_pro_priority_execution_20260405.md \
  docs/plan/archive/governance/project_roadmap_20260330.md \
  docs/superpowers/plans/2026-04-09-dual-machine-experiment-charter-implementation.md
git commit -m "docs: enforce dual machine experiment charter"
```
Expected: 本轮 charter 落地以单独 docs 提交收口。
