# Dual Machine Decoupled Scheduling Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the current “6000 waits for 6002” dual-machine wording with a throughput-first, fully decoupled scheduling policy in the live experiment docs.

**Architecture:** Treat the implementation as a docs-and-state migration, not a code change. First re-check live remote state so the edited text matches reality, then rewrite the charter as the authoritative rules source, rewrite `TRACKING.md` into parallel `6000`/`6002` live rows, and finally trim the handoff report so it remains a snapshot instead of a scheduler.

**Tech Stack:** Markdown, `rg`, `sed`, `git`, SSH to `6000`/`6002`

---

### Task 1: Freeze the live inputs before editing

**Files:**
- Modify: `TRACKING.md`
- Modify: `docs/plan/2026-04-09_dual_machine_experiment_charter.md`
- Modify: `reports/repository_status_handoff_20260409.md`

- [ ] **Step 1: Re-check `6000` GPU/process state**

Run:

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 \
  "nvidia-smi --query-gpu=index,name,utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits && \
   echo '---' && \
   ps -eo pid,etimes,cmd | grep -E 'train_ddp|alphagenome|teacher_v2|ntv2_' | grep -v grep"
```

Expected:

- If `6000` is still idle, both A6000 rows should show `0 util` and the process list should be empty.
- If a new job has already been started manually, carry that exact run name, GPU usage, and timestamp into the later doc edits instead of using stale assumptions.

- [ ] **Step 2: Re-check `6002` live run state**

Run:

```bash
ssh -i /home/zhengwei/.ssh/codex_6002_ed25519 -p 6002 zhengwei@127.0.0.1 \
  "nvidia-smi --query-gpu=index,name,utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits && \
   echo '---' && \
   ps -eo pid,etimes,cmd | grep 'teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r4' | grep -v grep && \
   echo '---LOG---' && \
   tail -n 40 /home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r4.log"
```

Expected:

- If `r4` is still running, keep it as the active `6002` row.
- If `r4` has already finished, carry the actual completion state and best validation evidence into the edited docs.

- [ ] **Step 3: Find all remaining coupling phrases before editing**

Run:

```bash
rg -n "等 `6002|等 6002|先把 `6002|先等 `6002|6002.*晋级候选|6000.*不再.*6000 run|继续等 `6002|唯一下一步" \
  TRACKING.md \
  docs/plan/2026-04-09_dual_machine_experiment_charter.md \
  reports/repository_status_handoff_20260409.md
```

Expected:

- Matches should show the exact rows and charter clauses that still encode the old coupling model.
- Use these hits as the required replacement list in Tasks 2-4.

- [ ] **Step 4: Do not edit yet; record the actual observed time anchors**

Record in a scratch note before touching the docs:

```text
6000_status_time=<absolute CST timestamp from this check>
6000_active_run=<idle or actual run name>
6002_status_time=<absolute CST timestamp from this check>
6002_active_run=<r4 running / r4 complete / other>
```

Expected:

- Later edits use concrete timestamps such as `2026-04-09 23:15 CST`, not vague relative words like “现在” or “刚才”.

### Task 2: Rewrite the charter as a decoupled rules source

**Files:**
- Modify: `docs/plan/2026-04-09_dual_machine_experiment_charter.md`

- [ ] **Step 1: Replace the machine-role section with throughput-first responsibilities**

Replace the current Section 3 logic with wording equivalent to:

```md
### 3.1 `6000 / A6000 x2`

- 默认承担 `full-budget`、正式 gate、长任务、或明显受益于更高吞吐的实验
- 可以跑 foundation，也可以跑 non-foundation
- 不做为了和 `6002` 保持同步而启动的镜像 run

### 3.2 `6002 / RTX 3080`

- 默认承担 `short10 cheap-screen`、确认性 rerun、轻量结构探索、工具/数据小验证
- 可以跑 foundation，也可以跑 non-foundation
- 不承担需要 A6000 吞吐才合理的正式 full-budget 判定
```

- [ ] **Step 2: Replace the current “双机协同门” with decoupling rules**

Insert language equivalent to:

```md
### 5.3 双机解耦门

- 两台机器各自维护 backlog，不等待对方 run 收口才决定自己的下一条
- 一台机器的结果只能提供建议，不能成为另一台机器发车的必需前提
- 只有共享前置物缺失时，才允许出现等待
- 明确禁止“为了同步而做镜像 run”与“共用一个唯一下一步”
```

- [ ] **Step 3: Rewrite the current queue section so `6000` and `6002` are parallel backlogs**

Ensure Section 5 says:

```md
### 5.1 `6000 / A6000 x2`

- `teacher-distill tutorial` 线已正式停表
- 当前不再为这条线分配后续 `6000` run
- `6000` 的默认动作是不等待 `6002`，而是启动自己的独立高价值任务
- 当前首选 backlog：`AlphaGenome matched raw-track slice`

### 5.2 `6002 / RTX 3080`

- `6002` 继续独立维护 `U-Net-lite r4 -> rigor closeout -> 下一个 cheap readout 变体` 这条队列
- 若暂无 ready 的新 readout 变体，再转轻量外部对照或小型工具验证
```

- [ ] **Step 4: Rewrite the “当前默认读法” so it no longer serializes the machines**

Make Section 10 read like:

```md
1. `6000` 不再等待 `6002`；当前应独立推进自己的高价值 backlog，首选 `AlphaGenome matched raw-track slice`
2. `6002` 继续把当前 `U-Net-lite r4` 跑完并完成 rigor closeout
3. 双机并行推进，不互相卡住；一台机器的负结果不自动改写另一台机器的默认下一步
```

- [ ] **Step 5: Verify the rewritten charter no longer contains the old coupling phrases**

Run:

```bash
rg -n "等 `6002|等 6002|先把 `6002|先等 `6002|晋级候选，不能自动|继续等 `6002" \
  docs/plan/2026-04-09_dual_machine_experiment_charter.md
```

Expected:

- No matches for the old “6000 waits for 6002” wording.
- It is acceptable to keep “结果只能提供建议” style wording if it explicitly denies blocking.

### Task 3: Rewrite `TRACKING.md` into parallel live rows

**Files:**
- Modify: `TRACKING.md`

- [ ] **Step 1: Rewrite the top “仓库任务盘点” row so `6000` and `6002` are parallel**

Edit the first live row so the conclusion and next-step text use this structure:

```md
当前 live 实验面按新解耦规则读取为：
- `6000`：<idle / active run / next backlog>
- `6002`：<active run or just-finished run>
```

The “下一步” cell must not say “等 `6002` 收口后再决定 `6000`”.

- [ ] **Step 2: Rewrite the `6000` row as an independent backlog owner**

Replace the current `6000 NT v2 teacher-distill short10` next-step cell with text equivalent to:

```md
`teacher-distill tutorial` 线已停表；`6000` 当前不再等待 `6002`。
下一次 A6000 占卡默认给自己的独立高价值 backlog，首选 `AlphaGenome matched raw-track slice`；
若该任务尚未 ready，则先完成其前置准备，不让 `6000` 空等。
```

- [ ] **Step 3: Keep the `6002` row focused on its own cheap backlog**

The `6002 RTX 3080 单卡实验侧线` row should keep:

```md
若 `r4` 仍落在 `r3` 附近，则正式收成可停表 verdict；
若明显改善，再重新评估是否保留在 3080 shortlist；
若 `v1` 判负，则继续下一个 cheap readout / structure 变体。
```

Do not mention `6000` as a waiting condition.

- [ ] **Step 4: Verify only the intended rows changed**

Run:

```bash
git diff -- TRACKING.md
```

Expected:

- Changes should be limited to the live experiment rows and the charter pointer text.
- Do not churn unrelated archive/resource sections.

### Task 4: Trim the handoff report so it remains a snapshot, not a scheduler

**Files:**
- Modify: `reports/repository_status_handoff_20260409.md`

- [ ] **Step 1: Rewrite the summary paragraph to say the machines are now decoupled**

Update the top section so it states:

```md
- `6000` 已结束 `teacher-distill`，当前按新调度规则独立维护自己的 backlog
- `6002` 继续自己的 `U-Net-lite r4` / cheap-screen 队列
- 本报告只记录现场快照，不再给出“双机必须串行推进”的默认动作
```

- [ ] **Step 2: Replace “当前最值得做的事” with parallel next actions**

Rewrite that section so the numbered list becomes:

```md
1. `6000` 侧：按解耦规则准备并启动独立 backlog 的下一条高价值任务，首选 `AlphaGenome matched raw-track slice`
2. `6002` 侧：继续完成 `U-Net-lite r4` 与 rigor closeout
3. 若要改双机优先级，以 charter 为准；本报告不再承担排程裁决
```

- [ ] **Step 3: Keep old run evidence intact**

Do not remove:

```md
- `teacher-distill r2` 的 run/log/gate evidence
- `U-Net-lite r1/r2/r3/r4` 的历史证据
- worktree / branch snapshot
```

Only remove the scheduling implication that `6000` must wait for `6002`.

- [ ] **Step 4: Verify the handoff still reads as a snapshot**

Run:

```bash
rg -n "本报告只负责快照|不负责裁决默认下一步|独立维护自己的 backlog|不再承担排程裁决" \
  reports/repository_status_handoff_20260409.md
```

Expected:

- The report should explicitly disclaim scheduler authority and point back to the charter.

### Task 5: Final verification and commit

**Files:**
- Modify: `docs/plan/2026-04-09_dual_machine_experiment_charter.md`
- Modify: `TRACKING.md`
- Modify: `reports/repository_status_handoff_20260409.md`

- [ ] **Step 1: Run coupling-regression searches across all three docs**

Run:

```bash
rg -n "等 `6002|等 6002|先把 `6002|继续等 `6002|唯一下一步|自动触发 `6000|自动改写另一台机器" \
  TRACKING.md \
  docs/plan/2026-04-09_dual_machine_experiment_charter.md \
  reports/repository_status_handoff_20260409.md
```

Expected:

- No remaining matches for the old serialized scheduling model.

- [ ] **Step 2: Run formatting checks**

Run:

```bash
git diff --check -- \
  TRACKING.md \
  docs/plan/2026-04-09_dual_machine_experiment_charter.md \
  reports/repository_status_handoff_20260409.md
```

Expected:

- Exit code `0`

- [ ] **Step 3: Review the final diff before committing**

Run:

```bash
git diff -- \
  TRACKING.md \
  docs/plan/2026-04-09_dual_machine_experiment_charter.md \
  reports/repository_status_handoff_20260409.md
```

Expected:

- The diff should only rewrite scheduling logic and live-state timestamps.
- No unrelated doc churn.

- [ ] **Step 4: Commit the rollout**

Run:

```bash
git add \
  TRACKING.md \
  docs/plan/2026-04-09_dual_machine_experiment_charter.md \
  reports/repository_status_handoff_20260409.md
git commit -m "docs: decouple dual machine experiment scheduling"
```

Expected:

- A single focused docs commit with only the three runtime-rule documents.

- [ ] **Step 5: Confirm clean status**

Run:

```bash
git status --short --branch
git log --oneline -3
```

Expected:

- Clean working tree
- Latest commit is `docs: decouple dual machine experiment scheduling`
