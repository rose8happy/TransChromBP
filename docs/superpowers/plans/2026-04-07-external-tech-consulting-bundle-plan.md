# External Tech-Consulting Bundle Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reframe the external handoff bundle so external models focus on technical route consultation, not reviewer-style story cleanup, while keeping the current evidence and stop-rule constraints intact.

**Architecture:** Keep the existing `2026-04-05` bundle as the evidence base, then change the bundle entrypoint, add one focused technical-questions file, and retune the open-questions / delta files so the package clearly asks for next-generation architecture and comparison advice. Keep old shortcut/bias-mismatch material only as downgraded background, not as the main prompt attractor.

**Tech Stack:** Markdown documentation, repository guidance files, ripgrep/sed/git diff verification

---

### Task 1: Rewrite The Bundle Entrypoint

**Files:**
- Modify: `reports/chatgpt_bundle_project_handoff_20260405/00_readme_upload_order_and_prompt.md`
- Verify: `reports/chatgpt_bundle_project_handoff_20260405/00_readme_upload_order_and_prompt.md`

- [ ] **Step 1: Replace the old reviewer-first framing with a technical-consulting-first framing**

Update the intro so it says the package is now used as a technical route consultation bundle, not primarily a reviewer cleanup bundle. The intro must explicitly say the `2026-04-05` base package now requires a later delta file.

```md
当前使用方式已更新到：`2026-04-07`。

基础主包仍然形成于 `2026-04-05`，但当前外发时应强制附带后续 delta 文件，
用于补充 `2026-04-06` 之后的关键结果和运行仓说明。
```

- [ ] **Step 2: Update the first-round upload order**

Change the upload order so the first round includes the new technical-questions file and the post-`2026-04-06` delta file.

```md
### 第一轮：先传 00-05 + 08 + 09

1. `00_readme_upload_order_and_prompt.md`
2. `01_project_state_claims_and_open_questions.md`
3. `02_architecture_code_map_and_best_model.md`
4. `03_experiment_history_and_evidence.md`
5. `04_foundation_attempts_and_current_restart.md`
6. `05_key_metrics_and_runs.csv`
7. `08_post_20260406_delta.md`
8. `09_technical_questions_for_external_models.md`
```

- [ ] **Step 3: Replace the prompt with a technical-consulting prompt**

Rewrite the prompt so it asks for:

- model-improvement routes
- AlphaGenome / U-Net style guidance
- ATAC metric recommendations
- why pretrained genomic FMs did not help yet
- how to compare against AlphaGenome

Also add a hard constraint that the old `Transformer/bias` mismatch suspicion is downgraded background, not the main issue.

- [ ] **Step 4: Verify the entrypoint text**

Run: `sed -n '1,260p' reports/chatgpt_bundle_project_handoff_20260405/00_readme_upload_order_and_prompt.md`

Expected:
- first-round upload list includes `08` and `09`
- prompt is technical-consulting-first
- prompt explicitly deprioritizes the old shortcut/bias suspicion

---

### Task 2: Add A Focused Technical-Questions File

**Files:**
- Create: `reports/chatgpt_bundle_project_handoff_20260405/09_technical_questions_for_external_models.md`
- Verify: `reports/chatgpt_bundle_project_handoff_20260405/09_technical_questions_for_external_models.md`

- [ ] **Step 1: Create the technical-questions document**

Write a new file that explicitly states this round wants the highest-potential route judgment, even if the path would require medium or large refactoring.

Required sections:

- why this round is different from the previous reviewer-oriented round
- the five technical question groups
- what counts as a useful answer
- what not to spend attention on

Include this kind of framing:

```md
本轮默认优先级不是“最小改动、最短周期”，而是：

> 在当前证据约束下，优先判断最有潜力的下一代方向，即使工程代价较大。
```

- [ ] **Step 2: Add the five technical question groups**

Include concrete subsections for:

```md
1. 我们的模型还能怎样提升
2. AlphaGenome / U-Net 风格是否可能提升
3. ATAC 领域常用哪些指标，我们该补哪些
4. 预训练基因组大模型官方如何做下游 ATAC/可及性任务；为什么我们当前接入没效果；接下来怎么做
5. 我们应该怎样与 AlphaGenome 做对比
```

- [ ] **Step 3: Add the “do not get distracted” section**

Add a short section saying the old `Transformer/bias` mismatch suspicion was an earlier misunderstanding and should not become the main analysis focus.

- [ ] **Step 4: Verify the new file**

Run: `sed -n '1,260p' reports/chatgpt_bundle_project_handoff_20260405/09_technical_questions_for_external_models.md`

Expected:
- the file clearly asks for route consultation
- all five question groups are present
- the downgraded old narrative is explicitly deprioritized

---

### Task 3: Retune The State And Delta Files

**Files:**
- Modify: `reports/chatgpt_bundle_project_handoff_20260405/01_project_state_claims_and_open_questions.md`
- Modify: `reports/chatgpt_bundle_project_handoff_20260405/08_post_20260406_delta.md`
- Verify: `reports/chatgpt_bundle_project_handoff_20260405/01_project_state_claims_and_open_questions.md`
- Verify: `reports/chatgpt_bundle_project_handoff_20260405/08_post_20260406_delta.md`

- [ ] **Step 1: Update the “open questions” section in `01`**

Retune the open-questions section so it no longer leads with reviewer cleanup or the old shortcut/bias concern. It should instead emphasize:

- next-generation improvement paths
- AlphaGenome/U-Net evaluation
- metric expansion
- FM downstream integration strategy
- AlphaGenome comparison design

Use wording like:

```md
当前最需要外部模型回答的，不再是“是否还要继续纠缠旧的 shortcut/bias 误会”，
而是“下一代提升路线应该朝哪里走”。
```

- [ ] **Step 2: Retune the delta file**

Keep the current bins16 negative result and stop-rule, but rewrite the landing point so it feeds the new purpose:

```md
这条新增结果的意义不是让我们继续沿当前 measured family 微调，
而是为“下一代不同方向该怎么选”提供更清楚的边界条件。
```

- [ ] **Step 3: Verify both files**

Run:

```bash
sed -n '1,240p' reports/chatgpt_bundle_project_handoff_20260405/01_project_state_claims_and_open_questions.md
sed -n '1,240p' reports/chatgpt_bundle_project_handoff_20260405/08_post_20260406_delta.md
```

Expected:
- `01` now asks technical route questions first
- `08` still preserves the stop-rule facts
- neither file foregrounds the old `Transformer/bias` suspicion

---

### Task 4: Audit For Leftover Narrative Drift And Sync Tracking

**Files:**
- Inspect: `reports/chatgpt_bundle_project_handoff_20260405/02_architecture_code_map_and_best_model.md`
- Inspect: `reports/chatgpt_bundle_project_handoff_20260405/03_experiment_history_and_evidence.md`
- Inspect: `reports/chatgpt_bundle_project_handoff_20260405/04_foundation_attempts_and_current_restart.md`
- Modify if needed: the specific file(s) above that still over-foreground the old narrative
- Modify: `TRACKING.md`

- [ ] **Step 1: Search for leftover old-narrative hotspots**

Run:

```bash
rg -n "shortcut|bias.*连接|Transformer-specific|stop-gradient|误会" reports/chatgpt_bundle_project_handoff_20260405
```

Expected:
- identify whether `02-04` still over-foreground the old concern

- [ ] **Step 2: Apply only minimal denoising edits if needed**

If a file still spends disproportionate space on the old issue, compress it to one short historical-background statement. Do not rewrite the full file unless the framing is obviously harmful.

Example replacement style:

```md
这条怀疑曾是中期排查重点，但当前已被降级为历史背景，不再是外部咨询的主问题。
```

- [ ] **Step 3: Update `TRACKING.md`**

Rewrite the ChatGPT/外发 row so it records:

- the bundle now targets technical route consultation
- `00/01/08/09` are the active bundle control points
- the old `Transformer/bias` mismatch narrative has been deprioritized

- [ ] **Step 4: Final verification**

Run:

```bash
sed -n '24,28p' TRACKING.md
sed -n '1,260p' reports/chatgpt_bundle_project_handoff_20260405/00_readme_upload_order_and_prompt.md
sed -n '1,260p' reports/chatgpt_bundle_project_handoff_20260405/09_technical_questions_for_external_models.md
git diff -- AGENTS.md TRACKING.md reports/chatgpt_bundle_project_handoff_20260405/00_readme_upload_order_and_prompt.md reports/chatgpt_bundle_project_handoff_20260405/01_project_state_claims_and_open_questions.md reports/chatgpt_bundle_project_handoff_20260405/08_post_20260406_delta.md reports/chatgpt_bundle_project_handoff_20260405/09_technical_questions_for_external_models.md
```

Expected:
- bundle entrypoint is now tech-consulting-first
- `09` exists and is focused
- `TRACKING` reflects the new purpose
- no unexpected file edits appear in the diff
