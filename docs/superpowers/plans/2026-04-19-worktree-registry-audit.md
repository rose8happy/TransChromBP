# Worktree Registry Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Audit mounted worktrees against registry, tags, and closeout evidence; remove redundant mounts and document any preserved snapshots.

**Architecture:** Treat `master` as the canonical archive, preserve any still-valuable dirty worktree state with explicit snapshot commits/tags, and keep `docs/experiments/registry.md` aligned with the actual mounted set. Avoid destructive cleanup unless the worktree is either already fully merged or first snapshotted.

**Tech Stack:** git worktrees, annotated git tags, Markdown registry/report docs, pytest doc-contract tests.

---

### Task 1: Audit Current Worktrees And Tags

**Files:**
- Modify: `docs/experiments/registry.md`
- Modify: `TRACKING.md`
- Create: `reports/governance/worktree_registry_audit_20260419.md`
- Test: `tests/test_factor_ladder_docs.py`

- [ ] **Step 1: Write failing test for the desired post-audit loss-balance registry state**

Expect `loss_balance_curriculum` to keep its final canonical report on `master`, but to reference an explicit snapshot tag for the preserved worktree state instead of relying on a live mounted dirty worktree.

- [ ] **Step 2: Verify the test fails**

Run: `pytest -q tests/test_factor_ladder_docs.py -k loss_balance`
Expected: FAIL because the registry row does not yet mention the preserved snapshot / closeout tags.

- [ ] **Step 3: Record the audit result in a reusable report**

Create `reports/governance/worktree_registry_audit_20260419.md` summarizing:
- mounted worktrees that are already clean and properly tagged
- redundant worktrees that can be removed immediately
- dirty worktrees that must be snapshotted before removal

- [ ] **Step 4: Update registry / tracking to match the audit result**

Keep `master` as canonical branch where appropriate, but add snapshot / closeout tag references and explicit notes when a worktree was removed after preservation.

### Task 2: Preserve Dirty Loss-Balance Worktree

**Files:**
- Preserve snapshot on branch: `loss-balance-20260417`
- Modify: `docs/experiments/registry.md`
- Modify: `reports/governance/worktree_registry_audit_20260419.md`

- [ ] **Step 1: Snapshot the dirty `loss-balance-20260417` worktree**

Run:

```bash
git -C /home/zhengwei/.config/superpowers/worktrees/chromBPNet/loss-balance-20260417 add -A
git -C /home/zhengwei/.config/superpowers/worktrees/chromBPNet/loss-balance-20260417 commit -m "wip: snapshot loss balance worktree"
```

Expected: a new archival commit preserving the exact local WIP state.

- [ ] **Step 2: Create annotated tags for the preserved milestone**

Run:

```bash
git tag -a snapshot/loss-balance-20260417/20260419 <snapshot-commit> -m "<message>"
git tag -a closeout/loss_balance_curriculum/20260418 <canonical-commit> -m "<message>"
```

Expected: the preserved code snapshot and the canonical closeout state both become addressable by tag.

- [ ] **Step 3: Remove the now-preserved mounted worktree**

Run:

```bash
git worktree remove /home/zhengwei/.config/superpowers/worktrees/chromBPNet/loss-balance-20260417
```

Expected: no mounted dirty worktree remains for that family.

### Task 3: Remove Redundant Bundle Worktree

**Files:**
- Modify: `reports/governance/worktree_registry_audit_20260419.md`

- [ ] **Step 1: Verify the bundle files already exist on `master`**

Run:

```bash
git ls-tree --name-only master -- reports/chatgpt_bundle_loss_balance_20260416
```

Expected: the bundle directory is already present on `master`.

- [ ] **Step 2: Remove the redundant worktree and branch**

Run:

```bash
git worktree remove --force /home/zhengwei/.config/superpowers/worktrees/chromBPNet/loss-balance-bundle-20260416
git branch -D loss-balance-bundle-20260416
```

Expected: the orphan bundle mount disappears without losing unique content.

### Task 4: Verify Final State

**Files:**
- Modify: `tests/test_factor_ladder_docs.py`
- Modify: `TRACKING.md`
- Modify: `docs/experiments/registry.md`

- [ ] **Step 1: Run doc-contract tests**

Run:

```bash
/home/zhengwei/project/python/TransChromBP/.venv/bin/python -m pytest -q tests/test_repository_governance_docs.py tests/test_factor_ladder_docs.py
```

Expected: PASS.

- [ ] **Step 2: Re-run sync script syntax check**

Run:

```bash
bash -n scripts/sync_project.sh
```

Expected: PASS.

- [ ] **Step 3: Verify the mounted worktree list**

Run:

```bash
git worktree list
```

Expected: redundant loss-balance bundle mount is gone, preserved loss-balance WIP is either snapshotted and removed or explicitly accounted for.
