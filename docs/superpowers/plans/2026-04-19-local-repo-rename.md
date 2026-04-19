# Local Repository Rename Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rename the local canonical repository root from `/home/zhengwei/project/python/chromBPNet` to `/home/zhengwei/project/python/TransChromBP`, rewrite repository-local absolute links to the new root, preserve historical worktree evidence correctly, and leave a compatibility symlink at the old path.

**Architecture:** Split the work into two phases. First, in an implementation branch/worktree, update tests, canonical docs, manifests, and historical report links so that repository-local absolute paths point at the new root while historical worktree paths stay explicit and non-fabricated. Second, after those text changes are merged and no linked implementation worktree remains, rename the canonical working tree on disk, add the compatibility symlink, and re-run the same governance/test checks from the new root.

**Tech Stack:** git, git worktrees, shell (`rg`, `perl`, `mv`, `ln -s`), Markdown/LaTeX docs, CSV manifests, pytest doc-contract tests.

---

## File Structure / Responsibility Map

- `AGENTS.md`
  - canonical repository instructions; must declare the new local canonical root.
- `README.md`
  - top-level operator-facing entrypoint; must stop pointing to the old local root.
- `TRACKING.md`
  - live status; must describe the new local canonical path and symlink transition accurately.
- `docs/env/repository_governance.md`
  - canonical topology doc; must change the local canonical root row to `TransChromBP`.
- `docs/experiments/registry.md`
  - canonical family/workstream index; must not imply unmounted worktree paths are still current local entrypoints.
- `docs/experiments/runs.csv`
  - run-level manifest; archived local worktree path fields that no longer represent accessible live paths should be normalized during the rename cleanup.
- `tests/test_repository_governance_docs.py`
  - doc-contract tests for canonical path, archive closeout tags, and governance invariants.
- `reports/`, `docs/plan/`, `docs/superpowers/`, `scripts/`, `workflows/`
  - historical and operator docs containing repository-local absolute links that must be rewritten to the new local root when they refer to the canonical repo, but must not fabricate renamed historical worktree paths.
- `reports/assets/local_repo_rename_20260419/local_repo_path_candidates.txt`
  - generated inventory of files that still mention the old local canonical root and therefore need review/rewrite.
- `reports/local_repo_rename_closeout_20260419.md`
  - reusable closeout describing the rename, symlink, verification, and any intentionally preserved historical exceptions.

### Task 0: Create The Implementation Worktree And Freeze The Baseline

**Files:**
- Modify: no repository files
- Operate on filesystem paths:
  - `/home/zhengwei/.config/superpowers/worktrees/chromBPNet/local-repo-rename-20260419`

- [ ] **Step 1: Create the dedicated implementation worktree**

Run:

```bash
git worktree add /home/zhengwei/.config/superpowers/worktrees/chromBPNet/local-repo-rename-20260419 -b local-repo-rename-20260419
```

Expected: a clean implementation branch/worktree exists at `/home/zhengwei/.config/superpowers/worktrees/chromBPNet/local-repo-rename-20260419`.

- [ ] **Step 2: Verify the baseline from inside the implementation worktree**

Run:

```bash
cd /home/zhengwei/.config/superpowers/worktrees/chromBPNet/local-repo-rename-20260419
bash -n scripts/sync_project.sh
/home/zhengwei/project/python/chromBPNet/.venv/bin/python -m pytest -q tests/test_sync_project_contract.py tests/test_repository_governance_docs.py tests/test_factor_ladder_docs.py
```

Expected: PASS. If this baseline fails, stop and fix the baseline before continuing, or explicitly re-plan around the failure.

- [ ] **Step 3: Commit nothing; keep the baseline worktree clean**

Run:

```bash
git status --short --branch
```

Expected: `## local-repo-rename-20260419` with no modified files before Task 1 begins.

### Task 1: Lock The Rename Contract In Tests And Build A Rewrite Inventory

**Files:**
- Modify: `tests/test_repository_governance_docs.py`
- Create: `reports/assets/local_repo_rename_20260419/local_repo_path_candidates.txt`
- Create: `reports/assets/local_repo_rename_20260419/local_repo_path_residuals_after_rewrite.txt`

- [ ] **Step 1: Add failing doc-contract assertions for the new local canonical root**

Add assertions in `tests/test_repository_governance_docs.py` for all of the following:

```python
def test_governance_doc_uses_transchrombp_local_root() -> None:
    text = read_text("docs/env/repository_governance.md")
    assert "/home/zhengwei/project/python/TransChromBP" in text
    assert "| `/home/zhengwei/project/python/chromBPNet` |" not in text


def test_no_fake_transchrombp_worktree_paths_exist() -> None:
    tree = subprocess.run(
        ["rg", "-n", "/home/zhengwei/\\.config/superpowers/worktrees/TransChromBP", "AGENTS.md", "README.md", "TRACKING.md", "docs", "reports", "scripts", "tests", "workflows"],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert tree.stdout.strip() == ""
```

- [ ] **Step 2: Run the focused governance tests and verify they fail**

Run:

```bash
/home/zhengwei/project/python/chromBPNet/.venv/bin/python -m pytest -q tests/test_repository_governance_docs.py -k 'transchrombp_local_root or fake_transchrombp_worktree'
```

Expected: FAIL, because the canonical docs still reference `/home/zhengwei/project/python/chromBPNet` as the active local root and the new assertions are not yet satisfied.

- [ ] **Step 3: Build a concrete candidate inventory for repository-local absolute path rewrites**

Run:

```bash
mkdir -p reports/assets/local_repo_rename_20260419
rg -l '/home/zhengwei/project/python/chromBPNet' \
  AGENTS.md README.md TRACKING.md docs reports scripts tests workflows \
  > reports/assets/local_repo_rename_20260419/local_repo_path_candidates.txt
```

Expected: `reports/assets/local_repo_rename_20260419/local_repo_path_candidates.txt` contains the exact file list to review in later tasks.

- [ ] **Step 4: Commit the failing-test baseline and candidate inventory**

Run:

```bash
git add tests/test_repository_governance_docs.py reports/assets/local_repo_rename_20260419/local_repo_path_candidates.txt
git commit -m "test: lock local repo rename contract"
```

Expected: a checkpoint commit that freezes the rename scope and failing contract before any doc rewrite starts.

### Task 2: Update Canonical Governance Docs And Run Manifests

**Files:**
- Modify: `AGENTS.md`
- Modify: `README.md`
- Modify: `TRACKING.md`
- Modify: `docs/env/repository_governance.md`
- Modify: `docs/experiments/registry.md`
- Modify: `docs/experiments/runs.csv`
- Modify: `tests/test_repository_governance_docs.py`

- [ ] **Step 1: Rewrite the active local canonical root in governance/live docs**

Apply the local-root replacement and wording adjustments in:

```text
AGENTS.md
README.md
TRACKING.md
docs/env/repository_governance.md
```

Required outcomes:
- the active local canonical root becomes `/home/zhengwei/project/python/TransChromBP`
- any text that previously said “current physical directory stays unchanged” is rewritten to the new reality
- any mention of the old local path is either removed or recast as the compatibility symlink / rollback path

- [ ] **Step 2: Normalize canonical manifests so they do not depend on unmounted local worktree paths**

Edit:

```text
docs/experiments/registry.md
docs/experiments/runs.csv
```

Required outcomes:
- `registry.md` continues to use `n/a` or branch/tag recovery for unmounted archival worktrees
- `runs.csv` stops presenting removed local worktree paths as current accessible local entrypoints
- archived rows that still need provenance should keep branch/tag/commit context in notes rather than obsolete filesystem paths

- [ ] **Step 3: Tighten the tests to match the canonical-doc outcomes**

Extend `tests/test_repository_governance_docs.py` so it also checks:

```python
def test_tracking_and_agents_reference_transchrombp_local_root() -> None:
    for path in ("AGENTS.md", "TRACKING.md", "README.md"):
        text = read_text(path)
        assert "/home/zhengwei/project/python/TransChromBP" in text


def test_runs_manifest_does_not_require_removed_local_worktrees() -> None:
    runs = read_text("docs/experiments/runs.csv")
    assert "/home/zhengwei/.config/superpowers/worktrees/chromBPNet/dual-track-20260409" not in runs
    assert "/home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-structure" not in runs
    assert "/home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-chrombpnet-externalization" not in runs
```

- [ ] **Step 4: Run the focused governance tests and make sure they pass**

Run:

```bash
/home/zhengwei/project/python/chromBPNet/.venv/bin/python -m pytest -q tests/test_repository_governance_docs.py -k 'transchrombp_local_root or fake_transchrombp_worktree or removed_local_worktrees'
```

Expected: PASS.

- [ ] **Step 5: Commit the canonical-doc and manifest rewrite**

Run:

```bash
git add AGENTS.md README.md TRACKING.md docs/env/repository_governance.md docs/experiments/registry.md docs/experiments/runs.csv tests/test_repository_governance_docs.py
git commit -m "docs: rename local canonical root references"
```

Expected: a clean checkpoint where live governance and manifests already speak the new local-root language.

### Task 3: Bulk Rewrite Repository-Local Absolute Links Across Historical Docs

**Files:**
- Modify: `reports/assets/local_repo_rename_20260419/local_repo_path_candidates.txt`
- Modify: files enumerated in `reports/assets/local_repo_rename_20260419/local_repo_path_candidates.txt`
- Modify: `reports/assets/local_repo_rename_20260419/local_repo_path_residuals_after_rewrite.txt`

- [ ] **Step 1: Apply the bulk local-root rewrite for repository-local absolute links**

Run:

```bash
rg -l -0 '/home/zhengwei/project/python/chromBPNet' \
  AGENTS.md README.md TRACKING.md docs reports scripts tests workflows \
| xargs -0 perl -0pi -e 's#/home/zhengwei/project/python/chromBPNet#/home/zhengwei/project/python/TransChromBP#g'
```

Expected: repository-local absolute links in Markdown, LaTeX, shell examples, and tests now point at the new local root.

- [ ] **Step 2: Manually revert any accidental rewrite that changed historical evidence instead of current navigation**

Review the diff for files that intentionally need old historical wording, especially:

```text
docs/superpowers/specs/2026-04-19-local-repo-rename-design.md
docs/superpowers/plans/2026-04-19-local-repo-rename.md
```

Required outcomes:
- keep old local path where it is part of the symlink, rollback, or historical comparison narrative
- do not touch remote official root `/data1/zhoujiazhen/bylw_atac/chromBPNet`
- do not invent `/home/zhengwei/.config/superpowers/worktrees/TransChromBP/...`

- [ ] **Step 3: Capture the residual old-path hits after the bulk rewrite**

Run:

```bash
rg -n '/home/zhengwei/project/python/chromBPNet' \
  AGENTS.md README.md TRACKING.md docs reports scripts tests workflows \
  > reports/assets/local_repo_rename_20260419/local_repo_path_residuals_after_rewrite.txt
```

Expected: the residual list contains only known exceptions that are explicitly justified by the spec/plan/closeout or compatibility-symlink wording.

- [ ] **Step 4: Commit the bulk historical-link rewrite**

Run:

```bash
git add reports/assets/local_repo_rename_20260419/local_repo_path_candidates.txt \
        reports/assets/local_repo_rename_20260419/local_repo_path_residuals_after_rewrite.txt \
        AGENTS.md README.md TRACKING.md docs reports scripts tests workflows
git commit -m "docs: rewrite local repo absolute links"
```

Expected: a single audit-friendly commit containing the bulk path rewrite and its residual inventory.

### Task 4: Verify Textual Invariants And Write Closeout Before Physical Rename

**Files:**
- Create: `reports/local_repo_rename_closeout_20260419.md`
- Modify: `TRACKING.md`
- Modify: `tests/test_repository_governance_docs.py`

- [ ] **Step 1: Run the full governance/doc verification suite before any filesystem rename**

Run:

```bash
bash -n scripts/sync_project.sh
/home/zhengwei/project/python/chromBPNet/.venv/bin/python -m pytest -q tests/test_sync_project_contract.py tests/test_repository_governance_docs.py tests/test_factor_ladder_docs.py
rg -n '/home/zhengwei/\.config/superpowers/worktrees/TransChromBP' AGENTS.md README.md TRACKING.md docs reports scripts tests workflows
```

Expected:
- `bash -n` PASS
- pytest PASS
- `rg` returns no fake renamed historical worktree paths

- [ ] **Step 2: Write the pre-rename closeout report**

Create `reports/local_repo_rename_closeout_20260419.md` with all of the following:

```text
- old local canonical root
- new local canonical root
- compatibility symlink plan
- verification commands and results before physical rename
- residual old-path exceptions and why they remain
- statement that remote official/runtime roots were intentionally left unchanged
```

- [ ] **Step 3: Update TRACKING to point at the closeout**

Add a concise line in `TRACKING.md` stating that:
- repository-local absolute links have been rewritten
- physical rename is next
- `reports/local_repo_rename_closeout_20260419.md` is the canonical audit trail

- [ ] **Step 4: Commit the verified pre-rename closeout state**

Run:

```bash
git add TRACKING.md reports/local_repo_rename_closeout_20260419.md tests/test_repository_governance_docs.py
git commit -m "docs: record local repo rename closeout"
```

Expected: a final pre-rename commit that can be merged to `master` before the on-disk rename.

### Task 5: Collapse Back To The Canonical Working Tree Before Renaming The Directory

**Files:**
- Modify: no repository files
- Operate on: the temporary implementation branch/worktree created for this plan

- [ ] **Step 1: Merge the implementation branch back to `master`**

Run from the canonical repo root:

```bash
git checkout master
git merge --ff-only local-repo-rename-20260419
```

Expected: `master` now contains all text/test/closeout changes required for the rename.

- [ ] **Step 2: Remove the temporary implementation worktree and branch**

Run:

```bash
git worktree remove /home/zhengwei/.config/superpowers/worktrees/chromBPNet/local-repo-rename-20260419
git branch -d local-repo-rename-20260419
```

Expected: no linked implementation worktree remains that would hold stale gitdir pointers when the main repository directory moves.

- [ ] **Step 3: Verify only the canonical working tree is mounted**

Run:

```bash
git worktree list
```

Expected: only `/home/zhengwei/project/python/chromBPNet` remains mounted before the physical rename.

### Task 6: Perform The Physical Rename And Add The Compatibility Symlink

**Files:**
- Modify: no repository files
- Operate on filesystem paths:
  - `/home/zhengwei/project/python/chromBPNet`
  - `/home/zhengwei/project/python/TransChromBP`

- [ ] **Step 1: Rename the repository directory in the parent folder**

Run:

```bash
cd /home/zhengwei/project/python
mv chromBPNet TransChromBP
```

Expected: the repository now physically lives at `/home/zhengwei/project/python/TransChromBP`.

- [ ] **Step 2: Re-create the old path as a compatibility symlink**

Run:

```bash
cd /home/zhengwei/project/python
ln -s TransChromBP chromBPNet
```

Expected: `/home/zhengwei/project/python/chromBPNet` resolves to the new canonical root and external stale callers continue to work.

- [ ] **Step 3: Enter the new root and verify git sees the moved repository correctly**

Run:

```bash
cd /home/zhengwei/project/python/TransChromBP
pwd
git rev-parse --show-toplevel
git status --short --branch
test -L /home/zhengwei/project/python/chromBPNet
```

Expected:
- `pwd` prints `/home/zhengwei/project/python/TransChromBP`
- `git rev-parse --show-toplevel` prints `/home/zhengwei/project/python/TransChromBP`
- `git status` remains clean on `master`
- `test -L ...` succeeds

### Task 7: Re-Verify From The New Root And Finalize The Rename

**Files:**
- Modify: `reports/local_repo_rename_closeout_20260419.md`
- Modify: `TRACKING.md`

- [ ] **Step 1: Re-run the full verification suite from the new root**

Run:

```bash
cd /home/zhengwei/project/python/TransChromBP
bash -n scripts/sync_project.sh
/home/zhengwei/project/python/TransChromBP/.venv/bin/python -m pytest -q tests/test_sync_project_contract.py tests/test_repository_governance_docs.py tests/test_factor_ladder_docs.py
rg -n '/home/zhengwei/project/python/chromBPNet' AGENTS.md README.md TRACKING.md docs reports scripts tests workflows
```

Expected:
- syntax check PASS
- pytest PASS
- residual `chromBPNet` hits are limited to intentional symlink/rollback/historical exceptions documented in the closeout

- [ ] **Step 2: Update the closeout with post-rename evidence**

Append to `reports/local_repo_rename_closeout_20260419.md`:

```text
- final `pwd`
- final `git rev-parse --show-toplevel`
- symlink verification
- post-rename test results
- final residual old-path list and why each hit is acceptable
```

- [ ] **Step 3: Refresh TRACKING to say the physical rename is complete**

Update `TRACKING.md` so the repository-governance cleanup row now states:
- the canonical physical root is already `/home/zhengwei/project/python/TransChromBP`
- the old `chromBPNet` path is a compatibility symlink
- the next small follow-up is only deciding when to remove that symlink

- [ ] **Step 4: Commit the post-rename finalization**

Run:

```bash
git add TRACKING.md reports/local_repo_rename_closeout_20260419.md
git commit -m "chore: finalize local repo rename"
```

Expected: a final commit that documents the completed rename and its verification evidence.
