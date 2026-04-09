# Worktree Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 把 `master` 的脏改动重新分流到正确分支，并让每条线都有可追溯的本地提交。

**Architecture:** 先做文件级归属映射，再按 `dual-track`、`structure`、`externalization`、`master docs` 四个桶顺序整理。代码线优先回到对应 worktree，主档案类 docs 留在 `master` 并单独提交。

**Tech Stack:** `git worktree`, `git status`, `git add`, `git commit`, shell 文件拷贝, Markdown docs

---

### Task 1: Build The Ownership Map

**Files:**
- Modify: `TRACKING.md`
- Create: `reports/repository_status_handoff_20260409.md`
- Create: `docs/superpowers/specs/2026-04-09-worktree-cleanup-design.md`
- Create: `docs/superpowers/plans/2026-04-09-worktree-cleanup.md`

- [ ] **Step 1: Capture dirty file lists for `master` and each worktree**

Run:
```bash
git status --porcelain=v1
git -C /home/zhengwei/.config/superpowers/worktrees/chromBPNet/dual-track-20260409 status --porcelain=v1
git -C /home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-structure status --porcelain=v1
git -C /home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-chrombpnet-externalization status --porcelain=v1
```
Expected: 得到每条线的 dirty 清单。

- [ ] **Step 2: Compare overlapping files between `master` and target worktrees**

Run:
```bash
sha256sum vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py \
  /home/zhengwei/.config/superpowers/worktrees/chromBPNet/dual-track-20260409/vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py \
  /home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-structure/vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py
```
Expected: 能判断 `master` 当前内容与哪条线一致。

- [ ] **Step 3: Record the agreed ownership rules**

Rules:
```text
master -> repo docs / paper / tracking / archive
dual-track-20260409 -> teacher-distill + U-Net-lite
autonomy/20260406-structure -> foundation cache contract
autonomy/20260406-chrombpnet-externalization -> official chrombpnet externalization
```

- [ ] **Step 4: Keep the latest live summary on `master`**

Run:
```bash
sed -n '1,80p' TRACKING.md
sed -n '1,220p' reports/repository_status_handoff_20260409.md
```
Expected: `master` 顶部状态与实际现场一致。

### Task 2: Rehome `dual-track` Code

**Files:**
- Modify: `vendor/transchrombp/transchrombp/models/transchrombp.py`
- Modify: `vendor/transchrombp/transchrombp/models/profile_decoder.py`
- Modify: `tests/test_multiscale_decoder_probe.py`
- Modify: `/home/zhengwei/.config/superpowers/worktrees/chromBPNet/dual-track-20260409/...`

- [ ] **Step 1: Copy `master` versions of dual-track-owned files into the dual-track worktree**

Run:
```bash
cp vendor/transchrombp/transchrombp/models/transchrombp.py \
  /home/zhengwei/.config/superpowers/worktrees/chromBPNet/dual-track-20260409/vendor/transchrombp/transchrombp/models/transchrombp.py
cp vendor/transchrombp/transchrombp/models/profile_decoder.py \
  /home/zhengwei/.config/superpowers/worktrees/chromBPNet/dual-track-20260409/vendor/transchrombp/transchrombp/models/profile_decoder.py
cp tests/test_multiscale_decoder_probe.py \
  /home/zhengwei/.config/superpowers/worktrees/chromBPNet/dual-track-20260409/tests/test_multiscale_decoder_probe.py
```
Expected: `dual-track` 拿到 `master` 上最新的该线代码。

- [ ] **Step 2: Commit the dual-track snapshot**

Run:
```bash
git -C /home/zhengwei/.config/superpowers/worktrees/chromBPNet/dual-track-20260409 add .
git -C /home/zhengwei/.config/superpowers/worktrees/chromBPNet/dual-track-20260409 commit -m "wip: snapshot dual track pivot"
```
Expected: `dual-track-20260409` 得到一个本地 `wip` 提交。

- [ ] **Step 3: Restore dual-track-owned files on `master` to HEAD**

Run:
```bash
git restore vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py \
  vendor/transchrombp/transchrombp/models/foundation_adapter.py \
  vendor/transchrombp/transchrombp/models/transchrombp.py \
  vendor/transchrombp/transchrombp/scripts/build_foundation_cache.py \
  vendor/transchrombp/transchrombp/scripts/run_ntv2_residual_short10.sh \
  vendor/transchrombp/transchrombp/training/train_ddp.py
git restore --staged --worktree tests/test_multiscale_decoder_probe.py 2>/dev/null || true
rm -f tests/test_multiscale_decoder_probe.py \
  tests/test_ntv2_bins16_residual_head.py \
  vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_msdec_v1.yaml \
  vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_msdec_v1_narrow.yaml \
  vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_msdec_v1_s2.yaml \
  vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_ntv2_bins16_residual.yaml \
  vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_skipprobe_v1.yaml \
  vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_skipprobe_v1_wide.yaml \
  vendor/transchrombp/transchrombp/configs/train/train_tutorial_foundation_short10_ntv2_bins16.yaml \
  vendor/transchrombp/transchrombp/models/profile_decoder.py \
  vendor/transchrombp/transchrombp/scripts/run_multiscale_decoder_probe.sh \
  vendor/transchrombp/transchrombp/scripts/run_short10_no_foundation_control.sh
```
Expected: `master` 不再混着 dual-track 代码。

### Task 3: Commit `structure` And `externalization`

**Files:**
- Modify: `/home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-structure/...`
- Modify: `/home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-chrombpnet-externalization/...`

- [ ] **Step 1: Commit the structure worktree**

Run:
```bash
git -C /home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-structure add .
git -C /home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-structure commit -m "wip: snapshot foundation cache contract"
```
Expected: `autonomy/20260406-structure` clean。

- [ ] **Step 2: Copy externalization-owned files from `master` into the externalization worktree**

Run:
```bash
cp scripts/run_remote_chrombpnet_dataset_prep.sh \
  /home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-chrombpnet-externalization/scripts/run_remote_chrombpnet_dataset_prep.sh
cp scripts/start_6000_chrombpnet_dataset_prep.sh \
  /home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-chrombpnet-externalization/scripts/start_6000_chrombpnet_dataset_prep.sh
cp scripts/start_6002_chrombpnet_dataset_prep.sh \
  /home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-chrombpnet-externalization/scripts/start_6002_chrombpnet_dataset_prep.sh
cp tests/test_chrombpnet_official_externalization.sh \
  /home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-chrombpnet-externalization/tests/test_chrombpnet_official_externalization.sh
cp workflows/tutorial/step3_get_background_regions.sh \
  /home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-chrombpnet-externalization/workflows/tutorial/step3_get_background_regions.sh
```
Expected: externalization 分支拿到 `master` 的最新实现。

- [ ] **Step 3: Commit the externalization worktree**

Run:
```bash
git -C /home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-chrombpnet-externalization add .
git -C /home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-chrombpnet-externalization commit -m "wip: snapshot chrombpnet externalization"
```
Expected: externalization 分支 clean。

- [ ] **Step 4: Restore externalization-owned files on `master` to HEAD**

Run:
```bash
git restore scripts/run_remote_chrombpnet_dataset_prep.sh \
  scripts/start_6000_chrombpnet_dataset_prep.sh \
  scripts/start_6002_chrombpnet_dataset_prep.sh \
  tests/test_chrombpnet_official_externalization.sh \
  workflows/tutorial/step3_get_background_regions.sh
```
Expected: `master` 不再混着 externalization 工具链改动。

### Task 4: Commit The Archive Docs On `master`

**Files:**
- Modify: `AGENTS.md`
- Modify: `DEVELOPMENT.md`
- Modify: `TRACKING.md`
- Modify: `TRACKING_archive.md`
- Modify: `docs/plan/*.md`
- Modify: `reports/*.md`
- Modify: `reports/*.tex`

- [ ] **Step 1: Review remaining `master` dirty files**

Run:
```bash
git status --short
```
Expected: 剩余内容应主要是 docs / paper / reports / tracking。

- [ ] **Step 2: Commit the archive/docs snapshot**

Run:
```bash
git add AGENTS.md DEVELOPMENT.md TRACKING.md TRACKING_archive.md docs reports
git commit -m "docs: snapshot repository handoff and archive updates"
```
Expected: `master` 的 docs 归档得到单独提交。

- [ ] **Step 3: Verify final boundaries**

Run:
```bash
git status --short --branch
git -C /home/zhengwei/.config/superpowers/worktrees/chromBPNet/dual-track-20260409 status --short --branch
git -C /home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-structure status --short --branch
git -C /home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-chrombpnet-externalization status --short --branch
git worktree list
```
Expected: 每条线都只剩自己的内容，且大部分工作树 clean。
