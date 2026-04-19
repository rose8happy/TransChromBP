# A6000 AlphaGenome Matched Slice v1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Launch and close out the `6000`-side `AlphaGenome matched raw-track slice v1` pilot without waiting for `6002`, while keeping `TRACKING.md`, the dual-machine charter, and a reusable closeout report in sync.

**Architecture:** Treat this as a short remote-ops plus docs workflow, not a model-code feature. First freeze the real remote prerequisites and copy the minimal pilot assets into the `6000` runtime repo, then launch one background pilot that runs the AlphaGenome fetch and merge step end-to-end, update live docs only after launch is verified, and finally write a single closeout report that classifies the run as `pass`, `fail`, or `launch-blocked`.

**Tech Stack:** Markdown, `git`, `ssh`, `scp`, `nohup`, Python 3.11 in `/data1/zhoujiazhen/bylw_atac/.mamba/envs/alphagenome`

---

### Task 1: Freeze live inputs and remote prerequisites

**Files:**
- Modify: `TRACKING.md`
- Modify: `docs/plan/dual_machine_experiment_charter_20260409.md`
- Create: `reports/closeout/alphagenome_matched_raw_track_slice_v1_closeout_20260410.md`
- Remote inspect only: `/data1/zhoujiazhen/bylw_atac/TransChromBP`

- [ ] **Step 1: Re-check `6000` host state and confirm the machine is available**

Run:

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 \
  "nvidia-smi --query-gpu=index,name,utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits && \
   echo '---' && \
   ps -eo pid,etimes,cmd | grep -E 'alphagenome|run_alphagenome_pilot|teacher_v2|ntv2_' | grep -v grep"
```

Expected:

- Both A6000 rows should still be idle or clearly available.
- No existing AlphaGenome pilot process should already be using the target output directory or log path.

- [ ] **Step 2: Re-check `6002` only to keep the top live row factual**

Run:

```bash
ssh -i /home/zhengwei/.ssh/codex_6002_ed25519 -p 6002 zhengwei@127.0.0.1 \
  "nvidia-smi --query-gpu=index,name,utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits && \
   echo '---' && \
   ps -eo pid,etimes,cmd | grep 'teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r4' | grep -v grep || true"
```

Expected:

- Record whether `r4` is still running or already complete.
- Do not let `6002` state change the `6000` launch decision; this check is only for later doc accuracy.

- [ ] **Step 3: Verify the `6000` AlphaGenome env and the comparison CSV are really available**

Run:

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 \
  "cd /data1/zhoujiazhen/bylw_atac/TransChromBP && \
   /data1/zhoujiazhen/bylw_atac/.mamba/envs/alphagenome/bin/python - <<'PY'
import sys
mods = {}
for name in ['alphagenome', 'numpy', 'pandas']:
    try:
        __import__(name)
        mods[name] = True
    except Exception as exc:
        mods[name] = f'ERR:{type(exc).__name__}:{exc}'
print(sys.version)
print(mods)
PY
   test -f outputs/reports/transchrombp_tutorial_test_20260315/selected_loci_prediction_totals.csv && \
   echo LOCAL_TOTALS_OK"
```

Expected:

- Python should report `3.11.x`.
- The module map should show `alphagenome`, `numpy`, and `pandas` as `True`.
- `LOCAL_TOTALS_OK` must appear.

- [ ] **Step 4: Confirm the runtime repo still lacks the local pilot scripts**

Run:

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 \
  "cd /data1/zhoujiazhen/bylw_atac/TransChromBP && \
   test -d scripts/alphagenome_pilot && echo REMOTE_PILOT_DIR_EXISTS || echo REMOTE_PILOT_DIR_MISSING"
```

Expected:

- `REMOTE_PILOT_DIR_MISSING` is acceptable and currently expected.
- This check determines whether Task 2 must create and sync the minimal runtime scripts.

### Task 2: Sync the minimal pilot assets into the `6000` runtime repo

**Files:**
- Local source: `scripts/alphagenome_pilot/run_alphagenome_pilot.py`
- Local source: `scripts/alphagenome_pilot/merge_locus_totals.py`
- Local source: `scripts/alphagenome_pilot/regions_k562_tutorial_selected_loci.csv`
- Remote create: `/data1/zhoujiazhen/bylw_atac/TransChromBP/scripts/alphagenome_pilot/run_alphagenome_pilot.py`
- Remote create: `/data1/zhoujiazhen/bylw_atac/TransChromBP/scripts/alphagenome_pilot/merge_locus_totals.py`
- Remote create: `/data1/zhoujiazhen/bylw_atac/TransChromBP/scripts/alphagenome_pilot/regions_k562_tutorial_selected_loci.csv`

- [ ] **Step 1: Stage a remote temp directory and upload the three required files**

Run:

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 \
  "rm -rf /tmp/alphagenome_pilot_sync_20260410 && mkdir -p /tmp/alphagenome_pilot_sync_20260410"

scp -P 6000 \
  scripts/alphagenome_pilot/run_alphagenome_pilot.py \
  scripts/alphagenome_pilot/merge_locus_totals.py \
  scripts/alphagenome_pilot/regions_k562_tutorial_selected_loci.csv \
  zhoujiazhen@127.0.0.1:/tmp/alphagenome_pilot_sync_20260410/
```

Expected:

- All three files should land in `/tmp/alphagenome_pilot_sync_20260410/` on `6000`.

- [ ] **Step 2: Atomically move the uploaded files into the runtime repo**

Run:

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 \
  "mkdir -p /data1/zhoujiazhen/bylw_atac/TransChromBP/scripts/alphagenome_pilot && \
   mv /tmp/alphagenome_pilot_sync_20260410/run_alphagenome_pilot.py \
      /data1/zhoujiazhen/bylw_atac/TransChromBP/scripts/alphagenome_pilot/run_alphagenome_pilot.py && \
   mv /tmp/alphagenome_pilot_sync_20260410/merge_locus_totals.py \
      /data1/zhoujiazhen/bylw_atac/TransChromBP/scripts/alphagenome_pilot/merge_locus_totals.py && \
   mv /tmp/alphagenome_pilot_sync_20260410/regions_k562_tutorial_selected_loci.csv \
      /data1/zhoujiazhen/bylw_atac/TransChromBP/scripts/alphagenome_pilot/regions_k562_tutorial_selected_loci.csv"
```

Expected:

- The runtime repo should now have a minimal `scripts/alphagenome_pilot/` directory with exactly three files.

- [ ] **Step 3: Verify the remote files are present and syntactically valid**

Run:

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 \
  "cd /data1/zhoujiazhen/bylw_atac/TransChromBP && \
   ls scripts/alphagenome_pilot && \
   /data1/zhoujiazhen/bylw_atac/.mamba/envs/alphagenome/bin/python -m py_compile \
     scripts/alphagenome_pilot/run_alphagenome_pilot.py \
     scripts/alphagenome_pilot/merge_locus_totals.py && \
   sed -n '1,10p' scripts/alphagenome_pilot/regions_k562_tutorial_selected_loci.csv"
```

Expected:

- `ls` should show the three uploaded files.
- `py_compile` should exit cleanly.
- The CSV preview should show the four fixed loci with `label,source,chrom,center,notes` header.

### Task 3: Launch the `v1` pilot and verify the background job

**Files:**
- Remote create: `/data1/zhoujiazhen/bylw_atac/TransChromBP/logs/alphagenome_matched_raw_track_slice_v1_20260410.log`
- Remote create: `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/alphagenome_pilot/alphagenome_matched_raw_track_slice_v1_20260410/`

- [ ] **Step 1: Start one background job that runs the pilot and merge end-to-end**

Run:

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 \
  "cd /data1/zhoujiazhen/bylw_atac/TransChromBP && \
   mkdir -p logs outputs/alphagenome_pilot/alphagenome_matched_raw_track_slice_v1_20260410 && \
   nohup bash -lc '
     /data1/zhoujiazhen/bylw_atac/.mamba/envs/alphagenome/bin/python \
       scripts/alphagenome_pilot/run_alphagenome_pilot.py \
       --regions-csv scripts/alphagenome_pilot/regions_k562_tutorial_selected_loci.csv \
       --output-dir outputs/alphagenome_pilot/alphagenome_matched_raw_track_slice_v1_20260410 && \
     /data1/zhoujiazhen/bylw_atac/.mamba/envs/alphagenome/bin/python \
       scripts/alphagenome_pilot/merge_locus_totals.py \
       --alpha-summary outputs/alphagenome_pilot/alphagenome_matched_raw_track_slice_v1_20260410/summary.csv \
       --local-totals outputs/reports/transchrombp_tutorial_test_20260315/selected_loci_prediction_totals.csv \
       --output-csv outputs/alphagenome_pilot/alphagenome_matched_raw_track_slice_v1_20260410/merged_locus_totals.csv
   ' > logs/alphagenome_matched_raw_track_slice_v1_20260410.log 2>&1 & \
   echo \$!"
```

Expected:

- The command should print one background PID.
- This is a `6000`-side independent pilot, but it is not expected to occupy both A6000 cards.

- [ ] **Step 2: Immediately verify that the process and log both exist**

Run:

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 \
  "cd /data1/zhoujiazhen/bylw_atac/TransChromBP && \
   ps -eo pid,etimes,cmd | grep 'alphagenome_matched_raw_track_slice_v1_20260410' | grep -v grep && \
   echo '---LOG---' && \
   tail -n 20 logs/alphagenome_matched_raw_track_slice_v1_20260410.log"
```

Expected:

- Either the process is still running, or the log already contains the AlphaGenome script output.
- If the process is absent and the log shows an immediate error, treat the run as `launch-blocked`.

- [ ] **Step 3: Record the absolute start time and derive the ETA window**

Run:

```bash
TZ=Asia/Shanghai date '+%Y-%m-%d %H:%M:%S CST'
```

Expected:

- Record one absolute launch timestamp.
- The ETA window to write into `TRACKING.md` is `start + 1h` to `start + 2h`.

- [ ] **Step 4: Tell the user how to monitor the run and then pause for completion confirmation**

Use exactly these monitor commands in the user-facing update:

```bash
tail -f /data1/zhoujiazhen/bylw_atac/TransChromBP/logs/alphagenome_matched_raw_track_slice_v1_20260410.log
grep -i "Wrote merged comparison" /data1/zhoujiazhen/bylw_atac/TransChromBP/logs/alphagenome_matched_raw_track_slice_v1_20260410.log
```

Expected:

- Do not start closeout tasks until the user confirms the run finished or the log clearly shows completion.

### Task 4: Update the live docs right after launch is verified

**Files:**
- Modify: `TRACKING.md`
- Modify: `docs/plan/dual_machine_experiment_charter_20260409.md`

- [ ] **Step 1: Update the top live row and `6000` row in `TRACKING.md`**

Edit `TRACKING.md` so it records:

```md
- `6000` 当前 active run：`alphagenome_matched_raw_track_slice_v1_20260410`
- 机器性质：`6000` 独立 pilot / API 推理，不是双卡训练
- 日志：`/data1/zhoujiazhen/bylw_atac/TransChromBP/logs/alphagenome_matched_raw_track_slice_v1_20260410.log`
- 输出目录：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/alphagenome_pilot/alphagenome_matched_raw_track_slice_v1_20260410`
- ETA：写成刚记录的绝对时间窗
- 下一步：等待 `v1` 收口并按 gate 写 closeout
```

Do not say that A6000 GPUs are actively occupied if `nvidia-smi` still shows them idle.

- [ ] **Step 2: Update the charter so `6000` backlog moves from candidate to active pilot**

Edit Section 5.1 / Section 10 so they state:

```md
- `6000` 当前 active 项已从“候选 backlog”变为 `AlphaGenome matched raw-track slice v1`
- 这是 `6000` 独立推进的 external-coordinate pilot
- 它与 `6002` 当前队列解耦
- 它不是双卡长训练，按 `1-2 h` pilot 管理
```

- [ ] **Step 3: Commit the launch-state doc update**

Run:

```bash
git add TRACKING.md docs/plan/dual_machine_experiment_charter_20260409.md
git commit -m "docs: track alphagenome pilot launch"
```

Expected:

- Commit succeeds and contains only the launch-state live doc changes.

### Task 5: Close out the pilot after completion

**Files:**
- Modify: `TRACKING.md`
- Modify: `docs/plan/dual_machine_experiment_charter_20260409.md`
- Create: `reports/closeout/alphagenome_matched_raw_track_slice_v1_closeout_20260410.md`

- [ ] **Step 1: Re-check the finished log and the expected output files**

Run:

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 \
  "cd /data1/zhoujiazhen/bylw_atac/TransChromBP && \
   tail -n 60 logs/alphagenome_matched_raw_track_slice_v1_20260410.log && \
   echo '---FILES---' && \
   ls outputs/alphagenome_pilot/alphagenome_matched_raw_track_slice_v1_20260410 && \
   echo '---META---' && \
   sed -n '1,40p' outputs/alphagenome_pilot/alphagenome_matched_raw_track_slice_v1_20260410/run_meta.json"
```

Expected:

- The output directory must contain `summary.csv`, `region_metadata.jsonl`, `profiles/`, `run_meta.json`, and `merged_locus_totals.csv` for a full `pass` candidate.
- If the run stopped earlier, capture the exact missing artifact or error instead of assuming success.

- [ ] **Step 2: Inspect summary/metadata and classify the result**

Run:

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 \
  "cd /data1/zhoujiazhen/bylw_atac/TransChromBP && \
   /data1/zhoujiazhen/bylw_atac/.mamba/envs/alphagenome/bin/python - <<'PY'
import json
from pathlib import Path
import pandas as pd

base = Path('outputs/alphagenome_pilot/alphagenome_matched_raw_track_slice_v1_20260410')
summary = pd.read_csv(base / 'summary.csv')
meta_rows = [json.loads(line) for line in (base / 'region_metadata.jsonl').read_text(encoding='utf-8').splitlines()]
print('num_summary_rows', len(summary))
print('labels', summary['label'].tolist())
print('num_tracks_used', summary['num_tracks_used'].tolist())
print('alpha_total_count', summary['alpha_total_count'].round(4).tolist())
print('metadata_rows', len(meta_rows))
print('track_names_preview', meta_rows[0].get('track_names', [])[:5] if meta_rows else [])
PY"
```

Expected:

- `num_summary_rows` should be `4` for a clean `pass` candidate.
- `num_tracks_used` should stay positive for every locus.
- Use this output plus the file presence from Step 1 to classify the run as `pass`, `fail`, or `launch-blocked`.

- [ ] **Step 3: Write the closeout report with a hard verdict**

Create `reports/closeout/alphagenome_matched_raw_track_slice_v1_closeout_20260410.md` with these sections:

```md
# AlphaGenome Matched Raw-Track Slice v1 Closeout (2026-04-10)

## Verdict
- `pass` / `fail` / `launch-blocked`

## Run Identity
- run name
- machine
- log path
- output dir
- start / end times

## Evidence
- whether all 4 loci succeeded
- whether `summary.csv`, `region_metadata.jsonl`, `profiles/*.npz`, `run_meta.json`, `merged_locus_totals.csv` exist
- how many ATAC tracks remained after ontology/filter
- whether the merged table formed cleanly

## Gate Decision
- if `pass`: next step is only `12-20 loci` expansion
- if `fail`: stop this pilot and do not silently rerun
- if `launch-blocked`: fix the blocker first; do not write this up as an experiment result
```

- [ ] **Step 4: Update `TRACKING.md` and the charter from active to closed**

Edit both docs so they reflect one of these states:

```md
- `pass`: `6000` completed `v1`; next allowed step is only `12-20 loci`
- `fail`: `6000` completed `v1`; pilot stopped; next move is the next genuinely new high-value task
- `launch-blocked`: `6000` launch failed before a valid run existed; blocker documented in the closeout report
```

Do not leave the run in `进行中` after the log and artifacts already show a terminal state.

- [ ] **Step 5: Commit the closeout docs**

Run:

```bash
git add TRACKING.md docs/plan/dual_machine_experiment_charter_20260409.md reports/closeout/alphagenome_matched_raw_track_slice_v1_closeout_20260410.md
git commit -m "docs: close out alphagenome pilot v1"
```

Expected:

- Commit succeeds and contains only the closeout verdict and report updates.
