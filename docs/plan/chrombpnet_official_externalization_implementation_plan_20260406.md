# Official ChromBPNet Externalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the local official `chrombpnet/` payload from this repository without breaking active official-compare or dataset-prep workflows, by bridging every live dependency to the 6000 official repository and rewriting canonical docs around the new repo role.

**Architecture:** Keep this repository as the `TransChromBP` main repo plus project archive, and treat `/data1/zhoujiazhen/bylw_atac/chrombpnet_official` on 6000 as the only canonical official ChromBPNet source tree. Active scripts in this repo must either call that external official root explicitly or fail with a clear message; only after all live dependencies are cut do we delete `chrombpnet/`, `chrombpnet.egg-info/`, `setup.py`, and `MANIFEST.in`.

**Tech Stack:** Bash, Python 3, SSH/SCP, ripgrep, Git, `bash -n`, `python3 -m py_compile`

## `2026-04-08` Execution Note

- Fresh 6000 verification this round showed `/data1/zhoujiazhen/bylw_atac/chrombpnet_official` does not currently exist; the real usable official source root is still `/data1/zhoujiazhen/bylw_atac/chromBPNet`.
- Therefore the `CHROMBPNET_OFFICIAL_ROOT=/data1/zhoujiazhen/bylw_atac/chrombpnet_official` strings below should currently be read as a target alias, not as an already-validated machine path.
- This round completed the helper-bridge part of the plan locally: [test_chrombpnet_official_externalization.sh](/home/zhengwei/project/python/chromBPNet/tests/test_chrombpnet_official_externalization.sh), [run_remote_chrombpnet_dataset_prep.sh](/home/zhengwei/project/python/chromBPNet/scripts/run_remote_chrombpnet_dataset_prep.sh), [start_6000_chrombpnet_dataset_prep.sh](/home/zhengwei/project/python/chromBPNet/scripts/start_6000_chrombpnet_dataset_prep.sh), [start_6002_chrombpnet_dataset_prep.sh](/home/zhengwei/project/python/chromBPNet/scripts/start_6002_chrombpnet_dataset_prep.sh), and [step3_get_background_regions.sh](/home/zhengwei/project/python/chromBPNet/workflows/tutorial/step3_get_background_regions.sh) now pass the local guard and `bash -n`.
- A real 6000 smoke was then launched as `chrombpnet_official_step3_bridge_smoke_20260408_110411`, explicitly exporting `CHROMBPNET_OFFICIAL_ROOT=/data1/zhoujiazhen/bylw_atac/chromBPNet` and running the patched tutorial `step3_get_background_regions.sh` in a fresh scratch output directory. This was chosen over the full dataset-prep launcher because existing `prep_v1` artifacts would cause most steps to `skip`, which would not truly validate the helper bridge.
- `2026-04-08 11:21 CST` re-check: the smoke passed. Evidence chain:
  - process exited cleanly and log ended with `Completed execution`
  - scratch output `negatives_with_summit.bed` exists at `/data1/zhoujiazhen/bylw_atac/.codex_jobs/chrombpnet_official_step3_bridge_smoke_20260408_110411/output/negatives_with_summit.bed`
  - `wc -l` on that file is `508429`
  - the actual executed helper path was `/data1/zhoujiazhen/bylw_atac/chromBPNet/chrombpnet/helpers/make_gc_matched_negatives/make_gc_matched_negatives.sh`
- Therefore the helper bridge is now validated on a real 6000 run. The remaining documentation action is no longer “wait for smoke”, but “choose one canonical root and make it consistent everywhere”: either promote `/data1/zhoujiazhen/bylw_atac/chromBPNet` to the validated canonical root, or create and verify a real `chrombpnet_official` alias first.

---

## File Map

- `scripts/paper_aligned_repro/select_best_epoch.py`
  - External best-epoch selector. Must stop importing local `chrombpnet` and always shell out to the official `predict.py` under `CHROMBPNET_OFFICIAL_ROOT` or `--official-root`.
- `scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh`
  - Paper-aligned training/eval driver. Must resolve the official root once, reuse it for bias/chrom prediction metrics, and stop pointing at `REPO_ROOT/chrombpnet/training/predict.py`.
- `scripts/run_remote_chrombpnet_dataset_prep.sh`
  - Remote dataset-prep driver. Must support two helper sources:
    - `--official-root` for in-place use on 6000
    - `--gc-helper-dir` for staged helper files on 6002
- `scripts/start_6000_chrombpnet_dataset_prep.sh`
  - 6000 launcher. Must stop uploading helper scripts from the local repo and pass `--official-root /data1/zhoujiazhen/bylw_atac/chrombpnet_official`.
- `scripts/start_6002_chrombpnet_dataset_prep.sh`
  - 6002 launcher. Must fetch the three GC helper python files from 6000’s official repo into a local temp dir, then upload them to the 6002 job dir as runtime staging only.
- `workflows/tutorial/step3_get_background_regions.sh`
  - Tutorial background-region workflow. Must resolve `make_gc_matched_negatives.sh` from `MAKE_GC_MATCHED_NEGATIVES_SCRIPT` or `CHROMBPNET_OFFICIAL_ROOT`, never `../../chrombpnet/...`.
- `tests/test_chrombpnet_official_externalization.sh`
  - Regression guard for active scripts/docs: no direct local official-path dependencies, required flags exist, and deleted payload paths stay gone.
- `reports/chrombpnet_official_patch_ledger_20260406.md`
  - Durable patch ledger for the official files we still depend on operationally after externalization (`predict.py`, `metrics.py`, GC helper scripts).
- `scripts/paper_aligned_repro/README.md`
  - Operator instructions for strict compare / best-epoch selection. Must tell users to set `CHROMBPNET_OFFICIAL_ROOT=/data1/zhoujiazhen/bylw_atac/chrombpnet_official`.
- `AGENTS.md`
  - Canonical repo guidance. Must describe the repo as `TransChromBP` main repo + archive, not as a local official ChromBPNet package.
- `DEVELOPMENT.md`
  - Dev workflow doc. Must remove “root repo contains `chrombpnet/` core code” and explain that official reproduction/code lookup happen on 6000’s external repo.
- `README.md`
  - Root repo landing page. Must be rewritten from the official ChromBPNet README to a project-specific overview.
- `TRACKING.md`
  - Live status. Must reflect the implementation progress and final verification state of this migration.
- Delete:
  - `chrombpnet/`
  - `chrombpnet.egg-info/`
  - `setup.py`
  - `MANIFEST.in`

### Task 1: Refactor Official Predict Entrypoints

**Files:**
- Create: `tests/test_chrombpnet_official_externalization.sh`
- Modify: `scripts/paper_aligned_repro/select_best_epoch.py`
- Modify: `scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh`
- Test: `tests/test_chrombpnet_official_externalization.sh`

- [ ] **Step 1: Write the failing regression guard for official predict consumers**

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

assert_no_match() {
  local pattern="$1"
  shift
  if rg -n "$pattern" "$@" >/dev/null; then
    echo "unexpected match for pattern: $pattern" >&2
    rg -n "$pattern" "$@" >&2
    exit 1
  fi
}

assert_contains_text() {
  local needle="$1"
  local file="$2"
  grep -Fq "$needle" "$file" || {
    echo "missing text '$needle' in $file" >&2
    exit 1
  }
}

assert_no_match 'REPO_ROOT/chrombpnet/training/predict.py|from chrombpnet.training import predict' \
  "$ROOT/scripts/paper_aligned_repro/select_best_epoch.py" \
  "$ROOT/scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh"

python3 "$ROOT/scripts/paper_aligned_repro/select_best_epoch.py" --help | grep -q -- '--official-root'
bash "$ROOT/scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh" --help | grep -q -- '--official-root'
assert_contains_text 'CHROMBPNET_OFFICIAL_ROOT' "$ROOT/scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh"
```

- [ ] **Step 2: Run the regression guard and verify it fails on current local-path usage**

Run: `bash tests/test_chrombpnet_official_externalization.sh`

Expected: FAIL with matches from `scripts/paper_aligned_repro/select_best_epoch.py` and `scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh` because they still reference local `chrombpnet/training/predict.py`.

- [ ] **Step 3: Replace local `predict.py` coupling with explicit official-root resolution**

```python
# scripts/paper_aligned_repro/select_best_epoch.py
def resolve_official_root(raw: str) -> Path:
    value = raw or os.environ.get("CHROMBPNET_OFFICIAL_ROOT", "")
    if not value:
        raise SystemExit(
            "missing official ChromBPNet root; pass --official-root or set CHROMBPNET_OFFICIAL_ROOT"
        )
    root = Path(value).expanduser().resolve()
    predict_py = root / "chrombpnet" / "training" / "predict.py"
    if not predict_py.is_file():
        raise SystemExit(f"official predict.py not found: {predict_py}")
    return root

parser.add_argument(
    "--official-root",
    default="",
    help="Path to the external official ChromBPNet repo; defaults to CHROMBPNET_OFFICIAL_ROOT",
)

def build_predict_cmd(model_path: Path, output_prefix: Path, args: argparse.Namespace, official_root: Path) -> list[str]:
    return [
        sys.executable,
        str(official_root / "chrombpnet" / "training" / "predict.py"),
        "-g", args.genome,
        "-b", args.bigwig,
        "-p", args.peaks,
        "-n", args.nonpeaks,
        "-o", str(output_prefix),
        "-fl", args.fold_json,
        "-m", str(model_path),
        "-bs", str(args.batch_size),
        "-s", str(args.seed),
        "-il", str(args.inputlen),
        "-ol", str(args.outputlen),
        "--split", args.split,
        "--metrics-only",
    ]
```

```bash
# scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh
OFFICIAL_ROOT="${CHROMBPNET_OFFICIAL_ROOT:-}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --official-root) OFFICIAL_ROOT="$2"; shift 2 ;;
    # keep existing cases below unchanged
  esac
done

if [[ -z "${OFFICIAL_ROOT}" ]]; then
  echo "ERROR: missing official ChromBPNet root; pass --official-root or set CHROMBPNET_OFFICIAL_ROOT" >&2
  exit 1
fi

PREDICT_PY="${OFFICIAL_ROOT}/chrombpnet/training/predict.py"
if [[ ! -f "${PREDICT_PY}" ]]; then
  echo "ERROR: official predict.py not found: ${PREDICT_PY}" >&2
  exit 1
fi

CUDA_VISIBLE_DEVICES="${predict_gpu}" CHROMBPNET_MULTI_GPU=0 python3 "${PREDICT_PY}" \
  -g "${GENOME}" \
  -b "${bias_bigwig}" \
  -p "${PEAKS}" \
  -n "${nonpeaks}" \
  -o "${bias_dir}/evaluation/bias" \
  -fl "${fold_json}" \
  -m "${bias_model}" \
  -bs "${PREDICT_BATCH_SIZE}" \
  -il 2114 \
  -ol 1000 \
  -s "${SEED}"
```

- [ ] **Step 4: Re-run the guard and syntax checks**

Run: `bash tests/test_chrombpnet_official_externalization.sh`

Expected: PASS

Run: `python3 -m py_compile scripts/paper_aligned_repro/select_best_epoch.py`

Expected: PASS with no output

Run: `bash -n scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh`

Expected: PASS with no output

- [ ] **Step 5: Commit**

```bash
git add tests/test_chrombpnet_official_externalization.sh \
  scripts/paper_aligned_repro/select_best_epoch.py \
  scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh
git commit -m "refactor official predict entrypoints"
```

### Task 2: Bridge Dataset Prep to the External Official Helper Source

**Files:**
- Modify: `tests/test_chrombpnet_official_externalization.sh`
- Modify: `scripts/run_remote_chrombpnet_dataset_prep.sh`
- Modify: `scripts/start_6000_chrombpnet_dataset_prep.sh`
- Modify: `scripts/start_6002_chrombpnet_dataset_prep.sh`
- Modify: `workflows/tutorial/step3_get_background_regions.sh`
- Test: `tests/test_chrombpnet_official_externalization.sh`

- [ ] **Step 1: Extend the regression guard for GC-helper path contracts**

```bash
assert_no_match 'REPO_ROOT/chrombpnet/helpers|../../chrombpnet/helpers' \
  "$ROOT/scripts/start_6000_chrombpnet_dataset_prep.sh" \
  "$ROOT/scripts/start_6002_chrombpnet_dataset_prep.sh" \
  "$ROOT/workflows/tutorial/step3_get_background_regions.sh"

bash "$ROOT/scripts/run_remote_chrombpnet_dataset_prep.sh" --help | grep -q -- '--official-root'
bash "$ROOT/scripts/run_remote_chrombpnet_dataset_prep.sh" --help | grep -q -- '--gc-helper-dir'
grep -Fq 'CHROMBPNET_OFFICIAL_ROOT' "$ROOT/workflows/tutorial/step3_get_background_regions.sh"
```

- [ ] **Step 2: Run the guard and verify it fails on current helper-script references**

Run: `bash tests/test_chrombpnet_official_externalization.sh`

Expected: FAIL with matches from `scripts/start_6000_chrombpnet_dataset_prep.sh`, `scripts/start_6002_chrombpnet_dataset_prep.sh`, or `workflows/tutorial/step3_get_background_regions.sh`.

- [ ] **Step 3: Implement explicit helper-source switching for 6000, 6002, and tutorial step 3**

```bash
# scripts/run_remote_chrombpnet_dataset_prep.sh
OFFICIAL_ROOT=""
GC_HELPER_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --official-root) OFFICIAL_ROOT="$2"; shift 2 ;;
    --gc-helper-dir) GC_HELPER_DIR="$2"; shift 2 ;;
    # keep existing cases below unchanged
  esac
done

if [[ -n "$OFFICIAL_ROOT" ]]; then
  GC_CONTENT_PY="$OFFICIAL_ROOT/chrombpnet/helpers/make_gc_matched_negatives/get_gc_content.py"
  GC_MATCH_PY="$OFFICIAL_ROOT/chrombpnet/helpers/make_gc_matched_negatives/get_gc_matched_negatives.py"
  GC_BINS_PY="$OFFICIAL_ROOT/chrombpnet/helpers/make_gc_matched_negatives/get_genomewide_gc_buckets/get_genomewide_gc_bins.py"
elif [[ -n "$GC_HELPER_DIR" ]]; then
  GC_CONTENT_PY="$GC_HELPER_DIR/get_gc_content.py"
  GC_MATCH_PY="$GC_HELPER_DIR/get_gc_matched_negatives.py"
  GC_BINS_PY="$GC_HELPER_DIR/get_genomewide_gc_bins.py"
else
  fail "missing helper source: pass --official-root or --gc-helper-dir"
fi
```

```bash
# scripts/start_6000_chrombpnet_dataset_prep.sh
OFFICIAL_ROOT="${CHROMBPNET_OFFICIAL_ROOT:-/data1/zhoujiazhen/bylw_atac/chrombpnet_official}"

"${scp_base[@]}" \
  "$REPO_ROOT/scripts/run_remote_chrombpnet_dataset_prep.sh" \
  "$REMOTE_HOST:$REMOTE_JOB_DIR/"

pid="$("${ssh_base[@]}" \
  "nohup bash '$REMOTE_JOB_DIR/run_remote_chrombpnet_dataset_prep.sh' \
    --root '$REMOTE_ROOT' \
    --env-dir '$REMOTE_ENV' \
    --python-bin '$REMOTE_PYTHON' \
    --datasets '$DATASETS' \
    --threads '$THREADS' \
    --nice-level '$NICE_LEVEL' \
    --official-root '$OFFICIAL_ROOT' \
    > '$REMOTE_LOG' 2>&1 < /dev/null & echo \$!")"
```

```bash
# scripts/start_6002_chrombpnet_dataset_prep.sh
OFFICIAL_HOST="${OFFICIAL_HOST:-zhoujiazhen@127.0.0.1}"
OFFICIAL_PORT="${OFFICIAL_PORT:-6000}"
OFFICIAL_ROOT="${CHROMBPNET_OFFICIAL_ROOT:-/data1/zhoujiazhen/bylw_atac/chrombpnet_official}"
LOCAL_STAGE="$(mktemp -d)"
trap 'rm -rf "$LOCAL_STAGE"' EXIT

scp -P "$OFFICIAL_PORT" \
  "$OFFICIAL_HOST:$OFFICIAL_ROOT/chrombpnet/helpers/make_gc_matched_negatives/get_gc_content.py" \
  "$OFFICIAL_HOST:$OFFICIAL_ROOT/chrombpnet/helpers/make_gc_matched_negatives/get_gc_matched_negatives.py" \
  "$OFFICIAL_HOST:$OFFICIAL_ROOT/chrombpnet/helpers/make_gc_matched_negatives/get_genomewide_gc_buckets/get_genomewide_gc_bins.py" \
  "$LOCAL_STAGE/"

"${scp_base[@]}" \
  "$REPO_ROOT/scripts/run_remote_chrombpnet_dataset_prep.sh" \
  "$LOCAL_STAGE/get_gc_content.py" \
  "$LOCAL_STAGE/get_gc_matched_negatives.py" \
  "$LOCAL_STAGE/get_genomewide_gc_bins.py" \
  "$REMOTE_HOST:$REMOTE_JOB_DIR/"

pid="$("${ssh_base[@]}" \
  "nohup bash '$REMOTE_JOB_DIR/run_remote_chrombpnet_dataset_prep.sh' \
    --root '$REMOTE_ROOT' \
    --env-dir '$REMOTE_ENV' \
    --python-bin '$REMOTE_PYTHON' \
    --datasets '$DATASETS' \
    --threads '$THREADS' \
    --nice-level '$NICE_LEVEL' \
    --gc-helper-dir '$REMOTE_JOB_DIR' \
    > '$REMOTE_LOG' 2>&1 < /dev/null & echo \$!")"
```

```bash
# workflows/tutorial/step3_get_background_regions.sh
make_gc_matched_negatives_script="${MAKE_GC_MATCHED_NEGATIVES_SCRIPT:-}"

if [[ -z "${make_gc_matched_negatives_script}" ]]; then
    official_root="${CHROMBPNET_OFFICIAL_ROOT:-}"
    if [[ -z "${official_root}" ]]; then
        echo "missing CHROMBPNET_OFFICIAL_ROOT or MAKE_GC_MATCHED_NEGATIVES_SCRIPT" >&2
        exit 1
    fi
    make_gc_matched_negatives_script="${official_root}/chrombpnet/helpers/make_gc_matched_negatives/make_gc_matched_negatives.sh"
fi

if [[ ! -x "${make_gc_matched_negatives_script}" ]]; then
    echo "missing helper script: ${make_gc_matched_negatives_script}" >&2
    exit 1
fi
```

- [ ] **Step 4: Run syntax checks and the updated guard**

Run: `bash tests/test_chrombpnet_official_externalization.sh`

Expected: PASS

Run: `bash -n scripts/run_remote_chrombpnet_dataset_prep.sh scripts/start_6000_chrombpnet_dataset_prep.sh scripts/start_6002_chrombpnet_dataset_prep.sh workflows/tutorial/step3_get_background_regions.sh`

Expected: PASS with no output

- [ ] **Step 5: Commit**

```bash
git add tests/test_chrombpnet_official_externalization.sh \
  scripts/run_remote_chrombpnet_dataset_prep.sh \
  scripts/start_6000_chrombpnet_dataset_prep.sh \
  scripts/start_6002_chrombpnet_dataset_prep.sh \
  workflows/tutorial/step3_get_background_regions.sh
git commit -m "bridge dataset prep to official helper source"
```

### Task 3: Record the Official Patch Ledger and Operator Contract

**Files:**
- Create: `reports/chrombpnet_official_patch_ledger_20260406.md`
- Modify: `scripts/paper_aligned_repro/README.md`
- Test: `reports/chrombpnet_official_patch_ledger_20260406.md`

- [ ] **Step 1: Write the failing operator-contract check**

```bash
rg -n 'CHROMBPNET_OFFICIAL_ROOT=/data1/zhoujiazhen/bylw_atac/chrombpnet_official|chrombpnet_official_patch_ledger_20260406.md' \
  scripts/paper_aligned_repro/README.md \
  reports/chrombpnet_official_patch_ledger_20260406.md
```

- [ ] **Step 2: Run the check and verify it fails before the doc changes**

Run: `rg -n 'CHROMBPNET_OFFICIAL_ROOT=/data1/zhoujiazhen/bylw_atac/chrombpnet_official|chrombpnet_official_patch_ledger_20260406.md' scripts/paper_aligned_repro/README.md reports/chrombpnet_official_patch_ledger_20260406.md`

Expected: FAIL because `reports/chrombpnet_official_patch_ledger_20260406.md` does not exist yet and the README still points operators at the local repo.

- [ ] **Step 3: Add the ledger and rewrite the operator README around the external official root**

```markdown
# reports/chrombpnet_official_patch_ledger_20260406.md
## Canonical official root

- 6000 official repo: `/data1/zhoujiazhen/bylw_atac/chrombpnet_official`
- This main repo no longer vendors official ChromBPNet source files.

## Active official files we still rely on

| Official file | Current use | Canonical runtime path | Why we still care |
| --- | --- | --- | --- |
| `chrombpnet/training/predict.py` | strict compare, external best-epoch selector, paper-aligned metrics | `${CHROMBPNET_OFFICIAL_ROOT}/chrombpnet/training/predict.py` | Official metrics path plus local classification-metric patch lineage |
| `chrombpnet/training/metrics.py` | predict-time metric generation | `${CHROMBPNET_OFFICIAL_ROOT}/chrombpnet/training/metrics.py` | Needed to explain why current official compare emits the extra metrics we cite |
| `chrombpnet/helpers/make_gc_matched_negatives/*` | dataset prep and tutorial background-region generation | `${CHROMBPNET_OFFICIAL_ROOT}/chrombpnet/helpers/make_gc_matched_negatives/...` | Needed by 6000 in place and by 6002 staged helper uploads |
```

````markdown
# scripts/paper_aligned_repro/README.md
## 1. 环境准备（6000）

```bash
export CHROMBPNET_ENV=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet
export CHROMBPNET_OFFICIAL_ROOT=/data1/zhoujiazhen/bylw_atac/chrombpnet_official
export PATH="$CHROMBPNET_ENV/bin:$PATH"
export LD_LIBRARY_PATH="$CHROMBPNET_ENV/lib:$LD_LIBRARY_PATH"
export PYTHONPATH=/data1/zhoujiazhen/bylw_atac/TransChromBP/vendor/transchrombp:$PYTHONPATH
```

- `select_best_epoch.py` 与 strict-compare launcher 默认依赖 `CHROMBPNET_OFFICIAL_ROOT`
- 官方 patch 台账见 `reports/chrombpnet_official_patch_ledger_20260406.md`
````

- [ ] **Step 4: Run the doc checks and placeholder scan**

Run: `rg -n 'CHROMBPNET_OFFICIAL_ROOT=/data1/zhoujiazhen/bylw_atac/chrombpnet_official|chrombpnet_official_patch_ledger_20260406.md' scripts/paper_aligned_repro/README.md reports/chrombpnet_official_patch_ledger_20260406.md`

Expected: PASS with matching lines in both files

Run: `PLACEHOLDER_PATTERNS='TO''DO|TB''D|之后再''说|稍''后|占''位'; rg -n "$PLACEHOLDER_PATTERNS" scripts/paper_aligned_repro/README.md reports/chrombpnet_official_patch_ledger_20260406.md`

Expected: no output

- [ ] **Step 5: Commit**

```bash
git add scripts/paper_aligned_repro/README.md \
  reports/chrombpnet_official_patch_ledger_20260406.md
git commit -m "document official chrombpnet contract"
```

### Task 4: Rewrite Canonical Repo Documentation

**Files:**
- Modify: `AGENTS.md`
- Modify: `DEVELOPMENT.md`
- Modify: `README.md`
- Test: `AGENTS.md`

- [ ] **Step 1: Write the failing stale-doc grep**

```bash
rg -n 'chrombpnet/ contains the Python package|pip install -e \\.|核心代码库|Bias factorized, base-resolution deep learning models' \
  AGENTS.md DEVELOPMENT.md README.md
```

- [ ] **Step 2: Run the grep and verify the old official-repo framing is still present**

Run: `rg -n 'chrombpnet/ contains the Python package|pip install -e \\.|核心代码库|Bias factorized, base-resolution deep learning models' AGENTS.md DEVELOPMENT.md README.md`

Expected: FAIL with hits from all three files.

- [ ] **Step 3: Replace the stale repo-role wording with the new main-repo contract**

```markdown
# README.md
# TransChromBP Main Repository

本仓库承担三类职责：
- `vendor/transchrombp/`：版本化的 `TransChromBP` snapshot
- `docs/`、`reports/`、`TRACKING.md`：项目文档、论文与证据链
- `scripts/`、`workflows/`：自研 launcher 与官方 bridge scripts

官方 `ChromBPNet` 源码不再保存在本仓。查阅官方实现、跑官方复现、或执行 official compare 时，请使用 6000 上的 `/data1/zhoujiazhen/bylw_atac/chrombpnet_official`。
```

```markdown
# AGENTS.md
## Project Structure & Module Organization
- `vendor/transchrombp/` is the versioned local snapshot of the TransChromBP codebase and helper scripts.
- `scripts/` holds project launchers, bridge scripts, benchmarking helpers, and maintenance utilities.
- 官方 `ChromBPNet` 不再常驻本仓；源码查阅与官方复现默认在 6000 `/data1/zhoujiazhen/bylw_atac/chrombpnet_official` 完成。

## Build, Test, and Development Commands
- `pip install -r requirements.txt` installs Python deps used by local docs/utility workflows.
- 不再使用 `pip install -e .` 把根仓当作官方 ChromBPNet 安装包。
- 自研 `transchrombp` 代码默认通过 `PYTHONPATH=vendor/transchrombp` 或远端 runtime 工作区使用。
```

```markdown
# DEVELOPMENT.md
以下内容默认应进入本地仓库并参与 Git 版本管理：

- `vendor/transchrombp/`、`scripts/`、`docs/`、`tests/` 等自研代码与文档目录
- `reports/*.tex` 与 `reports/assets/` 下的小型结果摘要
- `TRACKING.md`、`DEVELOPMENT.md`、实验安排文档

以下内容不再由本仓承载：

- 官方 `chrombpnet/` 源码树
- `setup.py`、`MANIFEST.in`、`chrombpnet.egg-info/` 这类官方 packaging 入口
- 官方复现与官方源码查阅；它们默认转到 6000 `/data1/zhoujiazhen/bylw_atac/chrombpnet_official`
```

- [ ] **Step 4: Re-run the stale-doc grep and confirm wrapper docs remain thin**

Run: `rg -n 'chrombpnet/ contains the Python package|pip install -e \\.|核心代码库|Bias factorized, base-resolution deep learning models' AGENTS.md DEVELOPMENT.md README.md`

Expected: no output

Run: `rg -n 'chrombpnet/' CLAUDE.md GEMINI.md`

Expected: no output

- [ ] **Step 5: Commit**

```bash
git add AGENTS.md DEVELOPMENT.md README.md
git commit -m "rewrite docs for external official chrombpnet"
```

### Task 5: Remove the Local Official Payload and Packaging

**Files:**
- Modify: `tests/test_chrombpnet_official_externalization.sh`
- Delete: `chrombpnet/`
- Delete: `chrombpnet.egg-info/`
- Delete: `setup.py`
- Delete: `MANIFEST.in`
- Test: `tests/test_chrombpnet_official_externalization.sh`

- [ ] **Step 1: Extend the regression guard to require the payload paths to be gone**

```bash
for path in \
  "$ROOT/chrombpnet" \
  "$ROOT/chrombpnet.egg-info" \
  "$ROOT/setup.py" \
  "$ROOT/MANIFEST.in"
do
  if [[ -e "$path" ]]; then
    echo "path should be removed after externalization: $path" >&2
    exit 1
  fi
done
```

- [ ] **Step 2: Run the guard and verify it fails before deletion**

Run: `bash tests/test_chrombpnet_official_externalization.sh`

Expected: FAIL because the local official payload still exists on disk.

- [ ] **Step 3: Delete the local official source tree and packaging files**

```bash
git rm -r chrombpnet chrombpnet.egg-info
git rm setup.py MANIFEST.in
```

- [ ] **Step 4: Re-run the guard and confirm no active packaging references remain**

Run: `bash tests/test_chrombpnet_official_externalization.sh`

Expected: PASS

Run: `rg -n 'chrombpnet\\.egg-info|pip install -e \\.|REPO_ROOT/chrombpnet' AGENTS.md DEVELOPMENT.md README.md scripts workflows tests`

Expected: no output

- [ ] **Step 5: Commit**

```bash
git add tests/test_chrombpnet_official_externalization.sh
git commit -m "remove local chrombpnet payload"
```

### Task 6: Sync Tracking and Run Final Verification

**Files:**
- Modify: `TRACKING.md`
- Test: `tests/test_chrombpnet_official_externalization.sh`

- [ ] **Step 1: Update `TRACKING.md` with the implementation result and next verification hook**

```markdown
| 官方 ChromBPNet 外置化实施 | 待验证 | design 与 implementation plan 已执行完成；本仓已切断 active local-official 依赖、补齐 patch ledger、删除本地 `chrombpnet/` 与官方 packaging 入口。 | 在 6000 用一次真实 official compare / dataset-prep smoke 验证桥接链路；若通过，再把该条移入 `TRACKING_archive.md`。 | `docs/plan/chrombpnet_official_externalization_design_20260406.md`、`docs/plan/chrombpnet_official_externalization_implementation_plan_20260406.md`、`reports/chrombpnet_official_patch_ledger_20260406.md` |
```

- [ ] **Step 2: Run the full verification bundle**

Run: `bash tests/test_chrombpnet_official_externalization.sh`

Expected: PASS

Run: `python3 -m py_compile scripts/paper_aligned_repro/select_best_epoch.py`

Expected: PASS with no output

Run: `bash -n scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh scripts/run_remote_chrombpnet_dataset_prep.sh scripts/start_6000_chrombpnet_dataset_prep.sh scripts/start_6002_chrombpnet_dataset_prep.sh workflows/tutorial/step3_get_background_regions.sh`

Expected: PASS with no output

Run: `find . -maxdepth 1 \\( -name 'chrombpnet' -o -name 'chrombpnet.egg-info' -o -name 'setup.py' -o -name 'MANIFEST.in' \\)`

Expected: no output

- [ ] **Step 3: Review the final tree state before closeout**

Run: `git status --short`

Expected: only the intended externalization changes appear in the worktree; no unrelated file reverts or accidental payload re-additions.

- [ ] **Step 4: Commit**

```bash
git add TRACKING.md
git commit -m "sync tracking for chrombpnet externalization"
```
