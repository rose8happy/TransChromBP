# Factor Ladder Unattended Matrix Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an unattended single-machine queue for `alphagenome_factor_ladder` that can attach to the live `E2 short10` run on `6000`, then continue through `E1 short10 -> E2 teacher30 -> teacher-cache export -> E3 distill short10` without manual intervention unless a technical fault occurs.

**Architecture:** Implement a small Python queue engine under the TransChromBP package to own queue config parsing, stage completion probes, checkpoint selection, state/event/summary writing, and narrow retry logic for port conflicts. Keep runtime orchestration in a thin bash wrapper plus a queue YAML file so the execution contract stays inspectable, while canonical docs/registers record the queue and downstream queued runs before the sidecar is launched.

**Tech Stack:** Python, bash, PyTorch runtime conventions, YAML configs, JSON/JSONL state files, pytest, existing `transchrombp.training.train_ddp`, existing `transchrombp.evaluation.model_teacher_cache_export`

**Execution Location:** Implement in the `6000` runtime repo `/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411` and backfill the same files into the local archive paths shown below. For local archive references, map `vendor/transchrombp/transchrombp/` to remote `src/transchrombp/`, and map `vendor/transchrombp/scripts/` / `vendor/transchrombp/configs/` to remote `scripts/` / `configs/`.

---

## File Structure

**Create:**

- `vendor/transchrombp/configs/queues/factor_ladder_unattended_20260411.yaml`
- `vendor/transchrombp/transchrombp/orchestration/__init__.py`
- `vendor/transchrombp/transchrombp/orchestration/factor_ladder_unattended_queue.py`
- `vendor/transchrombp/scripts/run_factor_ladder_unattended_matrix.sh`
- `tests/test_factor_ladder_unattended_queue.py`
- `tests/test_factor_ladder_unattended_launcher.py`

**Modify:**

- `tests/test_factor_ladder_docs.py`
- `TRACKING.md`
- `docs/experiments/registry.md`
- `docs/experiments/runs.csv`

**Responsibility Split:**

- `configs/queues/factor_ladder_unattended_20260411.yaml` is the explicit queue table: stage ids, run names, `model/train/data` config paths, dual-card launch parameters, completion probes, retry ports, and state-dir defaults.
- `transchrombp/orchestration/factor_ladder_unattended_queue.py` owns queue semantics only: load config, probe the live `E2`, launch subsequent stages, choose teacher checkpoints, validate exported manifests, and emit `queue_state.json`, `events.jsonl`, `summary.md`.
- `scripts/run_factor_ladder_unattended_matrix.sh` is a thin runtime wrapper that sets environment, points at the queue config/state dir, and execs the Python queue module.
- Root `tests/` cover queue config structure, helper semantics, wrapper wiring, and canonical doc registration of the unattended queue plus downstream queued runs.
- `TRACKING.md`, `registry.md`, and `runs.csv` stay canonical for queue registration and queued/run-state bookkeeping.

### Task 1: Scaffold the Queue Config and Its Structural Test

**Files:**

- Create: `vendor/transchrombp/configs/queues/factor_ladder_unattended_20260411.yaml`
- Create: `tests/test_factor_ladder_unattended_queue.py`

- [ ] **Step 1: Write the failing queue-config structure test**

```python
from __future__ import annotations

from pathlib import Path

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ARCHIVE_ROOT = PROJECT_ROOT / "vendor" / "transchrombp"
QUEUE_CONFIG = ARCHIVE_ROOT / "configs" / "queues" / "factor_ladder_unattended_20260411.yaml"


def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def test_factor_ladder_unattended_queue_config_declares_expected_stage_order() -> None:
    cfg = load_yaml(QUEUE_CONFIG)

    assert cfg["queue_id"] == "factor_ladder_unattended_20260411"
    assert cfg["runtime_root"] == "/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411"
    assert cfg["state_dir"] == "outputs/queue/factor_ladder_unattended_20260411"
    assert cfg["default_master_ports"] == [29982, 29983]

    stages = cfg["stages"]
    assert [stage["stage_id"] for stage in stages] == [
        "wait_existing_e2_short10",
        "e1_longctx4096_short10",
        "e2_hierdec4096_teacher30",
        "export_e2_teacher_cache",
        "e3_hierdec4096_distill_short10",
    ]
    assert [stage["stage_type"] for stage in stages] == [
        "wait_existing_run",
        "launch_run",
        "launch_run",
        "export_teacher_cache",
        "launch_run",
    ]

    wait_stage = stages[0]
    assert wait_stage["run_name"] == "teacher_v2_hierdec4096_short10_s42_6000_20260411_r1"
    assert wait_stage["completion"]["require_run_meta"] is True
    assert wait_stage["completion"]["checkpoint_policy"] == "best_or_latest_epoch"

    e1_stage = stages[1]
    assert e1_stage["run_name"] == "teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1"
    assert e1_stage["model_config"] == "configs/model/transchrombp_teacher_v2_center_pool_longctx4096.yaml"
    assert e1_stage["train_config"] == "configs/train/train_tutorial_teacher_v2_longctx4096_short10.yaml"
    assert e1_stage["data_config"] == "configs/data/data_tutorial_canonical_v1_longctx4096.yaml"
    assert e1_stage["train_gpu_ids"] == "0,1"
    assert e1_stage["nproc_per_node"] == 2

    e2_teacher_stage = stages[2]
    assert e2_teacher_stage["run_name"] == "teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1"
    assert e2_teacher_stage["model_config"] == "configs/model/transchrombp_teacher_v2_hierdec4096.yaml"
    assert e2_teacher_stage["train_config"] == "configs/train/train_tutorial_teacher_v2_hierdec4096_teacher30.yaml"
    assert e2_teacher_stage["train_gpu_ids"] == "0,1"
    assert e2_teacher_stage["nproc_per_node"] == 2

    export_stage = stages[3]
    assert export_stage["checkpoint_from_stage"] == "e2_hierdec4096_teacher30"
    assert export_stage["teacher_cache_dir"] == "outputs/teacher_cache/teacher_v2_hierdec4096_teacher30"

    distill_stage = stages[4]
    assert distill_stage["run_name"] == "teacher_v2_hierdec4096_distill_short10_s42_6000_20260411_r1"
    assert distill_stage["model_config"] == "configs/model/transchrombp_teacher_v2_hierdec4096.yaml"
    assert distill_stage["train_config"] == "configs/train/train_tutorial_teacher_v2_hierdec4096_distill_short10.yaml"
    assert distill_stage["data_config"] == "configs/data/data_tutorial_canonical_v1_longctx4096.yaml"
    assert distill_stage["train_gpu_ids"] == "0,1"
    assert distill_stage["nproc_per_node"] == 2
    assert distill_stage["env"]["TEACHER_CACHE_DIR"] == "outputs/teacher_cache/teacher_v2_hierdec4096_teacher30"
```

- [ ] **Step 2: Run the test to verify it fails**

Run:

```bash
pytest tests/test_factor_ladder_unattended_queue.py::test_factor_ladder_unattended_queue_config_declares_expected_stage_order -q
```

Expected:

```text
E   FileNotFoundError: [Errno 2] No such file or directory: '/repo/vendor/transchrombp/configs/queues/factor_ladder_unattended_20260411.yaml'
```

- [ ] **Step 3: Add the unattended queue YAML**

```yaml
# vendor/transchrombp/configs/queues/factor_ladder_unattended_20260411.yaml
queue_id: factor_ladder_unattended_20260411
runtime_root: /data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411
state_dir: outputs/queue/factor_ladder_unattended_20260411
default_master_ports: [29982, 29983]
poll_interval_sec: 30
gpu_wait_interval_sec: 60
gpu_idle_used_mem_threshold_mib: 4096

stages:
  - stage_id: wait_existing_e2_short10
    stage_type: wait_existing_run
    run_name: teacher_v2_hierdec4096_short10_s42_6000_20260411_r1
    log_path: logs/teacher_v2_hierdec4096_short10_s42_6000_20260411_r1.log
    completion:
      require_run_meta: true
      checkpoint_policy: best_or_latest_epoch
      completion_markers: ["[metrics] wrote", "[meta] wrote"]

  - stage_id: e1_longctx4096_short10
    stage_type: launch_run
    run_name: teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1
    model_config: configs/model/transchrombp_teacher_v2_center_pool_longctx4096.yaml
    train_config: configs/train/train_tutorial_teacher_v2_longctx4096_short10.yaml
    data_config: configs/data/data_tutorial_canonical_v1_longctx4096.yaml
    output_dir: outputs
    log_path: logs/teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1.log
    train_gpu_ids: "0,1"
    nproc_per_node: 2
    master_ports: [29982, 29983]
    env:
      RUN_NAME: teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1
    completion:
      require_run_meta: true
      checkpoint_policy: best_or_latest_epoch

  - stage_id: e2_hierdec4096_teacher30
    stage_type: launch_run
    run_name: teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1
    model_config: configs/model/transchrombp_teacher_v2_hierdec4096.yaml
    train_config: configs/train/train_tutorial_teacher_v2_hierdec4096_teacher30.yaml
    data_config: configs/data/data_tutorial_canonical_v1_longctx4096.yaml
    output_dir: outputs
    log_path: logs/teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1.log
    train_gpu_ids: "0,1"
    nproc_per_node: 2
    master_ports: [29982, 29983]
    env:
      RUN_NAME: teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1
    completion:
      require_run_meta: true
      checkpoint_policy: best_or_latest_epoch

  - stage_id: export_e2_teacher_cache
    stage_type: export_teacher_cache
    checkpoint_from_stage: e2_hierdec4096_teacher30
    teacher_cache_dir: outputs/teacher_cache/teacher_v2_hierdec4096_teacher30
    data_config: configs/data/data_tutorial_canonical_v1_longctx4096.yaml
    targets: ["profile16", "logcount"]
    splits: ["train", "valid"]

  - stage_id: e3_hierdec4096_distill_short10
    stage_type: launch_run
    run_name: teacher_v2_hierdec4096_distill_short10_s42_6000_20260411_r1
    model_config: configs/model/transchrombp_teacher_v2_hierdec4096.yaml
    train_config: configs/train/train_tutorial_teacher_v2_hierdec4096_distill_short10.yaml
    data_config: configs/data/data_tutorial_canonical_v1_longctx4096.yaml
    output_dir: outputs
    log_path: logs/teacher_v2_hierdec4096_distill_short10_s42_6000_20260411_r1.log
    train_gpu_ids: "0,1"
    nproc_per_node: 2
    master_ports: [29982, 29983]
    env:
      RUN_NAME: teacher_v2_hierdec4096_distill_short10_s42_6000_20260411_r1
      TEACHER_CACHE_DIR: outputs/teacher_cache/teacher_v2_hierdec4096_teacher30
    completion:
      require_run_meta: true
      checkpoint_policy: best_or_latest_epoch
```

- [ ] **Step 4: Re-run the structural test**

Run:

```bash
pytest tests/test_factor_ladder_unattended_queue.py::test_factor_ladder_unattended_queue_config_declares_expected_stage_order -q
```

Expected:

```text
1 passed
```

- [ ] **Step 5: Commit the queue-config scaffold**

```bash
git add vendor/transchrombp/configs/queues/factor_ladder_unattended_20260411.yaml tests/test_factor_ladder_unattended_queue.py
git commit -m "add factor ladder unattended queue config"
```

### Task 2: Implement Queue Parsing, Completion Probes, and State Writers

**Files:**

- Create: `vendor/transchrombp/transchrombp/orchestration/__init__.py`
- Create: `vendor/transchrombp/transchrombp/orchestration/factor_ladder_unattended_queue.py`
- Modify: `tests/test_factor_ladder_unattended_queue.py`

- [ ] **Step 1: Extend the queue test with failing helper-level expectations**

```python
from transchrombp.orchestration.factor_ladder_unattended_queue import (
    append_event,
    choose_checkpoint_path,
    load_queue_spec,
    render_summary_markdown,
    run_completion_probe,
    write_queue_state,
)


def test_choose_checkpoint_prefers_best_over_latest_epoch(tmp_path: Path) -> None:
    ckpt_dir = tmp_path / "checkpoints"
    ckpt_dir.mkdir()
    (ckpt_dir / "epoch_010.pt").write_text("epoch10", encoding="utf-8")
    (ckpt_dir / "epoch_030.pt").write_text("epoch30", encoding="utf-8")
    (ckpt_dir / "best.pt").write_text("best", encoding="utf-8")

    selected = choose_checkpoint_path(ckpt_dir, policy="best_or_latest_epoch")
    assert selected == ckpt_dir / "best.pt"


def test_run_completion_probe_requires_meta_and_checkpoint(tmp_path: Path) -> None:
    output_root = tmp_path / "outputs"
    run_name = "teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1"
    log_path = tmp_path / "logs" / f"{run_name}.log"
    log_path.parent.mkdir(parents=True)
    log_path.write_text("[metrics] wrote x\n[meta] wrote y\n", encoding="utf-8")

    meta_dir = output_root / "logs" / run_name
    meta_dir.mkdir(parents=True)
    (meta_dir / "run_meta.json").write_text("{}", encoding="utf-8")

    ckpt_dir = output_root / "checkpoints" / run_name
    ckpt_dir.mkdir(parents=True)
    (ckpt_dir / "epoch_010.pt").write_text("epoch", encoding="utf-8")

    probe = run_completion_probe(
        runtime_root=tmp_path,
        output_root=output_root,
        run_name=run_name,
        log_path=log_path,
        require_run_meta=True,
        checkpoint_policy="best_or_latest_epoch",
        completion_markers=["[metrics] wrote", "[meta] wrote"],
    )

    assert probe.status == "completed"
    assert probe.checkpoint_path == ckpt_dir / "epoch_010.pt"


def test_queue_state_and_summary_include_current_stage(tmp_path: Path) -> None:
    state_dir = tmp_path / "queue"
    state_path = write_queue_state(
        state_dir=state_dir,
        queue_id="factor_ladder_unattended_20260411",
        status="running",
        current_stage_id="e1_longctx4096_short10",
        current_run_name="teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1",
        current_log_path="logs/teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1.log",
    )
    append_event(
        state_dir=state_dir,
        event_type="stage_started",
        stage_id="e1_longctx4096_short10",
        payload={"run_name": "teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1"},
    )
    summary = render_summary_markdown(
        queue_id="factor_ladder_unattended_20260411",
        status="running",
        current_stage_id="e1_longctx4096_short10",
        events_path=state_dir / "events.jsonl",
    )

    assert state_path.exists()
    assert "e1_longctx4096_short10" in summary
    assert "stage_started" in summary
```

- [ ] **Step 2: Run the expanded queue test file and verify it fails on missing module imports**

Run:

```bash
PATH=./.venv/bin:$PATH PYTHONPATH=vendor/transchrombp:${PYTHONPATH:-} pytest tests/test_factor_ladder_unattended_queue.py -q
```

Expected:

```text
E   ModuleNotFoundError: No module named 'transchrombp.orchestration'
```

- [ ] **Step 3: Create the orchestration package and queue-spec dataclasses**

```python
# vendor/transchrombp/transchrombp/orchestration/__init__.py
from .factor_ladder_unattended_queue import (
    QueueSpec,
    StageCompletionProbe,
    StageSpec,
    append_event,
    choose_checkpoint_path,
    load_queue_spec,
    render_summary_markdown,
    run_completion_probe,
    write_queue_state,
)

__all__ = [
    "QueueSpec",
    "StageCompletionProbe",
    "StageSpec",
    "append_event",
    "choose_checkpoint_path",
    "load_queue_spec",
    "render_summary_markdown",
    "run_completion_probe",
    "write_queue_state",
]
```

```python
# vendor/transchrombp/transchrombp/orchestration/factor_ladder_unattended_queue.py
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import json
import yaml


@dataclass(frozen=True)
class StageSpec:
    stage_id: str
    stage_type: str
    run_name: str = ""
    model_config: str = ""
    train_config: str = ""
    data_config: str = ""
    output_dir: str = ""
    log_path: str = ""
    env: dict[str, str] = field(default_factory=dict)
    completion: dict[str, Any] = field(default_factory=dict)
    master_ports: list[int] = field(default_factory=list)
    train_gpu_ids: str = ""
    nproc_per_node: int = 1
    checkpoint_from_stage: str = ""
    teacher_cache_dir: str = ""
    splits: list[str] = field(default_factory=list)
    targets: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class QueueSpec:
    queue_id: str
    runtime_root: Path
    state_dir: Path
    default_master_ports: list[int]
    poll_interval_sec: int
    gpu_wait_interval_sec: int
    gpu_idle_used_mem_threshold_mib: int
    stages: list[StageSpec]


@dataclass(frozen=True)
class StageCompletionProbe:
    status: str
    checkpoint_path: Path | None
    reason: str


def load_queue_spec(path: Path) -> QueueSpec:
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    return QueueSpec(
        queue_id=str(raw["queue_id"]),
        runtime_root=Path(str(raw["runtime_root"])).resolve(),
        state_dir=Path(str(raw["state_dir"])),
        default_master_ports=[int(port) for port in raw.get("default_master_ports", [])],
        poll_interval_sec=int(raw.get("poll_interval_sec", 30)),
        gpu_wait_interval_sec=int(raw.get("gpu_wait_interval_sec", 60)),
        gpu_idle_used_mem_threshold_mib=int(raw.get("gpu_idle_used_mem_threshold_mib", 4096)),
        stages=[StageSpec(**stage) for stage in raw.get("stages", [])],
    )
```

- [ ] **Step 4: Add checkpoint selection, completion probes, and state/event writers**

```python
def choose_checkpoint_path(checkpoint_dir: Path, policy: str) -> Path | None:
    if policy != "best_or_latest_epoch":
        raise ValueError(f"Unsupported checkpoint policy: {policy}")
    best_path = checkpoint_dir / "best.pt"
    if best_path.exists():
        return best_path
    epoch_paths = sorted(checkpoint_dir.glob("epoch_*.pt"))
    return epoch_paths[-1] if epoch_paths else None


def run_completion_probe(
    *,
    runtime_root: Path,
    output_root: Path,
    run_name: str,
    log_path: Path,
    require_run_meta: bool,
    checkpoint_policy: str,
    completion_markers: list[str],
) -> StageCompletionProbe:
    if not log_path.exists():
        return StageCompletionProbe(status="running", checkpoint_path=None, reason="log_missing")

    log_text = log_path.read_text(encoding="utf-8", errors="replace")
    if not all(marker in log_text for marker in completion_markers):
        return StageCompletionProbe(status="running", checkpoint_path=None, reason="markers_missing")

    meta_path = output_root / "logs" / run_name / "run_meta.json"
    if require_run_meta and not meta_path.exists():
        return StageCompletionProbe(status="failed", checkpoint_path=None, reason="run_meta_missing")

    checkpoint_dir = output_root / "checkpoints" / run_name
    checkpoint_path = choose_checkpoint_path(checkpoint_dir, checkpoint_policy)
    if checkpoint_path is None:
        return StageCompletionProbe(status="failed", checkpoint_path=None, reason="checkpoint_missing")

    return StageCompletionProbe(status="completed", checkpoint_path=checkpoint_path, reason="ok")


def write_queue_state(
    *,
    state_dir: Path,
    queue_id: str,
    status: str,
    current_stage_id: str,
    current_run_name: str,
    current_log_path: str,
) -> Path:
    state_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "queue_id": queue_id,
        "status": status,
        "current_stage_id": current_stage_id,
        "current_run_name": current_run_name,
        "current_log_path": current_log_path,
    }
    out_path = state_dir / "queue_state.json"
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return out_path


def append_event(*, state_dir: Path, event_type: str, stage_id: str, payload: dict[str, Any]) -> Path:
    state_dir.mkdir(parents=True, exist_ok=True)
    event_path = state_dir / "events.jsonl"
    line = {"event_type": event_type, "stage_id": stage_id, **payload}
    with event_path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(line, sort_keys=True) + "\n")
    return event_path


def render_summary_markdown(*, queue_id: str, status: str, current_stage_id: str, events_path: Path) -> str:
    events = []
    if events_path.exists():
        events = [json.loads(line) for line in events_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    lines = [
        f"# {queue_id} summary",
        "",
        f"- status: `{status}`",
        f"- current_stage_id: `{current_stage_id}`",
        "",
        "## Recent Events",
    ]
    for event in events[-5:]:
        lines.append(f"- `{event['event_type']}` `{event['stage_id']}`")
    return "\n".join(lines) + "\n"
```

- [ ] **Step 5: Run queue helper tests and verify they pass**

Run:

```bash
PATH=./.venv/bin:$PATH PYTHONPATH=vendor/transchrombp:${PYTHONPATH:-} pytest tests/test_factor_ladder_unattended_queue.py -q
```

Expected:

```text
4 passed
```

- [ ] **Step 6: Commit the orchestration core**

```bash
git add \
  vendor/transchrombp/transchrombp/orchestration/__init__.py \
  vendor/transchrombp/transchrombp/orchestration/factor_ladder_unattended_queue.py \
  tests/test_factor_ladder_unattended_queue.py
git commit -m "add factor ladder unattended queue core"
```

### Task 3: Add the CLI Loop and the Thin Shell Wrapper

**Files:**

- Modify: `vendor/transchrombp/transchrombp/orchestration/factor_ladder_unattended_queue.py`
- Create: `vendor/transchrombp/scripts/run_factor_ladder_unattended_matrix.sh`
- Create: `tests/test_factor_ladder_unattended_launcher.py`

- [ ] **Step 1: Write the failing wrapper/dry-run launcher test**

```python
from __future__ import annotations

import json
import os
import stat
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ARCHIVE_ROOT = PROJECT_ROOT / "vendor" / "transchrombp"
LAUNCHER = ARCHIVE_ROOT / "scripts" / "run_factor_ladder_unattended_matrix.sh"


def write_stub_python(bin_dir: Path, capture_path: Path) -> None:
    stub = bin_dir / "python"
    stub.write_text(
        "\n".join(
            [
                "#!/usr/bin/env python3",
                "import json",
                "from pathlib import Path",
                "import sys",
                f'capture_path = r\"{capture_path}\"',
                "path = Path(capture_path)",
                "calls = []",
                "if path.exists():",
                "    calls = json.loads(path.read_text(encoding='utf-8'))",
                "calls.append(sys.argv[1:])",
                "path.write_text(json.dumps(calls), encoding='utf-8')",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    stub.chmod(stub.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def test_unattended_launcher_binds_queue_defaults_and_execs_python_module(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    capture_path = tmp_path / "argv.json"
    write_stub_python(bin_dir, capture_path)

    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}{os.pathsep}{env.get('PATH', '')}"
    env["VENV_DIR"] = str(tmp_path / "missing-venv")

    subprocess.run(["bash", "-n", str(LAUNCHER)], cwd=ARCHIVE_ROOT, check=True)
    subprocess.run(["bash", str(LAUNCHER), "--dry-run", "--once"], cwd=ARCHIVE_ROOT, env=env, check=True)

    calls = json.loads(capture_path.read_text(encoding="utf-8"))
    assert len(calls) == 1
    argv = calls[0]
    assert argv[:2] == ["-m", "transchrombp.orchestration.factor_ladder_unattended_queue"]
    assert argv[2:] == [
        "--queue-config",
        str(ARCHIVE_ROOT / "configs" / "queues" / "factor_ladder_unattended_20260411.yaml"),
        "--state-dir",
        str(ARCHIVE_ROOT / "outputs" / "queue" / "factor_ladder_unattended_20260411"),
        "--dry-run",
        "--once",
    ]
```

- [ ] **Step 2: Run the launcher test and verify it fails because the wrapper does not exist**

Run:

```bash
PATH=./.venv/bin:$PATH PYTHONPATH=vendor/transchrombp:${PYTHONPATH:-} pytest tests/test_factor_ladder_unattended_launcher.py -q
```

Expected:

```text
E   FileNotFoundError: [Errno 2] No such file or directory: '/repo/vendor/transchrombp/scripts/run_factor_ladder_unattended_matrix.sh'
```

- [ ] **Step 3: Add the queue CLI loop, GPU wait, and port-fallback logic**

```python
import argparse
import json
import os
import subprocess
import time


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the factor-ladder unattended matrix queue.")
    parser.add_argument("--queue-config", required=True)
    parser.add_argument("--state-dir", required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--once", action="store_true")
    return parser.parse_args()


def is_port_conflict(log_text: str) -> bool:
    needles = ("Address already in use", "EADDRINUSE", "errno 98")
    return any(needle in log_text for needle in needles)


def run_stage_command(command: list[str], stage_log_path: Path, env: dict[str, str]) -> tuple[int, str]:
    stage_log_path.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(command, env=env, text=True, capture_output=True, check=False)
    with stage_log_path.open("a", encoding="utf-8") as handle:
        handle.write(proc.stdout)
        handle.write(proc.stderr)
    combined = proc.stdout + proc.stderr
    return proc.returncode, combined


def wait_for_idle_gpus(spec: QueueSpec) -> None:
    while True:
        query = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=memory.used",
                "--format=csv,noheader,nounits",
            ],
            text=True,
            capture_output=True,
            check=True,
        )
        used = [int(line.strip()) for line in query.stdout.splitlines() if line.strip()]
        if used and all(value <= spec.gpu_idle_used_mem_threshold_mib for value in used):
            return
        time.sleep(spec.gpu_wait_interval_sec)


def wait_for_existing_stage_completion(spec: QueueSpec, stage: StageSpec) -> StageCompletionProbe:
    output_root = spec.runtime_root / "outputs"
    log_path = spec.runtime_root / stage.log_path
    completion_cfg = dict(stage.completion)
    while True:
        probe = run_completion_probe(
            runtime_root=spec.runtime_root,
            output_root=output_root,
            run_name=stage.run_name,
            log_path=log_path,
            require_run_meta=bool(completion_cfg.get("require_run_meta", True)),
            checkpoint_policy=str(completion_cfg.get("checkpoint_policy", "best_or_latest_epoch")),
            completion_markers=list(completion_cfg.get("completion_markers", [])),
        )
        if probe.status != "running":
            return probe
        time.sleep(spec.poll_interval_sec)


def build_launch_command(spec: QueueSpec, stage: StageSpec, master_port: int) -> tuple[list[str], dict[str, str]]:
    env = os.environ.copy()
    env.update(stage.env)
    env["CUDA_VISIBLE_DEVICES"] = stage.train_gpu_ids
    command = [
        "torchrun",
        "--standalone",
        "--nproc_per_node",
        str(stage.nproc_per_node),
        "--master_port",
        str(master_port),
        "-m",
        "transchrombp.training.train_ddp",
        "--model-config",
        str(spec.runtime_root / stage.model_config),
        "--train-config",
        str(spec.runtime_root / stage.train_config),
        "--data-config",
        str(spec.runtime_root / stage.data_config),
        "--output-dir",
        str(spec.runtime_root / stage.output_dir),
        "--run-name",
        stage.run_name,
    ]
    return command, env


def build_export_command(spec: QueueSpec, stage: StageSpec, checkpoint_path: Path) -> tuple[list[str], dict[str, str]]:
    env = os.environ.copy()
    command = [
        "python",
        "-m",
        "transchrombp.evaluation.model_teacher_cache_export",
        "--checkpoint",
        str(checkpoint_path),
        "--data-config",
        str(spec.runtime_root / stage.data_config),
        "--output-dir",
        str(spec.runtime_root / stage.teacher_cache_dir),
        "--splits",
        *stage.splits,
    ]
    return command, env


def run_queue(args: argparse.Namespace) -> int:
    spec = load_queue_spec(Path(args.queue_config))
    state_dir = Path(args.state_dir)
    output_root = spec.runtime_root / "outputs"

    for stage in spec.stages:
        write_queue_state(
            state_dir=state_dir,
            queue_id=spec.queue_id,
            status="running",
            current_stage_id=stage.stage_id,
            current_run_name=stage.run_name,
            current_log_path=stage.log_path,
        )
        append_event(state_dir=state_dir, event_type="stage_started", stage_id=stage.stage_id, payload={})

        if stage.stage_type == "wait_existing_run":
            if args.dry_run:
                append_event(
                    state_dir=state_dir,
                    event_type="stage_dry_run",
                    stage_id=stage.stage_id,
                    payload={"mode": "wait_existing_run"},
                )
                if args.once:
                    break
                continue
            probe = wait_for_existing_stage_completion(spec, stage)
            if probe.status != "completed":
                append_event(state_dir=state_dir, event_type="queue_halted", stage_id=stage.stage_id, payload={"reason": probe.reason})
                return 1
            append_event(
                state_dir=state_dir,
                event_type="stage_passed",
                stage_id=stage.stage_id,
                payload={"checkpoint_path": str(probe.checkpoint_path) if probe.checkpoint_path else ""},
            )
            continue

        if stage.stage_type == "launch_run":
            wait_for_idle_gpus(spec)
            stage_log_path = spec.runtime_root / stage.log_path
            ports = stage.master_ports or spec.default_master_ports
            last_log = ""
            for index, port in enumerate(ports):
                command, env = build_launch_command(spec, stage, port)
                if args.dry_run:
                    append_event(
                        state_dir=state_dir,
                        event_type="stage_dry_run",
                        stage_id=stage.stage_id,
                        payload={"command": command, "master_port": port},
                    )
                    last_log = json.dumps({"command": command, "master_port": port})
                    break
                return_code, last_log = run_stage_command(command, stage_log_path, env)
                if return_code == 0:
                    completion_cfg = dict(stage.completion)
                    probe = run_completion_probe(
                        runtime_root=spec.runtime_root,
                        output_root=output_root,
                        run_name=stage.run_name,
                        log_path=stage_log_path,
                        require_run_meta=bool(completion_cfg.get("require_run_meta", True)),
                        checkpoint_policy=str(completion_cfg.get("checkpoint_policy", "best_or_latest_epoch")),
                        completion_markers=list(completion_cfg.get("completion_markers", [])),
                    )
                    if probe.status != "completed":
                        append_event(
                            state_dir=state_dir,
                            event_type="queue_halted",
                            stage_id=stage.stage_id,
                            payload={"reason": probe.reason},
                        )
                        return 1
                    break
                if index == 0 and is_port_conflict(last_log) and len(ports) > 1:
                    append_event(
                        state_dir=state_dir,
                        event_type="stage_retry_port_conflict",
                        stage_id=stage.stage_id,
                        payload={"failed_port": port, "next_port": ports[1]},
                    )
                    continue
                append_event(
                    state_dir=state_dir,
                    event_type="queue_halted",
                    stage_id=stage.stage_id,
                    payload={"reason": "launch_failed", "log_excerpt": last_log[-400:]},
                )
                return 1
            append_event(state_dir=state_dir, event_type="stage_passed", stage_id=stage.stage_id, payload={"run_name": stage.run_name})

        if stage.stage_type == "export_teacher_cache":
            teacher_stage = next(item for item in spec.stages if item.stage_id == stage.checkpoint_from_stage)
            teacher_checkpoint = choose_checkpoint_path(
                output_root / "checkpoints" / teacher_stage.run_name,
                policy="best_or_latest_epoch",
            )
            if teacher_checkpoint is None:
                append_event(state_dir=state_dir, event_type="queue_halted", stage_id=stage.stage_id, payload={"reason": "teacher_checkpoint_missing"})
                return 1
            command, env = build_export_command(spec, stage, teacher_checkpoint)
            if not args.dry_run:
                return_code, combined_log = run_stage_command(
                    command,
                    spec.runtime_root / "logs" / f"{stage.stage_id}.log",
                    env,
                )
                if return_code != 0:
                    append_event(state_dir=state_dir, event_type="queue_halted", stage_id=stage.stage_id, payload={"reason": "export_failed", "log_excerpt": combined_log[-400:]})
                    return 1
                teacher_cache_root = spec.runtime_root / stage.teacher_cache_dir
                required_manifests = [
                    teacher_cache_root / "teacher_manifest_train.json",
                    teacher_cache_root / "teacher_manifest_valid.json",
                ]
                if not all(path.exists() for path in required_manifests):
                    append_event(state_dir=state_dir, event_type="queue_halted", stage_id=stage.stage_id, payload={"reason": "teacher_manifest_missing"})
                    return 1
            append_event(
                state_dir=state_dir,
                event_type="stage_passed",
                stage_id=stage.stage_id,
                payload={"teacher_cache_dir": str(spec.runtime_root / stage.teacher_cache_dir)},
            )

        if args.once:
            break

    write_queue_state(
        state_dir=state_dir,
        queue_id=spec.queue_id,
        status="completed",
        current_stage_id="",
        current_run_name="",
        current_log_path="",
    )
    return 0
```

- [ ] **Step 4: Add the thin shell wrapper**

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QUEUE_CONFIG="${QUEUE_CONFIG:-${ROOT_DIR}/configs/queues/factor_ladder_unattended_20260411.yaml}"
STATE_DIR="${STATE_DIR:-${ROOT_DIR}/outputs/queue/factor_ladder_unattended_20260411}"
VENV_DIR="${VENV_DIR:-/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b}"

export PATH="${VENV_DIR}/bin:${PATH:-}"
export LD_LIBRARY_PATH="${VENV_DIR}/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${ROOT_DIR}:${PYTHONPATH:-}"

exec python -m transchrombp.orchestration.factor_ladder_unattended_queue \
  --queue-config "${QUEUE_CONFIG}" \
  --state-dir "${STATE_DIR}" \
  "$@"
```

- [ ] **Step 5: Run the queue-core and launcher tests together**

Run:

```bash
PATH=./.venv/bin:$PATH PYTHONPATH=vendor/transchrombp:${PYTHONPATH:-} \
  pytest tests/test_factor_ladder_unattended_queue.py tests/test_factor_ladder_unattended_launcher.py -q
```

Expected:

```text
5 passed
```

- [ ] **Step 6: Commit the queue launcher**

```bash
git add \
  vendor/transchrombp/transchrombp/orchestration/factor_ladder_unattended_queue.py \
  vendor/transchrombp/scripts/run_factor_ladder_unattended_matrix.sh \
  tests/test_factor_ladder_unattended_launcher.py
git commit -m "add factor ladder unattended launcher"
```

### Task 4: Register the Queue and Downstream Queued Runs in Canonical Docs

**Files:**

- Modify: `TRACKING.md`
- Modify: `docs/experiments/registry.md`
- Modify: `docs/experiments/runs.csv`
- Modify: `tests/test_factor_ladder_docs.py`

- [ ] **Step 1: Extend the docs regression with failing queue expectations**

```python
import csv


def parse_runs_csv() -> list[dict[str, str]]:
    with (REPO_ROOT / "docs/experiments/runs.csv").open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def test_tracking_mentions_unattended_queue_registration() -> None:
    rows = parse_markdown_table("TRACKING.md")
    factor_row = next(row for row in rows if row["事项"] == "AlphaGenome-like factor ladder（E1/E2/E3）")

    assert "unattended matrix queue" in factor_row["当前结论 / 进度"]
    assert "outputs/queue/factor_ladder_unattended_20260411" in factor_row["当前结论 / 进度"]
    assert "teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1" in factor_row["下一步"]
    assert "teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1" in factor_row["下一步"]
    assert "teacher_v2_hierdec4096_distill_short10_s42_6000_20260411_r1" in factor_row["下一步"]


def test_runs_manifest_pre_registers_unattended_queue_downstream_runs() -> None:
    rows = parse_runs_csv()
    queued = {row["run_id"]: row for row in rows if row["family_id"] == "alphagenome_factor_ladder" and row["run_status"] == "queued"}

    assert queued["teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1"]["commit_sha"] == "079b603"
    assert queued["teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1"]["commit_sha"] == "079b603"
    assert queued["teacher_v2_hierdec4096_distill_short10_s42_6000_20260411_r1"]["commit_sha"] == "079b603"
```

- [ ] **Step 2: Run the docs regression and verify it fails**

Run:

```bash
PATH=./.venv/bin:$PATH PYTHONPATH=vendor/transchrombp:${PYTHONPATH:-} pytest tests/test_factor_ladder_docs.py -q
```

Expected:

```text
E   AssertionError: assert 'unattended matrix queue' in factor_row['当前结论 / 进度']
```

- [ ] **Step 3: Update canonical docs and queued-run rows**

```markdown
# TRACKING.md factor-ladder row additions
- 已登记 unattended matrix queue：`outputs/queue/factor_ladder_unattended_20260411`
- 当前 active stage 仍为 `teacher_v2_hierdec4096_short10_s42_6000_20260411_r1`
- downstream queued runs:
  - `teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1`
  - `teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1`
  - `teacher_v2_hierdec4096_distill_short10_s42_6000_20260411_r1`
```

```markdown
# docs/experiments/registry.md factor-ladder row changes
- status 保持 `running`
- next_allowed_action 改成：
  当前由 unattended matrix queue 接管后续 stage；等待现有 `E2 short10` 完成后，顺序触发 `E1 short10 -> E2 teacher30 -> export -> E3 distill short10`
- notes 增加：
  queue state dir=`outputs/queue/factor_ladder_unattended_20260411`
```

```csv
teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1,alphagenome_factor_ladder,6000,A6000x2,main,/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411,079b603,longctx4096_short10_formal;queued_by_unattended_matrix,/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411,/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411/logs/teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1.log,/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411/outputs,unknown,unknown,queued,n/a,n/a,n/a,queued by unattended matrix against canonical commit 079b603
teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1,alphagenome_factor_ladder,6000,A6000x2,main,/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411,079b603,hierdec4096_teacher30;queued_by_unattended_matrix,/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411,/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411/logs/teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1.log,/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411/outputs,unknown,unknown,queued,n/a,n/a,n/a,queued by unattended matrix against canonical commit 079b603
teacher_v2_hierdec4096_distill_short10_s42_6000_20260411_r1,alphagenome_factor_ladder,6000,A6000x2,main,/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411,079b603,hierdec4096_distill_short10;queued_by_unattended_matrix,/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411,/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411/logs/teacher_v2_hierdec4096_distill_short10_s42_6000_20260411_r1.log,/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411/outputs,unknown,unknown,queued,n/a,n/a,n/a,queued by unattended matrix against canonical commit 079b603
```

- [ ] **Step 4: Run docs regression after the updates**

Run:

```bash
PATH=./.venv/bin:$PATH PYTHONPATH=vendor/transchrombp:${PYTHONPATH:-} pytest tests/test_factor_ladder_docs.py -q
```

Expected:

```text
4 passed
```

- [ ] **Step 5: Commit the canonical registration**

```bash
git add TRACKING.md docs/experiments/registry.md docs/experiments/runs.csv tests/test_factor_ladder_docs.py
git commit -m "register factor ladder unattended queue"
```

### Task 5: Sync to 6000, Run a Queue Dry-Run Smoke, Then Launch the Sidecar

**Files:**

- Remote mirror of: `configs/queues/factor_ladder_unattended_20260411.yaml`
- Remote mirror of: `src/transchrombp/orchestration/__init__.py`
- Remote mirror of: `src/transchrombp/orchestration/factor_ladder_unattended_queue.py`
- Remote mirror of: `scripts/run_factor_ladder_unattended_matrix.sh`
- Modify after live launch: `TRACKING.md`
- Modify after live launch: `docs/experiments/registry.md`
- Modify after live launch: `docs/experiments/runs.csv`

- [ ] **Step 1: Sync the new queue files into the `6000` isolated runtime repo and syntax-check them**

Run:

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 'mkdir -p /data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411/configs/queues /data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411/src/transchrombp/orchestration'
scp -P 6000 vendor/transchrombp/configs/queues/factor_ladder_unattended_20260411.yaml zhoujiazhen@127.0.0.1:/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411/configs/queues/
scp -P 6000 vendor/transchrombp/transchrombp/orchestration/__init__.py zhoujiazhen@127.0.0.1:/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411/src/transchrombp/orchestration/
scp -P 6000 vendor/transchrombp/transchrombp/orchestration/factor_ladder_unattended_queue.py zhoujiazhen@127.0.0.1:/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411/src/transchrombp/orchestration/
scp -P 6000 vendor/transchrombp/scripts/run_factor_ladder_unattended_matrix.sh zhoujiazhen@127.0.0.1:/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411/scripts/
ssh -p 6000 zhoujiazhen@127.0.0.1 'cd /data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411 && bash -n scripts/run_factor_ladder_unattended_matrix.sh && python3 -m py_compile src/transchrombp/orchestration/factor_ladder_unattended_queue.py'
```

Expected:

```text
# no syntax errors
```

- [ ] **Step 2: Run local and remote queue-related tests**

Run:

```bash
PATH=./.venv/bin:$PATH PYTHONPATH=vendor/transchrombp:${PYTHONPATH:-} \
  pytest tests/test_factor_ladder_unattended_queue.py tests/test_factor_ladder_unattended_launcher.py tests/test_factor_ladder_docs.py -q

ssh -p 6000 zhoujiazhen@127.0.0.1 '
  cd /data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411 &&
  export PATH=/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b/bin:$PATH &&
  export PYTHONPATH=/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411/src:$PYTHONPATH &&
  python -m transchrombp.orchestration.factor_ladder_unattended_queue \
    --queue-config configs/queues/factor_ladder_unattended_20260411.yaml \
    --state-dir outputs/queue/factor_ladder_unattended_20260411_smoke \
    --dry-run --once
'
```

Expected:

```text
9 passed
# remote dry-run exits 0 and writes queue_state.json / events.jsonl / summary.md
```

- [ ] **Step 3: Verify the dry-run smoke artifacts**

Run:

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 '
  cd /data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411 &&
  wc -l outputs/queue/factor_ladder_unattended_20260411_smoke/events.jsonl &&
  sed -n "1,80p" outputs/queue/factor_ladder_unattended_20260411_smoke/queue_state.json &&
  sed -n "1,80p" outputs/queue/factor_ladder_unattended_20260411_smoke/summary.md
'
```

Expected:

```text
# summary points at wait_existing_e2_short10 or the next stage if current E2 already finished
```

- [ ] **Step 4: Launch the real unattended queue sidecar on `6000`**

Run:

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 '
  cd /data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411 &&
  nohup bash scripts/run_factor_ladder_unattended_matrix.sh > logs/factor_ladder_unattended_matrix_20260411.log 2>&1 < /dev/null &
  echo $!
'
```

Expected:

```text
<pid>
```

- [ ] **Step 5: Immediately update canonical docs with live queue status**

```markdown
# TRACKING.md / registry.md / runs.csv live-launch additions
- unattended queue log: `/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411/logs/factor_ladder_unattended_matrix_20260411.log`
- queue state dir: `/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411/outputs/queue/factor_ladder_unattended_20260411`
- 当前 active stage:
  - 若 `E2` 还在跑，则 `wait_existing_e2_short10`
  - 否则写实际已进入的下一段
```

- [ ] **Step 6: Re-run docs regression and commit the live-launch handoff**

Run:

```bash
PATH=./.venv/bin:$PATH PYTHONPATH=vendor/transchrombp:${PYTHONPATH:-} pytest tests/test_factor_ladder_docs.py -q
git add TRACKING.md docs/experiments/registry.md docs/experiments/runs.csv
git commit -m "launch factor ladder unattended queue"
```

Expected:

```text
4 passed
[master 1234abc] launch factor ladder unattended queue
```

## Self-Review Checklist

- Spec coverage:
  - Queue scope and stage order: Task 1
  - Scientific no-stop + technical fuse semantics: Tasks 2-3
  - State files and summary: Task 2
  - Thin wrapper and runtime launch: Task 3
  - Canonical queue registration and queued runs: Task 4
  - `6000` dry-run + live sidecar launch: Task 5
- Placeholder scan:
  - No red-flag placeholders remain
  - Each code-changing step includes concrete code or command blocks
- Type consistency:
  - Queue ids, stage ids, run names, queue state filenames, and config paths are consistent across tasks
