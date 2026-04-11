from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
from typing import Any

import yaml


FATAL_LOG_MARKERS = (
    "[bootstrap-fail]",
    "unbound variable",
)


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
    resolved_log_path = log_path if log_path.is_absolute() else runtime_root / log_path
    if not resolved_log_path.exists():
        return StageCompletionProbe(status="running", checkpoint_path=None, reason="log_missing")

    log_text = resolved_log_path.read_text(encoding="utf-8", errors="replace")
    if any(marker in log_text for marker in FATAL_LOG_MARKERS):
        return StageCompletionProbe(status="failed", checkpoint_path=None, reason="fatal_marker_detected")
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
