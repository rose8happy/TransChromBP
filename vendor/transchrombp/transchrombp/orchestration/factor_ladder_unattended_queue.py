from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import yaml


FATAL_LOG_MARKERS = (
    "[bootstrap-fail]",
    "unbound variable",
)

DEFAULT_QUEUE_CONFIG = Path("configs/queues/factor_ladder_unattended_20260411.yaml")
DEFAULT_STATE_DIR = Path("outputs/queue/factor_ladder_unattended_20260411")


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


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the factor ladder unattended queue")
    parser.add_argument("--queue-config", required=True)
    parser.add_argument("--state-dir", required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--once", action="store_true")
    return parser.parse_args(list(argv) if argv is not None else None)


def is_port_conflict(log_text: str) -> bool:
    lowered = log_text.lower()
    needles = ("address already in use", "eaddrinuse", "errno 98")
    return any(needle in lowered for needle in needles)


def _stage_log_path(spec: QueueSpec, stage: StageSpec) -> Path:
    return spec.runtime_root / (stage.log_path or f"logs/{stage.stage_id}.log")


def _stage_output_path(spec: QueueSpec, stage: StageSpec) -> Path:
    return spec.runtime_root / (stage.output_dir or "outputs")


def _stage_status(state_dir: Path, queue_id: str, status: str, stage: StageSpec) -> None:
    write_queue_state(
        state_dir=state_dir,
        queue_id=queue_id,
        status=status,
        current_stage_id=stage.stage_id,
        current_run_name=stage.run_name,
        current_log_path=stage.log_path or f"logs/{stage.stage_id}.log",
    )
    summary_path = state_dir / "summary.md"
    summary_path.write_text(
        render_summary_markdown(
            queue_id=queue_id,
            status=status,
            current_stage_id=stage.stage_id,
            events_path=state_dir / "events.jsonl",
        ),
        encoding="utf-8",
    )


def _finalize_queue(state_dir: Path, queue_id: str) -> None:
    write_queue_state(
        state_dir=state_dir,
        queue_id=queue_id,
        status="completed",
        current_stage_id="",
        current_run_name="",
        current_log_path="",
    )
    summary_path = state_dir / "summary.md"
    summary_path.write_text(
        render_summary_markdown(
            queue_id=queue_id,
            status="completed",
            current_stage_id="",
            events_path=state_dir / "events.jsonl",
        ),
        encoding="utf-8",
    )


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
        if not used or all(value <= spec.gpu_idle_used_mem_threshold_mib for value in used):
            return
        time.sleep(spec.gpu_wait_interval_sec)


def wait_for_existing_run_completion(spec: QueueSpec, stage: StageSpec) -> StageCompletionProbe:
    completion_cfg = dict(stage.completion)
    while True:
        probe = run_completion_probe(
            runtime_root=spec.runtime_root,
            output_root=_stage_output_path(spec, stage),
            run_name=stage.run_name,
            log_path=_stage_log_path(spec, stage),
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
    env["MASTER_ADDR"] = env.get("MASTER_ADDR", "127.0.0.1")
    env["MASTER_PORT"] = str(master_port)
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
        sys.executable,
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
    if stage.targets:
        command.extend(["--targets", *stage.targets])
    return command, env


def run_stage_command(command: list[str], stage_log_path: Path, env: dict[str, str]) -> tuple[int, str]:
    stage_log_path.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(command, env=env, text=True, capture_output=True, check=False)
    combined = (proc.stdout or "") + (proc.stderr or "")
    with stage_log_path.open("a", encoding="utf-8") as handle:
        if proc.stdout:
            handle.write(proc.stdout)
        if proc.stderr:
            handle.write(proc.stderr)
    return proc.returncode, combined


def _record_event(
    *,
    state_dir: Path,
    queue_id: str,
    status: str,
    stage: StageSpec,
    event_type: str,
    payload: dict[str, Any],
) -> None:
    append_event(state_dir=state_dir, event_type=event_type, stage_id=stage.stage_id, payload=payload)
    _stage_status(state_dir, queue_id, status, stage)


def _find_stage(spec: QueueSpec, stage_id: str) -> StageSpec:
    for stage in spec.stages:
        if stage.stage_id == stage_id:
            return stage
    raise KeyError(f"Unknown stage_id: {stage_id}")


def run_queue(args: argparse.Namespace) -> int:
    spec = load_queue_spec(Path(args.queue_config))
    state_dir = Path(args.state_dir)
    dry_run = bool(getattr(args, "dry_run", False))
    once = bool(getattr(args, "once", False))
    last_stage = StageSpec(stage_id="", stage_type="")

    for stage in spec.stages:
        last_stage = stage
        _record_event(
            state_dir=state_dir,
            queue_id=spec.queue_id,
            status="running",
            stage=stage,
            event_type="stage_started",
            payload={
                "stage_type": stage.stage_type,
                "run_name": stage.run_name,
                "log_path": stage.log_path or f"logs/{stage.stage_id}.log",
            },
        )

        if stage.stage_type == "wait_existing_run":
            if dry_run:
                _record_event(
                    state_dir=state_dir,
                    queue_id=spec.queue_id,
                    status="running",
                    stage=stage,
                    event_type="stage_dry_run",
                    payload={"mode": "wait_existing_run"},
                )
                if once:
                    break
                continue

            probe = wait_for_existing_run_completion(spec, stage)
            if probe.status != "completed":
                _record_event(
                    state_dir=state_dir,
                    queue_id=spec.queue_id,
                    status="halted",
                    stage=stage,
                    event_type="queue_halted",
                    payload={"reason": probe.reason},
                )
                return 1
            _record_event(
                state_dir=state_dir,
                queue_id=spec.queue_id,
                status="running",
                stage=stage,
                event_type="stage_passed",
                payload={
                    "checkpoint_path": str(probe.checkpoint_path) if probe.checkpoint_path else "",
                },
            )
            if once:
                break
            continue

        if stage.stage_type == "launch_run":
            if dry_run:
                ports = stage.master_ports or spec.default_master_ports
                command, _ = build_launch_command(spec, stage, ports[0] if ports else 0)
                _record_event(
                    state_dir=state_dir,
                    queue_id=spec.queue_id,
                    status="running",
                    stage=stage,
                    event_type="stage_dry_run",
                    payload={"command": command, "master_port": ports[0] if ports else 0},
                )
                if once:
                    break
                continue

            wait_for_idle_gpus(spec)
            stage_log_path = _stage_log_path(spec, stage)
            ports = stage.master_ports or spec.default_master_ports
            if not ports:
                ports = [0]

            last_log = ""
            for index, port in enumerate(ports):
                command, env = build_launch_command(spec, stage, port)
                return_code, last_log = run_stage_command(command, stage_log_path, env)
                if return_code == 0:
                    completion_cfg = dict(stage.completion)
                    probe = run_completion_probe(
                        runtime_root=spec.runtime_root,
                        output_root=_stage_output_path(spec, stage),
                        run_name=stage.run_name,
                        log_path=stage_log_path,
                        require_run_meta=bool(completion_cfg.get("require_run_meta", True)),
                        checkpoint_policy=str(completion_cfg.get("checkpoint_policy", "best_or_latest_epoch")),
                        completion_markers=list(completion_cfg.get("completion_markers", [])),
                    )
                    if probe.status != "completed":
                        _record_event(
                            state_dir=state_dir,
                            queue_id=spec.queue_id,
                            status="halted",
                            stage=stage,
                            event_type="queue_halted",
                            payload={"reason": probe.reason},
                        )
                        return 1
                    break

                if index == 0 and len(ports) > 1 and is_port_conflict(last_log):
                    _record_event(
                        state_dir=state_dir,
                        queue_id=spec.queue_id,
                        status="running",
                        stage=stage,
                        event_type="stage_retry_port_conflict",
                        payload={"failed_port": port, "next_port": ports[1]},
                    )
                    continue

                _record_event(
                    state_dir=state_dir,
                    queue_id=spec.queue_id,
                    status="halted",
                    stage=stage,
                    event_type="queue_halted",
                    payload={"reason": "launch_failed", "log_excerpt": last_log[-400:]},
                )
                return 1

            _record_event(
                state_dir=state_dir,
                queue_id=spec.queue_id,
                status="running",
                stage=stage,
                event_type="stage_passed",
                payload={"run_name": stage.run_name},
            )
            if once:
                break
            continue

        if stage.stage_type == "export_teacher_cache":
            if dry_run:
                checkpoint_stub = stage.checkpoint_from_stage or stage.stage_id
                try:
                    teacher_stage = _find_stage(spec, stage.checkpoint_from_stage)
                except KeyError:
                    teacher_run_name = checkpoint_stub
                else:
                    teacher_run_name = teacher_stage.run_name or checkpoint_stub
                checkpoint_path = spec.runtime_root / "outputs" / "checkpoints" / teacher_run_name / "best.pt"
                command, _ = build_export_command(spec, stage, checkpoint_path)
                _record_event(
                    state_dir=state_dir,
                    queue_id=spec.queue_id,
                    status="running",
                    stage=stage,
                    event_type="stage_dry_run",
                    payload={"command": command, "checkpoint_from_stage": stage.checkpoint_from_stage},
                )
                if once:
                    break
                continue

            teacher_stage = _find_stage(spec, stage.checkpoint_from_stage)
            teacher_checkpoint = choose_checkpoint_path(
                _stage_output_path(spec, teacher_stage) / "checkpoints" / teacher_stage.run_name,
                policy="best_or_latest_epoch",
            )
            if teacher_checkpoint is None:
                _record_event(
                    state_dir=state_dir,
                    queue_id=spec.queue_id,
                    status="halted",
                    stage=stage,
                    event_type="queue_halted",
                    payload={"reason": "teacher_checkpoint_missing"},
                )
                return 1

            command, env = build_export_command(spec, stage, teacher_checkpoint)
            export_log_path = _stage_log_path(spec, stage)
            return_code, combined_log = run_stage_command(command, export_log_path, env)
            if return_code != 0:
                _record_event(
                    state_dir=state_dir,
                    queue_id=spec.queue_id,
                    status="halted",
                    stage=stage,
                    event_type="queue_halted",
                    payload={"reason": "export_failed", "log_excerpt": combined_log[-400:]},
                )
                return 1

            cache_root = spec.runtime_root / stage.teacher_cache_dir
            required_manifests = [
                cache_root / "teacher_manifest_train.json",
                cache_root / "teacher_manifest_valid.json",
            ]
            if not all(path.exists() for path in required_manifests):
                _record_event(
                    state_dir=state_dir,
                    queue_id=spec.queue_id,
                    status="halted",
                    stage=stage,
                    event_type="queue_halted",
                    payload={"reason": "teacher_manifest_missing"},
                )
                return 1

            _record_event(
                state_dir=state_dir,
                queue_id=spec.queue_id,
                status="running",
                stage=stage,
                event_type="stage_passed",
                payload={"teacher_cache_dir": str(cache_root)},
            )
            if once:
                break
            continue

        raise ValueError(f"Unsupported stage_type: {stage.stage_type}")

    _record_event(
        state_dir=state_dir,
        queue_id=spec.queue_id,
        status="completed",
        stage=last_stage,
        event_type="queue_completed",
        payload={},
    )
    write_queue_state(
        state_dir=state_dir,
        queue_id=spec.queue_id,
        status="completed",
        current_stage_id="",
        current_run_name="",
        current_log_path="",
    )
    summary_path = state_dir / "summary.md"
    summary_path.write_text(
        render_summary_markdown(
            queue_id=spec.queue_id,
            status="completed",
            current_stage_id="",
            events_path=state_dir / "events.jsonl",
        ),
        encoding="utf-8",
    )
    return 0


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    return run_queue(args)


if __name__ == "__main__":
    raise SystemExit(main())
