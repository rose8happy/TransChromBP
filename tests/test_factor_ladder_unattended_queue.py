from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from transchrombp.orchestration.factor_ladder_unattended_queue import (
    append_event,
    choose_checkpoint_path,
    load_queue_spec,
    render_summary_markdown,
    parse_args,
    run_completion_probe,
    run_queue,
    write_queue_state,
)


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
    assert cfg["poll_interval_sec"] == 30
    assert cfg["gpu_wait_interval_sec"] == 60
    assert cfg["gpu_idle_used_mem_threshold_mib"] == 4096

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
    assert wait_stage["log_path"] == "logs/teacher_v2_hierdec4096_short10_s42_6000_20260411_r1.log"
    assert wait_stage["completion"]["require_run_meta"] is True
    assert wait_stage["completion"]["checkpoint_policy"] == "best_or_latest_epoch"

    e1_stage = stages[1]
    assert e1_stage["run_name"] == "teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1"
    assert e1_stage["model_config"] == "configs/model/transchrombp_teacher_v2_center_pool_longctx4096.yaml"
    assert e1_stage["train_config"] == "configs/train/train_tutorial_teacher_v2_longctx4096_short10.yaml"
    assert e1_stage["data_config"] == "configs/data/data_tutorial_canonical_v1_longctx4096.yaml"
    assert e1_stage["output_dir"] == "outputs"
    assert e1_stage["log_path"] == "logs/teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1.log"
    assert e1_stage["train_gpu_ids"] == "0,1"
    assert e1_stage["nproc_per_node"] == 2
    assert e1_stage["master_ports"] == [29982, 29983]
    assert e1_stage["env"]["RUN_NAME"] == "teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1"
    assert e1_stage["completion"]["require_run_meta"] is True
    assert e1_stage["completion"]["checkpoint_policy"] == "best_or_latest_epoch"

    e2_teacher_stage = stages[2]
    assert e2_teacher_stage["run_name"] == "teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1"
    assert e2_teacher_stage["model_config"] == "configs/model/transchrombp_teacher_v2_hierdec4096.yaml"
    assert e2_teacher_stage["train_config"] == "configs/train/train_tutorial_teacher_v2_hierdec4096_teacher30.yaml"
    assert e2_teacher_stage["data_config"] == "configs/data/data_tutorial_canonical_v1_longctx4096.yaml"
    assert e2_teacher_stage["output_dir"] == "outputs"
    assert e2_teacher_stage["log_path"] == "logs/teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1.log"
    assert e2_teacher_stage["train_gpu_ids"] == "0,1"
    assert e2_teacher_stage["nproc_per_node"] == 2
    assert e2_teacher_stage["master_ports"] == [29982, 29983]
    assert e2_teacher_stage["env"]["RUN_NAME"] == "teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1"
    assert e2_teacher_stage["completion"]["require_run_meta"] is True
    assert e2_teacher_stage["completion"]["checkpoint_policy"] == "best_or_latest_epoch"

    export_stage = stages[3]
    assert export_stage["checkpoint_from_stage"] == "e2_hierdec4096_teacher30"
    assert export_stage["teacher_cache_dir"] == "outputs/teacher_cache/teacher_v2_hierdec4096_teacher30"
    assert export_stage["data_config"] == "configs/data/data_tutorial_canonical_v1_longctx4096.yaml"
    assert export_stage["targets"] == ["profile16", "logcount"]
    assert export_stage["splits"] == ["train", "valid"]
    assert "env" not in export_stage

    distill_stage = stages[4]
    assert distill_stage["run_name"] == "teacher_v2_hierdec4096_distill_short10_s42_6000_20260411_r1"
    assert distill_stage["model_config"] == "configs/model/transchrombp_teacher_v2_hierdec4096.yaml"
    assert distill_stage["train_config"] == "configs/train/train_tutorial_teacher_v2_hierdec4096_distill_short10.yaml"
    assert distill_stage["data_config"] == "configs/data/data_tutorial_canonical_v1_longctx4096.yaml"
    assert distill_stage["output_dir"] == "outputs"
    assert distill_stage["log_path"] == "logs/teacher_v2_hierdec4096_distill_short10_s42_6000_20260411_r1.log"
    assert distill_stage["train_gpu_ids"] == "0,1"
    assert distill_stage["nproc_per_node"] == 2
    assert distill_stage["master_ports"] == [29982, 29983]
    assert distill_stage["env"]["RUN_NAME"] == "teacher_v2_hierdec4096_distill_short10_s42_6000_20260411_r1"
    assert distill_stage["env"]["TEACHER_CACHE_DIR"] == "outputs/teacher_cache/teacher_v2_hierdec4096_teacher30"
    assert distill_stage["completion"]["require_run_meta"] is True
    assert distill_stage["completion"]["checkpoint_policy"] == "best_or_latest_epoch"


def test_load_queue_spec_parses_runtime_defaults_and_stage_specs() -> None:
    spec = load_queue_spec(QUEUE_CONFIG)

    assert spec.queue_id == "factor_ladder_unattended_20260411"
    assert spec.runtime_root == Path("/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411")
    assert spec.state_dir == Path("outputs/queue/factor_ladder_unattended_20260411")
    assert spec.default_master_ports == [29982, 29983]
    assert spec.poll_interval_sec == 30
    assert spec.gpu_wait_interval_sec == 60
    assert spec.gpu_idle_used_mem_threshold_mib == 4096
    assert [stage.stage_id for stage in spec.stages] == [
        "wait_existing_e2_short10",
        "e1_longctx4096_short10",
        "e2_hierdec4096_teacher30",
        "export_e2_teacher_cache",
        "e3_hierdec4096_distill_short10",
    ]


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


def test_run_completion_probe_treats_explicit_fatal_marker_as_failed_when_completion_markers_are_missing(
    tmp_path: Path,
) -> None:
    output_root = tmp_path / "outputs"
    run_name = "teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1"
    log_path = tmp_path / "logs" / f"{run_name}.log"
    log_path.parent.mkdir(parents=True)
    log_path.write_text(
        "[bootstrap-fail] /bin/bash: line 12: FOO: unbound variable\n",
        encoding="utf-8",
    )

    probe = run_completion_probe(
        runtime_root=tmp_path,
        output_root=output_root,
        run_name=run_name,
        log_path=log_path,
        require_run_meta=True,
        checkpoint_policy="best_or_latest_epoch",
        completion_markers=["[metrics] wrote", "[meta] wrote"],
    )

    assert probe.status == "failed"
    assert probe.reason == "fatal_marker_detected"


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
    state_payload = json.loads(state_path.read_text(encoding="utf-8"))
    assert state_payload["current_stage_id"] == "e1_longctx4096_short10"
    assert "e1_longctx4096_short10" in summary
    assert "stage_started" in summary


def test_parse_args_binds_queue_cli_flags() -> None:
    args = parse_args(
        [
            "--queue-config",
            "configs/queues/factor_ladder_unattended_20260411.yaml",
            "--state-dir",
            "outputs/queue/factor_ladder_unattended_20260411",
            "--dry-run",
            "--once",
        ]
    )

    assert args.queue_config == "configs/queues/factor_ladder_unattended_20260411.yaml"
    assert args.state_dir == "outputs/queue/factor_ladder_unattended_20260411"
    assert args.dry_run is True
    assert args.once is True


def test_run_queue_dry_run_writes_queue_state_events_and_summary(tmp_path: Path) -> None:
    args = argparse.Namespace(
        queue_config=QUEUE_CONFIG,
        state_dir=tmp_path / "queue-state",
        dry_run=True,
        once=True,
    )

    exit_code = run_queue(args)

    assert exit_code == 0
    state_path = args.state_dir / "queue_state.json"
    events_path = args.state_dir / "events.jsonl"
    summary_path = args.state_dir / "summary.md"

    assert state_path.exists()
    assert events_path.exists()
    assert summary_path.exists()

    state_payload = json.loads(state_path.read_text(encoding="utf-8"))
    assert state_payload["queue_id"] == "factor_ladder_unattended_20260411"
    assert state_payload["status"] == "completed"
    assert state_payload["current_stage_id"] == ""
    events = [json.loads(line) for line in events_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    assert [event["event_type"] for event in events] == ["stage_started", "stage_dry_run", "queue_completed"]
    assert "stage_dry_run" in summary_path.read_text(encoding="utf-8")


def test_run_queue_launch_stage_dry_run_records_command_without_executing(
    tmp_path: Path,
    monkeypatch,
) -> None:
    queue_config = tmp_path / "queue.yaml"
    queue_config.write_text(
        "\n".join(
            [
                "queue_id: launch_dry_run_queue",
                f"runtime_root: {tmp_path}",
                "state_dir: outputs/queue/launch_dry_run_queue",
                "default_master_ports: [29982, 29983]",
                "poll_interval_sec: 1",
                "gpu_wait_interval_sec: 1",
                "gpu_idle_used_mem_threshold_mib: 1",
                "stages:",
                "  - stage_id: launch_stage",
                "    stage_type: launch_run",
                "    run_name: run_one",
                "    model_config: configs/model/transchrombp_teacher_v2_hierdec4096.yaml",
                "    train_config: configs/train/train_tutorial_teacher_v2_hierdec4096_teacher30.yaml",
                "    data_config: configs/data/data_tutorial_canonical_v1_longctx4096.yaml",
                "    output_dir: outputs",
                "    log_path: logs/run_one.log",
                '    train_gpu_ids: "0,1"',
                "    nproc_per_node: 2",
                "    master_ports: [29982, 29983]",
                "    env:",
                "      RUN_NAME: run_one",
                "    completion:",
                "      require_run_meta: true",
                "      checkpoint_policy: best_or_latest_epoch",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    def fail_if_called(*args, **kwargs):  # type: ignore[no-untyped-def]
        raise AssertionError("launch command should not run during dry-run")

    monkeypatch.setattr(
        "transchrombp.orchestration.factor_ladder_unattended_queue.run_stage_command",
        fail_if_called,
    )
    monkeypatch.setattr(
        "transchrombp.orchestration.factor_ladder_unattended_queue.wait_for_idle_gpus",
        fail_if_called,
    )

    args = argparse.Namespace(
        queue_config=queue_config,
        state_dir=tmp_path / "state",
        dry_run=True,
        once=True,
    )

    exit_code = run_queue(args)

    assert exit_code == 0
    events = [json.loads(line) for line in (args.state_dir / "events.jsonl").read_text(encoding="utf-8").splitlines() if line.strip()]
    assert events[1]["event_type"] == "stage_dry_run"
    assert events[1]["stage_id"] == "launch_stage"
    assert events[1]["command"][0] == "torchrun"
    assert events[1]["command"][1] == "--standalone"


def test_run_queue_export_stage_dry_run_skips_checkpoint_lookup(
    tmp_path: Path,
    monkeypatch,
) -> None:
    queue_config = tmp_path / "queue_export.yaml"
    queue_config.write_text(
        "\n".join(
            [
                "queue_id: export_dry_run_queue",
                f"runtime_root: {tmp_path}",
                "state_dir: outputs/queue/export_dry_run_queue",
                "default_master_ports: [29982, 29983]",
                "poll_interval_sec: 1",
                "gpu_wait_interval_sec: 1",
                "gpu_idle_used_mem_threshold_mib: 1",
                "stages:",
                "  - stage_id: export_stage",
                "    stage_type: export_teacher_cache",
                "    checkpoint_from_stage: teacher_stage",
                "    teacher_cache_dir: outputs/teacher_cache/export_dry_run_queue",
                "    data_config: configs/data/data_tutorial_canonical_v1_longctx4096.yaml",
                '    targets: ["profile16", "logcount"]',
                '    splits: ["train", "valid"]',
                "  - stage_id: teacher_stage",
                "    stage_type: launch_run",
                "    run_name: teacher_run",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    def fail_if_called(*args, **kwargs):  # type: ignore[no-untyped-def]
        raise AssertionError("checkpoint lookup should not run during dry-run")

    monkeypatch.setattr(
        "transchrombp.orchestration.factor_ladder_unattended_queue.choose_checkpoint_path",
        fail_if_called,
    )

    args = argparse.Namespace(
        queue_config=queue_config,
        state_dir=tmp_path / "export-state",
        dry_run=True,
        once=True,
    )

    exit_code = run_queue(args)

    assert exit_code == 0
    events = [json.loads(line) for line in (args.state_dir / "events.jsonl").read_text(encoding="utf-8").splitlines() if line.strip()]
    assert events[1]["event_type"] == "stage_dry_run"
    assert events[1]["stage_id"] == "export_stage"
    assert "--targets" in events[1]["command"]
    target_index = events[1]["command"].index("--targets")
    assert events[1]["command"][target_index + 1 : target_index + 3] == ["profile16", "logcount"]
    checkpoint_index = events[1]["command"].index("--checkpoint")
    assert events[1]["command"][checkpoint_index + 1] == str(
        tmp_path / "outputs" / "checkpoints" / "teacher_run" / "best.pt"
    )


def test_run_queue_export_stage_uses_teacher_checkpoint_under_outputs_root(
    tmp_path: Path,
    monkeypatch,
) -> None:
    queue_config = tmp_path / "queue_export_run.yaml"
    queue_config.write_text(
        "\n".join(
            [
                "queue_id: export_run_queue",
                f"runtime_root: {tmp_path}",
                "state_dir: outputs/queue/export_run_queue",
                "default_master_ports: [29982, 29983]",
                "poll_interval_sec: 1",
                "gpu_wait_interval_sec: 1",
                "gpu_idle_used_mem_threshold_mib: 1",
                "stages:",
                "  - stage_id: export_stage",
                "    stage_type: export_teacher_cache",
                "    checkpoint_from_stage: teacher_stage",
                "    teacher_cache_dir: outputs/teacher_cache/export_run_queue",
                "    data_config: configs/data/data_tutorial_canonical_v1_longctx4096.yaml",
                '    targets: ["profile16", "logcount"]',
                '    splits: ["train", "valid"]',
                "  - stage_id: teacher_stage",
                "    stage_type: launch_run",
                "    run_name: teacher_run",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    observed: dict[str, object] = {}

    def fake_choose_checkpoint_path(checkpoint_dir: Path, policy: str) -> Path:
        observed["checkpoint_dir"] = checkpoint_dir
        observed["policy"] = policy
        return checkpoint_dir / "best.pt"

    def fake_run_stage_command(command, stage_log_path, env):  # type: ignore[no-untyped-def]
        observed["command"] = list(command)
        observed["stage_log_path"] = stage_log_path
        cache_root = tmp_path / "outputs" / "teacher_cache" / "export_run_queue"
        cache_root.mkdir(parents=True, exist_ok=True)
        (cache_root / "teacher_manifest_train.json").write_text("{}", encoding="utf-8")
        (cache_root / "teacher_manifest_valid.json").write_text("{}", encoding="utf-8")
        return 0, ""

    monkeypatch.setattr(
        "transchrombp.orchestration.factor_ladder_unattended_queue.choose_checkpoint_path",
        fake_choose_checkpoint_path,
    )
    monkeypatch.setattr(
        "transchrombp.orchestration.factor_ladder_unattended_queue.run_stage_command",
        fake_run_stage_command,
    )

    args = argparse.Namespace(
        queue_config=queue_config,
        state_dir=tmp_path / "export-run-state",
        dry_run=False,
        once=True,
    )

    exit_code = run_queue(args)

    assert exit_code == 0
    assert observed["checkpoint_dir"] == tmp_path / "outputs" / "checkpoints" / "teacher_run"
    assert observed["policy"] == "best_or_latest_epoch"
    assert "--targets" in observed["command"]


def test_run_queue_launch_stage_retries_once_on_port_conflict(
    tmp_path: Path,
    monkeypatch,
) -> None:
    queue_config = tmp_path / "queue_retry.yaml"
    queue_config.write_text(
        "\n".join(
            [
                "queue_id: retry_queue",
                f"runtime_root: {tmp_path}",
                "state_dir: outputs/queue/retry_queue",
                "default_master_ports: [29982, 29983]",
                "poll_interval_sec: 1",
                "gpu_wait_interval_sec: 1",
                "gpu_idle_used_mem_threshold_mib: 1",
                "stages:",
                "  - stage_id: launch_stage",
                "    stage_type: launch_run",
                "    run_name: run_retry",
                "    model_config: configs/model/transchrombp_teacher_v2_hierdec4096.yaml",
                "    train_config: configs/train/train_tutorial_teacher_v2_hierdec4096_teacher30.yaml",
                "    data_config: configs/data/data_tutorial_canonical_v1_longctx4096.yaml",
                "    output_dir: outputs",
                "    log_path: logs/run_retry.log",
                '    train_gpu_ids: "0,1"',
                "    nproc_per_node: 2",
                "    master_ports: [12345, 12346]",
                "    env:",
                "      RUN_NAME: run_retry",
                "    completion:",
                "      require_run_meta: true",
                "      checkpoint_policy: best_or_latest_epoch",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(
        "transchrombp.orchestration.factor_ladder_unattended_queue.wait_for_idle_gpus",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        "transchrombp.orchestration.factor_ladder_unattended_queue.run_completion_probe",
        lambda *args, **kwargs: type("Probe", (), {"status": "completed", "checkpoint_path": tmp_path / "checkpoints" / "run_retry" / "best.pt", "reason": "ok"})(),
    )

    calls: list[tuple[list[str], dict[str, str]]] = []

    def fake_run_stage_command(command, stage_log_path, env):  # type: ignore[no-untyped-def]
        calls.append((list(command), dict(env)))
        if len(calls) == 1:
            return 1, "Address already in use"
        return 0, "done"

    monkeypatch.setattr(
        "transchrombp.orchestration.factor_ladder_unattended_queue.run_stage_command",
        fake_run_stage_command,
    )

    args = argparse.Namespace(
        queue_config=queue_config,
        state_dir=tmp_path / "retry-state",
        dry_run=False,
        once=True,
    )

    exit_code = run_queue(args)

    assert exit_code == 0
    assert len(calls) == 2
    assert calls[0][1]["MASTER_PORT"] == "12345"
    assert calls[1][1]["MASTER_PORT"] == "12346"
    events = [json.loads(line) for line in (args.state_dir / "events.jsonl").read_text(encoding="utf-8").splitlines() if line.strip()]
    assert any(event["event_type"] == "stage_retry_port_conflict" for event in events)
