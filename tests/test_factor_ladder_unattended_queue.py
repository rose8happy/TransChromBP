from __future__ import annotations

import json
from pathlib import Path

import yaml

from transchrombp.orchestration.factor_ladder_unattended_queue import (
    append_event,
    choose_checkpoint_path,
    load_queue_spec,
    render_summary_markdown,
    run_completion_probe,
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
