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
