"""Tests for NT v2 teacher-distill configs."""

from __future__ import annotations

import sys
from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]
IMPORT_ROOTS = [
    REPO_ROOT / "vendor" / "transchrombp",
    REPO_ROOT / "src",
    REPO_ROOT,
]
for import_root in IMPORT_ROOTS:
    if (import_root / "transchrombp").exists() and str(import_root) not in sys.path:
        sys.path.insert(0, str(import_root))
        break

from transchrombp.models.transchrombp import build_transchrombp_from_config


def test_ntv2_teacher_distill_configs_build_and_expose_teacher_targets() -> None:
    model_config_dir = REPO_ROOT / "vendor" / "transchrombp" / "transchrombp" / "configs" / "model"
    train_config_dir = REPO_ROOT / "vendor" / "transchrombp" / "transchrombp" / "configs" / "train"

    model_cfg = yaml.safe_load(
        (model_config_dir / "transchrombp_teacher_v2_center_pool_ntv2_distill.yaml").read_text(
            encoding="utf-8"
        )
    )
    model = build_transchrombp_from_config(model_cfg)
    assert model_cfg["foundation_model"]["enabled"] is True
    assert model_cfg["foundation_model"]["mode"] == "distill_only"
    assert model.foundation_cross_adapter is None
    assert model.foundation_residual_head is None

    short_cfg = yaml.safe_load(
        (train_config_dir / "train_tutorial_teacher_v2_ntv2_distill_short10.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert short_cfg["loss"]["distill_profile_weight"] > 0.0
    assert short_cfg["loss"]["distill_count_weight"] > 0.0
    assert short_cfg["loss"]["distill_rank_weight"] == 0.0
    assert short_cfg["data"]["teacher_target_names"] == ["profile16", "logcount"]
    assert short_cfg["data"]["foundation_cache_features"] == ["layer_07__bins16_mean"]

    full_cfg = yaml.safe_load(
        (train_config_dir / "train_tutorial_teacher_v2_ntv2_distill_full.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert full_cfg["max_epochs"] >= short_cfg["max_epochs"]
    assert full_cfg["data"]["teacher_target_names"] == ["profile16", "logcount"]
    assert full_cfg["loss"]["distill_profile_weight"] == short_cfg["loss"]["distill_profile_weight"]
    assert full_cfg["loss"]["distill_count_weight"] == short_cfg["loss"]["distill_count_weight"]
