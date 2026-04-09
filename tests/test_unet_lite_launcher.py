"""Smoke checks for the U-Net-lite decoder launcher."""

from __future__ import annotations

import subprocess
from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_unet_lite_launcher_and_train_config_exist() -> None:
    train_config_path = (
        REPO_ROOT
        / "vendor"
        / "transchrombp"
        / "transchrombp"
        / "configs"
        / "train"
        / "train_tutorial_teacher_v2_readout_short10.yaml"
    )
    launcher_path = (
        REPO_ROOT
        / "vendor"
        / "transchrombp"
        / "transchrombp"
        / "scripts"
        / "run_unet_lite_decoder_probe.sh"
    )
    data_config_6002_path = (
        REPO_ROOT
        / "vendor"
        / "transchrombp"
        / "transchrombp"
        / "configs"
        / "data"
        / "data_tutorial_canonical_v1_6002.yaml"
    )

    assert train_config_path.is_file()
    train_cfg = yaml.safe_load(train_config_path.read_text(encoding="utf-8"))
    assert train_cfg["loss"]["distill_profile_weight"] == 0.0
    assert train_cfg["loss"]["distill_count_weight"] == 0.0
    assert train_cfg["data"]["teacher_target_names"] == []
    assert data_config_6002_path.is_file()
    data_cfg_6002 = yaml.safe_load(data_config_6002_path.read_text(encoding="utf-8"))
    assert data_cfg_6002["sampling"]["nonpeak_ratio"] == 0.1

    assert launcher_path.is_file()
    content = launcher_path.read_text(encoding="utf-8")
    assert "train_tutorial_teacher_v2_readout_short10.yaml" in content
    assert "transchrombp_teacher_v2_center_pool_unet_lite_v1.yaml" in content
    assert "resolve_venv_dir" in content
    assert "resolve_data_config" in content
    assert "data_tutorial_canonical_v1_6002.yaml" in content
    assert "/home/zhengwei/bylw_atac/.mamba/envs/transchrombp" in content

    result = subprocess.run(["bash", "-n", str(launcher_path)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
