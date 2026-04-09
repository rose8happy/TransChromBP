"""Smoke checks for the NT v2 teacher-distill launcher."""

from __future__ import annotations

import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_ntv2_teacher_distill_launcher_exists_and_is_bash_valid() -> None:
    script_path = (
        REPO_ROOT
        / "vendor"
        / "transchrombp"
        / "transchrombp"
        / "scripts"
        / "run_ntv2_teacher_distill.sh"
    )
    assert script_path.is_file()

    content = script_path.read_text(encoding="utf-8")
    assert "teacher_cache_export" in content
    assert "train_tutorial_teacher_v2_ntv2_distill_short10.yaml" in content
    assert "train_tutorial_teacher_v2_ntv2_distill_full.yaml" in content

    result = subprocess.run(["bash", "-n", str(script_path)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
