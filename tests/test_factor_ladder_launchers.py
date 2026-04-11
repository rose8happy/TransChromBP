from __future__ import annotations

import json
import os
import stat
import subprocess
from pathlib import Path

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ARCHIVE_ROOT = PROJECT_ROOT / "vendor/transchrombp"
BASELINE_LONGCTX_SHORT10_CONFIG = ARCHIVE_ROOT / "configs/train/train_tutorial_teacher_v2_longctx4096_short10.yaml"
SHORT10_CONFIG = ARCHIVE_ROOT / "configs/train/train_tutorial_teacher_v2_hierdec4096_short10.yaml"
TEACHER30_CONFIG = ARCHIVE_ROOT / "configs/train/train_tutorial_teacher_v2_hierdec4096_teacher30.yaml"
DISTILL_SHORT10_CONFIG = ARCHIVE_ROOT / "configs/train/train_tutorial_teacher_v2_hierdec4096_distill_short10.yaml"
MODEL_CONFIG = ARCHIVE_ROOT / "configs/model/transchrombp_teacher_v2_hierdec4096.yaml"
DATA_CONFIG = ARCHIVE_ROOT / "configs/data/data_tutorial_canonical_v1_longctx4096.yaml"
LAUNCHER = ARCHIVE_ROOT / "scripts/run_teacher_v2_hierdec4096_probe.sh"
DISTILL_LAUNCHER = ARCHIVE_ROOT / "scripts/run_teacher_v2_hierdec4096_distill.sh"
EXPECTED_OUTPUT_DIR = "outputs"
DISTILL_TEACHER_CACHE_SUFFIX = Path("teacher_cache/teacher_v2_hierdec4096_teacher30")


def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def flatten_dict(tree: dict, prefix: str = "") -> dict[str, object]:
    flat: dict[str, object] = {}
    for key, value in tree.items():
        path = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            flat.update(flatten_dict(value, path))
        else:
            flat[path] = value
    return flat


def diff_paths(lhs: dict, rhs: dict) -> set[str]:
    lhs_flat = flatten_dict(lhs)
    rhs_flat = flatten_dict(rhs)
    return {key for key in lhs_flat.keys() | rhs_flat.keys() if lhs_flat.get(key) != rhs_flat.get(key)}


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


def capture_launcher_calls(
    tmp_path: Path,
    launcher: Path,
    extra_env: dict[str, str] | None = None,
) -> list[list[str]]:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    capture_path = tmp_path / "argv.json"
    write_stub_python(bin_dir, capture_path)

    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}{os.pathsep}{env.get('PATH', '')}"
    env["VENV_DIR"] = str(tmp_path / "missing-venv")
    if extra_env:
        env.update(extra_env)
    subprocess.run(["bash", "-n", str(launcher)], cwd=ARCHIVE_ROOT, env=env, check=True)
    subprocess.run(["bash", str(launcher)], cwd=ARCHIVE_ROOT, env=env, check=True)

    return json.loads(capture_path.read_text(encoding="utf-8"))


def capture_launcher_argv(tmp_path: Path) -> list[str]:
    calls = capture_launcher_calls(tmp_path, LAUNCHER)
    assert len(calls) == 1
    return calls[0]


def write_teacher_manifest(cache_dir: Path, split: str) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = cache_dir / f"teacher_manifest_{split}.json"
    manifest_path.write_text(json.dumps({"split": split}) + "\n", encoding="utf-8")


def test_hierdec_short10_config_diff_is_limited_to_task3_contract() -> None:
    baseline_cfg = load_yaml(BASELINE_LONGCTX_SHORT10_CONFIG)
    train_cfg = load_yaml(SHORT10_CONFIG)

    assert diff_paths(baseline_cfg, train_cfg) == {"logging.output_dir", "logging.run_name"}
    assert train_cfg["logging"]["output_dir"] == EXPECTED_OUTPUT_DIR
    assert train_cfg["logging"]["run_name"] == "teacher_v2_hierdec4096_short10"


def test_hierdec_short10_config_matches_task3_contract() -> None:
    train_cfg = load_yaml(SHORT10_CONFIG)

    assert train_cfg["seed"] == 42
    assert train_cfg["max_epochs"] == 10
    assert train_cfg["data"]["config_path"] == "configs/data/data_tutorial_canonical_v1_longctx4096.yaml"
    assert train_cfg["data"]["max_seq_len"] == 4096
    assert train_cfg["data"]["input_len"] == 4096
    assert train_cfg["data"]["output_len"] == 1000
    assert train_cfg["data"]["supervised_bp"] == 1000
    assert train_cfg["data"]["teacher_target_names"] == []
    assert train_cfg["loss"]["distill_profile_weight"] == 0.0
    assert train_cfg["loss"]["distill_count_weight"] == 0.0
    assert train_cfg["loss"]["distill_rank_weight"] == 0.0
    assert train_cfg["logging"]["output_dir"] == EXPECTED_OUTPUT_DIR
    assert train_cfg["logging"]["run_name"] == "teacher_v2_hierdec4096_short10"


def test_hierdec_teacher30_config_only_extends_budget() -> None:
    short10_cfg = load_yaml(SHORT10_CONFIG)
    teacher30_cfg = load_yaml(TEACHER30_CONFIG)

    assert diff_paths(short10_cfg, teacher30_cfg) == {"logging.run_name", "max_epochs"}
    assert teacher30_cfg["logging"]["output_dir"] == EXPECTED_OUTPUT_DIR
    assert teacher30_cfg["max_epochs"] == 30
    assert teacher30_cfg["logging"]["run_name"] == "teacher_v2_hierdec4096_teacher30"


def test_hierdec_distill_short10_config_matches_task5_contract() -> None:
    train_cfg = load_yaml(DISTILL_SHORT10_CONFIG)

    assert train_cfg["seed"] == 42
    assert train_cfg["max_epochs"] == 10
    assert train_cfg["loss"]["distill_profile_weight"] == 0.5
    assert train_cfg["loss"]["distill_count_weight"] == 0.05
    assert train_cfg["data"]["teacher_target_names"] == ["profile16", "logcount"]
    assert train_cfg["data"]["teacher_cache_dir"] == "outputs/teacher_cache/teacher_v2_hierdec4096_teacher30"
    assert train_cfg["logging"]["output_dir"] == EXPECTED_OUTPUT_DIR
    assert train_cfg["logging"]["run_name"] == "teacher_v2_hierdec4096_distill_short10"


def test_hierdec_distill_short10_config_diff_is_limited_to_distill_contract() -> None:
    short10_cfg = load_yaml(SHORT10_CONFIG)
    distill_cfg = load_yaml(DISTILL_SHORT10_CONFIG)

    assert diff_paths(short10_cfg, distill_cfg) == {
        "data.teacher_cache_dir",
        "data.teacher_target_names",
        "logging.run_name",
        "loss.distill_count_weight",
        "loss.distill_profile_weight",
    }


def test_hierdec_probe_launcher_executes_expected_default_contract(tmp_path: Path) -> None:
    assert MODEL_CONFIG.exists()
    assert DATA_CONFIG.exists()

    argv = capture_launcher_argv(tmp_path)

    assert argv[:2] == ["-m", "transchrombp.training.train_ddp"]
    assert argv[2:] == [
        "--model-config",
        str(MODEL_CONFIG),
        "--train-config",
        str(SHORT10_CONFIG),
        "--data-config",
        str(DATA_CONFIG),
        "--output-dir",
        str(ARCHIVE_ROOT / "outputs"),
        "--run-name",
        "teacher_v2_hierdec4096_short10",
    ]


def test_hierdec_distill_launcher_exports_teacher_cache_then_trains(tmp_path: Path) -> None:
    assert MODEL_CONFIG.exists()
    assert DATA_CONFIG.exists()
    teacher_ckpt = tmp_path / "teacher30.pt"
    teacher_ckpt.write_text("stub checkpoint", encoding="utf-8")
    output_base = tmp_path / "isolated_outputs"
    teacher_cache_dir = output_base / DISTILL_TEACHER_CACHE_SUFFIX
    write_teacher_manifest(teacher_cache_dir, "valid")

    calls = capture_launcher_calls(
        tmp_path,
        DISTILL_LAUNCHER,
        extra_env={
            "TEACHER_CKPT": str(teacher_ckpt),
            "OUTPUT_BASE": str(output_base),
        },
    )

    assert len(calls) == 2

    exporter_argv, train_argv = calls
    expected_teacher_cache_dir = teacher_cache_dir

    assert exporter_argv[:2] == ["-m", "transchrombp.evaluation.model_teacher_cache_export"]
    assert exporter_argv[2:] == [
        "--checkpoint",
        str(teacher_ckpt),
        "--data-config",
        str(DATA_CONFIG),
        "--output-dir",
        str(expected_teacher_cache_dir),
        "--splits",
        "train",
        "valid",
    ]

    assert train_argv[:2] == ["-m", "transchrombp.training.train_ddp"]
    assert train_argv[2:] == [
        "--model-config",
        str(MODEL_CONFIG),
        "--train-config",
        str(DISTILL_SHORT10_CONFIG),
        "--data-config",
        str(DATA_CONFIG),
        "--output-dir",
        str(output_base),
        "--run-name",
        "teacher_v2_hierdec4096_distill_short10",
        "--teacher-cache-dir",
        str(expected_teacher_cache_dir),
    ]


def test_hierdec_distill_launcher_skips_export_when_train_and_valid_manifests_exist(tmp_path: Path) -> None:
    assert MODEL_CONFIG.exists()
    assert DATA_CONFIG.exists()
    output_base = tmp_path / "isolated_outputs"
    teacher_cache_dir = output_base / DISTILL_TEACHER_CACHE_SUFFIX
    write_teacher_manifest(teacher_cache_dir, "train")
    write_teacher_manifest(teacher_cache_dir, "valid")

    calls = capture_launcher_calls(
        tmp_path,
        DISTILL_LAUNCHER,
        extra_env={
            "OUTPUT_BASE": str(output_base),
        },
    )

    assert len(calls) == 1

    train_argv = calls[0]
    assert train_argv[:2] == ["-m", "transchrombp.training.train_ddp"]
    assert train_argv[2:] == [
        "--model-config",
        str(MODEL_CONFIG),
        "--train-config",
        str(DISTILL_SHORT10_CONFIG),
        "--data-config",
        str(DATA_CONFIG),
        "--output-dir",
        str(output_base),
        "--run-name",
        "teacher_v2_hierdec4096_distill_short10",
        "--teacher-cache-dir",
        str(teacher_cache_dir),
    ]
