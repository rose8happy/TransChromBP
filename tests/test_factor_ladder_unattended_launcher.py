from __future__ import annotations

import json
import os
import stat
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ARCHIVE_ROOT = PROJECT_ROOT / "vendor" / "transchrombp"
LAUNCHER = ARCHIVE_ROOT / "scripts" / "run_factor_ladder_unattended_matrix.sh"
QUEUE_CONFIG = ARCHIVE_ROOT / "configs" / "queues" / "factor_ladder_unattended_20260411.yaml"
STATE_DIR = ARCHIVE_ROOT / "outputs" / "queue" / "factor_ladder_unattended_20260411"


def write_stub_python(bin_dir: Path, capture_path: Path) -> None:
    stub = bin_dir / "python"
    stub.write_text(
        "\n".join(
            [
                "#!/usr/bin/env python3",
                "import json",
                "import os",
                "from pathlib import Path",
                "import sys",
                f'capture_path = r"{capture_path}"',
                "path = Path(capture_path)",
                "calls = []",
                "if path.exists():",
                "    calls = json.loads(path.read_text(encoding='utf-8'))",
                "calls.append({'argv': sys.argv[1:], 'pythonpath': os.environ.get('PYTHONPATH', '')})",
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

    subprocess.run(["bash", "-n", str(LAUNCHER)], cwd=ARCHIVE_ROOT, env=env, check=True)
    subprocess.run(["bash", str(LAUNCHER), "--dry-run", "--once"], cwd=ARCHIVE_ROOT, env=env, check=True)

    calls = json.loads(capture_path.read_text(encoding="utf-8"))
    assert len(calls) == 1
    call = calls[0]
    argv = call["argv"]
    assert argv[:2] == ["-m", "transchrombp.orchestration.factor_ladder_unattended_queue"]
    assert argv[2:] == [
        "--queue-config",
        str(QUEUE_CONFIG),
        "--state-dir",
        str(STATE_DIR),
        "--dry-run",
        "--once",
    ]
    assert call["pythonpath"].split(os.pathsep)[0] == str(ARCHIVE_ROOT)


def test_unattended_launcher_prefers_src_import_root_when_src_layout_exists(tmp_path: Path) -> None:
    fake_root = tmp_path / "runtime_repo"
    scripts_dir = fake_root / "scripts"
    scripts_dir.mkdir(parents=True)
    (fake_root / "src" / "transchrombp").mkdir(parents=True)
    launcher_copy = scripts_dir / LAUNCHER.name
    launcher_copy.write_text(LAUNCHER.read_text(encoding="utf-8"), encoding="utf-8")
    launcher_copy.chmod(launcher_copy.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    capture_path = tmp_path / "argv_src.json"
    write_stub_python(bin_dir, capture_path)

    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}{os.pathsep}{env.get('PATH', '')}"
    env["VENV_DIR"] = str(tmp_path / "missing-venv")

    subprocess.run(["bash", str(launcher_copy), "--dry-run", "--once"], cwd=fake_root, env=env, check=True)

    calls = json.loads(capture_path.read_text(encoding="utf-8"))
    assert len(calls) == 1
    assert calls[0]["pythonpath"].split(os.pathsep)[0] == str(fake_root / "src")
