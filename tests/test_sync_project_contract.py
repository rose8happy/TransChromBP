from __future__ import annotations

import json
import os
import stat
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts/sync_project.sh"


def write_rsync_stub(bin_dir: Path, capture_path: Path) -> None:
    stub = bin_dir / "rsync"
    stub.write_text(
        "\n".join(
            [
                "#!/usr/bin/env python3",
                "import json",
                "from pathlib import Path",
                "import sys",
                f'capture_path = Path(r"{capture_path}")',
                "calls = []",
                "if capture_path.exists():",
                "    calls = json.loads(capture_path.read_text(encoding='utf-8'))",
                "calls.append(sys.argv[1:])",
                "capture_path.write_text(json.dumps(calls), encoding='utf-8')",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    stub.chmod(stub.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def run_script(args: list[str], env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)
    return subprocess.run(
        ["bash", str(SCRIPT), *args],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        env=merged_env,
        check=False,
    )


def test_help_lists_explicit_runtime_commands() -> None:
    result = run_script(["help"])

    assert result.returncode == 0
    help_text = result.stdout + result.stderr
    assert "publish-runtime-6000" in help_text
    assert "publish-runtime-6002" in help_text
    assert "pull-results-6000" in help_text
    assert "pull-results-6002" in help_text
    assert "status-all" in help_text


def test_legacy_deploy_command_exits_with_migration_hint() -> None:
    result = run_script(["deploy"])

    assert result.returncode != 0
    hint_text = result.stdout + result.stderr
    assert "deprecated" in hint_text.lower()
    assert "publish-runtime-6000" in hint_text


def test_publish_runtime_6000_dry_run_targets_runtime_repo(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    capture_path = tmp_path / "rsync_calls.json"
    write_rsync_stub(bin_dir, capture_path)

    result = run_script(
        ["publish-runtime-6000", "--dry-run"],
        env={"PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}"},
    )

    assert result.returncode == 0
    calls = json.loads(capture_path.read_text(encoding="utf-8"))
    assert len(calls) == 1

    argv = calls[0]
    argv_text = " ".join(argv)
    assert "--dry-run" in argv
    assert "/data1/zhoujiazhen/bylw_atac/TransChromBP/" in argv_text
    assert "/data1/zhoujiazhen/bylw_atac/chromBPNet/" not in argv_text
    assert "--exclude=logs/" in argv
    assert "--exclude=outputs/" in argv
    assert "--exclude=.git/" in argv


def test_publish_runtime_excludes_local_only_and_staging_paths(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    capture_path = tmp_path / "rsync_calls.json"
    write_rsync_stub(bin_dir, capture_path)

    result = run_script(
        ["publish-runtime-6000", "--dry-run"],
        env={"PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}"},
    )

    assert result.returncode == 0
    calls = json.loads(capture_path.read_text(encoding="utf-8"))
    assert len(calls) == 1

    argv = calls[0]
    for expected in (
        "--exclude=.agents/",
        "--exclude=.claude/",
        "--exclude=.codex",
        "--exclude=.codex_remote_edit/",
        "--exclude=.idea/",
        "--exclude=tmp_remote_edit/",
        "--exclude=references/local-only/",
    ):
        assert expected in argv


def test_publish_runtime_deletes_stale_repo_tracked_paths(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    capture_path = tmp_path / "rsync_calls.json"
    write_rsync_stub(bin_dir, capture_path)

    result = run_script(
        ["publish-runtime-6000", "--dry-run"],
        env={"PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}"},
    )

    assert result.returncode == 0
    calls = json.loads(capture_path.read_text(encoding="utf-8"))
    assert len(calls) == 1

    argv = calls[0]
    assert "--delete" in argv


def test_pull_results_6000_dry_run_only_requests_logs_and_reports(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    capture_path = tmp_path / "rsync_calls.json"
    write_rsync_stub(bin_dir, capture_path)

    result = run_script(
        ["pull-results-6000", "--dry-run"],
        env={"PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}"},
    )

    assert result.returncode == 0
    calls = json.loads(capture_path.read_text(encoding="utf-8"))
    assert len(calls) == 2

    call_texts = [" ".join(call) for call in calls]
    assert any("/data1/zhoujiazhen/bylw_atac/TransChromBP/logs/" in text for text in call_texts)
    assert any("/data1/zhoujiazhen/bylw_atac/TransChromBP/reports/" in text for text in call_texts)
    assert not any("/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/" in text for text in call_texts)
