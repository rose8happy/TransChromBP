import csv
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def parse_markdown_table(path: str) -> list[dict[str, str]]:
    text = (REPO_ROOT / path).read_text(encoding="utf-8")
    rows = [line.strip() for line in text.splitlines() if line.strip().startswith("|")]
    header = [cell.strip() for cell in rows[0].strip("|").split("|")]
    data_rows = rows[2:]

    parsed_rows = []
    for row in data_rows:
        cells = [cell.strip() for cell in row.strip("|").split("|")]
        if len(cells) != len(header):
            continue
        parsed_rows.append(dict(zip(header, cells)))

    return parsed_rows


def parse_runs_csv() -> list[dict[str, str]]:
    with (REPO_ROOT / "docs/experiments/runs.csv").open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def test_tracking_registers_factor_ladder_as_live_row() -> None:
    rows = parse_markdown_table("TRACKING.md")
    factor_row = next(
        row for row in rows if row["事项"] == "AlphaGenome-like factor ladder（E1/E2/E3）"
    )

    assert factor_row["当前状态"] == "进行中"
    assert "docs/spec/plan" in factor_row["当前结论 / 进度"]
    assert "remote isolated workdir" in factor_row["当前结论 / 进度"]
    assert "single-GPU" in factor_row["当前结论 / 进度"]
    assert "2110" in factor_row["当前结论 / 进度"]
    assert "5877077" in factor_row["当前结论 / 进度"]
    assert "079b603" in factor_row["当前结论 / 进度"]
    assert "已回填" in factor_row["当前结论 / 进度"]
    assert "teacher_v2_hierdec4096_short10_s42_6000_20260411_r1" in factor_row["当前结论 / 进度"]
    assert "已于 `2026-04-11 21:46:19 CST` 启动" in factor_row["当前结论 / 进度"]
    assert "23:20-23:45" in factor_row["当前结论 / 进度"]
    assert "unattended matrix queue" in factor_row["当前结论 / 进度"]
    assert "outputs/queue/factor_ladder_unattended_20260411" in factor_row["当前结论 / 进度"]
    assert "E2" in factor_row["下一步"]
    assert "不并行再开 `E3`" in factor_row["下一步"]
    assert "当前 run log" in factor_row["下一步"]
    assert "GPU" in factor_row["下一步"]
    assert "teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1" in factor_row["下一步"]
    assert "teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1" in factor_row["下一步"]
    assert "teacher_v2_hierdec4096_distill_short10_s42_6000_20260411_r1" in factor_row["下一步"]


def test_registry_registers_factor_ladder_family_row() -> None:
    rows = parse_markdown_table("docs/experiments/registry.md")
    factor_row = next(
        row for row in rows if row["family_id"] == "`alphagenome_factor_ladder`"
    )

    assert factor_row["status"] == "`running`"
    assert factor_row["latest_terminal_run"] == "`hierdec4096_exporter_distill_smoke_20260411_2110`"
    assert "teacher_v2_hierdec4096_short10_s42_6000_20260411_r1" in factor_row["next_allowed_action"]
    assert "不并行再开 `E3`" in factor_row["next_allowed_action"]
    assert "remote runtime" in factor_row["notes"]
    assert "079b603" in factor_row["notes"]
    assert "5877077" in factor_row["notes"]
    assert "single-GPU" in factor_row["notes"]
    assert "genos_summary" in factor_row["notes"]
    assert "n_records" in factor_row["notes"]
    assert "2026-04-11 21:46:19 CST" in factor_row["notes"]
    assert "23:20-23:45" in factor_row["notes"]
    assert "unattended matrix queue" in factor_row["next_allowed_action"]
    assert "outputs/queue/factor_ladder_unattended_20260411" in factor_row["notes"]


def test_runs_manifest_pre_registers_unattended_queue_downstream_runs() -> None:
    rows = parse_runs_csv()
    queued = {
        row["run_id"]: row
        for row in rows
        if row["family_id"] == "alphagenome_factor_ladder" and row["run_status"] == "queued"
    }

    assert queued["teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1"]["commit_sha"] == "079b603"
    assert queued["teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1"]["commit_sha"] == "079b603"
    assert queued["teacher_v2_hierdec4096_distill_short10_s42_6000_20260411_r1"]["commit_sha"] == "079b603"
