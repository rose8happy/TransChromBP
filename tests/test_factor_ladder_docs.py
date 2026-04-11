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
    assert "尚未回填" in factor_row["当前结论 / 进度"]
    assert "正式 short10 / formal run" in factor_row["当前结论 / 进度"]
    assert "回填" in factor_row["下一步"]
    assert "E2/E3 short10" in factor_row["下一步"]
    assert "technical smoke" in factor_row["下一步"]


def test_registry_registers_factor_ladder_family_row() -> None:
    rows = parse_markdown_table("docs/experiments/registry.md")
    factor_row = next(
        row for row in rows if row["family_id"] == "`alphagenome_factor_ladder`"
    )

    assert factor_row["status"] == "`staged`"
    assert factor_row["latest_terminal_run"] == "`hierdec4096_exporter_distill_smoke_20260411_2110`"
    assert "回填 canonical archive" in factor_row["next_allowed_action"]
    assert "E2/E3 short10" in factor_row["next_allowed_action"]
    assert "technical smoke" in factor_row["next_allowed_action"]
    assert "remote isolated workdir" in factor_row["notes"]
    assert "5877077" in factor_row["notes"]
    assert "single-GPU" in factor_row["notes"]
    assert "尚未回填" in factor_row["notes"]
    assert "genos_summary" in factor_row["notes"]
    assert "n_records" in factor_row["notes"]
