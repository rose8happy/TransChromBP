import csv
from pathlib import Path
import subprocess


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

    assert factor_row["当前状态"] == "待处理"
    assert "2026-04-12 08:46 CST" in factor_row["当前结论 / 进度"]
    assert "export_e2_teacher_cache" in factor_row["当前结论 / 进度"]
    assert "teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1" in factor_row["当前结论 / 进度"]
    assert "teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1" in factor_row["当前结论 / 进度"]
    assert "export_e2_teacher_cache.log" in factor_row["当前结论 / 进度"]
    assert "--targets profile16 logcount" in factor_row["当前结论 / 进度"]
    assert "E3 distill" in factor_row["当前结论 / 进度"]
    assert "先修 exporter 参数漂移" in factor_row["下一步"]
    assert "model_teacher_cache_export.py" in factor_row["下一步"]
    assert "最小 smoke" in factor_row["下一步"]
    assert "不要继续推进 E3 distill" in factor_row["下一步"]
    assert "不要重新启动 factor-ladder unattended queue" in factor_row["下一步"]


def test_tracking_registers_loss_balance_closeout_row() -> None:
    rows = parse_markdown_table("TRACKING.md")
    loss_row = next(
        row for row in rows if row["事项"] == "loss balance curriculum（2026-04-17）"
    )

    assert loss_row["当前状态"] == "进行中"
    assert "teacher_v2_center_pool_lossbal_e0_selector_s42_6000_20260417_r1" in loss_row["当前结论 / 进度"]
    assert "selector_jsd_gap=0.00114" in loss_row["当前结论 / 进度"]
    assert "selector-only `s1234`" in loss_row["当前结论 / 进度"]
    assert "2026-04-17 14:55:46 CST" in loss_row["当前结论 / 进度"]
    assert "teacher_v2_center_pool_lossbal_e2_dynamic_count_s42_6000_20260417_r1" in loss_row["当前结论 / 进度"]
    assert "2026-04-17 21:45:24 CST" in loss_row["当前结论 / 进度"]
    assert "count_weight` 从 `0.10` 增到 `0.23186`" in loss_row["当前结论 / 进度"]
    assert "dynamic-count 没拿到更好的 JSD，也没有更好的 count" in loss_row["当前结论 / 进度"]
    assert "没有给出 clear gain" in loss_row["下一步"]
    assert "静态 safe-envelope" in loss_row["下一步"]
    assert "暂不继续扩 dynamic-count 的 sweep" in loss_row["下一步"]


def test_registry_registers_factor_ladder_family_row() -> None:
    rows = parse_markdown_table("docs/experiments/registry.md")
    factor_row = next(
        row for row in rows if row["family_id"] == "`alphagenome_factor_ladder`"
    )

    assert factor_row["status"] == "`halted`"
    assert factor_row["latest_terminal_run"] == "`teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1`"
    assert "export_e2_teacher_cache" in factor_row["next_allowed_action"]
    assert "最小 smoke" in factor_row["next_allowed_action"]
    assert "不继续启动新的 factor-ladder run" in factor_row["next_allowed_action"]
    assert "不推进 `teacher_v2_hierdec4096_distill_short10_s42_6000_20260411_r1`" in factor_row["next_allowed_action"]
    assert "remote runtime" in factor_row["notes"]
    assert "079b603" in factor_row["notes"]
    assert "ac55153" in factor_row["notes"]
    assert "2026-04-12 08:46 CST" in factor_row["notes"]
    assert "queue_halted" in factor_row["notes"]
    assert "model_teacher_cache_export.py" in factor_row["notes"]
    assert "--targets profile16 logcount" in factor_row["notes"]
    assert "outputs/queue/factor_ladder_unattended_20260411" in factor_row["notes"]


def test_registry_registers_loss_balance_family_row() -> None:
    rows = parse_markdown_table("docs/experiments/registry.md")
    loss_row = next(
        row for row in rows if row["family_id"] == "`loss_balance_curriculum`"
    )

    assert loss_row["status"] == "`running`"
    assert (
        loss_row["latest_terminal_run"]
        == "`teacher_v2_center_pool_lossbal_e2_dynamic_count_s42_6000_20260417_r1`"
    )
    assert "当前没有 active run" in loss_row["next_allowed_action"]
    assert "静态 safe-envelope" in loss_row["next_allowed_action"]
    assert "不要重开 selector-only `s1234`" in loss_row["next_allowed_action"]
    assert "不要直接扩 dynamic-count sweep" in loss_row["next_allowed_action"]
    assert loss_row["mounted_worktree"] == "`n/a`"
    assert "closeout/loss_balance_curriculum/20260418" in loss_row["closeout_tags"]
    assert "snapshot/loss-balance-20260417/20260419" in loss_row["closeout_tags"]
    assert "selector_jsd_gap=0.00114" in loss_row["notes"]
    assert "2026-04-17 14:55:46 CST" in loss_row["notes"]
    assert "2026-04-17 21:45:24 CST" in loss_row["notes"]
    assert "count_weight` 从 `0.10` 增到 `0.23186`" in loss_row["notes"]
    assert "clear gain" in loss_row["notes"]
    assert "exact local dynamic-count worktree state" in loss_row["notes"]


def test_loss_balance_snapshot_and_closeout_tags_exist() -> None:
    tag_list = subprocess.run(
        ["git", "tag", "--list", "closeout/loss_balance_curriculum/20260418", "snapshot/loss-balance-20260417/20260419"],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=True,
    ).stdout

    assert "closeout/loss_balance_curriculum/20260418" in tag_list
    assert "snapshot/loss-balance-20260417/20260419" in tag_list


def test_runs_manifest_records_factor_ladder_halt_and_loss_balance_closeout() -> None:
    rows = parse_runs_csv()
    indexed = {row["run_id"]: row for row in rows}

    assert indexed["teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1"]["run_status"] == "completed"
    assert indexed["teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1"]["end_time"] == "2026-04-12 05:35:01 CST"
    assert indexed["teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1"]["run_status"] == "completed"
    assert indexed["teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1"]["end_time"] == "2026-04-12 08:46:08 CST"
    assert indexed["export_e2_teacher_cache_20260412"]["run_status"] == "halted"
    assert indexed["export_e2_teacher_cache_20260412"]["gate_verdict"] == "technical_fail"
    assert "--targets profile16 logcount" in indexed["export_e2_teacher_cache_20260412"]["notes"]
    assert indexed["teacher_v2_hierdec4096_distill_short10_s42_6000_20260411_r1"]["run_status"] == "blocked"
    assert "queue_halted" in indexed["teacher_v2_hierdec4096_distill_short10_s42_6000_20260411_r1"]["notes"]

    assert indexed["teacher_v2_center_pool_lossbal_e0_selector_s42_6000_20260417_r1"]["run_status"] == "completed"
    assert (
        indexed["teacher_v2_center_pool_lossbal_e0_selector_s42_6000_20260417_r1"]["report_path"]
        == "reports/loss_balance_selector_wave1_seed42_20260417.md"
    )
    assert indexed["teacher_v2_center_pool_lossbal_e0_selector_s1234_6000_20260417_r1"]["run_status"] == "stopped"
    assert indexed["teacher_v2_center_pool_lossbal_e0_selector_s1234_6000_20260417_r1"]["gate_verdict"] == "reprioritized"
    assert indexed["teacher_v2_center_pool_lossbal_e2_dynamic_count_smoke_s42_6000_20260417_r1"]["gate_verdict"] == "technical_pass"
    assert indexed["teacher_v2_center_pool_lossbal_e2_dynamic_count_s42_6000_20260417_r1"]["run_status"] == "completed"
    assert indexed["teacher_v2_center_pool_lossbal_e2_dynamic_count_s42_6000_20260417_r1"]["gate_verdict"] == "no_clear_gain"
    assert "count_weight rose from 0.10 to 0.23186" in indexed["teacher_v2_center_pool_lossbal_e2_dynamic_count_s42_6000_20260417_r1"]["notes"]
