#!/usr/bin/env python3
"""Compare best-by-JSD and best-by-loss checkpoints across loss-balance runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def get(dct: dict[str, Any], path: str) -> Any:
    cur: Any = dct
    for part in path.split("."):
        cur = cur[part]
    return cur


def format_value(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.6f}"
    return str(value)


def collect_rows(root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for metrics_path in sorted(root.glob("*/epoch_metrics.jsonl")):
        hist = [json.loads(line) for line in metrics_path.read_text(encoding="utf-8").splitlines() if line.strip()]
        hist = [entry for entry in hist if entry.get("val")]
        if not hist:
            continue

        best_jsd = min(hist, key=lambda h: get(h, "val.peak.profile_target_jsd_full_mean"))
        best_loss = min(hist, key=lambda h: get(h, "val.peak.loss_total"))

        rows.append(
            {
                "run": metrics_path.parent.name,
                "best_jsd_epoch": best_jsd["epoch"],
                "best_jsd": get(best_jsd, "val.peak.profile_target_jsd_full_mean"),
                "best_jsd_count_r": get(best_jsd, "val.peak.count_pearson_full"),
                "best_jsd_gap": get(best_jsd, "val.peak.profile_full_debiased_jsd"),
                "best_jsd_count_shift": get(best_jsd, "val.peak.count_full_debiased_abs"),
                "best_loss_epoch": best_loss["epoch"],
                "loss_sel_jsd": get(best_loss, "val.peak.profile_target_jsd_full_mean"),
                "loss_sel_count_r": get(best_loss, "val.peak.count_pearson_full"),
                "loss_sel_gap": get(best_loss, "val.peak.profile_full_debiased_jsd"),
                "loss_sel_count_shift": get(best_loss, "val.peak.count_full_debiased_abs"),
                "selector_epoch_gap": abs(best_jsd["epoch"] - best_loss["epoch"]),
                "selector_jsd_gap": get(best_loss, "val.peak.profile_target_jsd_full_mean")
                - get(best_jsd, "val.peak.profile_target_jsd_full_mean"),
            }
        )

    rows.sort(key=lambda row: (row["best_jsd"], -row["best_jsd_count_r"]))
    return rows


def render_table(rows: list[dict[str, Any]]) -> str:
    if not rows:
        return "No runs found."

    columns = list(rows[0].keys())
    widths = {
        col: max(len(col), *(len(format_value(row[col])) for row in rows))
        for col in columns
    }

    header = "  ".join(col.ljust(widths[col]) for col in columns)
    divider = "  ".join("-" * widths[col] for col in columns)
    body = [
        "  ".join(format_value(row[col]).ljust(widths[col]) for col in columns)
        for row in rows
    ]
    return "\n".join([header, divider, *body])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs_root", type=Path, help="Directory containing per-run output subdirs")
    args = parser.parse_args()

    rows = collect_rows(args.logs_root)
    print(render_table(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
