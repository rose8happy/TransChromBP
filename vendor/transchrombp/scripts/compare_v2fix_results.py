#!/usr/bin/env python3
"""Compare V2 fix ablation results across experiments.

Reads epoch_metrics.jsonl from each run, extracts best-epoch metrics,
and prints a comparison table.

Usage:
  python scripts/compare_v2fix_results.py v2fix_20260321

  # Also include baseline runs for comparison
  python scripts/compare_v2fix_results.py v2fix_20260321 --include-baseline
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional


@dataclass
class RunResult:
    run_name: str
    group: str
    seed: int
    best_epoch: int
    metrics: Dict[str, Any] = field(default_factory=dict)


def extract_best_metrics(metrics_jsonl: Path) -> Optional[Dict[str, Any]]:
    """Read epoch_metrics.jsonl and return the best epoch's full record."""
    if not metrics_jsonl.exists():
        return None

    records = []
    with open(metrics_jsonl) as f:
        for line in f:
            line = line.strip()
            if line:
                records.append(json.loads(line))

    if not records:
        return None

    # 最后一条记录通常包含 best_epoch 信息
    last = records[-1]
    best_epoch = last.get("best_epoch")

    if best_epoch is not None:
        # 找到 best_epoch 对应的记录
        for rec in records:
            if rec.get("epoch") == best_epoch:
                return rec
    # Fallback: 返回最后一条
    return last


def get_metric(record: Dict[str, Any], region: str, key: str) -> Optional[float]:
    """Extract a nested metric like peak.profile_target_jsd_full_mean."""
    val = record.get("val", {})
    region_data = val.get(region, {})
    return region_data.get(key)


def discover_runs(output_dir: Path, tag: str) -> List[RunResult]:
    """Discover all runs matching the tag pattern."""
    logs_dir = output_dir / "logs"
    results = []

    if not logs_dir.exists():
        return results

    for run_dir in sorted(logs_dir.iterdir()):
        if not run_dir.is_dir():
            continue
        name = run_dir.name
        if not name.startswith(tag):
            continue

        metrics_file = run_dir / "epoch_metrics.jsonl"
        rec = extract_best_metrics(metrics_file)
        if rec is None:
            continue

        # Parse group and seed from run name: v2fix_20260321_freeze_s42
        parts = name.replace(tag + "_", "").rsplit("_s", 1)
        group = parts[0] if parts else name
        seed = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0

        results.append(RunResult(
            run_name=name,
            group=group,
            seed=seed,
            best_epoch=rec.get("best_epoch", rec.get("epoch", -1)),
            metrics=rec,
        ))

    return results


def discover_baseline_runs(output_dir: Path) -> List[RunResult]:
    """Discover baseline runs from ablation_tf_20260318."""
    logs_dir = output_dir / "logs"
    results = []
    baseline_prefix = "ablation_tf_20260318_full_s"

    if not logs_dir.exists():
        return results

    for run_dir in sorted(logs_dir.iterdir()):
        if not run_dir.is_dir():
            continue
        name = run_dir.name
        if not name.startswith(baseline_prefix):
            continue

        metrics_file = run_dir / "epoch_metrics.jsonl"
        rec = extract_best_metrics(metrics_file)
        if rec is None:
            continue

        seed_str = name.replace("ablation_tf_20260318_full_s", "")
        seed = int(seed_str) if seed_str.isdigit() else 0

        results.append(RunResult(
            run_name=name,
            group="baseline",
            seed=seed,
            best_epoch=rec.get("best_epoch", rec.get("epoch", -1)),
            metrics=rec,
        ))

    return results


METRIC_COLUMNS = [
    ("peak JSD full", "peak", "profile_target_jsd_full_mean", "↓"),
    ("peak JSD debiased", "peak", "profile_target_jsd_debiased_mean", "↓"),
    ("full-debiased gap", "peak", "profile_full_debiased_jsd", "↓"),
    ("count_pearson debiased", "peak", "count_pearson_debiased", "↑"),
    ("count_pearson full", "peak", "count_pearson_full", "↑"),
    ("count_mae debiased", "peak", "count_mae_debiased", "↓"),
    ("nonpeak JSD", "nonpeak", "profile_target_jsd_full_mean", "↓"),
]


def print_table(runs: List[RunResult]) -> None:
    """Print a formatted comparison table."""
    # Header
    header = f"{'group':>14s} {'seed':>5s} {'epoch':>5s}"
    for col_name, _, _, direction in METRIC_COLUMNS:
        header += f" {col_name + direction:>18s}"
    print(header)
    print("-" * len(header))

    # Group runs
    groups: Dict[str, List[RunResult]] = {}
    for r in runs:
        groups.setdefault(r.group, []).append(r)

    for group_name, group_runs in groups.items():
        values_by_col: Dict[str, List[float]] = {col[0]: [] for col in METRIC_COLUMNS}

        for r in sorted(group_runs, key=lambda x: x.seed):
            row = f"{r.group:>14s} {r.seed:>5d} {r.best_epoch:>5d}"
            for col_name, region, key, _ in METRIC_COLUMNS:
                val = get_metric(r.metrics, region, key)
                if val is not None:
                    row += f" {val:>18.5f}"
                    values_by_col[col_name].append(val)
                else:
                    row += f" {'N/A':>18s}"
            print(row)

        # Print mean ± std if multiple seeds
        if len(group_runs) > 1:
            import statistics
            row = f"{'  mean±std':>14s} {'':>5s} {'':>5s}"
            for col_name, _, _, _ in METRIC_COLUMNS:
                vals = values_by_col[col_name]
                if len(vals) >= 2:
                    mean = statistics.mean(vals)
                    std = statistics.stdev(vals)
                    row += f" {mean:>10.5f}±{std:<6.5f}"
                elif len(vals) == 1:
                    row += f" {vals[0]:>18.5f}"
                else:
                    row += f" {'N/A':>18s}"
            print(row)
        print()


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare V2 fix ablation results")
    parser.add_argument("tag", help="Experiment tag, e.g. v2fix_20260321")
    parser.add_argument("--output-dir", default=None,
                        help="TransChromBP outputs directory (default: auto-detect)")
    parser.add_argument("--include-baseline", action="store_true",
                        help="Include baseline (ablation_tf_20260318_full) for comparison")
    parser.add_argument("--json", default=None, help="Output results as JSON to file")
    args = parser.parse_args()

    # Auto-detect output dir
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        # 尝试几个常见位置
        candidates = [
            Path(__file__).resolve().parent.parent / "outputs",
            Path("/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs"),
        ]
        output_dir = None
        for c in candidates:
            if c.exists():
                output_dir = c
                break
        if output_dir is None:
            print("[error] Cannot find outputs directory. Use --output-dir.", file=sys.stderr)
            sys.exit(1)

    print(f"Output dir: {output_dir}")
    print(f"Tag: {args.tag}")
    print()

    runs = discover_runs(output_dir, args.tag)
    if args.include_baseline:
        runs.extend(discover_baseline_runs(output_dir))

    if not runs:
        print("[warn] No completed runs found.")
        sys.exit(0)

    print(f"Found {len(runs)} completed runs\n")
    print_table(runs)

    # Also check for soup eval results
    eval_dir = output_dir / "eval"
    if eval_dir.exists():
        eval_files = sorted(eval_dir.glob("*.json"))
        if eval_files:
            print("=== Soup / Checkpoint Evaluation Results ===\n")
            for ef in eval_files:
                with open(ef) as f:
                    data = json.load(f)
                ckpt = data.get("checkpoint", "?")
                metrics = data.get("metrics", {})
                peak = metrics.get("peak", {})
                print(f"  {ef.name}:")
                jsd_f = peak.get("profile_target_jsd_full_mean", "N/A")
                jsd_d = peak.get("profile_target_jsd_debiased_mean", "N/A")
                count_r = peak.get("count_pearson_full", "N/A")
                if isinstance(jsd_f, float):
                    print(f"    peak JSD_full={jsd_f:.5f}  JSD_debiased={jsd_d:.5f}  count_r={count_r:.4f}")
                else:
                    print(f"    peak JSD_full={jsd_f}  JSD_debiased={jsd_d}  count_r={count_r}")
            print()

    if args.json:
        json_data = []
        for r in runs:
            entry = {"run_name": r.run_name, "group": r.group, "seed": r.seed, "best_epoch": r.best_epoch}
            for col_name, region, key, _ in METRIC_COLUMNS:
                entry[col_name] = get_metric(r.metrics, region, key)
            json_data.append(entry)
        with open(args.json, "w") as f:
            json.dump(json_data, f, indent=2)
        print(f"JSON results saved to {args.json}")


if __name__ == "__main__":
    main()
