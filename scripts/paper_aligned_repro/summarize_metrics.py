#!/usr/bin/env python3
"""Summarize paper-aligned ChromBPNet reproduction metrics."""

import argparse
import csv
import json
import math
import statistics
from pathlib import Path


METRIC_KEYS = [
    "chrom_counts_pearson_peaks",
    "chrom_counts_spearman_peaks",
    "chrom_counts_mse_peaks",
    "chrom_profile_median_jsd_peaks",
    "chrom_profile_median_norm_jsd_peaks",
    "bias_counts_pearson_nonpeaks",
    "bias_counts_pearson_peaks",
    "bias_profile_median_jsd_nonpeaks",
    "bias_profile_median_jsd_peaks",
    "nobias_max_bias_response_mean",
]


def safe_get(d, path):
    cur = d
    for key in path:
        if key not in cur:
            return math.nan
        cur = cur[key]
    return cur


def parse_max_bias_response(path: Path):
    # Format example: corrected_0.001_0.001/0.001/0.001/0.001/0.001
    text = path.read_text(encoding="utf-8").strip()
    parts = text.split("_")
    if len(parts) < 2:
        return math.nan
    try:
        return float(parts[1])
    except ValueError:
        return math.nan


def row_from_run(seed_dir: Path):
    chrom_metrics_path = seed_dir / "chrombpnet" / "evaluation" / "chrombpnet_metrics.json"
    bias_metrics_path = seed_dir / "bias" / "evaluation" / "bias_metrics.json"
    if not bias_metrics_path.exists():
        bias_metrics_path = seed_dir / "bias" / "bias_metrics.json"
    max_bias_path = seed_dir / "chrombpnet" / "evaluation" / "chrombpnet_nobias_max_bias_response.txt"

    if not chrom_metrics_path.exists() or not bias_metrics_path.exists():
        return None

    chrom = json.loads(chrom_metrics_path.read_text(encoding="utf-8"))
    bias = json.loads(bias_metrics_path.read_text(encoding="utf-8"))

    row = {
        "fold": seed_dir.parent.name,
        "seed": seed_dir.name.replace("seed_", ""),
        "seed_dir": str(seed_dir),
        "chrom_counts_pearson_peaks": safe_get(chrom, ["counts_metrics", "peaks", "pearsonr"]),
        "chrom_counts_spearman_peaks": safe_get(chrom, ["counts_metrics", "peaks", "spearmanr"]),
        "chrom_counts_mse_peaks": safe_get(chrom, ["counts_metrics", "peaks", "mse"]),
        "chrom_profile_median_jsd_peaks": safe_get(chrom, ["profile_metrics", "peaks", "median_jsd"]),
        "chrom_profile_median_norm_jsd_peaks": safe_get(chrom, ["profile_metrics", "peaks", "median_norm_jsd"]),
        "bias_counts_pearson_nonpeaks": safe_get(bias, ["counts_metrics", "nonpeaks", "pearsonr"]),
        "bias_counts_pearson_peaks": safe_get(bias, ["counts_metrics", "peaks", "pearsonr"]),
        "bias_profile_median_jsd_nonpeaks": safe_get(bias, ["profile_metrics", "nonpeaks", "median_jsd"]),
        "bias_profile_median_jsd_peaks": safe_get(bias, ["profile_metrics", "peaks", "median_jsd"]),
        "nobias_max_bias_response_mean": parse_max_bias_response(max_bias_path) if max_bias_path.exists() else math.nan,
    }
    return row


def group_summary(rows, group_key):
    groups = {}
    for row in rows:
        groups.setdefault(row[group_key], []).append(row)

    summary_rows = []
    for group_name, grow in sorted(groups.items(), key=lambda x: x[0]):
        out = {group_key: group_name, "n_runs": len(grow)}
        for m in METRIC_KEYS:
            vals = [r[m] for r in grow if isinstance(r[m], (int, float)) and not math.isnan(r[m])]
            if vals:
                out[f"{m}_mean"] = statistics.mean(vals)
                out[f"{m}_std"] = statistics.stdev(vals) if len(vals) > 1 else 0.0
            else:
                out[f"{m}_mean"] = math.nan
                out[f"{m}_std"] = math.nan
        summary_rows.append(out)
    return summary_rows


def write_csv(rows, out_path: Path):
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def overall_summary(rows):
    out = {"n_runs": len(rows)}
    for m in METRIC_KEYS:
        vals = [r[m] for r in rows if isinstance(r[m], (int, float)) and not math.isnan(r[m])]
        if vals:
            out[f"{m}_mean"] = statistics.mean(vals)
            out[f"{m}_std"] = statistics.stdev(vals) if len(vals) > 1 else 0.0
        else:
            out[f"{m}_mean"] = math.nan
            out[f"{m}_std"] = math.nan
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, help="Path like .../runs/K562_ATAC")
    parser.add_argument("--output-dir", required=True, help="Summary output directory")
    args = parser.parse_args()

    run_root = Path(args.run_root).resolve()
    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for fold_dir in sorted(run_root.glob("fold_*"), key=lambda p: p.name):
        if not fold_dir.is_dir():
            continue
        for seed_dir in sorted(fold_dir.glob("seed_*"), key=lambda p: p.name):
            if not seed_dir.is_dir():
                continue
            row = row_from_run(seed_dir)
            if row is not None:
                rows.append(row)

    if not rows:
        raise SystemExit(f"No completed runs found under {run_root}")

    by_seed = group_summary(rows, "seed")
    by_fold = group_summary(rows, "fold")
    overall = overall_summary(rows)

    write_csv(rows, out_dir / "run_metrics.csv")
    write_csv(by_seed, out_dir / "summary_by_seed.csv")
    write_csv(by_fold, out_dir / "summary_by_fold.csv")
    (out_dir / "summary_overall.json").write_text(json.dumps(overall, indent=2), encoding="utf-8")

    print(f"Wrote {len(rows)} runs to {out_dir / 'run_metrics.csv'}")
    print(f"Wrote seed summary to {out_dir / 'summary_by_seed.csv'}")
    print(f"Wrote fold summary to {out_dir / 'summary_by_fold.csv'}")
    print(f"Wrote overall summary to {out_dir / 'summary_overall.json'}")


if __name__ == "__main__":
    main()
