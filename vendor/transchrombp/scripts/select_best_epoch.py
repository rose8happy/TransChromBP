#!/usr/bin/env python3
"""Evaluate TransChromBP epoch checkpoints and select the external best checkpoint."""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any


EPOCH_RE = re.compile(r"epoch_(\d+)\.pt$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate TransChromBP epoch checkpoints with a unified external metric and select the best checkpoint."
    )
    parser.add_argument("--checkpoint-glob", required=True, help="Glob for epoch checkpoints, e.g. outputs/checkpoints/run/epoch_*.pt")
    parser.add_argument("--output-dir", required=True, help="Directory for per-epoch eval JSONs and summary files")
    parser.add_argument("--metric", default="results.peak.profile_target_jsd_full_median", help="Metric path used to select the best checkpoint")
    parser.add_argument("--mode", choices=["min", "max"], default="min", help="Whether lower or higher metric is better")
    parser.add_argument("--split", default="valid", help="Evaluation split passed to evaluate_checkpoint")
    parser.add_argument("--data-config", default="", help="Optional override data config yaml")
    parser.add_argument("--batch-size", type=int, default=0, help="Optional batch size override")
    parser.add_argument("--num-workers", type=int, default=-1, help="Optional num_workers override")
    parser.add_argument("--max-regions", type=int, default=0, help="Optional max_regions override")
    parser.add_argument("--region-source", default="", help="Optional region source override")
    parser.add_argument("--nonpeak-ratio", type=float, default=-1.0, help="Optional nonpeak ratio override")
    parser.add_argument("--device", default="auto", help="Device passed to evaluate_checkpoint")
    parser.add_argument("--devices", default="", help="Comma-separated devices used to shard checkpoint evaluation, e.g. cuda:0,cuda:1")
    parser.add_argument("--max-workers", type=int, default=0, help="Optional worker cap for --devices; 0 means use all listed devices")
    parser.add_argument("--force", action="store_true", help="Recompute eval JSON even if it already exists")
    return parser.parse_args()


def extract_metric(payload: dict[str, Any], path: str) -> float:
    current: Any = payload
    for part in path.split("."):
        if not isinstance(current, dict) or part not in current:
            return float("nan")
        current = current[part]
    try:
        return float(current)
    except (TypeError, ValueError):
        return float("nan")


def epoch_sort_key(path_str: str) -> tuple[int, str]:
    path = Path(path_str)
    match = EPOCH_RE.search(path.name)
    if match:
        return (int(match.group(1)), path.name)
    return (10**9, path.name)


def is_better(candidate: float, current_best: float, mode: str) -> bool:
    if math.isnan(current_best):
        return not math.isnan(candidate)
    if math.isnan(candidate):
        return False
    if mode == "max":
        return candidate > current_best
    return candidate < current_best


def write_csv(rows: list[dict[str, Any]], out_path: Path) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_eval_cmd(args: argparse.Namespace, checkpoint_path: Path, output_json: Path, device: str) -> list[str]:
    cmd = [
        sys.executable,
        "-m",
        "transchrombp.evaluation.evaluate_checkpoint",
        "--checkpoint",
        str(checkpoint_path),
        "--split",
        args.split,
        "--device",
        device,
        "--output",
        str(output_json),
    ]
    if args.data_config:
        cmd += ["--data-config", args.data_config]
    if args.batch_size > 0:
        cmd += ["--batch-size", str(args.batch_size)]
    if args.num_workers >= 0:
        cmd += ["--num-workers", str(args.num_workers)]
    if args.max_regions > 0:
        cmd += ["--max-regions", str(args.max_regions)]
    if args.region_source:
        cmd += ["--region-source", args.region_source]
    if args.nonpeak_ratio >= 0.0:
        cmd += ["--nonpeak-ratio", str(args.nonpeak_ratio)]
    return cmd


def parse_device_list(raw: str) -> list[str]:
    return [item.strip() for item in raw.split(",") if item.strip()]


def output_json_for(checkpoint_path: Path, out_dir: Path, split: str) -> Path:
    return out_dir / f"{checkpoint_path.stem}__{split}.json"


def evaluate_assigned_checkpoints(
    checkpoint_paths: list[Path],
    out_dir: Path,
    args: argparse.Namespace,
    device: str,
) -> None:
    for checkpoint_path in checkpoint_paths:
        output_json = output_json_for(checkpoint_path, out_dir, args.split)
        if args.force or not output_json.exists():
            print(f"[selector] device={device} evaluate {checkpoint_path.name}")
            subprocess.run(build_eval_cmd(args, checkpoint_path, output_json, device), check=True)


def evaluate_checkpoints(
    checkpoint_paths: list[Path],
    out_dir: Path,
    args: argparse.Namespace,
) -> None:
    device_list = parse_device_list(args.devices)
    if len(device_list) <= 1:
        evaluate_assigned_checkpoints(checkpoint_paths, out_dir, args, args.device)
        return

    worker_count = len(device_list) if args.max_workers <= 0 else min(args.max_workers, len(device_list))
    worker_devices = device_list[:worker_count]
    assignments: dict[str, list[Path]] = {device: [] for device in worker_devices}
    for idx, checkpoint_path in enumerate(checkpoint_paths):
        device = worker_devices[idx % worker_count]
        assignments[device].append(checkpoint_path)

    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = [
            executor.submit(evaluate_assigned_checkpoints, assigned_paths, out_dir, args, device)
            for device, assigned_paths in assignments.items()
            if assigned_paths
        ]
        for future in futures:
            future.result()


def main() -> None:
    args = parse_args()
    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    checkpoint_path_strs = sorted(glob.glob(args.checkpoint_glob), key=epoch_sort_key)
    if not checkpoint_path_strs:
        raise SystemExit(f"No checkpoints matched: {args.checkpoint_glob}")
    checkpoint_paths = [Path(checkpoint_path_str).resolve() for checkpoint_path_str in checkpoint_path_strs]

    rows: list[dict[str, Any]] = []
    best_row: dict[str, Any] | None = None
    best_value = float("nan")

    evaluate_checkpoints(checkpoint_paths, out_dir, args)

    for checkpoint_path in checkpoint_paths:
        output_json = output_json_for(checkpoint_path, out_dir, args.split)
        payload = json.loads(output_json.read_text(encoding="utf-8"))
        metric_value = extract_metric(payload, args.metric)
        mean_jsd = extract_metric(payload, "results.peak.profile_target_jsd_full_mean")
        median_jsd = extract_metric(payload, "results.peak.profile_target_jsd_full_median")
        count_r = extract_metric(payload, "results.peak.count_pearson_full")
        cls_auroc = extract_metric(payload, "results.overall.peak_auroc_logcount_full")
        cls_auprc = extract_metric(payload, "results.overall.peak_auprc_logcount_full")
        cls_f1 = extract_metric(payload, "results.overall.peak_best_f1_logcount_full")

        match = EPOCH_RE.search(checkpoint_path.name)
        row = {
            "epoch": int(match.group(1)) if match else "",
            "checkpoint_path": str(checkpoint_path),
            "eval_json": str(output_json),
            "selection_metric": args.metric,
            "selection_mode": args.mode,
            "selection_value": metric_value,
            "peak_profile_mean_jsd": mean_jsd,
            "peak_profile_median_jsd": median_jsd,
            "peak_count_pearson_full": count_r,
            "peak_vs_nonpeak_auroc_logcount_full": cls_auroc,
            "peak_vs_nonpeak_auprc_logcount_full": cls_auprc,
            "peak_vs_nonpeak_best_f1_logcount_full": cls_f1,
        }
        rows.append(row)

        if is_better(metric_value, best_value, args.mode):
            best_row = row
            best_value = metric_value

    if best_row is None:
        raise SystemExit("No valid metric values were produced; cannot select a best checkpoint")

    write_csv(rows, out_dir / "epoch_metrics.csv")
    (out_dir / "best_epoch.json").write_text(
        json.dumps(
            {
                "metric": args.metric,
                "mode": args.mode,
                "split": args.split,
                "n_checkpoints": len(rows),
                "best": best_row,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"Wrote epoch summary to {out_dir / 'epoch_metrics.csv'}")
    print(f"Wrote best checkpoint payload to {out_dir / 'best_epoch.json'}")
    print(
        "Best checkpoint: "
        f"epoch={best_row['epoch']} value={best_row['selection_value']:.6f} "
        f"ckpt={best_row['checkpoint_path']}"
    )


if __name__ == "__main__":
    main()
