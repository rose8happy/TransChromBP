#!/usr/bin/env python3
"""Evaluate per-epoch ChromBPNet checkpoints and select the external best epoch."""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

EPOCH_RE = re.compile(r"\.epoch_(\d+)\.h5$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate ChromBPNet epoch checkpoints with a unified external metric and select the best one."
    )
    parser.add_argument("--model-glob", required=True, help="Glob for epoch checkpoints, e.g. /path/chrombpnet.epoch_*.h5")
    parser.add_argument("--genome", required=True, help="Genome fasta")
    parser.add_argument("--bigwig", required=True, help="Observed shifted/unstranded bigWig")
    parser.add_argument("--peaks", default="None", help="Peak BED for evaluation")
    parser.add_argument("--nonpeaks", default="None", help="Nonpeak BED for evaluation")
    parser.add_argument("--fold-json", required=True, help="Fold json used for evaluation")
    parser.add_argument("--output-dir", required=True, help="Directory for per-epoch metrics and summary files")
    parser.add_argument("--metric", default="profile_metrics.peaks.median_jsd", help="Metric path used to select best epoch")
    parser.add_argument("--mode", choices=["min", "max"], default="min", help="Whether lower or higher metric is better")
    parser.add_argument("--split", choices=["train", "valid", "test"], default="valid", help="Fold split used for external checkpoint selection")
    parser.add_argument(
        "--official-root",
        default="",
        help="Path to the official ChromBPNet repo root; defaults to $CHROMBPNET_OFFICIAL_ROOT",
    )
    parser.add_argument("--batch-size", type=int, default=512, help="Predict batch size")
    parser.add_argument("--seed", type=int, default=1234, help="Predict seed")
    parser.add_argument("--inputlen", type=int, default=2114, help="Input sequence length")
    parser.add_argument("--outputlen", type=int, default=1000, help="Prediction output length")
    parser.add_argument("--gpus", default="", help="Comma-separated GPU ids used to shard checkpoint evaluation, e.g. 0,1")
    parser.add_argument("--max-workers", type=int, default=0, help="Optional worker cap for --gpus; 0 means use all listed GPUs")
    parser.add_argument("--force", action="store_true", help="Recompute epoch metrics even if JSON already exists")
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


def resolve_official_root(args: argparse.Namespace) -> tuple[Path, Path]:
    raw_root = args.official_root or os.environ.get("CHROMBPNET_OFFICIAL_ROOT", "")
    if not raw_root:
        raise SystemExit("ERROR: missing official ChromBPNet root; pass --official-root or set CHROMBPNET_OFFICIAL_ROOT")

    official_root = Path(raw_root).expanduser().resolve()
    predict_py = official_root / "chrombpnet" / "training" / "predict.py"
    if not predict_py.is_file():
        raise SystemExit(
            f"ERROR: official predict.py not found: {predict_py} "
            "(set --official-root or CHROMBPNET_OFFICIAL_ROOT to the official ChromBPNet repo root)"
        )
    return official_root, predict_py


def build_predict_cmd(
    model_path: Path,
    output_prefix: Path,
    args: argparse.Namespace,
    predict_py: Path,
) -> list[str]:
    return [
        sys.executable,
        str(predict_py),
        "-g",
        args.genome,
        "-b",
        args.bigwig,
        "-p",
        args.peaks,
        "-n",
        args.nonpeaks,
        "-o",
        str(output_prefix),
        "-fl",
        args.fold_json,
        "-m",
        str(model_path),
        "-bs",
        str(args.batch_size),
        "-s",
        str(args.seed),
        "-il",
        str(args.inputlen),
        "-ol",
        str(args.outputlen),
        "--split",
        args.split,
        "--metrics-only",
    ]


def parse_gpu_list(raw: str) -> list[str]:
    return [item.strip() for item in raw.split(",") if item.strip()]


def metrics_path_for(model_path: Path, out_dir: Path) -> Path:
    return Path(f"{out_dir / model_path.stem}_metrics.json")


def run_predict_subprocess(model_path: Path, out_dir: Path, args: argparse.Namespace, gpu: str, predict_py: Path) -> None:
    output_prefix = out_dir / model_path.stem
    env = os.environ.copy()
    if gpu:
        env["CUDA_VISIBLE_DEVICES"] = gpu
    env["CHROMBPNET_MULTI_GPU"] = "0"
    env["PYTHONPATH"] = (
        f"{args.official_root}:{env['PYTHONPATH']}" if env.get("PYTHONPATH") else str(args.official_root)
    )
    print(f"[selector] gpu={gpu} evaluate {model_path.name}")
    subprocess.run(
        build_predict_cmd(model_path, output_prefix, args, predict_py),
        check=True,
        env=env,
    )


def evaluate_assigned_models(
    model_paths: list[Path],
    out_dir: Path,
    args: argparse.Namespace,
    gpu: str,
    predict_py: Path,
) -> None:
    for model_path in model_paths:
        metrics_path = metrics_path_for(model_path, out_dir)
        if args.force or not metrics_path.exists():
            run_predict_subprocess(model_path, out_dir, args, gpu, predict_py)


def evaluate_models(
    model_paths: list[Path],
    out_dir: Path,
    args: argparse.Namespace,
    predict_py: Path,
) -> None:
    gpu_list = parse_gpu_list(args.gpus)
    gpu = gpu_list[0] if gpu_list else ""
    if len(gpu_list) <= 1:
        for model_path in model_paths:
            metrics_path = metrics_path_for(model_path, out_dir)
            if args.force or not metrics_path.exists():
                run_predict_subprocess(model_path, out_dir, args, gpu, predict_py)
        return

    worker_count = len(gpu_list) if args.max_workers <= 0 else min(args.max_workers, len(gpu_list))
    worker_gpus = gpu_list[:worker_count]
    assignments: dict[str, list[Path]] = {gpu: [] for gpu in worker_gpus}
    for idx, model_path in enumerate(model_paths):
        gpu = worker_gpus[idx % worker_count]
        assignments[gpu].append(model_path)

    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = [
            executor.submit(evaluate_assigned_models, assigned_paths, out_dir, args, gpu, predict_py)
            for gpu, assigned_paths in assignments.items()
            if assigned_paths
        ]
        for future in futures:
            future.result()


def write_csv(rows: list[dict[str, Any]], out_path: Path) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def is_better(candidate: float, current_best: float, mode: str) -> bool:
    if math.isnan(current_best):
        return not math.isnan(candidate)
    if math.isnan(candidate):
        return False
    if mode == "max":
        return candidate > current_best
    return candidate < current_best


def main() -> None:
    args = parse_args()
    official_root, predict_py = resolve_official_root(args)
    args.official_root = str(official_root)

    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    model_path_strs = sorted(glob.glob(args.model_glob), key=epoch_sort_key)
    if not model_path_strs:
        raise SystemExit(f"No checkpoints matched: {args.model_glob}")
    model_paths = [Path(model_path_str).resolve() for model_path_str in model_path_strs]

    rows: list[dict[str, Any]] = []
    best_row: dict[str, Any] | None = None
    best_value = float("nan")

    evaluate_models(model_paths, out_dir, args, predict_py)

    for model_path in model_paths:
        metrics_path = metrics_path_for(model_path, out_dir)
        payload = json.loads(metrics_path.read_text(encoding="utf-8"))
        metric_value = extract_metric(payload, args.metric)
        mean_jsd = extract_metric(payload, "profile_metrics.peaks.mean_jsd")
        median_jsd = extract_metric(payload, "profile_metrics.peaks.median_jsd")
        count_r = extract_metric(payload, "counts_metrics.peaks.pearsonr")
        cls_auroc = extract_metric(payload, "classification_metrics.peaks_vs_nonpeaks.logcounts.auroc")
        cls_auprc = extract_metric(payload, "classification_metrics.peaks_vs_nonpeaks.logcounts.auprc")
        cls_f1 = extract_metric(payload, "classification_metrics.peaks_vs_nonpeaks.logcounts.best_f1")

        epoch_match = EPOCH_RE.search(model_path.name)
        row = {
            "epoch": int(epoch_match.group(1)) if epoch_match else "",
            "model_path": str(model_path),
            "metrics_path": str(metrics_path),
            "selection_metric": args.metric,
            "selection_mode": args.mode,
            "selection_value": metric_value,
            "peak_profile_mean_jsd": mean_jsd,
            "peak_profile_median_jsd": median_jsd,
            "peak_count_pearsonr": count_r,
            "peak_vs_nonpeak_auroc_logcounts": cls_auroc,
            "peak_vs_nonpeak_auprc_logcounts": cls_auprc,
            "peak_vs_nonpeak_best_f1_logcounts": cls_f1,
        }
        rows.append(row)

        if is_better(metric_value, best_value, args.mode):
            best_row = row
            best_value = metric_value

    if best_row is None:
        raise SystemExit("No valid metric values were produced; cannot select a best epoch")

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
    print(f"Wrote best epoch payload to {out_dir / 'best_epoch.json'}")
    print(
        "Best epoch: "
        f"epoch={best_row['epoch']} value={best_row['selection_value']:.6f} "
        f"model={best_row['model_path']}"
    )


if __name__ == "__main__":
    main()
