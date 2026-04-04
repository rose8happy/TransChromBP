#!/usr/bin/env python3
"""Build CSV assets for the 2026-03-20 training performance benchmark."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parent
SUMMARY_PATH = ROOT / "benchmark_summary_20260320_004701.txt"


@dataclass
class RunRecord:
    run_name: str
    nproc: int
    bs_per_gpu: int
    world_size: int
    global_batch_size: int
    avg_data_ms: float
    avg_fwd_bwd_ms: float
    avg_opt_ms: float
    avg_step_ms: float
    samples_per_sec: float
    pct_data: float
    pct_fwd_bwd: float
    pct_opt: float
    peak_allocated_gb: float
    sample_peak_mem_mib: float | None = None
    sample_mean_util_pct: float | None = None
    sample_peak_util_pct: float | None = None
    sample_mean_power_w: float | None = None
    sample_peak_power_w: float | None = None
    active_sample_count: int = 0


def _extract_float(block: str, label: str) -> float:
    pattern = rf"{re.escape(label)}:\s+([0-9.]+)"
    match = re.search(pattern, block)
    if not match:
        raise ValueError(f"missing {label!r}")
    return float(match.group(1))


def parse_summary(path: Path) -> list[RunRecord]:
    text = path.read_text(encoding="utf-8")
    pattern = re.compile(
        r"--- (?P<run_name>\S+) \(nproc=(?P<nproc>\d+), bs=(?P<bs>\d+)\) ---\n(?P<body>.*?)(?=\n--- |\Z)",
        re.DOTALL,
    )
    runs: list[RunRecord] = []
    for match in pattern.finditer(text):
        run_name = match.group("run_name")
        body = match.group("body")
        gpu_mem_match = re.search(r"GPU 0: peak_allocated=([0-9.]+) GB", body)
        if not gpu_mem_match:
            raise ValueError(f"missing GPU memory summary for {run_name}")
        runs.append(
            RunRecord(
                run_name=run_name,
                nproc=int(match.group("nproc")),
                bs_per_gpu=int(match.group("bs")),
                world_size=int(_extract_float(body, "world_size")),
                global_batch_size=int(_extract_float(body, "global_batch_size")),
                avg_data_ms=_extract_float(body, "avg_data_time"),
                avg_fwd_bwd_ms=_extract_float(body, "avg_fwd_bwd_time"),
                avg_opt_ms=_extract_float(body, "avg_opt_time"),
                avg_step_ms=_extract_float(body, "avg_step_time"),
                samples_per_sec=_extract_float(body, "avg_samples_per_sec"),
                pct_data=float(re.search(r"avg_data_time:\s+[0-9.]+ ms\s+\(([0-9.]+)%\)", body).group(1)),
                pct_fwd_bwd=float(re.search(r"avg_fwd_bwd_time:\s+[0-9.]+ ms\s+\(([0-9.]+)%\)", body).group(1)),
                pct_opt=float(re.search(r"avg_opt_time:\s+[0-9.]+ ms\s+\(([0-9.]+)%\)", body).group(1)),
                peak_allocated_gb=float(gpu_mem_match.group(1)),
            )
        )
    return runs


def enrich_with_gpu_samples(runs: list[RunRecord]) -> None:
    for run in runs:
        sample_path = ROOT / f"{run.run_name}_gpu_samples_20260320_004701.csv"
        if not sample_path.exists():
            continue
        rows = []
        with sample_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                rows.append(
                    {
                        "mem_used_mib": float(row["mem_used_MiB"].strip()),
                        "gpu_util_pct": float(row["gpu_util_pct"].strip()),
                        "power_w": float(row["power_W"].strip()),
                    }
                )
        active_rows = [
            row
            for row in rows
            if row["mem_used_mib"] > 1000 or row["gpu_util_pct"] >= 10 or row["power_w"] >= 100
        ]
        if not active_rows:
            continue
        run.active_sample_count = len(active_rows)
        run.sample_peak_mem_mib = max(row["mem_used_mib"] for row in active_rows)
        run.sample_mean_util_pct = mean(row["gpu_util_pct"] for row in active_rows)
        run.sample_peak_util_pct = max(row["gpu_util_pct"] for row in active_rows)
        run.sample_mean_power_w = mean(row["power_w"] for row in active_rows)
        run.sample_peak_power_w = max(row["power_w"] for row in active_rows)


def write_csv(path: Path, runs: list[RunRecord]) -> None:
    fieldnames = [
        "run_name",
        "nproc",
        "bs_per_gpu",
        "world_size",
        "global_batch_size",
        "avg_step_ms",
        "samples_per_sec",
        "pct_data",
        "pct_fwd_bwd",
        "pct_opt",
        "peak_allocated_gb",
        "sample_peak_mem_mib",
        "sample_mean_util_pct",
        "sample_peak_util_pct",
        "sample_mean_power_w",
        "sample_peak_power_w",
        "active_sample_count",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for run in runs:
            writer.writerow({name: getattr(run, name) for name in fieldnames})


def main() -> None:
    runs = parse_summary(SUMMARY_PATH)
    enrich_with_gpu_samples(runs)
    write_csv(ROOT / "run_summary.csv", runs)
    write_csv(
        ROOT / "batch_sweep_2gpu.csv",
        [run for run in runs if run.world_size == 2 and "profiler" not in run.run_name.lower()],
    )
    write_csv(
        ROOT / "scaling_baseline.csv",
        [
            run
            for run in runs
            if (run.world_size == 1 and run.bs_per_gpu == 16)
            or (run.world_size == 2 and run.bs_per_gpu == 16 and "profiler" not in run.run_name.lower())
        ],
    )


if __name__ == "__main__":
    main()
