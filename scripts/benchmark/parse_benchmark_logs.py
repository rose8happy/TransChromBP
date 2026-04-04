#!/usr/bin/env python3
"""Parse benchmark logs and produce a summary table.

Usage:
    python3 scripts/benchmark/parse_benchmark_logs.py outputs/logs/benchmark/bench_*.log
    python3 scripts/benchmark/parse_benchmark_logs.py outputs/logs/benchmark/  # parse all logs in dir
"""

import re
import sys
from pathlib import Path
from typing import Dict, List, Optional


def parse_timing_summary(text: str) -> Optional[Dict[str, str]]:
    """Extract TIMING SUMMARY block from log text."""
    m = re.search(r"=== TIMING SUMMARY ===(.+?)={10,}", text, re.DOTALL)
    if not m:
        return None
    block = m.group(1)
    result = {}
    for line in block.strip().splitlines():
        line = line.strip()
        # e.g. "avg_step_time:       12.34 ms"
        match = re.match(r"([A-Za-z0-9_ ]+):\s+(.+)", line)
        if match:
            key = match.group(1).strip().replace(" ", "_")
            value = match.group(2).strip()
            pct_match = re.search(r"\(([\d.]+)%\)", value)
            if pct_match:
                pct_key = {
                    "avg_data_time": "pct_data",
                    "avg_fwd_bwd_time": "pct_fwd_bwd",
                    "avg_opt_time": "pct_opt",
                }.get(key)
                if pct_key:
                    result[pct_key] = f"{pct_match.group(1)}%"
                value = re.sub(r"\s+\([\d.]+%\)$", "", value)
            result[key] = value
    return result


def parse_gpu_memory(text: str) -> Optional[str]:
    """Extract peak GPU memory from GPU MEMORY SUMMARY block."""
    m = re.search(r"=== GPU MEMORY SUMMARY ===(.+?)(?:={10,}|\Z)", text, re.DOTALL)
    if not m:
        return None
    block = m.group(1).strip()
    # Find max peak_allocated across GPUs
    peaks = re.findall(r"peak_allocated=(\d+\.\d+)\s*GB", block)
    if peaks:
        max_peak = max(float(p) for p in peaks)
        return f"{max_peak:.2f} GB"
    return block.split("\n")[0] if block else None


def parse_setup_info(text: str) -> Dict[str, str]:
    """Extract setup info from [bench-setup] line."""
    m = re.search(r"\[bench-setup\]\s+(.+)", text)
    if not m:
        return {}
    info = {}
    for kv in re.findall(r"(\w+)=(\S+)", m.group(1)):
        info[kv[0]] = kv[1]
    return info


def check_oom(text: str) -> bool:
    """Check if log indicates OOM."""
    oom_patterns = [
        "CUDA out of memory",
        "OutOfMemoryError",
        "RuntimeError: CUDA error",
        "torch.cuda.OutOfMemoryError",
    ]
    return any(p in text for p in oom_patterns)


def parse_log_file(path: Path) -> Dict[str, str]:
    """Parse a single benchmark log file."""
    text = path.read_text(errors="replace")
    result = {"file": path.name}

    # Check OOM
    if check_oom(text):
        result["OOM"] = "YES"
        return result
    result["OOM"] = "no"

    # Setup info
    setup = parse_setup_info(text)
    result.update(setup)

    # Timing summary
    timing = parse_timing_summary(text)
    if timing:
        result.update(timing)

    # GPU memory
    mem = parse_gpu_memory(text)
    if mem:
        result["peak_mem"] = mem

    return result


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <log_file_or_dir> [...]")
        sys.exit(1)

    log_files: List[Path] = []
    for arg in sys.argv[1:]:
        p = Path(arg)
        if p.is_dir():
            log_files.extend(sorted(p.glob("bench_*.log")))
        elif p.is_file():
            log_files.append(p)
        else:
            print(f"Warning: {arg} not found, skipping")

    if not log_files:
        print("No log files found.")
        sys.exit(1)

    results = [parse_log_file(f) for f in log_files]

    # Print summary table
    print("\n" + "=" * 100)
    print("BENCHMARK SUMMARY TABLE")
    print("=" * 100)

    header = f"{'run_name':<30} {'OOM':<5} {'world':<6} {'bs/gpu':<7} {'avg_step(ms)':<13} {'samples/s':<10} {'peak_mem':<12} {'pct_data':<9} {'pct_fwd':<9} {'pct_opt':<9}"
    print(header)
    print("-" * len(header))

    for r in results:
        name = r.get("file", "?")
        # Simplify name
        name = re.sub(r"_\d{8}_\d{6}\.log$", "", name)
        oom = r.get("OOM", "?")
        ws = r.get("world_size", "?")
        bs = r.get("batch_size_per_gpu", r.get("bs_per_gpu", "?"))
        step = r.get("avg_step_time", "?")
        sps = r.get("avg_samples_per_sec", "?")
        mem = r.get("peak_mem", "?")

        # Extract percentages
        pct_d = r.get("pct_data", "?")
        pct_f = r.get("pct_fwd_bwd", "?")
        pct_o = r.get("pct_opt", "?")

        # Clean up values (remove units for alignment)
        if isinstance(step, str) and "ms" in step:
            step = step.replace(" ms", "").strip()
            step = f"{step} ms"

        print(f"{name:<30} {oom:<5} {ws:<6} {bs:<7} {step:<13} {sps:<10} {mem:<12} {pct_d:<9} {pct_f:<9} {pct_o:<9}")

    # Scaling analysis
    print("\n" + "=" * 100)
    print("SCALING ANALYSIS")
    print("=" * 100)

    # Find 1gpu and 2gpu runs with same bs
    by_bs: Dict[str, Dict[str, float]] = {}
    for r in results:
        if "profiler" in r.get("file", "").lower():
            continue
        if r.get("OOM") == "YES":
            continue
        sps = r.get("avg_samples_per_sec")
        if not sps:
            continue
        try:
            sps_val = float(sps)
        except (ValueError, TypeError):
            continue
        ws = r.get("world_size", "1")
        bs = r.get("batch_size_per_gpu", r.get("bs_per_gpu", "?"))
        key = f"bs{bs}"
        by_bs.setdefault(key, {})[f"{ws}gpu"] = sps_val

    for bs_key, gpu_data in sorted(by_bs.items()):
        if "1gpu" in gpu_data and "2gpu" in gpu_data:
            s1 = gpu_data["1gpu"]
            s2 = gpu_data["2gpu"]
            speedup = s2 / s1 if s1 > 0 else 0
            efficiency = speedup / 2
            print(f"{bs_key}: 1GPU={s1:.1f} sps, 2GPU={s2:.1f} sps, speedup={speedup:.2f}x, efficiency={efficiency:.1%}")
            if efficiency > 0.85:
                print(f"  → Scaling is GOOD (efficiency > 85%)")
            elif efficiency > 0.70:
                print(f"  → Scaling is ACCEPTABLE but has room for improvement")
            else:
                print(f"  → Scaling is POOR — investigate communication/sync overhead")

    # Batch sweep throughput comparison
    print("\n" + "=" * 100)
    print("BATCH SIZE SWEEP (2 GPU)")
    print("=" * 100)
    print(f"{'bs_per_gpu':<12} {'samples/s':<12} {'step_time(ms)':<15} {'peak_mem':<12} {'improvement':<12}")
    print("-" * 63)

    baseline_sps = None
    for r in sorted(results, key=lambda x: int(x.get("batch_size_per_gpu", x.get("bs_per_gpu", "0")))):
        if "profiler" in r.get("file", "").lower():
            continue
        ws = r.get("world_size", "1")
        if ws != "2":
            continue
        bs = r.get("batch_size_per_gpu", r.get("bs_per_gpu", "?"))
        sps = r.get("avg_samples_per_sec", "?")
        step = r.get("avg_step_time", "?")
        mem = r.get("peak_mem", "?")
        oom = r.get("OOM", "no")

        if oom == "YES":
            print(f"{bs:<12} {'OOM':<12} {'---':<15} {'---':<12} {'---':<12}")
            continue

        try:
            sps_val = float(sps)
        except (ValueError, TypeError):
            sps_val = None

        improvement = "baseline"
        if baseline_sps is None and sps_val is not None:
            baseline_sps = sps_val
        elif baseline_sps is not None and sps_val is not None:
            pct = (sps_val / baseline_sps - 1) * 100
            improvement = f"{pct:+.1f}%"

        print(f"{bs:<12} {sps:<12} {step:<15} {mem:<12} {improvement:<12}")

    print("\n" + "=" * 100)


if __name__ == "__main__":
    main()
