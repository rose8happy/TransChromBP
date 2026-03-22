"""Benchmark wrapper for train_ddp — adds step-level timing and optional profiler.

Usage:
    torchrun ... -m transchrombp.training.train_ddp_bench [same args as train_ddp]

Env vars:
    BENCH_TIMING=1       Enable step-level timing breakdown (default: 1)
    ENABLE_PROFILER=1    Enable PyTorch profiler for first N steps
    PROFILER_OUTPUT_DIR  Where to save profiler traces (default: ./profiler_traces)

This module patches train_ddp.main() to inject timing around the training loop.
It does NOT modify the original train_ddp.py.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from collections import defaultdict
from contextlib import nullcontext
from pathlib import Path
from typing import Any, Dict, List, Optional

import torch
import torch.distributed as dist
import torch.nn.functional as F
from torch import nn
from torch.nn.parallel import DistributedDataParallel as DDP

# Re-use everything from original train_ddp
from transchrombp.training.train_ddp import (
    BIAS_DIAGNOSTIC_KEYS,
    LOSS_METRIC_KEYS,
    auto_configure_loss_weights,
    build_dataloaders,
    build_model,
    build_optimizer,
    build_scheduler,
    check_profile_resolution_compatibility,
    check_sequence_length_compatibility,
    compute_bias_diagnostics,
    compute_losses,
    format_bias_diag_suffix,
    init_distributed,
    is_main_process,
    load_yaml,
    maybe_compile,
    reduce_weighted_metrics,
    resolve_amp_dtype,
    resolve_data_config_path,
    apply_external_data_config_defaults,
    apply_training_mode_defaults,
    save_checkpoint,
    set_seed,
    to_scalar,
    validate_model_config,
)


def _print_env_info():
    """Print runtime environment for verification."""
    print("=" * 60)
    print("=== Runtime Environment ===")
    print(f"PyTorch version: {torch.__version__}")
    print(f"CUDA version: {torch.version.cuda}")
    try:
        print(f"cuDNN version: {torch.backends.cudnn.version()}")
    except Exception:
        print("cuDNN version: unavailable")
    print(f"GPU count: {torch.cuda.device_count()}")
    for i in range(torch.cuda.device_count()):
        props = torch.cuda.get_device_properties(i)
        print(f"  GPU {i}: {props.name} ({props.total_memory / 1e9:.1f} GB)")
    print(f"torch.backends.cuda.matmul.allow_tf32: {torch.backends.cuda.matmul.allow_tf32}")
    print(f"torch.backends.cudnn.allow_tf32: {torch.backends.cudnn.allow_tf32}")
    try:
        nccl_ver = torch.cuda.nccl.version()
        print(f"NCCL version: {nccl_ver}")
    except Exception:
        print("NCCL version: unavailable")
    print("=" * 60)


class TimingStats:
    """Lightweight step-level timing accumulator."""

    def __init__(self):
        self.data_times: List[float] = []
        self.fwd_bwd_times: List[float] = []
        self.opt_times: List[float] = []
        self.step_times: List[float] = []
        self.batch_sizes: List[int] = []

    def record(self, data_t: float, fwd_bwd_t: float, opt_t: float, step_t: float, bs: int):
        self.data_times.append(data_t)
        self.fwd_bwd_times.append(fwd_bwd_t)
        self.opt_times.append(opt_t)
        self.step_times.append(step_t)
        self.batch_sizes.append(bs)

    def rolling_avg(self, window: int = 10):
        """Return rolling average of last `window` steps."""
        w = min(window, len(self.step_times))
        if w == 0:
            return {}
        return {
            "data": sum(self.data_times[-w:]) / w,
            "fwd_bwd": sum(self.fwd_bwd_times[-w:]) / w,
            "opt": sum(self.opt_times[-w:]) / w,
            "step": sum(self.step_times[-w:]) / w,
            "samples_per_sec": sum(self.batch_sizes[-w:]) / sum(self.step_times[-w:]),
        }

    def summary(self, skip_first: int = 10):
        """Return overall averages, skipping warmup steps."""
        n = len(self.step_times)
        s = min(skip_first, n // 2)
        if n - s <= 0:
            s = 0
        data = self.data_times[s:]
        fwd_bwd = self.fwd_bwd_times[s:]
        opt = self.opt_times[s:]
        step = self.step_times[s:]
        bs = self.batch_sizes[s:]
        count = len(step)
        if count == 0:
            return {}
        return {
            "n_steps": count,
            "skipped_warmup": s,
            "avg_data": sum(data) / count,
            "avg_fwd_bwd": sum(fwd_bwd) / count,
            "avg_opt": sum(opt) / count,
            "avg_step": sum(step) / count,
            "avg_samples_per_sec": sum(bs) / sum(step),
            "pct_data": sum(data) / sum(step) * 100,
            "pct_fwd_bwd": sum(fwd_bwd) / sum(step) * 100,
            "pct_opt": sum(opt) / sum(step) * 100,
        }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark DDP trainer for TransChromBP")
    parser.add_argument("--model-config", type=str, default="configs/model/transchrombp_base.yaml")
    parser.add_argument("--train-config", type=str, default="configs/train/train_benchmark_base.yaml")
    parser.add_argument("--data-config", type=str, default="")
    parser.add_argument("--output-dir", type=str, default="")
    parser.add_argument("--run-name", type=str, default="")
    parser.add_argument("--max-epochs", type=int, default=0)
    parser.add_argument("--dry-run-steps", type=int, default=300)
    parser.add_argument("--batch-size-per-gpu", type=int, default=0)
    parser.add_argument("--num-workers", type=int, default=-1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    model_cfg = load_yaml(args.model_config)
    train_cfg = load_yaml(args.train_config)
    data_config_path = resolve_data_config_path(train_cfg, args.data_config)
    data_source_cfg = load_yaml(data_config_path) if data_config_path else {}
    if data_source_cfg:
        apply_external_data_config_defaults(train_cfg, data_source_cfg)
    apply_training_mode_defaults(train_cfg)

    if args.max_epochs > 0:
        train_cfg["max_epochs"] = args.max_epochs
    if args.batch_size_per_gpu > 0:
        train_cfg.setdefault("data", {})["batch_size_per_gpu"] = args.batch_size_per_gpu
    if args.num_workers >= 0:
        train_cfg.setdefault("data", {})["num_workers"] = args.num_workers

    trainer_cfg = train_cfg.get("trainer", {})
    backend = str(trainer_cfg.get("backend", "nccl"))
    dist_env = init_distributed(backend=backend)

    enable_timing = os.environ.get("BENCH_TIMING", "1") == "1"
    enable_profiler = os.environ.get("ENABLE_PROFILER", "0") == "1"
    profiler_output_dir = os.environ.get("PROFILER_OUTPUT_DIR", "./profiler_traces")

    if is_main_process(dist_env.rank):
        _print_env_info()

    try:
        set_seed(int(train_cfg.get("seed", 1234)), dist_env.rank)
        validate_model_config(model_cfg)
        check_sequence_length_compatibility(model_cfg, train_cfg, dist_env.rank)
        check_profile_resolution_compatibility(model_cfg, train_cfg, dist_env.rank)

        output_dir = Path(args.output_dir or train_cfg.get("logging", {}).get("output_dir", "outputs"))
        run_name = args.run_name or train_cfg.get("logging", {}).get("run_name", "benchmark")

        model, training_mode = build_model(model_cfg, train_cfg, dist_env.rank)
        model = model.to(dist_env.device)
        optimizer = build_optimizer(model, train_cfg, dist_env.rank)
        model = maybe_compile(model, train_cfg, dist_env.rank)

        if dist_env.distributed:
            if dist_env.device.type == "cuda":
                model = DDP(
                    model,
                    device_ids=[dist_env.local_rank],
                    output_device=dist_env.local_rank,
                    find_unused_parameters=bool(trainer_cfg.get("find_unused_parameters", False)),
                )
            else:
                model = DDP(model, find_unused_parameters=bool(trainer_cfg.get("find_unused_parameters", False)))

        train_loader, val_loader, train_sampler = build_dataloaders(train_cfg, dist_env, data_source_cfg)
        scheduler = build_scheduler(optimizer, train_cfg, steps_per_epoch=len(train_loader))

        precision = str(trainer_cfg.get("precision", "bf16"))
        amp_dtype = resolve_amp_dtype(precision)
        use_amp = (dist_env.device.type == "cuda") and (amp_dtype != torch.float32)
        use_scaler = use_amp and (amp_dtype == torch.float16)
        scaler = torch.amp.GradScaler(device=dist_env.device.type, enabled=use_scaler)

        grad_accum_steps = int(trainer_cfg.get("grad_accum_steps", 1))
        grad_clip = float(trainer_cfg.get("clip_grad_norm", 0.0))
        log_every = int(trainer_cfg.get("log_every_steps", 10))
        max_epochs = int(train_cfg.get("max_epochs", 999))
        dry_limit = int(args.dry_run_steps)

        batch_size_per_gpu = int(train_cfg.get("data", {}).get("batch_size_per_gpu", 16))
        world_size = max(1, dist_env.world_size)
        global_batch = batch_size_per_gpu * world_size * max(1, grad_accum_steps)

        if is_main_process(dist_env.rank):
            print(
                f"[bench-setup] device={dist_env.device} world_size={world_size} "
                f"precision={precision} use_amp={use_amp} amp_dtype={amp_dtype} "
                f"bs_per_gpu={batch_size_per_gpu} global_batch={global_batch} "
                f"dry_run_steps={dry_limit} timing={enable_timing} profiler={enable_profiler}"
            )

        # ── Profiler setup ──────────────────────────────────────
        profiler_ctx = nullcontext()
        prof = None
        if enable_profiler and dist_env.device.type == "cuda":
            from torch.profiler import ProfilerActivity, profile, schedule

            os.makedirs(profiler_output_dir, exist_ok=True)
            prof_schedule = schedule(wait=2, warmup=3, active=10, repeat=1)
            prof = profile(
                activities=[ProfilerActivity.CPU, ProfilerActivity.CUDA],
                schedule=prof_schedule,
                on_trace_ready=torch.profiler.tensorboard_trace_handler(profiler_output_dir),
                record_shapes=True,
                profile_memory=True,
                with_stack=False,
            )
            profiler_ctx = prof

        # ── Timing setup ────────────────────────────────────────
        timing = TimingStats() if enable_timing else None

        global_step = 0
        loss_calibration = auto_configure_loss_weights(train_cfg, train_loader.dataset, dist_env)
        optimizer.zero_grad(set_to_none=True)

        # ── Training loop (benchmark version) ───────────────────
        with profiler_ctx:
            for epoch in range(1, max_epochs + 1):
                model.train()
                if train_sampler is not None:
                    train_sampler.set_epoch(epoch)

                t_epoch = time.time()
                running = {key: 0.0 for key in LOSS_METRIC_KEYS + BIAS_DIAGNOSTIC_KEYS}
                running_batches = 0
                running_examples = 0.0

                t_data_start = time.perf_counter()

                for step_idx, batch in enumerate(train_loader, start=1):
                    # ── Data time ──
                    if enable_timing and dist_env.device.type == "cuda":
                        torch.cuda.synchronize()
                    t_data_end = time.perf_counter()
                    data_time = t_data_end - t_data_start

                    seq = batch["seq"].to(dist_env.device, non_blocking=True)
                    profile_counts = batch["profile_counts"].to(dist_env.device, non_blocking=True)

                    should_step = (step_idx % grad_accum_steps == 0)
                    no_sync_ctx = model.no_sync() if (isinstance(model, DDP) and not should_step) else nullcontext()

                    # ── Forward + Backward time ──
                    if enable_timing and dist_env.device.type == "cuda":
                        torch.cuda.synchronize()
                    t_fwd_start = time.perf_counter()

                    with no_sync_ctx:
                        with torch.autocast(
                            device_type=dist_env.device.type,
                            dtype=amp_dtype,
                            enabled=use_amp,
                        ):
                            outputs = model(seq)
                            loss, metrics = compute_losses(outputs, profile_counts, train_cfg.get("loss", {}))
                            diagnostics = compute_bias_diagnostics(model, outputs)
                            loss = loss / grad_accum_steps

                        if use_scaler:
                            scaler.scale(loss).backward()
                        else:
                            loss.backward()

                    if enable_timing and dist_env.device.type == "cuda":
                        torch.cuda.synchronize()
                    t_fwd_end = time.perf_counter()
                    fwd_bwd_time = t_fwd_end - t_fwd_start

                    # ── Optimizer time ──
                    t_opt_start = time.perf_counter()
                    if should_step:
                        if grad_clip > 0:
                            if use_scaler:
                                scaler.unscale_(optimizer)
                            nn.utils.clip_grad_norm_(model.parameters(), grad_clip)

                        if use_scaler:
                            scaler.step(optimizer)
                            scaler.update()
                        else:
                            optimizer.step()

                        optimizer.zero_grad(set_to_none=True)

                        if scheduler is not None:
                            scheduler.step()

                        global_step += 1

                    if enable_timing and dist_env.device.type == "cuda":
                        torch.cuda.synchronize()
                    t_opt_end = time.perf_counter()
                    opt_time = t_opt_end - t_opt_start

                    step_time = data_time + fwd_bwd_time + opt_time  # total step time

                    # ── Record timing ──
                    actual_bs = int(seq.shape[0]) * world_size
                    if timing is not None:
                        timing.record(data_time, fwd_bwd_time, opt_time, step_time, actual_bs)

                    # ── Running metrics ──
                    batch_weight = float(seq.shape[0])
                    running_batches += 1
                    running_examples += batch_weight
                    combined_metrics = {**metrics, **diagnostics}
                    for k in running:
                        running[k] += to_scalar(combined_metrics[k]) * batch_weight

                    # ── Logging ──
                    if (step_idx % log_every == 0) and is_main_process(dist_env.rank):
                        elapsed = time.time() - t_epoch
                        avg_loss = running["loss_total"] / max(1.0, running_examples)
                        lr = optimizer.param_groups[0]["lr"]

                        msg = (
                            f"[bench] epoch={epoch} step={step_idx}/{len(train_loader)} "
                            f"global_step={global_step} loss={avg_loss:.5f} lr={lr:.3e} "
                            f"time={elapsed:.1f}s"
                        )
                        if timing is not None:
                            ra = timing.rolling_avg(log_every)
                            msg += (
                                f" | data={ra['data']*1000:.1f}ms fwd_bwd={ra['fwd_bwd']*1000:.1f}ms "
                                f"opt={ra['opt']*1000:.1f}ms step={ra['step']*1000:.1f}ms "
                                f"samples/s={ra['samples_per_sec']:.1f}"
                            )
                        print(msg)

                    # ── Profiler step ──
                    if prof is not None:
                        prof.step()

                    # ── Dry run limit ──
                    if dry_limit > 0 and global_step >= dry_limit:
                        break

                    # ── Next data timing starts here ──
                    t_data_start = time.perf_counter()

                if dry_limit > 0 and global_step >= dry_limit:
                    if is_main_process(dist_env.rank):
                        print(f"[done] dry run limit reached: {dry_limit} steps")
                    break

        # ── Print timing summary ────────────────────────────────
        if timing is not None and is_main_process(dist_env.rank):
            summary = timing.summary(skip_first=10)
            if summary:
                print("")
                print("=" * 60)
                print("=== TIMING SUMMARY ===")
                print(f"Total steps measured: {summary['n_steps']} (skipped first {summary['skipped_warmup']} warmup)")
                print(f"avg_data_time:       {summary['avg_data']*1000:.2f} ms  ({summary['pct_data']:.1f}%)")
                print(f"avg_fwd_bwd_time:    {summary['avg_fwd_bwd']*1000:.2f} ms  ({summary['pct_fwd_bwd']:.1f}%)")
                print(f"avg_opt_time:        {summary['avg_opt']*1000:.2f} ms  ({summary['pct_opt']:.1f}%)")
                print(f"avg_step_time:       {summary['avg_step']*1000:.2f} ms")
                print(f"avg_samples_per_sec: {summary['avg_samples_per_sec']:.1f}")
                print(f"world_size:          {world_size}")
                print(f"batch_size_per_gpu:  {batch_size_per_gpu}")
                print(f"global_batch_size:   {global_batch}")
                print("=" * 60)

        # ── Print profiler summary ──────────────────────────────
        if prof is not None and is_main_process(dist_env.rank):
            print("")
            print("=" * 60)
            print("=== PROFILER SUMMARY (top 20 by CUDA time) ===")
            print(prof.key_averages().table(sort_by="cuda_time_total", row_limit=20))
            print("")
            print("=== PROFILER SUMMARY (top 10 by CPU time) ===")
            print(prof.key_averages().table(sort_by="cpu_time_total", row_limit=10))
            print("=" * 60)

        # ── GPU memory summary ──────────────────────────────────
        if dist_env.device.type == "cuda" and is_main_process(dist_env.rank):
            print("")
            print("=== GPU MEMORY SUMMARY ===")
            for i in range(torch.cuda.device_count()):
                allocated = torch.cuda.max_memory_allocated(i) / 1e9
                reserved = torch.cuda.max_memory_reserved(i) / 1e9
                print(f"GPU {i}: peak_allocated={allocated:.2f} GB, peak_reserved={reserved:.2f} GB")

    except Exception as e:
        print(f"[bench-error] rank={dist_env.rank} {type(e).__name__}: {e}")
        raise
    finally:
        if dist_env.distributed:
            dist.destroy_process_group()


if __name__ == "__main__":
    main()
