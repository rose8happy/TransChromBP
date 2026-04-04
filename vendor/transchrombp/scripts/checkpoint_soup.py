#!/usr/bin/env python3
"""Checkpoint soup: average weights from multiple checkpoints.

Usage:
  # 最后 5 个 epoch 的 weight average
  python scripts/checkpoint_soup.py \
    --ckpt-dir outputs/checkpoints/run_name/ \
    --mode last_k --k 5 \
    --output outputs/checkpoints/run_name/soup_last5.pt

  # validation JSD 最好的 3 个 epoch
  python scripts/checkpoint_soup.py \
    --ckpt-dir outputs/checkpoints/run_name/ \
    --mode top_k_jsd --k 3 \
    --metrics-jsonl outputs/logs/run_name/epoch_metrics.jsonl \
    --output outputs/checkpoints/run_name/soup_top3.pt
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import torch


def load_model_state(ckpt_path: Path) -> Dict[str, torch.Tensor]:
    """Load model_state from a checkpoint file."""
    ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    if "model_state" in ckpt:
        return ckpt["model_state"]
    # 兼容只有 state_dict 的情况
    return ckpt


def average_states(state_dicts: List[Dict[str, torch.Tensor]]) -> Dict[str, torch.Tensor]:
    """Element-wise average of multiple state dicts."""
    avg: Dict[str, torch.Tensor] = {}
    n = len(state_dicts)
    for key in state_dicts[0]:
        tensors = [sd[key].float() for sd in state_dicts]
        avg[key] = torch.stack(tensors).mean(dim=0)
    return avg


def select_last_k(ckpt_dir: Path, k: int) -> List[Path]:
    """Select the last K checkpoints by epoch number."""
    all_ckpts = sorted(ckpt_dir.glob("epoch_*.pt"))
    if not all_ckpts:
        raise FileNotFoundError(f"No epoch_*.pt checkpoints found in {ckpt_dir}")
    selected = all_ckpts[-k:]
    return selected


def select_top_k_jsd(
    ckpt_dir: Path,
    metrics_jsonl: Path,
    k: int,
    metric_key: str = "peak.profile_target_jsd_debiased_mean",
) -> List[Path]:
    """Select K checkpoints with best (lowest) validation JSD."""
    if not metrics_jsonl.exists():
        raise FileNotFoundError(f"Metrics file not found: {metrics_jsonl}")

    records: List[Dict[str, Any]] = []
    with open(metrics_jsonl) as f:
        for line in f:
            line = line.strip()
            if line:
                records.append(json.loads(line))

    # 解析 metric_key（支持 "peak.xxx" 格式）
    scored = []
    for rec in records:
        epoch = rec.get("epoch")
        val = rec.get("val", {})
        parts = metric_key.split(".", 1)
        if len(parts) == 2:
            value = val.get(parts[0], {}).get(parts[1])
        else:
            value = val.get(parts[0])

        if value is not None and epoch is not None:
            scored.append((float(value), int(epoch)))

    if not scored:
        raise ValueError(f"No valid values found for metric '{metric_key}' in {metrics_jsonl}")

    # 按 JSD 升序（越小越好）
    scored.sort(key=lambda x: x[0])
    top_epochs = [epoch for _, epoch in scored[:k]]

    selected = []
    for epoch in top_epochs:
        ckpt_path = ckpt_dir / f"epoch_{epoch:03d}.pt"
        if ckpt_path.exists():
            selected.append(ckpt_path)
        else:
            print(f"[warn] checkpoint for epoch {epoch} not found: {ckpt_path}", file=sys.stderr)

    if not selected:
        raise FileNotFoundError(f"No matching checkpoint files found for top epochs: {top_epochs}")

    return selected


def main() -> None:
    parser = argparse.ArgumentParser(description="Checkpoint soup: average weights from multiple checkpoints")
    parser.add_argument("--ckpt-dir", required=True, help="Directory containing epoch_NNN.pt checkpoints")
    parser.add_argument("--mode", choices=["last_k", "top_k_jsd"], default="last_k",
                        help="Selection mode: last_k (last K epochs) or top_k_jsd (best K by validation JSD)")
    parser.add_argument("--k", type=int, default=5, help="Number of checkpoints to average")
    parser.add_argument("--metrics-jsonl", default=None,
                        help="Path to epoch_metrics.jsonl (required for top_k_jsd mode)")
    parser.add_argument("--metric-key", default="peak.profile_target_jsd_debiased_mean",
                        help="Metric key for top_k_jsd selection (default: peak.profile_target_jsd_debiased_mean)")
    parser.add_argument("--output", required=True, help="Output path for the soup checkpoint")
    args = parser.parse_args()

    ckpt_dir = Path(args.ckpt_dir)
    output_path = Path(args.output)

    # Select checkpoints
    if args.mode == "last_k":
        selected = select_last_k(ckpt_dir, args.k)
    elif args.mode == "top_k_jsd":
        if not args.metrics_jsonl:
            parser.error("--metrics-jsonl is required for top_k_jsd mode")
        selected = select_top_k_jsd(ckpt_dir, Path(args.metrics_jsonl), args.k, args.metric_key)

    print(f"Averaging {len(selected)} checkpoints:")
    for p in selected:
        print(f"  {p.name}")

    # Load and average
    states = [load_model_state(p) for p in selected]
    avg_state = average_states(states)

    # Save: reuse the last checkpoint as template, replace model_state
    base_ckpt = torch.load(selected[-1], map_location="cpu", weights_only=False)
    if isinstance(base_ckpt, dict):
        base_ckpt["model_state"] = avg_state
        base_ckpt["soup_sources"] = [str(p) for p in selected]
        base_ckpt["soup_mode"] = args.mode
        base_ckpt["soup_k"] = args.k
    else:
        base_ckpt = avg_state

    output_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(base_ckpt, output_path)
    print(f"Saved soup checkpoint to {output_path}")


if __name__ == "__main__":
    main()
