#!/usr/bin/env python3
"""Evaluate a single checkpoint on validation/test data.

Reuses the validation pipeline from train_ddp.py but runs on a single GPU
without DDP. Supports evaluating soup checkpoints or any epoch checkpoint.

Usage:
  python scripts/evaluate_checkpoint.py \
    --model-config configs/model/transchrombp_teacher_v2.yaml \
    --train-config configs/train/train_ablation_v2_main.yaml \
    --checkpoint outputs/checkpoints/run_name/best.pt \
    --output-json outputs/eval/run_name_best.json

  # Evaluate a soup checkpoint
  python scripts/evaluate_checkpoint.py \
    --model-config configs/model/transchrombp_teacher_v2.yaml \
    --train-config configs/train/train_ablation_v2_main.yaml \
    --checkpoint outputs/checkpoints/run_name/soup_last5.pt \
    --output-json outputs/eval/run_name_soup_last5.json
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict

import torch
import yaml
from torch.utils.data import DataLoader

# ---- 插入项目 src 到 path（如果直接运行脚本） ----
_SCRIPT_DIR = Path(__file__).resolve().parent
_SRC_DIR = _SCRIPT_DIR.parent / "src"
if _SRC_DIR.is_dir() and str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

from transchrombp.data import ChromBPNetBigWigDataset
from transchrombp.models import build_transchrombp_from_config


@dataclass
class FakeDistEnv:
    """Minimal DistEnv-compatible object for single-GPU eval."""
    rank: int = 0
    local_rank: int = 0
    world_size: int = 1
    distributed: bool = False
    device: torch.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def load_yaml(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def build_model_for_checkpoint(
    cli_model_cfg: Dict[str, Any],
    ckpt: Any,
) -> tuple[torch.nn.Module, str]:
    """Build a model and load checkpoint weights, falling back to checkpoint config if needed."""
    ckpt_model_cfg = ckpt.get("model_config") if isinstance(ckpt, dict) else None
    state_dict = ckpt["model_state"] if isinstance(ckpt, dict) and "model_state" in ckpt else ckpt

    def load_with_compat(model: torch.nn.Module) -> bool:
        try:
            model.load_state_dict(state_dict)
            return True
        except RuntimeError as exc:
            missing_allowed = {
                "count_pool_proj.weight",
                "count_pool_proj.bias",
            }
            model_keys = set(model.state_dict().keys())
            state_keys = set(state_dict.keys())
            missing = sorted(model_keys - state_keys)
            unexpected = sorted(state_keys - model_keys)
            if any(key not in missing_allowed for key in missing):
                raise exc
            if any(not key.startswith("count_pool_proj.") for key in unexpected):
                raise exc
            model.load_state_dict(state_dict, strict=False)
            print(
                "[warn] Loaded checkpoint with count_pool_proj backward-compat mode; "
                "old attention-pooling checkpoints were trained before the pooling layer was implemented."
            )
            return True

    model = build_transchrombp_from_config(cli_model_cfg)
    try:
        load_with_compat(model)
        return model, "cli"
    except RuntimeError as exc:
        if not isinstance(ckpt_model_cfg, dict):
            raise
        print(
            "[warn] CLI model config failed to load checkpoint strictly; "
            "retrying with model_config embedded in checkpoint."
        )
        print(f"[warn] Original load_state_dict error: {exc}")
        model = build_transchrombp_from_config(ckpt_model_cfg)
        load_with_compat(model)
        return model, "checkpoint"


def build_eval_loader(
    train_cfg: Dict[str, Any],
    data_source_cfg: Dict[str, Any],
    split: str,
) -> DataLoader:
    """Build eval loader with the same real-data API used by train_ddp."""
    from transchrombp.training.train_ddp import SyntheticChromDataset, _resolve_real_data_value

    data_cfg = train_cfg.get("data", {})
    source = str(data_cfg.get("source", "synthetic")).lower()
    input_len = int(data_cfg.get("input_len", 2114))
    output_len = int(data_cfg.get("output_len", 1000))
    supervised_bp = int(data_cfg.get("supervised_bp", output_len))
    profile_bin_size = int(data_cfg.get("profile_bin_size", 1))
    batch_size = int(data_cfg.get("batch_size_per_gpu", 16))
    num_workers = int(data_cfg.get("num_workers", 2))
    pin_memory = bool(data_cfg.get("pin_memory", True))
    prefetch_factor = int(data_cfg.get("prefetch_factor", 2))
    persistent_workers = bool(data_cfg.get("persistent_workers", True))
    seed = int(train_cfg.get("seed", 1234))

    if source == "synthetic":
        sample_key = "val_samples" if split == "valid" else "test_samples"
        n_samples = int(data_cfg.get(sample_key, data_cfg.get("val_samples", 1024)))
        dataset = SyntheticChromDataset(
            n_samples=n_samples,
            input_len=input_len,
            output_len=output_len,
            seed=seed + (10_000 if split == "valid" else 20_000),
            count_scale=float(data_cfg.get("synthetic_count_scale", 20.0)),
        )
    elif source in {"chrombpnet_bigwig", "bigwig"}:
        nonpeak_ratio = float(_resolve_real_data_value(data_cfg, data_source_cfg, "nonpeak_ratio", 1.0))
        revcomp_prob = float(_resolve_real_data_value(data_cfg, data_source_cfg, "revcomp_prob", 0.5))
        track_total_count_target = float(
            _resolve_real_data_value(data_cfg, data_source_cfg, "track_total_count_target", 0.0)
        )
        default_region_source = str(data_cfg.get("region_source", "both"))
        if split == "valid":
            eval_region_source = str(data_cfg.get("val_region_source", default_region_source))
            max_records = int(data_cfg.get("max_val_regions", 0))
            split_seed = seed + 10_000
        else:
            eval_region_source = str(
                data_cfg.get("test_region_source", data_cfg.get("val_region_source", default_region_source))
            )
            max_records = int(data_cfg.get("max_test_regions", 0))
            split_seed = seed + 20_000

        genome_fasta = str(_resolve_real_data_value(data_cfg, data_source_cfg, "genome_fasta", ""))
        folds_json = str(_resolve_real_data_value(data_cfg, data_source_cfg, "folds_json", ""))
        peaks_bed = str(_resolve_real_data_value(data_cfg, data_source_cfg, "peaks_bed", ""))
        nonpeaks_bed = str(_resolve_real_data_value(data_cfg, data_source_cfg, "nonpeaks_bed", ""))
        bigwig_path = str(_resolve_real_data_value(data_cfg, data_source_cfg, "bigwig", ""))

        required = {
            "genome_fasta": genome_fasta,
            "folds_json": folds_json,
            "peaks_bed": peaks_bed,
            "bigwig": bigwig_path,
        }
        missing = [key for key, value in required.items() if not value]
        if missing:
            raise ValueError(
                "Real-data backend requires the following paths: " + ", ".join(missing)
            )

        dataset = ChromBPNetBigWigDataset(
            genome_fasta=genome_fasta,
            bigwig_path=bigwig_path,
            peaks_bed=peaks_bed,
            nonpeaks_bed=nonpeaks_bed,
            folds_json=folds_json,
            split=split,
            input_len=input_len,
            supervised_bp=supervised_bp,
            profile_bin_size=profile_bin_size,
            max_jitter=0,
            peak_max_jitter=0,
            nonpeak_max_jitter=0,
            seed=split_seed,
            nonpeak_ratio=nonpeak_ratio,
            max_records=max_records,
            region_source=eval_region_source,
            random_revcomp=False,
            revcomp_prob=revcomp_prob,
            track_total_count_target=track_total_count_target,
        )
    else:
        raise ValueError(f"Unsupported data.source: {source}")

    loader_kwargs: Dict[str, Any] = dict(
        batch_size=batch_size,
        shuffle=False,
        num_workers=num_workers,
        pin_memory=pin_memory,
        drop_last=False,
    )
    if num_workers > 0:
        loader_kwargs["prefetch_factor"] = prefetch_factor
        loader_kwargs["persistent_workers"] = persistent_workers
    return DataLoader(dataset, **loader_kwargs)


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate a checkpoint on validation data")
    parser.add_argument("--model-config", required=True, help="Model config YAML")
    parser.add_argument("--train-config", required=True, help="Train config YAML (for data and loss settings)")
    parser.add_argument("--checkpoint", required=True, help="Checkpoint .pt file to evaluate")
    parser.add_argument("--output-json", default=None, help="Output JSON file for metrics (default: print to stdout)")
    parser.add_argument("--device", default=None, help="Device (default: auto-detect)")
    parser.add_argument("--data-config", default="", help="Optional explicit data config YAML")
    parser.add_argument("--split", choices=["valid", "test"], default="valid", help="Eval split")
    parser.add_argument("--no-amp", action="store_true", help="Disable AMP")
    args = parser.parse_args()

    # Device
    if args.device:
        device = torch.device(args.device)
    else:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    from transchrombp.training.train_ddp import (
        apply_external_data_config_defaults,
        apply_training_mode_defaults,
        resolve_amp_dtype,
        resolve_data_config_path,
        run_validation,
        validate_semantics_profile,
    )

    # Load configs
    model_cfg = load_yaml(args.model_config)
    train_cfg = load_yaml(args.train_config)
    loss_cfg = train_cfg.get("loss", {})
    data_config_path = resolve_data_config_path(train_cfg, args.data_config)
    data_source_cfg = load_yaml(data_config_path) if data_config_path else {}
    if data_source_cfg:
        apply_external_data_config_defaults(train_cfg, data_source_cfg)
    apply_training_mode_defaults(train_cfg)
    validate_semantics_profile(train_cfg, data_source_cfg, rank=0)

    # Load checkpoint
    print(f"Loading checkpoint: {args.checkpoint}")
    ckpt = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    print("Building model...")
    model, model_config_source = build_model_for_checkpoint(model_cfg, ckpt)
    if isinstance(ckpt, dict) and "model_state" in ckpt:
        ckpt_info = {
            "epoch": ckpt.get("epoch"),
            "global_step": ckpt.get("global_step"),
            "soup_mode": ckpt.get("soup_mode"),
            "soup_k": ckpt.get("soup_k"),
        }
    else:
        ckpt_info = {}

    model = model.to(device)
    model.eval()

    # Build data loader
    print(f"Loading {args.split} data...")
    eval_loader = build_eval_loader(train_cfg, data_source_cfg, split=args.split)
    print(f"{args.split} samples: {len(eval_loader.dataset)}")

    # Run evaluation
    print("Running evaluation...")
    trainer_cfg = train_cfg.get("trainer", {})
    amp_dtype = resolve_amp_dtype(str(trainer_cfg.get("precision", "bf16")))
    use_amp = (device.type == "cuda") and (not args.no_amp) and (amp_dtype != torch.float32)
    metrics = run_validation(
        model=model,
        val_loader=eval_loader,
        dist_env=FakeDistEnv(device=device),
        loss_cfg=loss_cfg,
        use_amp=use_amp,
        amp_dtype=amp_dtype,
    )

    # Output
    result = {
        "checkpoint": str(args.checkpoint),
        "model_config": str(args.model_config),
        "model_config_source": model_config_source,
        "train_config": str(args.train_config),
        "data_config": str(data_config_path),
        "split": args.split,
        **ckpt_info,
        "metrics": metrics,
    }

    if args.output_json:
        output_path = Path(args.output_json)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(result, f, indent=2)
        print(f"Metrics saved to {output_path}")
    else:
        print(json.dumps(result, indent=2))

    # Print summary
    print("\n=== Summary ===")
    for region in ["overall", "peak", "nonpeak"]:
        if region in metrics:
            m = metrics[region]
            jsd_full = m.get("profile_target_jsd_full_mean", float("nan"))
            jsd_debiased = m.get("profile_target_jsd_debiased_mean", float("nan"))
            count_r = m.get("count_pearson_full", float("nan"))
            count_mae = m.get("count_mae_full", float("nan"))
            fd_gap = m.get("profile_full_debiased_jsd", float("nan"))
            print(f"  {region:>8s}: JSD_full={jsd_full:.5f}  JSD_debiased={jsd_debiased:.5f}  "
                  f"gap={fd_gap:.5f}  count_r={count_r:.4f}  count_mae={count_mae:.2f}")


if __name__ == "__main__":
    main()
