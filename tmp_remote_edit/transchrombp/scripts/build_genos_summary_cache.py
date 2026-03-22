#!/usr/bin/env python3
"""Build Genos cached summary features (global_mean, optional bins4_mean).

Offline one-time pre-computation: runs frozen Genos-1.2B over all dataset
records and saves features as fp16 numpy memmap + manifest.json.

Usage:
    python build_genos_summary_cache.py \
        --data_config configs/data/data_tutorial_canonical_v1.yaml \
        --train_config configs/train/train_genos_cached_short10.yaml \
        --genos_model_path /path/to/Genos-1.2B \
        --output_dir /path/to/genos_cache \
        --splits train valid \
        --features global_mean bins4_mean \
        --batch_size 8 \
        --bins4_budget_gib 4.0
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import torch
import yaml


_HERE = Path(__file__).resolve()
_IMPORT_ROOTS = [
    _HERE.parents[1] / "src",  # deployed repo layout: <repo>/scripts/*.py + <repo>/src/transchrombp
    _HERE.parents[2],          # local snapshot layout: tmp_remote_edit/transchrombp/scripts/*.py
]
for _root in _IMPORT_ROOTS:
    if (_root / "transchrombp").exists():
        if str(_root) not in sys.path:
            sys.path.insert(0, str(_root))
        break


def load_yaml(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def compute_record_sha1(records) -> str:
    h = hashlib.sha1()
    for r in records:
        h.update(f"{r.chrom}:{r.center}:{r.source}\n".encode("utf-8"))
    return h.hexdigest()


def resolve_data_value(train_data_cfg: dict, data_source_cfg: dict, key: str, default=None):
    if key in train_data_cfg and train_data_cfg[key] not in ("", None):
        return train_data_cfg[key]
    if key in data_source_cfg and data_source_cfg[key] not in ("", None):
        return data_source_cfg[key]
    for section in ("input", "window", "sampling"):
        section_cfg = data_source_cfg.get(section, {})
        if key in section_cfg and section_cfg[key] not in ("", None):
            return section_cfg[key]
    return default


def pool_bins4(feat: torch.Tensor, n_bins: int = 4) -> torch.Tensor:
    """Pool [B, L, H] into [B, n_bins, H] by splitting L into n_bins."""
    B, L, H = feat.shape
    bin_size = L // n_bins
    # Trim to exact multiple
    feat = feat[:, :bin_size * n_bins, :]
    return feat.reshape(B, n_bins, bin_size, H).mean(dim=2)


def build_dataset_for_cache(split: str, train_config: dict, data_source_config: dict):
    """Build dataset with augmentation disabled for cache view."""
    from transchrombp.data import ChromBPNetBigWigDataset

    data_cfg = train_config.get("data", {})
    seed = int(train_config.get("seed", 1234))

    input_len = int(data_cfg.get("input_len", 2114))
    supervised_bp = int(data_cfg.get("supervised_bp", 1000))
    profile_bin_size = int(data_cfg.get("profile_bin_size", 1))
    nonpeak_ratio = float(resolve_data_value(data_cfg, data_source_config, "nonpeak_ratio", 1.0))
    region_source = str(data_cfg.get("region_source", "both"))

    if split == "train":
        region_source = str(data_cfg.get("train_region_source", region_source))
    else:
        region_source = str(data_cfg.get("val_region_source", region_source))

    genome_fasta = str(resolve_data_value(data_cfg, data_source_config, "genome_fasta", ""))
    folds_json = str(resolve_data_value(data_cfg, data_source_config, "folds_json", ""))
    peaks_bed = str(resolve_data_value(data_cfg, data_source_config, "peaks_bed", ""))
    nonpeaks_bed = str(resolve_data_value(data_cfg, data_source_config, "nonpeaks_bed", ""))
    bigwig_path = str(resolve_data_value(data_cfg, data_source_config, "bigwig", ""))

    ds = ChromBPNetBigWigDataset(
        genome_fasta=genome_fasta,
        bigwig_path=bigwig_path,
        peaks_bed=peaks_bed,
        nonpeaks_bed=nonpeaks_bed,
        folds_json=folds_json,
        split=split if split != "valid" else "valid",
        input_len=input_len,
        supervised_bp=supervised_bp,
        profile_bin_size=profile_bin_size,
        max_jitter=0,
        peak_max_jitter=0,
        nonpeak_max_jitter=0,
        seed=seed if split == "train" else seed + 10_000,
        nonpeak_ratio=nonpeak_ratio,
        max_records=0,
        region_source=region_source,
        random_revcomp=False,
        revcomp_prob=0.0,
    )
    return ds


def main():
    parser = argparse.ArgumentParser(description="Build Genos summary cache")
    parser.add_argument("--data_config", type=str, required=True,
                        help="Path to data source YAML config")
    parser.add_argument("--train_config", type=str, required=True,
                        help="Path to training YAML config")
    parser.add_argument("--genos_model_path", type=str, required=True,
                        help="Path to Genos-1.2B model directory")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Output directory for cache files")
    parser.add_argument("--splits", nargs="+", default=["train", "valid"],
                        help="Splits to process")
    parser.add_argument("--features", nargs="+", default=["global_mean"],
                        help="Features to extract: global_mean, bins4_mean")
    parser.add_argument("--batch_size", type=int, default=8)
    parser.add_argument("--layer", type=int, default=6)
    parser.add_argument("--bidirectional", action="store_true", default=True)
    parser.add_argument("--no_bidirectional", dest="bidirectional", action="store_false")
    parser.add_argument("--bins4_budget_gib", type=float, default=4.0,
                        help="Max GiB for bins4_mean; skip if exceeded")
    parser.add_argument("--device", type=str, default="cuda")
    parser.add_argument("--dry_run", action="store_true",
                        help="Only print budget, do not extract")
    args = parser.parse_args()

    train_config = load_yaml(args.train_config)
    data_source_config = load_yaml(args.data_config)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for split in args.splits:
        print(f"\n{'='*60}")
        print(f"Processing split: {split}")
        print(f"{'='*60}")

        ds = build_dataset_for_cache(split, train_config, data_source_config)
        n_records = len(ds)
        if getattr(ds, "supports_epoch_resampling", False):
            epoch_regions = int(getattr(ds, "peak_record_count", 0)) + min(
                int(getattr(ds, "nonpeak_record_count", 0)),
                int(getattr(ds, "target_nonpeak_count", 0)),
            )
        else:
            epoch_regions = n_records
        record_sha1 = compute_record_sha1(ds.records)

        print(f"  {split}_regions = {n_records}")
        if split == "train":
            print(f"  train_epoch_regions ≈ {epoch_regions}")
        print(f"  record_sha1 = {record_sha1}")

        # Budget estimation
        genos_dim = 1024
        est_global_bytes = n_records * genos_dim * 2  # float16
        est_bins4_bytes = n_records * 4 * genos_dim * 2
        est_global_gib = est_global_bytes / (1024**3)
        est_bins4_gib = est_bins4_bytes / (1024**3)

        print(f"  global_mean estimated: {est_global_gib:.3f} GiB")
        print(f"  bins4_mean estimated:  {est_bins4_gib:.3f} GiB")

        # Check disk space
        stat = os.statvfs(str(output_dir))
        free_gib = (stat.f_bavail * stat.f_frsize) / (1024**3)
        print(f"  output_dir free space: {free_gib:.1f} GiB")

        features = list(args.features)
        if "bins4_mean" in features and est_bins4_gib > args.bins4_budget_gib:
            print(f"  WARNING: bins4_mean ({est_bins4_gib:.3f} GiB) exceeds budget "
                  f"({args.bins4_budget_gib} GiB). Skipping bins4_mean.")
            features = [f for f in features if f != "bins4_mean"]

        if args.dry_run:
            print("  [dry_run] Skipping extraction.")
            continue

        # Build extractor
        from transchrombp.models.genos_adapter import GenosFeatureExtractor

        extractor = GenosFeatureExtractor(
            model_path=args.genos_model_path,
            layer=args.layer,
            device=args.device,
            bidirectional=args.bidirectional,
        )
        print(f"  Genos extractor loaded: layer={extractor.layer} "
              f"hidden_size={extractor.hidden_size} bidir={extractor.bidirectional}")

        # Open memmap files
        memmaps = {}
        if "global_mean" in features:
            path = output_dir / f"{split}_global_mean.f16.npy"
            arr = np.lib.format.open_memmap(
                str(path), mode="w+", dtype=np.float16, shape=(n_records, genos_dim),
            )
            memmaps["global_mean"] = arr

        if "bins4_mean" in features:
            path = output_dir / f"{split}_bins4_mean.f16.npy"
            arr = np.lib.format.open_memmap(
                str(path), mode="w+", dtype=np.float16, shape=(n_records, 4, genos_dim),
            )
            memmaps["bins4_mean"] = arr

        # Extract features
        t0 = time.time()
        batch_size = args.batch_size
        n_batches = (n_records + batch_size - 1) // batch_size

        for batch_idx in range(n_batches):
            start = batch_idx * batch_size
            end = min(start + batch_size, n_records)

            # Fetch canonical one-hot sequences (no jitter, no revcomp)
            seqs = []
            for i in range(start, end):
                record = ds.records[i]
                seq_start = int(record.center) - (ds.input_len // 2)
                seq = ds._fetch_onehot(record.chrom, seq_start, ds.input_len)
                seqs.append(seq)

            seq_batch = torch.from_numpy(np.stack(seqs)).to(args.device)
            feat = extractor.extract(seq_batch)  # [B, L, H]

            if "global_mean" in memmaps:
                global_mean = feat.mean(dim=1).cpu().numpy().astype(np.float16)
                memmaps["global_mean"][start:end] = global_mean

            if "bins4_mean" in memmaps:
                bins4 = pool_bins4(feat).cpu().numpy().astype(np.float16)
                memmaps["bins4_mean"][start:end] = bins4

            if (batch_idx + 1) % 50 == 0 or batch_idx == n_batches - 1:
                elapsed = time.time() - t0
                rate = (end) / elapsed
                eta = (n_records - end) / max(rate, 1e-6)
                print(f"  [{split}] {end}/{n_records} records "
                      f"({elapsed:.0f}s elapsed, ETA {eta:.0f}s)")

        elapsed = time.time() - t0
        print(f"  [{split}] Done in {elapsed:.1f}s")

        # Flush
        for arr in memmaps.values():
            del arr

        # Write manifest
        manifest = {
            "data_config_path": os.path.abspath(args.data_config),
            "train_config_path": os.path.abspath(args.train_config),
            "split": split,
            "input_len": ds.input_len,
            "n_records": n_records,
            "genos_model_path": args.genos_model_path,
            "layer": args.layer,
            "bidirectional": args.bidirectional,
            "features": features,
            "dtype": "float16",
            "record_sha1": record_sha1,
        }
        if split == "train":
            manifest["train_epoch_regions"] = epoch_regions
            manifest["supports_epoch_resampling"] = getattr(ds, "supports_epoch_resampling", False)

        manifest_path = output_dir / f"manifest_{split}.json"
        with open(manifest_path, "w", encoding="utf-8") as f:
            json.dump(manifest, f, indent=2, ensure_ascii=False)
        print(f"  Manifest written: {manifest_path}")

    print("\nAll done.")


if __name__ == "__main__":
    main()
