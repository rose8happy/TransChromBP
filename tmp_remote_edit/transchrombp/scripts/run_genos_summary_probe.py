#!/usr/bin/env python3
"""Genos summary probes: count residual ridge + complement logistic.

Probe A (count residual ridge):
  - Uses G0 best.pt to predict logcount on validation set
  - Fits ridge regression: genos_global → (true_logcount - pred_logcount)
  - Reports R² and Pearson r

Probe B (complement logistic):
  - Uses G0 pooled encoded and genos_global
  - Fits logistic regression on peak/nonpeak labels
  - Reports AUC(encoded_only), AUC(genos_only), AUC(concat)

Usage:
    python run_genos_summary_probe.py \
        --g0_checkpoint /path/to/G0/best.pt \
        --model_config /path/to/v2fix_baseline.yaml \
        --train_config /path/to/train_genos_cached_short10.yaml \
        --data_config /path/to/data_tutorial_canonical_v1.yaml \
        --genos_cache_dir /path/to/genos_cache \
        --output_dir /path/to/probe_results
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import torch
import yaml
from scipy import stats
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import r2_score, roc_auc_score
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler


_PKG_PARENT = Path(__file__).resolve().parents[2]
if str(_PKG_PARENT) not in sys.path:
    sys.path.insert(0, str(_PKG_PARENT))


def load_yaml(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def load_checkpoint(path: str) -> dict:
    try:
        return torch.load(path, map_location="cpu", weights_only=False)
    except TypeError:
        return torch.load(path, map_location="cpu")


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


def ridge_oof_metrics(X: np.ndarray, y: np.ndarray, n_splits: int = 5, alpha: float = 1.0) -> tuple[float, float, float]:
    n_samples = int(X.shape[0])
    n_splits = max(2, min(n_splits, n_samples))
    splitter = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    oof_pred = np.zeros(n_samples, dtype=np.float64)

    for train_idx, val_idx in splitter.split(X):
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X[train_idx])
        X_val = scaler.transform(X[val_idx])
        model = Ridge(alpha=alpha)
        model.fit(X_train, y[train_idx])
        oof_pred[val_idx] = model.predict(X_val)

    r2 = r2_score(y, oof_pred)
    pearson_r, _ = stats.pearsonr(y, oof_pred)
    return float(r2), float(pearson_r), oof_pred


def logistic_oof_auc(X: np.ndarray, y: np.ndarray, n_splits: int = 5) -> float:
    pos = int(y.sum())
    neg = int((1 - y).sum())
    n_splits = min(n_splits, pos, neg)
    if n_splits < 2:
        raise ValueError("Need at least 2 examples per class for stratified OOF AUC")

    splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    oof_prob = np.zeros(X.shape[0], dtype=np.float64)
    for train_idx, val_idx in splitter.split(X, y):
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X[train_idx])
        X_val = scaler.transform(X[val_idx])
        model = LogisticRegression(max_iter=2000, C=1.0)
        model.fit(X_train, y[train_idx])
        oof_prob[val_idx] = model.predict_proba(X_val)[:, 1]
    return float(roc_auc_score(y, oof_prob))


def main():
    parser = argparse.ArgumentParser(description="Genos summary probes")
    parser.add_argument("--g0_checkpoint", type=str, required=True)
    parser.add_argument("--model_config", type=str, required=True)
    parser.add_argument("--train_config", type=str, required=True)
    parser.add_argument("--data_config", type=str, required=True)
    parser.add_argument("--genos_cache_dir", type=str, required=True)
    parser.add_argument("--output_dir", type=str, required=True)
    parser.add_argument("--device", type=str, default="cuda")
    parser.add_argument("--batch_size", type=int, default=64)
    args = parser.parse_args()

    from transchrombp.data import ChromBPNetBigWigDataset
    from transchrombp.models import build_transchrombp_from_config

    model_cfg = load_yaml(args.model_config)
    train_cfg = load_yaml(args.train_config)
    data_source_cfg = load_yaml(args.data_config)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build model and load G0 checkpoint
    model = build_transchrombp_from_config(model_cfg)
    ckpt = load_checkpoint(args.g0_checkpoint)
    state = ckpt.get("model_state", ckpt)
    # Strip DDP module. prefix if present
    state = {k.replace("module.", ""): v for k, v in state.items()}
    model.load_state_dict(state, strict=False)
    model = model.to(args.device).eval()
    print(f"[probe] Loaded G0 from {args.g0_checkpoint}")

    # Build validation dataset with genos cache
    data_cfg = train_cfg.get("data", {})

    val_ds = ChromBPNetBigWigDataset(
        genome_fasta=str(resolve_data_value(data_cfg, data_source_cfg, "genome_fasta", "")),
        bigwig_path=str(resolve_data_value(data_cfg, data_source_cfg, "bigwig", "")),
        peaks_bed=str(resolve_data_value(data_cfg, data_source_cfg, "peaks_bed", "")),
        nonpeaks_bed=str(resolve_data_value(data_cfg, data_source_cfg, "nonpeaks_bed", "")),
        folds_json=str(resolve_data_value(data_cfg, data_source_cfg, "folds_json", "")),
        split="valid",
        input_len=int(data_cfg.get("input_len", 2114)),
        supervised_bp=int(data_cfg.get("supervised_bp", 1000)),
        profile_bin_size=int(data_cfg.get("profile_bin_size", 1)),
        max_jitter=0,
        seed=int(train_cfg.get("seed", 1234)) + 10_000,
        nonpeak_ratio=float(resolve_data_value(data_cfg, data_source_cfg, "nonpeak_ratio", 1.0)),
        region_source=str(data_cfg.get("val_region_source", data_cfg.get("region_source", "both"))),
        genos_cache_dir=args.genos_cache_dir,
        genos_cache_features=["global_mean"],
    )

    val_loader = torch.utils.data.DataLoader(
        val_ds, batch_size=args.batch_size, shuffle=False,
        num_workers=2, pin_memory=True,
    )
    print(f"[probe] Validation set: {len(val_ds)} records")

    # Forward pass to collect features
    all_pred_logcount = []
    all_true_logcount = []
    all_pooled = []
    all_genos_global = []
    all_sources = []

    with torch.no_grad():
        for batch in val_loader:
            seq = batch["seq"].to(args.device, non_blocking=True)
            profile_counts = batch["profile_counts"]

            outputs = model(seq)

            pred_logcount = outputs.logcount_full.squeeze(-1).cpu().numpy()
            true_total = profile_counts.sum(dim=-1).float()
            true_logcount = torch.log1p(true_total).numpy()

            # Get pooled encoded
            x = model.conv_stem(_to_channel_first(seq))
            local_feat = model.local_tower(x)
            tokens = local_feat.transpose(1, 2)
            if hasattr(model, 'genos_adapter') and model.genos_adapter is not None:
                pass  # baseline has no genos_adapter
            encoded = model.transformer(tokens)
            pooled = encoded.mean(dim=1).cpu().numpy()

            genos_global = batch["genos_global_mean"].numpy()

            all_pred_logcount.append(pred_logcount)
            all_true_logcount.append(true_logcount)
            all_pooled.append(pooled)
            all_genos_global.append(genos_global)
            all_sources.extend(batch.get("region_source", ["unknown"] * seq.shape[0]))

    pred_logcount = np.concatenate(all_pred_logcount)
    true_logcount = np.concatenate(all_true_logcount)
    pooled = np.concatenate(all_pooled)
    genos_global = np.concatenate(all_genos_global)
    sources = np.array(all_sources)

    # ── Probe A: count residual ridge ──
    print("\n" + "="*60)
    print("Probe A: Count Residual Ridge Regression")
    print("="*60)

    residual = true_logcount - pred_logcount

    r2, pearson_r, pred_residual = ridge_oof_metrics(genos_global, residual, n_splits=5, alpha=1.0)
    pearson_p = stats.pearsonr(residual, pred_residual).pvalue

    probe_a_results = {
        "R2": float(r2),
        "pearson_r": float(pearson_r),
        "pearson_p": float(pearson_p),
        "evaluation": "5-fold OOF",
        "residual_mean": float(residual.mean()),
        "residual_std": float(residual.std()),
        "n_samples": int(len(residual)),
    }
    print(f"  R² = {r2:.4f}")
    print(f"  Pearson r = {pearson_r:.4f} (p={pearson_p:.2e})")
    print(f"  Residual: mean={residual.mean():.4f}, std={residual.std():.4f}")

    # ── Probe B: complement logistic ──
    print("\n" + "="*60)
    print("Probe B: Complement Logistic Regression")
    print("="*60)

    labels = (sources == "peak").astype(int)
    if len(np.unique(labels)) < 2:
        print("  SKIP: only one class present in validation set")
        probe_b_results = {"skipped": True, "reason": "single class"}
    else:
        auc_a = logistic_oof_auc(pooled, labels, n_splits=5)
        auc_b = logistic_oof_auc(genos_global, labels, n_splits=5)
        auc_c = logistic_oof_auc(np.concatenate([pooled, genos_global], axis=1), labels, n_splits=5)

        probe_b_results = {
            "AUC_encoded_only": float(auc_a),
            "AUC_genos_only": float(auc_b),
            "AUC_concat": float(auc_c),
            "evaluation": "5-fold OOF",
            "n_peak": int(labels.sum()),
            "n_nonpeak": int((1 - labels).sum()),
        }
        print(f"  AUC(encoded_only) = {auc_a:.4f}")
        print(f"  AUC(genos_only)   = {auc_b:.4f}")
        print(f"  AUC(concat)       = {auc_c:.4f}")
        print(f"  peak={labels.sum()}, nonpeak={(1-labels).sum()}")

    # Save results
    results = {
        "probe_a_count_residual_ridge": probe_a_results,
        "probe_b_complement_logistic": probe_b_results,
        "g0_checkpoint": args.g0_checkpoint,
        "genos_cache_dir": args.genos_cache_dir,
    }
    results_path = output_dir / "probe_results.json"
    with open(results_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\nResults saved to {results_path}")


def _to_channel_first(seq_onehot):
    if seq_onehot.dim() == 3 and seq_onehot.size(2) == 4:
        return seq_onehot.transpose(1, 2)
    return seq_onehot


if __name__ == "__main__":
    main()
