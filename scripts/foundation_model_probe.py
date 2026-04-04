#!/usr/bin/env python3
"""Unified foundation-model probe on cached features."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy import stats
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import r2_score, roc_auc_score
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records-jsonl", required=True)
    parser.add_argument("--feature-npz", default="")
    parser.add_argument("--cache-dir", default="")
    parser.add_argument("--split", default="valid")
    parser.add_argument("--baseline-checkpoint", required=True)
    parser.add_argument("--baseline-model-config", required=True)
    parser.add_argument("--transchrombp-root", required=True)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--profile-bins", type=int, default=16)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def add_import_root(root: Path) -> None:
    if (root / "src" / "transchrombp").exists():
        sys.path.insert(0, str(root / "src"))
        return
    if (root / "transchrombp").exists():
        sys.path.insert(0, str(root))
        return
    raise FileNotFoundError(f"Cannot find transchrombp package under {root}")


def load_metadata_records(path: Path) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    sequences: list[str] = []
    labels: list[int] = []
    logcounts: list[float] = []
    profile16: list[list[float]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            sequences.append(str(row["sequence"]))
            labels.append(1 if row["source"] == "peak" else 0)
            logcounts.append(float(row["true_logcount"]))
            bins = row.get("true_profile16", [])
            if not bins:
                raise ValueError(
                    f"{path} is missing true_profile16; regenerate the extract bundle with the upgraded extractor"
                )
            profile16.append([float(x) for x in bins])
    return (
        sequences,
        np.asarray(labels, dtype=np.int64),
        np.asarray(logcounts, dtype=np.float32),
        np.asarray(profile16, dtype=np.float32),
    )


def load_feature_bundle(feature_npz: str, cache_dir: str, split: str) -> dict[str, np.ndarray]:
    if feature_npz:
        bundle = np.load(feature_npz)
        return {key: bundle[key].astype(np.float32) for key in bundle.files}
    cache_root = Path(cache_dir).resolve()
    manifest_path = cache_root / f"manifest_{split}.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing cache manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    out: dict[str, np.ndarray] = {}
    for feature_name in manifest.get("features", []):
        feature_path = cache_root / f"{split}_{feature_name}.f32.npy"
        if not feature_path.is_file():
            feature_path = cache_root / f"{split}_{feature_name}.npy"
        if not feature_path.is_file():
            raise FileNotFoundError(f"Missing feature array for {feature_name!r} under {cache_root}")
        out[str(feature_name)] = np.load(feature_path).astype(np.float32)
    return out


def onehot_encode_batch(sequences: list[str]) -> np.ndarray:
    mapping = {
        "A": [1.0, 0.0, 0.0, 0.0],
        "C": [0.0, 1.0, 0.0, 0.0],
        "G": [0.0, 0.0, 1.0, 0.0],
        "T": [0.0, 0.0, 0.0, 1.0],
    }
    out = np.zeros((len(sequences), len(sequences[0]), 4), dtype=np.float32)
    for i, seq in enumerate(sequences):
        for j, base in enumerate(seq):
            out[i, j] = mapping.get(base.upper(), [0.0, 0.0, 0.0, 0.0])
    return out


def sum_into_bins(values: np.ndarray, n_bins: int) -> np.ndarray:
    chunks = np.array_split(np.asarray(values, dtype=np.float32), int(n_bins))
    return np.asarray([float(chunk.sum(dtype=np.float32)) for chunk in chunks], dtype=np.float32)


def load_checkpoint(path: str) -> dict[str, Any]:
    import torch

    try:
        return torch.load(path, map_location="cpu", weights_only=False)
    except TypeError:
        return torch.load(path, map_location="cpu")


def logistic_oof_auc(X: np.ndarray, y: np.ndarray, n_splits: int = 5) -> float:
    pos = int(y.sum())
    neg = int((1 - y).sum())
    n_splits = min(n_splits, pos, neg)
    if n_splits < 2:
        raise ValueError("Need at least 2 samples per class for stratified OOF AUC")
    splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    oof_prob = np.zeros(len(X), dtype=np.float64)
    for train_idx, valid_idx in splitter.split(X, y):
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X[train_idx])
        X_valid = scaler.transform(X[valid_idx])
        model = LogisticRegression(max_iter=2000, C=1.0)
        model.fit(X_train, y[train_idx])
        oof_prob[valid_idx] = model.predict_proba(X_valid)[:, 1]
    return float(roc_auc_score(y, oof_prob))


def ridge_oof_metrics(X: np.ndarray, y: np.ndarray, n_splits: int = 5, alpha: float = 1.0) -> dict[str, float]:
    n_splits = max(2, min(n_splits, len(X)))
    splitter = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    oof = np.zeros_like(y, dtype=np.float64)
    for train_idx, valid_idx in splitter.split(X):
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X[train_idx])
        X_valid = scaler.transform(X[valid_idx])
        model = Ridge(alpha=alpha)
        model.fit(X_train, y[train_idx])
        oof[valid_idx] = model.predict(X_valid)
    metrics = {"r2": float(r2_score(y, oof))}
    if y.ndim == 1:
        pearson_r, pearson_p = stats.pearsonr(y, oof)
        metrics["pearson_r"] = float(pearson_r)
        metrics["pearson_p"] = float(pearson_p)
    return metrics


def extract_baseline_features(
    sequences: list[str],
    model_config_path: str,
    checkpoint_path: str,
    transchrombp_root: str,
    batch_size: int,
    device: str,
    profile_bins: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    import torch
    import yaml

    root = Path(transchrombp_root).resolve()
    add_import_root(root)
    from transchrombp.models import build_transchrombp_from_config
    from transchrombp.models.transchrombp import _to_channel_first

    with open(model_config_path, "r", encoding="utf-8") as handle:
        model_cfg = yaml.safe_load(handle) or {}

    model = build_transchrombp_from_config(model_cfg)
    ckpt = load_checkpoint(checkpoint_path)
    state = ckpt.get("model_state", ckpt)
    state = {k.replace("module.", ""): v for k, v in state.items()}
    model.load_state_dict(state, strict=False)
    model = model.to(device).eval()

    pooled_chunks: list[np.ndarray] = []
    pred_logcount_chunks: list[np.ndarray] = []
    pred_profile16_chunks: list[np.ndarray] = []
    for start_idx in range(0, len(sequences), batch_size):
        batch_sequences = sequences[start_idx : start_idx + batch_size]
        seq_np = onehot_encode_batch(batch_sequences)
        seq = torch.from_numpy(seq_np).to(device)
        with torch.no_grad():
            outputs = model(seq)
            x = model.conv_stem(_to_channel_first(seq))
            local_feat = model.local_tower(x)
            tokens = local_feat.transpose(1, 2)
            encoded = model.transformer(tokens)
            pooled = model._pool_for_count(encoded)
            pred_logcount = outputs.logcount_full.squeeze(-1)
            pred_profile = torch.softmax(outputs.profile_logits_full, dim=-1).detach().cpu().numpy()

        pooled_chunks.append(pooled.detach().cpu().numpy().astype(np.float32))
        pred_logcount_chunks.append(pred_logcount.detach().cpu().numpy().astype(np.float32))
        pred_profile16_chunks.append(
            np.stack([sum_into_bins(row, profile_bins) for row in pred_profile], axis=0).astype(np.float32)
        )
    return (
        np.concatenate(pooled_chunks, axis=0),
        np.concatenate(pred_logcount_chunks, axis=0),
        np.concatenate(pred_profile16_chunks, axis=0),
    )


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    sequences, labels, true_logcount, true_profile16 = load_metadata_records(Path(args.records_jsonl).resolve())
    feature_bundle = load_feature_bundle(args.feature_npz, args.cache_dir, args.split)
    encoded, pred_logcount, pred_profile16 = extract_baseline_features(
        sequences=sequences,
        model_config_path=args.baseline_model_config,
        checkpoint_path=args.baseline_checkpoint,
        transchrombp_root=args.transchrombp_root,
        batch_size=args.batch_size,
        device=args.device,
        profile_bins=args.profile_bins,
    )

    count_residual = true_logcount - pred_logcount
    profile_residual = true_profile16 - pred_profile16
    baseline = {
        "encoded_only_auc": logistic_oof_auc(encoded, labels),
        "pred_logcount_r2": float(r2_score(true_logcount, pred_logcount)),
        "pred_logcount_pearson_r": float(stats.pearsonr(true_logcount, pred_logcount)[0]),
        "pred_profile16_r2": float(r2_score(true_profile16, pred_profile16)),
    }

    rows: list[dict[str, Any]] = []
    for feature_key in sorted(feature_bundle):
        feat = feature_bundle[feature_key].astype(np.float32)
        flat_feat = feat.reshape(feat.shape[0], -1)
        fm_only_auc = logistic_oof_auc(flat_feat, labels)
        concat_auc = logistic_oof_auc(np.concatenate([encoded, flat_feat], axis=1), labels)
        count_direct = ridge_oof_metrics(flat_feat, true_logcount)
        count_resid = ridge_oof_metrics(flat_feat, count_residual)
        profile_direct = ridge_oof_metrics(flat_feat, true_profile16)
        profile_resid = ridge_oof_metrics(flat_feat, profile_residual)
        rows.append(
            {
                "feature_key": feature_key,
                "feature_dim": int(flat_feat.shape[1]),
                "fm_only_auc": fm_only_auc,
                "concat_auc": concat_auc,
                "concat_delta_vs_encoded": concat_auc - baseline["encoded_only_auc"],
                "count_direct_r2": count_direct["r2"],
                "count_direct_pearson_r": count_direct.get("pearson_r", float("nan")),
                "count_residual_r2": count_resid["r2"],
                "count_residual_pearson_r": count_resid.get("pearson_r", float("nan")),
                "profile16_direct_r2": profile_direct["r2"],
                "profile16_residual_r2": profile_resid["r2"],
            }
        )

    summary = {
        "n_samples": len(sequences),
        "n_peak": int(labels.sum()),
        "n_nonpeak": int((1 - labels).sum()),
        "profile_bins": int(args.profile_bins),
        "baseline": baseline,
        "best_fm_only_auc": max(rows, key=lambda row: row["fm_only_auc"]),
        "best_concat_auc": max(rows, key=lambda row: row["concat_auc"]),
        "best_count_residual_r2": max(rows, key=lambda row: row["count_residual_r2"]),
        "best_profile16_residual_r2": max(rows, key=lambda row: row["profile16_residual_r2"]),
        "results": rows,
    }

    promotion = {
        "concat_gate_pass": bool(summary["best_concat_auc"]["concat_delta_vs_encoded"] >= -0.003),
        "count_residual_gate_pass": bool(summary["best_count_residual_r2"]["count_residual_r2"] > 0.0),
        "profile16_residual_gate_pass": bool(summary["best_profile16_residual_r2"]["profile16_residual_r2"] > 0.0),
    }
    promotion["promotion_recommended"] = bool(
        promotion["concat_gate_pass"]
        and (promotion["count_residual_gate_pass"] or promotion["profile16_residual_gate_pass"])
    )

    (output_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (output_dir / "promotion_decision.json").write_text(
        json.dumps(promotion, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with (output_dir / "analysis_results.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"summary": summary, "promotion": promotion}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
