#!/usr/bin/env python3
"""NT v2 feature extraction and probe analysis."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy import stats
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import r2_score, roc_auc_score
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler


@dataclass(frozen=True)
class RegionRecord:
    chrom: str
    center: int
    source: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    extract = subparsers.add_parser("extract", help="Extract NT v2 features on sampled tutorial-valid loci.")
    extract.add_argument("--model-dir", required=True)
    extract.add_argument("--genome-fasta", required=True)
    extract.add_argument("--bigwig", required=True)
    extract.add_argument("--peaks-bed", required=True)
    extract.add_argument("--nonpeaks-bed", required=True)
    extract.add_argument("--folds-json", required=True)
    extract.add_argument("--split", default="valid")
    extract.add_argument("--input-len", type=int, default=2114)
    extract.add_argument("--supervised-bp", type=int, default=1000)
    extract.add_argument("--sample-size-per-class", type=int, default=500)
    extract.add_argument("--layers", default="auto_quartiles")
    extract.add_argument("--feature-types", default="global_mean,bins4_mean")
    extract.add_argument("--bin-count", type=int, default=4)
    extract.add_argument("--batch-size", type=int, default=8)
    extract.add_argument("--device", default="cuda")
    extract.add_argument("--seed", type=int, default=42)
    extract.add_argument("--profile-probe-bins", type=int, default=16)
    extract.add_argument("--output-dir", required=True)

    analyze = subparsers.add_parser("analyze", help="Run probe metrics with NT features and baseline encoded features.")
    analyze.add_argument("--probe-dir", required=True)
    analyze.add_argument("--baseline-checkpoint", required=True)
    analyze.add_argument("--baseline-model-config", required=True)
    analyze.add_argument("--transchrombp-root", required=True)
    analyze.add_argument("--batch-size", type=int, default=64)
    analyze.add_argument("--device", default="cuda")
    analyze.add_argument("--output-dir", default="")
    return parser.parse_args()


def load_fold_chroms(folds_json: str, split: str) -> list[str]:
    with open(folds_json, "r", encoding="utf-8") as f:
        folds = json.load(f)
    chroms = folds.get(split) or folds.get("val" if split == "valid" else split)
    if not isinstance(chroms, list) or not chroms:
        raise ValueError(f"No chromosomes found for split={split!r} in {folds_json}")
    return [str(chrom) for chrom in chroms]


def _maybe_int(value: str) -> int | None:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def load_regions_from_bed(path: str, allowed_chroms: set[str], source: str) -> list[RegionRecord]:
    records: list[RegionRecord] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            if chrom not in allowed_chroms:
                continue
            start = int(fields[1])
            end = int(fields[2])
            summit = _maybe_int(fields[9]) if len(fields) > 9 else None
            center = start + summit if summit is not None else (start + end) // 2
            records.append(RegionRecord(chrom=chrom, center=center, source=source))
    if not records:
        raise ValueError(f"No usable records found in {path}")
    return records


def sample_records(records: list[RegionRecord], n: int, seed: int) -> list[RegionRecord]:
    if len(records) <= n:
        return list(records)
    rng = np.random.default_rng(seed)
    indices = np.sort(rng.choice(len(records), size=n, replace=False))
    return [records[int(i)] for i in indices]


def sum_into_bins(values: np.ndarray, n_bins: int) -> np.ndarray:
    chunks = np.array_split(np.asarray(values, dtype=np.float32), int(n_bins))
    return np.asarray([float(chunk.sum(dtype=np.float32)) for chunk in chunks], dtype=np.float32)


def fetch_sequences_and_targets(
    records: list[RegionRecord],
    genome_fasta: str,
    bigwig_path: str,
    input_len: int,
    supervised_bp: int,
    profile_probe_bins: int,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    from pyBigWig import open as open_bigwig
    from pyfaidx import Fasta

    fasta = Fasta(genome_fasta)
    bigwig = open_bigwig(bigwig_path)
    if bigwig is None:
        raise ValueError(f"Failed to open bigwig: {bigwig_path}")

    sequences: list[str] = []
    logcounts = np.zeros(len(records), dtype=np.float32)
    profile16 = np.zeros((len(records), profile_probe_bins), dtype=np.float32)
    try:
        for idx, record in enumerate(records):
            start = max(0, record.center - input_len // 2)
            end = start + input_len
            seq = str(fasta[record.chrom][start:end]).upper()
            if len(seq) < input_len:
                seq = seq + ("N" * (input_len - len(seq)))
            sequences.append(seq)

            prof_start = max(0, record.center - supervised_bp // 2)
            prof_end = prof_start + supervised_bp
            values = np.asarray(bigwig.values(record.chrom, prof_start, prof_end, numpy=True), dtype=np.float32)
            values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
            logcounts[idx] = float(np.log1p(values.sum()))
            profile16[idx] = sum_into_bins(values, profile_probe_bins)
    finally:
        bigwig.close()
    return sequences, logcounts, profile16


def resolve_layers(layer_spec: str, num_hidden_layers: int) -> list[int]:
    if layer_spec == "auto_quartiles":
        raw = [
            max(1, round(num_hidden_layers * 0.25)),
            max(1, round(num_hidden_layers * 0.5)),
            max(1, round(num_hidden_layers * 0.75)),
            num_hidden_layers,
        ]
        ordered: list[int] = []
        for value in raw:
            if value not in ordered:
                ordered.append(value)
        return ordered
    layers = [int(part.strip()) for part in layer_spec.split(",") if part.strip()]
    if not layers:
        raise ValueError(f"Invalid --layers spec: {layer_spec!r}")
    return layers


def feature_key(layer: int, feature_type: str) -> str:
    return f"layer_{layer:02d}__{feature_type}"


def masked_global_mean(hidden: torch.Tensor, mask: torch.Tensor) -> np.ndarray:
    masked = hidden * mask.unsqueeze(-1)
    denom = mask.sum(dim=1, keepdim=True).clamp_min(1)
    return (masked.sum(dim=1) / denom).detach().cpu().numpy().astype(np.float32)


def masked_bins_mean(hidden: torch.Tensor, mask: torch.Tensor, bin_count: int) -> np.ndarray:
    outputs: list[np.ndarray] = []
    hidden_np = hidden.detach().cpu().numpy()
    mask_np = mask.detach().cpu().numpy()

    for sample_hidden, sample_mask in zip(hidden_np, mask_np):
        valid_len = int(sample_mask.sum())
        valid_hidden = sample_hidden[:valid_len]
        if valid_len == 0:
            outputs.append(np.zeros((bin_count, sample_hidden.shape[-1]), dtype=np.float32))
            continue
        bins: list[np.ndarray] = []
        for index in range(bin_count):
            start = math.floor(index * valid_len / bin_count)
            end = math.floor((index + 1) * valid_len / bin_count)
            chunk = valid_hidden[start:end]
            if len(chunk) == 0:
                chunk = valid_hidden[max(0, min(start, valid_len - 1)) : max(1, min(start + 1, valid_len))]
            bins.append(chunk.mean(axis=0, dtype=np.float32))
        outputs.append(np.stack(bins, axis=0))
    return np.stack(outputs, axis=0).reshape(len(outputs), -1).astype(np.float32)


def compute_sha1(records: list[RegionRecord]) -> str:
    h = hashlib.sha1()
    for record in records:
        h.update(f"{record.chrom}:{record.center}:{record.source}\n".encode("utf-8"))
    return h.hexdigest()


def extract_command(args: argparse.Namespace) -> None:
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    allowed = set(load_fold_chroms(args.folds_json, args.split))
    peaks = sample_records(load_regions_from_bed(args.peaks_bed, allowed, "peak"), args.sample_size_per_class, args.seed)
    nonpeaks = sample_records(
        load_regions_from_bed(args.nonpeaks_bed, allowed, "nonpeak"),
        args.sample_size_per_class,
        args.seed + 1,
    )
    records = peaks + nonpeaks
    sequences, logcounts, profile16 = fetch_sequences_and_targets(
        records=records,
        genome_fasta=args.genome_fasta,
        bigwig_path=args.bigwig,
        input_len=args.input_len,
        supervised_bp=args.supervised_bp,
        profile_probe_bins=args.profile_probe_bins,
    )

    tokenizer = AutoTokenizer.from_pretrained(args.model_dir, trust_remote_code=True, local_files_only=True)
    model = AutoModelForMaskedLM.from_pretrained(args.model_dir, trust_remote_code=True, local_files_only=True)
    model = model.to(args.device).eval()

    layers = resolve_layers(args.layers, int(getattr(model.config, "num_hidden_layers", 0)))
    feature_types = [item.strip() for item in args.feature_types.split(",") if item.strip()]
    features: dict[str, list[np.ndarray]] = {feature_key(layer, feature): [] for layer in layers for feature in feature_types}
    token_lengths: list[int] = []
    timings: list[float] = []

    for start_idx in range(0, len(sequences), args.batch_size):
        batch_sequences = sequences[start_idx : start_idx + args.batch_size]
        batch = tokenizer.batch_encode_plus(batch_sequences, return_tensors="pt", padding=True, truncation=True)
        mask = batch["attention_mask"].bool()
        batch = {k: v.to(args.device) for k, v in batch.items()}
        mask = mask.to(args.device)

        started = time.perf_counter()
        with torch.no_grad():
            outputs = model(
                batch["input_ids"],
                attention_mask=mask,
                encoder_attention_mask=mask,
                output_hidden_states=True,
            )
        timings.append(time.perf_counter() - started)
        token_lengths.extend(mask.sum(dim=1).detach().cpu().tolist())

        for layer in layers:
            hidden = outputs.hidden_states[layer]
            for feature_type in feature_types:
                key = feature_key(layer, feature_type)
                if feature_type == "global_mean":
                    features[key].append(masked_global_mean(hidden, mask))
                elif feature_type == "bins4_mean":
                    features[key].append(masked_bins_mean(hidden, mask, args.bin_count))
                else:
                    raise ValueError(f"Unsupported feature type: {feature_type}")

    arrays = {key: np.concatenate(chunks, axis=0) for key, chunks in features.items()}
    np.savez_compressed(output_dir / "nt_features.npz", **arrays)
    for key, value in arrays.items():
        np.save(output_dir / f"{args.split}_{key}.f32.npy", value.astype(np.float32, copy=False))

    metadata_path = output_dir / "records.jsonl"
    with metadata_path.open("w", encoding="utf-8") as handle:
        for index, (record, sequence, logcount, profile_row) in enumerate(zip(records, sequences, logcounts, profile16)):
            handle.write(
                json.dumps(
                    {
                        "index": index,
                        "chrom": record.chrom,
                        "center": record.center,
                        "source": record.source,
                        "sequence": sequence,
                        "true_logcount": float(logcount),
                        "true_profile16": [float(x) for x in profile_row.tolist()],
                    },
                    ensure_ascii=True,
                )
                + "\n"
            )

    meta = {
        "model_dir": str(Path(args.model_dir).resolve()),
        "n_samples": len(records),
        "n_peak": len(peaks),
        "n_nonpeak": len(nonpeaks),
        "split": args.split,
        "input_len": args.input_len,
        "supervised_bp": args.supervised_bp,
        "layers": layers,
        "feature_types": feature_types,
        "feature_shapes": {key: list(value.shape) for key, value in arrays.items()},
        "record_sha1": compute_sha1(records),
        "profile_probe_bins": int(args.profile_probe_bins),
        "token_length_min": int(min(token_lengths)),
        "token_length_max": int(max(token_lengths)),
        "token_length_mean": float(np.mean(token_lengths)),
        "batch_size": args.batch_size,
        "device": args.device,
        "elapsed_sec_total": float(sum(timings)),
        "elapsed_sec_mean_batch": float(np.mean(timings)),
    }
    (output_dir / "extract_meta.json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "split": args.split,
        "n_records": len(records),
        "input_len": int(args.input_len),
        "record_sha1": compute_sha1(records),
        "features": sorted(arrays.keys()),
        "profile_probe_bins": int(args.profile_probe_bins),
    }
    (output_dir / f"manifest_{args.split}.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(meta, indent=2, sort_keys=True))


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
            profile16.append([float(x) for x in row.get("true_profile16", [])])
    return (
        sequences,
        np.asarray(labels, dtype=np.int64),
        np.asarray(logcounts, dtype=np.float32),
        np.asarray(profile16, dtype=np.float32),
    )


def ridge_oof_metrics(X: np.ndarray, y: np.ndarray, n_splits: int = 5, alpha: float = 1.0) -> dict[str, float]:
    n_splits = max(2, min(n_splits, len(X)))
    splitter = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    oof = np.zeros(len(X), dtype=np.float64)
    for train_idx, valid_idx in splitter.split(X):
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X[train_idx])
        X_valid = scaler.transform(X[valid_idx])
        model = Ridge(alpha=alpha)
        model.fit(X_train, y[train_idx])
        oof[valid_idx] = model.predict(X_valid)

    pearson_r, pearson_p = stats.pearsonr(y, oof)
    return {
        "r2": float(r2_score(y, oof)),
        "pearson_r": float(pearson_r),
        "pearson_p": float(pearson_p),
    }


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


def load_checkpoint(path: str) -> dict[str, Any]:
    import torch

    try:
        return torch.load(path, map_location="cpu", weights_only=False)
    except TypeError:
        return torch.load(path, map_location="cpu")


def extract_baseline_features(
    sequences: list[str],
    model_config_path: str,
    checkpoint_path: str,
    transchrombp_root: str,
    batch_size: int,
    device: str,
) -> tuple[np.ndarray, np.ndarray]:
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
    pred_chunks: list[np.ndarray] = []
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
            pooled = encoded.mean(dim=1)
            pred = outputs.logcount_full.squeeze(-1)

        pooled_chunks.append(pooled.detach().cpu().numpy().astype(np.float32))
        pred_chunks.append(pred.detach().cpu().numpy().astype(np.float32))
    return np.concatenate(pooled_chunks, axis=0), np.concatenate(pred_chunks, axis=0)


def analyze_command(args: argparse.Namespace) -> None:
    probe_dir = Path(args.probe_dir).resolve()
    output_dir = Path(args.output_dir).resolve() if args.output_dir else probe_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    feature_bundle = np.load(probe_dir / "nt_features.npz")
    sequences, labels, true_logcount, _ = load_metadata_records(probe_dir / "records.jsonl")
    encoded, pred_logcount = extract_baseline_features(
        sequences=sequences,
        model_config_path=args.baseline_model_config,
        checkpoint_path=args.baseline_checkpoint,
        transchrombp_root=args.transchrombp_root,
        batch_size=args.batch_size,
        device=args.device,
    )
    residual = true_logcount - pred_logcount

    baseline = {
        "encoded_only_auc": logistic_oof_auc(encoded, labels),
        "pred_logcount_r2": float(r2_score(true_logcount, pred_logcount)),
        "pred_logcount_pearson_r": float(stats.pearsonr(true_logcount, pred_logcount)[0]),
    }

    results: list[dict[str, Any]] = []
    for key in sorted(feature_bundle.files):
        nt_feat = feature_bundle[key].astype(np.float32)
        nt_only_auc = logistic_oof_auc(nt_feat, labels)
        concat_auc = logistic_oof_auc(np.concatenate([encoded, nt_feat], axis=1), labels)
        nt_logcount = ridge_oof_metrics(nt_feat, true_logcount)
        nt_residual = ridge_oof_metrics(nt_feat, residual)
        result = {
            "feature_key": key,
            "feature_dim": int(nt_feat.shape[1]),
            "nt_only_auc": nt_only_auc,
            "concat_auc": concat_auc,
            "concat_delta_vs_encoded": concat_auc - baseline["encoded_only_auc"],
            "nt_logcount_r2": nt_logcount["r2"],
            "nt_logcount_pearson_r": nt_logcount["pearson_r"],
            "nt_residual_r2": nt_residual["r2"],
            "nt_residual_pearson_r": nt_residual["pearson_r"],
        }
        results.append(result)

    summary = {
        "n_samples": len(sequences),
        "n_peak": int(labels.sum()),
        "n_nonpeak": int((1 - labels).sum()),
        "baseline": baseline,
        "results": results,
        "best_nt_only_auc": max(results, key=lambda row: row["nt_only_auc"]),
        "best_concat_auc": max(results, key=lambda row: row["concat_auc"]),
        "best_residual_r2": max(results, key=lambda row: row["nt_residual_r2"]),
    }

    (output_dir / "analysis_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with (output_dir / "analysis_results.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)
    print(json.dumps(summary, indent=2, sort_keys=True))


def main() -> None:
    args = parse_args()
    if args.command == "extract":
        extract_command(args)
    elif args.command == "analyze":
        analyze_command(args)
    else:  # pragma: no cover
        raise ValueError(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    main()
