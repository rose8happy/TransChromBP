"""Export record-aligned teacher targets from cached foundation features."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler


def _prepare_feature_matrix(features: np.ndarray, *, feature_name: str) -> np.ndarray:
    """Convert cached foundation features into a 2D matrix for sklearn."""

    matrix = np.asarray(features, dtype=np.float32)
    if matrix.ndim == 1:
        return matrix.reshape(-1, 1)
    if matrix.ndim == 2:
        return matrix
    if matrix.ndim == 3:
        # Cached bins*_mean features are [records, bins, hidden]. We reduce the
        # hidden axis so the resulting matrix stays tractable while preserving the
        # coarse positional bins that align with profile16 targets.
        return matrix.mean(axis=-1, dtype=np.float32)
    raise ValueError(
        f"Unsupported cached feature rank for {feature_name!r}: expected 1D/2D/3D, got shape={matrix.shape!r}"
    )


def _load_manifest(cache_dir: Path, split: str) -> dict:
    manifest_path = cache_dir / f"manifest_{split}.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing foundation cache manifest: {manifest_path}")
    return json.loads(manifest_path.read_text(encoding="utf-8"))


def _resolve_records_jsonl(cache_dir: Path, manifest: dict, split: str) -> Path:
    records_jsonl = str(manifest.get("records_jsonl", "")).strip()
    if not records_jsonl:
        raise ValueError(
            f"Foundation cache manifest for split={split!r} is missing records_jsonl; "
            "teacher export requires record-level metadata."
        )
    path = Path(records_jsonl)
    if not path.is_absolute():
        path = (cache_dir / path).resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Missing records_jsonl for split={split!r}: {path}")
    return path


def _resolve_feature_array(cache_dir: Path, split: str, feature_name: str) -> Path:
    candidates = [
        cache_dir / f"{split}_{feature_name}.f32.npy",
        cache_dir / f"{split}_{feature_name}.f16.npy",
        cache_dir / f"{split}_{feature_name}.npy",
    ]
    for path in candidates:
        if path.is_file():
            return path
    raise FileNotFoundError(
        f"Missing cached feature array for split={split!r}, feature={feature_name!r}; tried {candidates!r}"
    )


def _load_records_targets(records_jsonl: Path) -> tuple[np.ndarray, np.ndarray]:
    logcounts: list[float] = []
    profile16: list[list[float]] = []
    with records_jsonl.open("r", encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            logcounts.append(float(row["true_logcount"]))
            bins = row.get("true_profile16", [])
            if not bins:
                raise ValueError(f"{records_jsonl} is missing true_profile16 rows")
            profile16.append([float(x) for x in bins])
    return (
        np.asarray(logcounts, dtype=np.float32),
        np.asarray(profile16, dtype=np.float32),
    )


def _sanitize_profile_predictions(pred: np.ndarray) -> np.ndarray:
    pred = np.asarray(pred, dtype=np.float32)
    pred = np.clip(pred, a_min=0.0, a_max=None)
    row_sums = pred.sum(axis=1, keepdims=True)
    zero_rows = row_sums.squeeze(1) <= 0.0
    if np.any(zero_rows):
        pred[zero_rows] = 1.0 / max(pred.shape[1], 1)
    return pred


def export_teacher_cache(
    *,
    cache_dir: Path | str,
    output_dir: Path | str,
    feature_name: str,
    train_split: str = "train",
    predict_splits: Sequence[str] = ("train", "valid"),
    alpha: float = 1.0,
) -> None:
    cache_root = Path(cache_dir).resolve()
    output_root = Path(output_dir).resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    train_manifest = _load_manifest(cache_root, train_split)
    train_records_jsonl = _resolve_records_jsonl(cache_root, train_manifest, train_split)
    train_features = _prepare_feature_matrix(
        np.load(_resolve_feature_array(cache_root, train_split, feature_name), mmap_mode="r"),
        feature_name=feature_name,
    )
    train_logcount, train_profile16 = _load_records_targets(train_records_jsonl)
    if train_features.shape[0] != train_logcount.shape[0]:
        raise ValueError(
            f"Train feature rows ({train_features.shape[0]}) do not match teacher targets ({train_logcount.shape[0]})"
        )

    scaler = StandardScaler()
    train_features_scaled = scaler.fit_transform(train_features)

    count_model = Ridge(alpha=float(alpha))
    count_model.fit(train_features_scaled, train_logcount)

    profile_model = Ridge(alpha=float(alpha))
    profile_model.fit(train_features_scaled, train_profile16)

    for split in predict_splits:
        manifest = _load_manifest(cache_root, split)
        records_jsonl = _resolve_records_jsonl(cache_root, manifest, split)
        features = _prepare_feature_matrix(
            np.load(_resolve_feature_array(cache_root, split, feature_name), mmap_mode="r"),
            feature_name=feature_name,
        )
        logcount_targets, profile_targets = _load_records_targets(records_jsonl)
        if features.shape[0] != int(manifest.get("n_records", 0)):
            raise ValueError(
                f"Feature row count mismatch for split={split!r}: "
                f"{features.shape[0]} vs manifest {manifest.get('n_records', 0)}"
            )
        if features.shape[0] != logcount_targets.shape[0]:
            raise ValueError(
                f"records_jsonl row count mismatch for split={split!r}: "
                f"features={features.shape[0]} targets={logcount_targets.shape[0]}"
            )

        features_scaled = scaler.transform(features)
        pred_logcount = count_model.predict(features_scaled).astype(np.float32, copy=False)
        pred_profile16 = _sanitize_profile_predictions(profile_model.predict(features_scaled))

        np.save(output_root / f"{split}_logcount.f32.npy", pred_logcount)
        np.save(output_root / f"{split}_profile16.f32.npy", pred_profile16.astype(np.float32, copy=False))

        teacher_manifest = {
            "split": split,
            "n_records": int(features.shape[0]),
            "record_sha1": manifest.get("record_sha1", ""),
            "source_cache_dir": str(cache_root),
            "feature_name": feature_name,
            "targets": ["profile16", "logcount"],
            "records_jsonl": str(records_jsonl),
            "train_split": train_split,
            "ridge_alpha": float(alpha),
            "profile_bins": int(profile_targets.shape[1]),
        }
        (output_root / f"teacher_manifest_{split}.json").write_text(
            json.dumps(teacher_manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--feature-name", required=True)
    parser.add_argument("--train-split", default="train")
    parser.add_argument("--predict-splits", nargs="+", default=["train", "valid"])
    parser.add_argument("--alpha", type=float, default=1.0)
    return parser.parse_args(list(argv) if argv is not None else None)


def main(argv: Iterable[str] | None = None) -> None:
    args = parse_args(argv)
    export_teacher_cache(
        cache_dir=args.cache_dir,
        output_dir=args.output_dir,
        feature_name=args.feature_name,
        train_split=args.train_split,
        predict_splits=args.predict_splits,
        alpha=args.alpha,
    )


if __name__ == "__main__":
    main()
