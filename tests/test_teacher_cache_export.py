"""Tests for exporting record-aligned teacher caches from foundation features."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
IMPORT_ROOTS = [
    REPO_ROOT / "vendor" / "transchrombp",
    REPO_ROOT / "src",
    REPO_ROOT,
]
for import_root in IMPORT_ROOTS:
    if (import_root / "transchrombp").exists() and str(import_root) not in sys.path:
        sys.path.insert(0, str(import_root))
        break

from transchrombp.evaluation.teacher_cache_export import export_teacher_cache


def _write_records_jsonl(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=True) + "\n")


def _write_manifest(
    cache_dir: Path,
    *,
    split: str,
    feature_name: str,
    n_records: int,
    record_sha1: str,
    records_jsonl: Path,
) -> None:
    manifest = {
        "split": split,
        "n_records": n_records,
        "record_sha1": record_sha1,
        "features": [feature_name],
        "records_jsonl": str(records_jsonl),
    }
    (cache_dir / f"manifest_{split}.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _make_profile_row(scale: float) -> list[float]:
    base = np.linspace(1.0, 16.0, 16, dtype=np.float32)
    return (base * np.float32(scale)).tolist()


def test_export_teacher_cache_writes_manifest_and_arrays(tmp_path: Path) -> None:
    cache_dir = tmp_path / "foundation_cache"
    output_dir = tmp_path / "teacher_cache"
    cache_dir.mkdir()

    feature_name = "layer_07__bins16_mean"
    train_features = np.asarray(
        [
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
        ],
        dtype=np.float32,
    )
    valid_features = np.asarray(
        [
            [2.0, 1.0],
            [1.0, 2.0],
        ],
        dtype=np.float32,
    )
    np.save(cache_dir / f"train_{feature_name}.f32.npy", train_features)
    np.save(cache_dir / f"valid_{feature_name}.f32.npy", valid_features)

    train_rows = [
        {
            "chrom": "chr1",
            "center": 100,
            "source": "peak",
            "true_logcount": 1.0,
            "true_profile16": _make_profile_row(1.0),
        },
        {
            "chrom": "chr1",
            "center": 200,
            "source": "nonpeak",
            "true_logcount": 2.0,
            "true_profile16": _make_profile_row(2.0),
        },
        {
            "chrom": "chr1",
            "center": 300,
            "source": "peak",
            "true_logcount": 3.0,
            "true_profile16": _make_profile_row(3.0),
        },
    ]
    valid_rows = [
        {
            "chrom": "chr2",
            "center": 400,
            "source": "peak",
            "true_logcount": 4.0,
            "true_profile16": _make_profile_row(4.0),
        },
        {
            "chrom": "chr2",
            "center": 500,
            "source": "nonpeak",
            "true_logcount": 5.0,
            "true_profile16": _make_profile_row(5.0),
        },
    ]

    train_records_path = cache_dir / "records_train.jsonl"
    valid_records_path = cache_dir / "records_valid.jsonl"
    _write_records_jsonl(train_records_path, train_rows)
    _write_records_jsonl(valid_records_path, valid_rows)
    record_sha1 = "teacher-cache-export-test"
    _write_manifest(
        cache_dir,
        split="train",
        feature_name=feature_name,
        n_records=len(train_rows),
        record_sha1=record_sha1,
        records_jsonl=train_records_path,
    )
    _write_manifest(
        cache_dir,
        split="valid",
        feature_name=feature_name,
        n_records=len(valid_rows),
        record_sha1=record_sha1,
        records_jsonl=valid_records_path,
    )

    export_teacher_cache(
        cache_dir=cache_dir,
        output_dir=output_dir,
        feature_name=feature_name,
        train_split="train",
        predict_splits=("train", "valid"),
        alpha=1.0e-6,
    )

    train_manifest = json.loads((output_dir / "teacher_manifest_train.json").read_text(encoding="utf-8"))
    valid_manifest = json.loads((output_dir / "teacher_manifest_valid.json").read_text(encoding="utf-8"))
    assert train_manifest["targets"] == ["profile16", "logcount"]
    assert valid_manifest["targets"] == ["profile16", "logcount"]
    assert train_manifest["record_sha1"] == record_sha1
    assert valid_manifest["record_sha1"] == record_sha1

    train_profile16 = np.load(output_dir / "train_profile16.f32.npy")
    train_logcount = np.load(output_dir / "train_logcount.f32.npy")
    valid_profile16 = np.load(output_dir / "valid_profile16.f32.npy")
    valid_logcount = np.load(output_dir / "valid_logcount.f32.npy")

    assert train_profile16.shape == (3, 16)
    assert train_logcount.shape == (3,)
    assert valid_profile16.shape == (2, 16)
    assert valid_logcount.shape == (2,)
    assert np.isfinite(train_profile16).all()
    assert np.isfinite(valid_profile16).all()
    assert np.isfinite(train_logcount).all()
    assert np.isfinite(valid_logcount).all()
    assert (train_profile16 >= 0.0).all()
    assert (valid_profile16 >= 0.0).all()
    assert (train_profile16.sum(axis=1) > 0.0).all()
    assert (valid_profile16.sum(axis=1) > 0.0).all()


def test_export_teacher_cache_accepts_bins16_feature_tensors(tmp_path: Path) -> None:
    cache_dir = tmp_path / "foundation_cache"
    output_dir = tmp_path / "teacher_cache"
    cache_dir.mkdir()

    feature_name = "layer_07__bins16_mean"
    train_features = np.asarray(
        [
            [[1.0, 0.5], [0.0, 1.0], [1.0, 1.5]],
            [[0.5, 1.0], [1.0, 0.5], [1.5, 1.0]],
            [[1.5, 1.5], [0.5, 0.5], [1.0, 1.0]],
        ],
        dtype=np.float32,
    )
    valid_features = np.asarray(
        [
            [[2.0, 1.0], [1.0, 2.0], [1.5, 2.5]],
            [[1.0, 2.0], [2.0, 1.0], [2.5, 1.5]],
        ],
        dtype=np.float32,
    )
    np.save(cache_dir / f"train_{feature_name}.f32.npy", train_features)
    np.save(cache_dir / f"valid_{feature_name}.f32.npy", valid_features)

    train_rows = [
        {
            "chrom": "chr1",
            "center": 100,
            "source": "peak",
            "true_logcount": 1.0,
            "true_profile16": _make_profile_row(1.0),
        },
        {
            "chrom": "chr1",
            "center": 200,
            "source": "nonpeak",
            "true_logcount": 2.0,
            "true_profile16": _make_profile_row(2.0),
        },
        {
            "chrom": "chr1",
            "center": 300,
            "source": "peak",
            "true_logcount": 3.0,
            "true_profile16": _make_profile_row(3.0),
        },
    ]
    valid_rows = [
        {
            "chrom": "chr2",
            "center": 400,
            "source": "peak",
            "true_logcount": 4.0,
            "true_profile16": _make_profile_row(4.0),
        },
        {
            "chrom": "chr2",
            "center": 500,
            "source": "nonpeak",
            "true_logcount": 5.0,
            "true_profile16": _make_profile_row(5.0),
        },
    ]

    train_records_path = cache_dir / "records_train.jsonl"
    valid_records_path = cache_dir / "records_valid.jsonl"
    _write_records_jsonl(train_records_path, train_rows)
    _write_records_jsonl(valid_records_path, valid_rows)
    record_sha1 = "teacher-cache-export-bins16"
    _write_manifest(
        cache_dir,
        split="train",
        feature_name=feature_name,
        n_records=len(train_rows),
        record_sha1=record_sha1,
        records_jsonl=train_records_path,
    )
    _write_manifest(
        cache_dir,
        split="valid",
        feature_name=feature_name,
        n_records=len(valid_rows),
        record_sha1=record_sha1,
        records_jsonl=valid_records_path,
    )

    export_teacher_cache(
        cache_dir=cache_dir,
        output_dir=output_dir,
        feature_name=feature_name,
        train_split="train",
        predict_splits=("train", "valid"),
        alpha=1.0e-6,
    )

    train_profile16 = np.load(output_dir / "train_profile16.f32.npy")
    train_logcount = np.load(output_dir / "train_logcount.f32.npy")
    valid_profile16 = np.load(output_dir / "valid_profile16.f32.npy")
    valid_logcount = np.load(output_dir / "valid_logcount.f32.npy")

    assert train_profile16.shape == (3, 16)
    assert train_logcount.shape == (3,)
    assert valid_profile16.shape == (2, 16)
    assert valid_logcount.shape == (2,)
    assert np.isfinite(train_profile16).all()
    assert np.isfinite(valid_profile16).all()
    assert np.isfinite(train_logcount).all()
    assert np.isfinite(valid_logcount).all()


def test_export_teacher_cache_requires_records_jsonl_metadata(tmp_path: Path) -> None:
    cache_dir = tmp_path / "foundation_cache"
    output_dir = tmp_path / "teacher_cache"
    cache_dir.mkdir()

    feature_name = "layer_07__bins16_mean"
    np.save(cache_dir / f"train_{feature_name}.f32.npy", np.ones((2, 2), dtype=np.float32))
    (cache_dir / "manifest_train.json").write_text(
        json.dumps(
            {
                "split": "train",
                "n_records": 2,
                "record_sha1": "missing-records",
                "features": [feature_name],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="records_jsonl"):
        export_teacher_cache(
            cache_dir=cache_dir,
            output_dir=output_dir,
            feature_name=feature_name,
            train_split="train",
            predict_splits=("train",),
            alpha=1.0,
        )
