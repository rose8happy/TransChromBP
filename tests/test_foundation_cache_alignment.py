"""Regression tests for foundation-cache/dataset alignment."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

pytest.importorskip("numpy")
pytest.importorskip("pyBigWig")
pytest.importorskip("pyfaidx")

import pyBigWig


REPO_ROOT = Path(__file__).resolve().parents[1]
TRANSCHROMBP_ROOT = REPO_ROOT / "vendor" / "transchrombp"
if not TRANSCHROMBP_ROOT.exists():
    TRANSCHROMBP_ROOT = REPO_ROOT
if str(TRANSCHROMBP_ROOT) not in sys.path:
    sys.path.insert(0, str(TRANSCHROMBP_ROOT))

from transchrombp.data.real_data import ChromBPNetBigWigDataset, compute_record_sha1
from transchrombp.evaluation.evaluate_checkpoint import build_eval_dataset


def _load_build_foundation_cache_module():
    module_path = TRANSCHROMBP_ROOT / "transchrombp" / "scripts" / "build_foundation_cache.py"
    if not module_path.exists():
        module_path = TRANSCHROMBP_ROOT / "scripts" / "build_foundation_cache.py"
    spec = importlib.util.spec_from_file_location("build_foundation_cache_module", module_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_fasta(path: Path) -> None:
    seq = ("AACCGGTTACGTAAGTCCGA" * 120)[:2400]
    path.write_text(f">chr1\n{seq}\n", encoding="utf-8")


def _write_bigwig(path: Path) -> None:
    chrom_sizes = [("chr1", 2400)]
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader(chrom_sizes)
    starts = list(range(2400))
    ends = [start + 1 for start in starts]
    values = [float((start % 7) + 1) for start in starts]
    bw.addEntries(["chr1"] * 2400, starts, ends=ends, values=values)
    bw.close()


def _write_regions(path: Path, centers: list[int]) -> None:
    lines = []
    for idx, center in enumerate(centers, start=1):
        start = center - 57
        end = center + 57
        lines.append(f"chr1\t{start}\t{end}\tregion{idx}\t1\t.\t0\t0\t0\t57")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _toy_train_and_data_configs(tmp_path: Path) -> tuple[dict, dict]:
    fasta = tmp_path / "toy.fa"
    bigwig = tmp_path / "toy.bw"
    peaks = tmp_path / "peaks.bed"
    nonpeaks = tmp_path / "nonpeaks.bed"
    folds = tmp_path / "folds.json"

    _write_fasta(fasta)
    _write_bigwig(bigwig)
    _write_regions(peaks, [200, 320, 440, 560, 680, 800, 920, 1040])
    _write_regions(nonpeaks, [1200, 1320, 1440, 1560, 1680, 1800, 1920, 2040])
    folds.write_text(json.dumps({"train": ["chr1"], "valid": ["chr1"], "test": ["chr1"]}), encoding="utf-8")

    train_cfg = {
        "seed": 1234,
        "data": {
            "input_len": 128,
            "supervised_bp": 64,
            "profile_bin_size": 1,
            "nonpeak_ratio": 1.0,
            "test_region_source": "both",
            "max_test_regions": 5,
        },
    }
    data_source_cfg = {
        "genome_fasta": str(fasta),
        "folds_json": str(folds),
        "input": {
            "bigwig": str(bigwig),
            "peaks_bed": str(peaks),
            "nonpeaks_bed": str(nonpeaks),
        },
    }
    return train_cfg, data_source_cfg


def test_test_split_cache_builder_matches_eval_records(tmp_path: Path) -> None:
    module = _load_build_foundation_cache_module()
    train_cfg, data_source_cfg = _toy_train_and_data_configs(tmp_path)

    cache_ds = module.build_dataset_for_cache("test", train_cfg, data_source_cfg)
    eval_ds = build_eval_dataset(train_cfg, data_source_cfg, split="test")

    assert len(cache_ds) == len(eval_ds) == 5
    assert compute_record_sha1(cache_ds.records) == compute_record_sha1(eval_ds.records)


def test_real_data_dataset_exposes_max_records_for_manifest_metadata(tmp_path: Path) -> None:
    train_cfg, data_source_cfg = _toy_train_and_data_configs(tmp_path)

    dataset = ChromBPNetBigWigDataset(
        genome_fasta=data_source_cfg["genome_fasta"],
        bigwig_path=data_source_cfg["input"]["bigwig"],
        peaks_bed=data_source_cfg["input"]["peaks_bed"],
        nonpeaks_bed=data_source_cfg["input"]["nonpeaks_bed"],
        folds_json=data_source_cfg["folds_json"],
        split="test",
        input_len=train_cfg["data"]["input_len"],
        supervised_bp=train_cfg["data"]["supervised_bp"],
        profile_bin_size=train_cfg["data"]["profile_bin_size"],
        max_jitter=0,
        seed=train_cfg["seed"] + 10_000,
        nonpeak_ratio=train_cfg["data"]["nonpeak_ratio"],
        max_records=train_cfg["data"]["max_test_regions"],
        region_source=train_cfg["data"]["test_region_source"],
    )

    assert dataset.max_records == 5
