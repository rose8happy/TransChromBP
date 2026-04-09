"""Real-data dataset utilities for TransChromBP."""

from __future__ import annotations

import hashlib
import json
import os
import warnings
from collections import Counter
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Optional, Sequence

import numpy as np
import pyBigWig
import torch
from torch import Tensor
from torch.utils.data import Dataset

if TYPE_CHECKING:
    from pyfaidx import Fasta


_BASE_TO_INDEX = np.full(256, -1, dtype=np.int8)
for _base, _idx in (("A", 0), ("C", 1), ("G", 2), ("T", 3), ("a", 0), ("c", 1), ("g", 2), ("t", 3)):
    _BASE_TO_INDEX[ord(_base)] = _idx


@dataclass(frozen=True)
class RegionRecord:
    """Single centered genomic region."""

    chrom: str
    center: int
    source: str


def load_fold_chroms(folds_json: str, split: str) -> list[str]:
    with open(folds_json, "r", encoding="utf-8") as f:
        folds = json.load(f)
    chroms = folds.get(split, [])
    if not isinstance(chroms, list) or not chroms:
        raise ValueError(f"No chromosomes found for split={split!r} in {folds_json}")
    return [str(chrom) for chrom in chroms]


def _maybe_int(value: str) -> Optional[int]:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def load_regions_from_bed(path: str, allowed_chroms: Optional[Sequence[str]] = None, source: str = "regions") -> list[RegionRecord]:
    allowed = set(allowed_chroms or [])
    records: list[RegionRecord] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"BED line in {path} has fewer than 3 columns: {line}")

            chrom = fields[0]
            if allowed and chrom not in allowed:
                continue

            start = int(fields[1])
            end = int(fields[2])
            summit = _maybe_int(fields[9]) if len(fields) > 9 else None
            center = start + summit if summit is not None else (start + end) // 2
            records.append(RegionRecord(chrom=chrom, center=center, source=source))

    if not records:
        raise ValueError(f"No usable regions found in {path}")
    return records


def load_bigwig_chrom_sizes(path: str) -> dict[str, int]:
    bigwig = pyBigWig.open(path)
    if bigwig is None:
        raise ValueError(f"Failed to open bigWig file: {path}")
    try:
        chrom_sizes = bigwig.chroms() or {}
        return {str(chrom): int(size) for chrom, size in chrom_sizes.items()}
    finally:
        bigwig.close()


def filter_records_by_chroms(
    records: Sequence[RegionRecord],
    available_chroms: Sequence[str],
    path: str,
    source: str,
) -> list[RegionRecord]:
    available = set(str(chrom) for chrom in available_chroms)
    if not available:
        raise ValueError(f"No chromosomes found in bigWig while filtering {source} records from {path}")

    filtered: list[RegionRecord] = []
    dropped: Counter[str] = Counter()
    for record in records:
        if record.chrom in available:
            filtered.append(record)
        else:
            dropped[record.chrom] += 1

    if dropped:
        dropped_desc = ", ".join(f"{chrom}x{count}" for chrom, count in sorted(dropped.items()))
        warnings.warn(
            f"Filtered {sum(dropped.values())} {source} records from {path} because the bigWig lacks "
            f"matching chromosomes: {dropped_desc}",
            RuntimeWarning,
            stacklevel=2,
        )
    return filtered


def _subsample_records(records: Sequence[RegionRecord], max_records: int, seed: int) -> list[RegionRecord]:
    if max_records <= 0 or len(records) <= max_records:
        return list(records)
    rng = np.random.default_rng(seed)
    indices = np.sort(rng.choice(len(records), size=max_records, replace=False))
    return [records[int(idx)] for idx in indices]


def _reverse_complement_onehot(seq: np.ndarray) -> np.ndarray:
    if seq.ndim != 2 or seq.shape[1] != 4:
        raise ValueError(f"Expected one-hot sequence with shape (L, 4), got {tuple(seq.shape)}")
    return np.ascontiguousarray(seq[::-1, ::-1])


def _reverse_profile_counts(profile_counts: np.ndarray) -> np.ndarray:
    if profile_counts.ndim != 1:
        raise ValueError(f"Expected 1D profile counts, got shape {tuple(profile_counts.shape)}")
    return np.ascontiguousarray(profile_counts[::-1])


def compute_record_sha1(records: Sequence[RegionRecord]) -> str:
    """Compute SHA1 hash of the ordered record list for cache validation."""
    h = hashlib.sha1()
    for r in records:
        h.update(f"{r.chrom}:{r.center}:{r.source}\n".encode("utf-8"))
    return h.hexdigest()


def resolve_dataset_seed(base_seed: int, split: str) -> int:
    """Keep non-train splits on a shared deterministic seed schedule."""
    normalized = str(split).strip().lower()
    return int(base_seed) if normalized == "train" else int(base_seed) + 10_000


def _load_fasta_class() -> type["Fasta"]:
    try:
        from pyfaidx import Fasta
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "pyfaidx is required to read genome FASTA files; install it before using "
            "ChromBPNetBigWigDataset sequence-loading paths."
        ) from exc
    return Fasta


class ChromBPNetBigWigDataset(Dataset):
    """Real-data dataset backed by FASTA + bigWig + BED regions."""

    def __init__(
        self,
        genome_fasta: str,
        bigwig_path: str,
        peaks_bed: str,
        nonpeaks_bed: str,
        folds_json: str,
        split: str,
        input_len: int,
        supervised_bp: int,
        profile_bin_size: int = 1,
        max_jitter: int = 0,
        peak_max_jitter: Optional[int] = None,
        nonpeak_max_jitter: Optional[int] = None,
        seed: int = 1234,
        nonpeak_ratio: float = 1.0,
        max_records: int = 0,
        region_source: str = "both",
        random_revcomp: bool = False,
        revcomp_prob: float = 0.5,
        track_total_count_target: float = 0.0,
        genos_cache_dir: str = "",
        genos_cache_features: Sequence[str] = (),
        foundation_cache_dir: str = "",
        foundation_cache_features: Sequence[str] = (),
        teacher_cache_dir: str = "",
        teacher_target_names: Sequence[str] = (),
    ) -> None:
        if input_len <= 0:
            raise ValueError(f"input_len must be positive, got {input_len}")
        if supervised_bp <= 0:
            raise ValueError(f"supervised_bp must be positive, got {supervised_bp}")
        if profile_bin_size <= 0:
            raise ValueError(f"profile_bin_size must be positive, got {profile_bin_size}")
        if supervised_bp % profile_bin_size != 0:
            raise ValueError(
                f"supervised_bp ({supervised_bp}) must be divisible by profile_bin_size ({profile_bin_size})"
            )

        self.genome_fasta = str(genome_fasta)
        self.bigwig_path = str(bigwig_path)
        self.split = str(split)
        self.input_len = int(input_len)
        self.supervised_bp = int(supervised_bp)
        self.profile_bin_size = int(profile_bin_size)
        self.output_len = self.supervised_bp // self.profile_bin_size
        self.max_jitter = int(max_jitter) if self.split == "train" else 0
        self.peak_max_jitter = self.max_jitter if peak_max_jitter is None else int(peak_max_jitter)
        self.nonpeak_max_jitter = self.max_jitter if nonpeak_max_jitter is None else int(nonpeak_max_jitter)
        if self.peak_max_jitter < 0 or self.nonpeak_max_jitter < 0:
            raise ValueError(
                f"peak/nonpeak jitter must be non-negative, got {self.peak_max_jitter} and {self.nonpeak_max_jitter}"
            )
        if self.split != "train":
            self.peak_max_jitter = 0
            self.nonpeak_max_jitter = 0
        self.seed = int(seed)
        self.nonpeak_ratio = float(nonpeak_ratio)
        self.max_records = int(max_records)
        self.region_source = str(region_source).lower()
        self.random_revcomp = bool(random_revcomp) if self.split == "train" else False
        self.revcomp_prob = float(revcomp_prob)
        if not 0.0 <= self.revcomp_prob <= 1.0:
            raise ValueError(f"revcomp_prob must lie in [0, 1], got {self.revcomp_prob}")
        self.track_total_count_target = float(track_total_count_target)
        if self.track_total_count_target < 0.0:
            raise ValueError(
                f"track_total_count_target must be non-negative, got {self.track_total_count_target}"
            )
        if self.region_source not in {"both", "peak", "nonpeak"}:
            raise ValueError(
                f"region_source must be one of 'both', 'peak', 'nonpeak'; got {region_source!r}"
            )
        self._bigwig_chrom_sizes = load_bigwig_chrom_sizes(self.bigwig_path)

        split_chroms = load_fold_chroms(folds_json, split)
        peaks: list[RegionRecord] = []
        nonpeaks: list[RegionRecord] = []
        records: list[RegionRecord] = []
        target_nonpeaks = 0

        if self.region_source in {"both", "peak"}:
            peaks = load_regions_from_bed(peaks_bed, split_chroms, source="peak")
            peaks = filter_records_by_chroms(peaks, self._bigwig_chrom_sizes.keys(), peaks_bed, source="peak")
            records.extend(peaks)

        nonpeaks_bed = str(nonpeaks_bed or "").strip()
        if self.region_source in {"both", "nonpeak"}:
            if not nonpeaks_bed:
                raise ValueError("nonpeaks_bed is required when region_source includes nonpeak regions")
            nonpeaks = load_regions_from_bed(nonpeaks_bed, split_chroms, source="nonpeak")
            nonpeaks = filter_records_by_chroms(
                nonpeaks,
                self._bigwig_chrom_sizes.keys(),
                nonpeaks_bed,
                source="nonpeak",
            )

            if self.region_source == "both" and peaks:
                if self.nonpeak_ratio > 0:
                    target_nonpeaks = min(len(nonpeaks), max(1, int(round(len(peaks) * self.nonpeak_ratio))))
                    if self.split == "train" and max_records <= 0:
                        # Keep the full nonpeak pool so the sampler can draw a fresh subset every epoch.
                        records.extend(nonpeaks)
                    else:
                        records.extend(_subsample_records(nonpeaks, target_nonpeaks, self.seed + 97))
            else:
                if self.nonpeak_ratio <= 0:
                    raise ValueError("nonpeak_ratio must be positive when training on nonpeak-only data")
                if self.nonpeak_ratio >= 1.0:
                    target_nonpeaks = len(nonpeaks)
                else:
                    target_nonpeaks = max(1, int(round(len(nonpeaks) * self.nonpeak_ratio)))
                records.extend(_subsample_records(nonpeaks, target_nonpeaks, self.seed + 97))

        if self.max_records > 0:
            records = _subsample_records(records, self.max_records, self.seed + 131)
        if not records:
            raise ValueError(f"No records available for split={split!r}")

        self.peak_record_count = len(peaks)
        self.nonpeak_record_count = len(nonpeaks)
        self.target_nonpeak_count = target_nonpeaks
        self.supports_epoch_resampling = bool(
            self.split == "train"
            and self.region_source == "both"
            and self.peak_record_count > 0
            and self.nonpeak_record_count > 0
            and self.target_nonpeak_count > 0
            and self.max_records <= 0
        )
        self.records = records
        self._fasta: Optional["Fasta"] = None
        self._fasta_pid: Optional[int] = None
        self._bigwig: Optional[pyBigWig.pyBigWig] = None
        self._bigwig_pid: Optional[int] = None
        self._track_scale_factor: Optional[float] = None
        self._track_total_count_observed: Optional[float] = None

        # ── Genos cache ──
        self._genos_caches: dict[str, np.ndarray] = {}
        self._foundation_caches: dict[str, np.ndarray] = {}
        self._teacher_targets: dict[str, np.ndarray] = {}
        genos_cache_dir = str(genos_cache_dir or "").strip()
        genos_cache_features = list(genos_cache_features or [])
        if genos_cache_dir and genos_cache_features:
            self._load_genos_cache(genos_cache_dir, genos_cache_features)

        foundation_cache_dir = str(foundation_cache_dir or "").strip()
        foundation_cache_features = list(foundation_cache_features or [])
        if foundation_cache_dir and foundation_cache_features:
            self._load_foundation_cache(foundation_cache_dir, foundation_cache_features)

        teacher_cache_dir = str(teacher_cache_dir or "").strip()
        teacher_target_names = list(teacher_target_names or [])
        if teacher_cache_dir and teacher_target_names:
            self._load_teacher_cache(teacher_cache_dir, teacher_target_names)

    def _find_cache_array_path(self, cache_dir: str, name: str) -> str:
        candidates = [
            os.path.join(cache_dir, f"{self.split}_{name}.f16.npy"),
            os.path.join(cache_dir, f"{self.split}_{name}.f32.npy"),
            os.path.join(cache_dir, f"{self.split}_{name}.npy"),
        ]
        for path in candidates:
            if os.path.isfile(path):
                return path
        raise FileNotFoundError(
            f"Cache array for {name!r} not found under {cache_dir}; tried {candidates!r}"
        )

    def _load_genos_cache(self, cache_dir: str, features: list[str]) -> None:
        """Load and validate Genos cached summary features from memmap files."""
        manifest_path = os.path.join(cache_dir, f"manifest_{self.split}.json")
        if not os.path.isfile(manifest_path):
            raise FileNotFoundError(f"Genos cache manifest not found: {manifest_path}")

        with open(manifest_path, "r", encoding="utf-8") as f:
            manifest = json.load(f)

        # Validate record count
        if int(manifest.get("n_records", 0)) != len(self.records):
            raise ValueError(
                f"Genos cache n_records mismatch: manifest says {manifest['n_records']}, "
                f"dataset has {len(self.records)} records"
            )
        if str(manifest.get("split", "")) != self.split:
            raise ValueError(
                f"Genos cache split mismatch: manifest says {manifest.get('split', '')!r}, "
                f"dataset split is {self.split!r}"
            )
        if int(manifest.get("input_len", 0)) != self.input_len:
            raise ValueError(
                f"Genos cache input_len mismatch: manifest says {manifest.get('input_len', 0)}, "
                f"dataset input_len is {self.input_len}"
            )
        # Validate record_sha1
        expected_sha1 = compute_record_sha1(self.records)
        if manifest.get("record_sha1", "") != expected_sha1:
            raise ValueError(
                f"Genos cache record_sha1 mismatch: manifest={manifest.get('record_sha1', '')!r}, "
                f"dataset={expected_sha1!r}. Records order may have changed."
            )
        manifest_features = set(str(x) for x in manifest.get("features", []))

        for feat_name in features:
            if feat_name not in manifest_features:
                raise ValueError(
                    f"Requested Genos cache feature {feat_name!r} is missing from manifest features "
                    f"{sorted(manifest_features)!r}"
                )
            npy_path = self._find_cache_array_path(cache_dir, feat_name)
            mmap = np.load(npy_path, mmap_mode="r")
            if mmap.shape[0] != len(self.records):
                raise ValueError(
                    f"Genos cache {feat_name} has {mmap.shape[0]} rows, "
                    f"expected {len(self.records)}"
                )
            self._genos_caches[feat_name] = mmap

    def _load_foundation_cache(self, cache_dir: str, features: list[str]) -> None:
        """Load generic foundation feature cache with the same record validation as Genos."""
        manifest_path = os.path.join(cache_dir, f"manifest_{self.split}.json")
        if not os.path.isfile(manifest_path):
            raise FileNotFoundError(f"Foundation cache manifest not found: {manifest_path}")
        with open(manifest_path, "r", encoding="utf-8") as f:
            manifest = json.load(f)

        if int(manifest.get("n_records", 0)) != len(self.records):
            raise ValueError(
                f"Foundation cache n_records mismatch: manifest says {manifest.get('n_records', 0)}, "
                f"dataset has {len(self.records)} records"
            )
        if str(manifest.get("split", "")) != self.split:
            raise ValueError(
                f"Foundation cache split mismatch: manifest says {manifest.get('split', '')!r}, "
                f"dataset split is {self.split!r}"
            )
        expected_sha1 = compute_record_sha1(self.records)
        if manifest.get("record_sha1", "") != expected_sha1:
            raise ValueError(
                f"Foundation cache record_sha1 mismatch: manifest={manifest.get('record_sha1', '')!r}, "
                f"dataset={expected_sha1!r}"
            )
        manifest_features = set(str(x) for x in manifest.get("features", []))
        for feat_name in features:
            if feat_name not in manifest_features:
                raise ValueError(
                    f"Requested foundation cache feature {feat_name!r} is missing from manifest "
                    f"features {sorted(manifest_features)!r}"
                )
            mmap = np.load(self._find_cache_array_path(cache_dir, feat_name), mmap_mode="r")
            if mmap.shape[0] != len(self.records):
                raise ValueError(
                    f"Foundation cache {feat_name} has {mmap.shape[0]} rows, expected {len(self.records)}"
                )
            self._foundation_caches[feat_name] = mmap

    def _load_teacher_cache(self, cache_dir: str, targets: list[str]) -> None:
        """Load teacher targets aligned to the dataset records."""
        manifest_path = os.path.join(cache_dir, f"teacher_manifest_{self.split}.json")
        if not os.path.isfile(manifest_path):
            fallback = os.path.join(cache_dir, f"manifest_{self.split}.json")
            if os.path.isfile(fallback):
                manifest_path = fallback
        if not os.path.isfile(manifest_path):
            raise FileNotFoundError(
                f"Teacher cache manifest not found under {cache_dir} for split={self.split!r}"
            )
        with open(manifest_path, "r", encoding="utf-8") as f:
            manifest = json.load(f)

        if int(manifest.get("n_records", 0)) != len(self.records):
            raise ValueError(
                f"Teacher cache n_records mismatch: manifest says {manifest.get('n_records', 0)}, "
                f"dataset has {len(self.records)} records"
            )
        if str(manifest.get("split", "")) != self.split:
            raise ValueError(
                f"Teacher cache split mismatch: manifest says {manifest.get('split', '')!r}, "
                f"dataset split is {self.split!r}"
            )
        expected_sha1 = compute_record_sha1(self.records)
        if manifest.get("record_sha1", "") != expected_sha1:
            raise ValueError(
                f"Teacher cache record_sha1 mismatch: manifest={manifest.get('record_sha1', '')!r}, "
                f"dataset={expected_sha1!r}"
            )
        manifest_targets = set(str(x) for x in manifest.get("targets", manifest.get("features", [])))
        for target_name in targets:
            if target_name not in manifest_targets:
                raise ValueError(
                    f"Requested teacher target {target_name!r} is missing from manifest entries "
                    f"{sorted(manifest_targets)!r}"
                )
            mmap = np.load(self._find_cache_array_path(cache_dir, target_name), mmap_mode="r")
            if mmap.shape[0] != len(self.records):
                raise ValueError(
                    f"Teacher cache {target_name} has {mmap.shape[0]} rows, expected {len(self.records)}"
                )
            self._teacher_targets[target_name] = mmap

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(self, idx: int) -> dict[str, Any]:
        record = self.records[idx]

        jitter = 0
        jitter_limit = self.peak_max_jitter if record.source == "peak" else self.nonpeak_max_jitter
        if jitter_limit > 0:
            jitter = int(np.random.randint(-jitter_limit, jitter_limit + 1))
        center = int(record.center + jitter)

        seq_start = center - (self.input_len // 2)
        profile_start = center - (self.supervised_bp // 2)

        seq = self._fetch_onehot(record.chrom, seq_start, self.input_len)
        profile_counts = self._fetch_profile(record.chrom, profile_start, self.supervised_bp)
        if self.random_revcomp and np.random.random() < self.revcomp_prob:
            seq = _reverse_complement_onehot(seq)
            profile_counts = _reverse_profile_counts(profile_counts)

        sample = {
            "seq": torch.from_numpy(seq),
            "profile_counts": torch.from_numpy(profile_counts),
            "region_source": record.source,
        }
        for feat_name, mmap in self._genos_caches.items():
            sample[f"genos_{feat_name}"] = torch.from_numpy(
                np.array(mmap[idx], dtype=np.float32)
            )
        for feat_name, mmap in self._foundation_caches.items():
            sample[f"foundation_{feat_name}"] = torch.from_numpy(
                np.array(mmap[idx], dtype=np.float32)
            )
        for target_name, mmap in self._teacher_targets.items():
            target = np.array(mmap[idx], dtype=np.float32)
            if self.random_revcomp and target.ndim == 1 and "profile" in target_name:
                target = np.ascontiguousarray(target[::-1])
            sample[f"teacher_{target_name}"] = torch.from_numpy(target)
        return sample

    def _get_fasta(self) -> "Fasta":
        current_pid = os.getpid()
        if self._fasta is None or self._fasta_pid != current_pid:
            fasta_cls = _load_fasta_class()
            self._fasta = fasta_cls(self.genome_fasta, as_raw=True, sequence_always_upper=True)
            self._fasta_pid = current_pid
        return self._fasta

    def _get_bigwig(self) -> pyBigWig.pyBigWig:
        current_pid = os.getpid()
        if self._bigwig is None or self._bigwig_pid != current_pid:
            if self._bigwig is not None:
                try:
                    self._bigwig.close()
                except Exception:
                    pass
            self._bigwig = pyBigWig.open(self.bigwig_path)
            self._bigwig_pid = current_pid
            if self._track_scale_factor is None:
                self._track_scale_factor = self._compute_track_scale_factor(self._bigwig)
        return self._bigwig

    def _compute_track_scale_factor(self, bigwig: pyBigWig.pyBigWig) -> float:
        if self.track_total_count_target <= 0.0:
            self._track_total_count_observed = None
            return 1.0

        header = bigwig.header() or {}
        observed_total = float(header.get("sumData", 0.0))
        self._track_total_count_observed = observed_total
        if observed_total <= 0.0:
            return 1.0
        return self.track_total_count_target / observed_total

    def _chrom_len(self, chrom: str) -> int:
        return len(self._get_fasta()[chrom])

    def _fetch_onehot(self, chrom: str, start: int, length: int) -> np.ndarray:
        chrom_len = self._chrom_len(chrom)
        end = start + length
        out = np.zeros((length, 4), dtype=np.float32)

        clipped_start = max(0, start)
        clipped_end = min(end, chrom_len)
        if clipped_start >= clipped_end:
            return out

        seq = self._get_fasta()[chrom][clipped_start:clipped_end]
        if not isinstance(seq, str):
            seq = str(seq)
        encoded = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
        base_indices = _BASE_TO_INDEX[encoded]
        valid = base_indices >= 0

        if np.any(valid):
            rows = np.arange(encoded.shape[0], dtype=np.int64)[valid] + (clipped_start - start)
            out[rows, base_indices[valid]] = 1.0
        return out

    def _fetch_profile(self, chrom: str, start: int, length: int) -> np.ndarray:
        bigwig_chrom_len = int(self._bigwig_chrom_sizes.get(chrom, 0))
        if bigwig_chrom_len <= 0:
            return np.zeros((self.output_len,), dtype=np.float32) if self.profile_bin_size > 1 else np.zeros((length,), dtype=np.float32)

        chrom_len = min(self._chrom_len(chrom), bigwig_chrom_len)
        end = start + length
        out = np.zeros((length,), dtype=np.float32)

        clipped_start = max(0, start)
        clipped_end = min(end, chrom_len)
        if clipped_start < clipped_end:
            values = self._get_bigwig().values(chrom, clipped_start, clipped_end, numpy=True)
            values = np.asarray(values, dtype=np.float32)
            values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
            offset = clipped_start - start
            out[offset : offset + values.shape[0]] = values

        if self.profile_bin_size > 1:
            out = out.reshape(self.output_len, self.profile_bin_size).sum(axis=1)
        if self._track_scale_factor is None:
            self._track_scale_factor = self._compute_track_scale_factor(self._get_bigwig())
        if self._track_scale_factor != 1.0:
            out = out * self._track_scale_factor
        return out.astype(np.float32, copy=False)

    def estimate_total_count_statistics(
        self,
        sample_size: int = 0,
        seed: Optional[int] = None,
        quantile_low: float = 0.0,
        quantile_high: float = 1.0,
    ) -> dict[str, float]:
        if not 0.0 <= quantile_low <= quantile_high <= 1.0:
            raise ValueError(
                "quantile_low/quantile_high must satisfy 0 <= low <= high <= 1, "
                f"got {quantile_low} and {quantile_high}"
            )

        n_records = len(self.records)
        if n_records <= 0:
            raise ValueError("Cannot estimate count statistics on an empty dataset")

        rng = np.random.default_rng(self.seed + 4096 if seed is None else int(seed))
        if sample_size <= 0 or sample_size >= n_records:
            indices = np.arange(n_records, dtype=np.int64)
        else:
            indices = np.sort(rng.choice(n_records, size=int(sample_size), replace=False))

        totals = np.empty(indices.shape[0], dtype=np.float32)
        profile_half_width = self.supervised_bp // 2
        for pos, idx in enumerate(indices):
            record = self.records[int(idx)]
            profile_start = int(record.center) - profile_half_width
            totals[pos] = float(self._fetch_profile(record.chrom, profile_start, self.supervised_bp).sum())

        lower_threshold = float(np.quantile(totals, quantile_low))
        upper_threshold = float(np.quantile(totals, quantile_high))
        keep_mask = (totals >= lower_threshold) & (totals <= upper_threshold)
        trimmed = totals[keep_mask]
        if trimmed.size == 0:
            trimmed = totals

        return {
            "sample_size": float(indices.shape[0]),
            "retained_size": float(trimmed.shape[0]),
            "median_total_count": float(np.median(trimmed)),
            "mean_total_count": float(np.mean(trimmed)),
            "lower_threshold": lower_threshold,
            "upper_threshold": upper_threshold,
            "track_scale_factor": float(self._track_scale_factor or 1.0),
            "track_total_count_observed": float(self._track_total_count_observed or 0.0),
        }

    def __del__(self) -> None:
        if self._bigwig is not None:
            try:
                self._bigwig.close()
            except Exception:
                pass
