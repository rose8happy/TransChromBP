#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import warnings
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np
import yaml


PANEL_QUANTILES = [0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9, 0.98]
WINDOW_BP = 1000
OUTPUT_COLUMNS = [
    "label",
    "source",
    "quantile",
    "chrom",
    "center",
    "observed_total",
    "notes",
]


@dataclass(frozen=True)
class _FallbackRegionRecord:
    chrom: str
    center: int
    source: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a deterministic matched AlphaGenome tutorial panel."
    )
    parser.add_argument("--data-config", type=Path, required=True)
    parser.add_argument("--split", default="test")
    parser.add_argument("--output-csv", type=Path, required=True)
    return parser.parse_args()


def _import_real_data_helpers() -> tuple[Any, Any, Any, Any]:
    try:
        from transchrombp.data.real_data import (
            filter_records_by_chroms,
            load_bigwig_chrom_sizes,
            load_fold_chroms,
            load_regions_from_bed,
        )
    except ModuleNotFoundError as exc:
        if exc.name not in {"transchrombp", "transchrombp.data", "transchrombp.data.real_data"}:
            raise
        return (
            _fallback_load_fold_chroms,
            _fallback_load_regions_from_bed,
            _fallback_load_bigwig_chrom_sizes,
            _fallback_filter_records_by_chroms,
        )

    return (
        load_fold_chroms,
        load_regions_from_bed,
        load_bigwig_chrom_sizes,
        filter_records_by_chroms,
    )


def _fallback_load_fold_chroms(folds_json: str, split: str) -> list[str]:
    with open(folds_json, "r", encoding="utf-8") as handle:
        folds = json.load(handle)
    chroms = folds.get(split, [])
    if not isinstance(chroms, list) or not chroms:
        raise ValueError(f"No chromosomes found for split={split!r} in {folds_json}")
    return [str(chrom) for chrom in chroms]


def _fallback_maybe_int(value: str) -> Optional[int]:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _fallback_load_regions_from_bed(
    path: str,
    allowed_chroms: Optional[Sequence[str]] = None,
    source: str = "regions",
) -> list[_FallbackRegionRecord]:
    allowed = set(allowed_chroms or [])
    records: list[_FallbackRegionRecord] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
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
            summit = _fallback_maybe_int(fields[9]) if len(fields) > 9 else None
            center = start + summit if summit is not None else (start + end) // 2
            records.append(_FallbackRegionRecord(chrom=chrom, center=center, source=source))

    if not records:
        raise ValueError(f"No usable regions found in {path}")
    return records


def _fallback_load_bigwig_chrom_sizes(path: str) -> dict[str, int]:
    try:
        import pyBigWig
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "pyBigWig is required to build the matched panel; use the TransChromBP runtime environment."
        ) from exc

    bigwig = pyBigWig.open(path)
    if bigwig is None:
        raise ValueError(f"Failed to open bigWig file: {path}")
    try:
        chrom_sizes = bigwig.chroms() or {}
        return {str(chrom): int(size) for chrom, size in chrom_sizes.items()}
    finally:
        bigwig.close()


def _fallback_filter_records_by_chroms(
    records: Sequence[_FallbackRegionRecord],
    available_chroms: Sequence[str],
    path: str,
    source: str,
) -> list[_FallbackRegionRecord]:
    available = set(str(chrom) for chrom in available_chroms)
    if not available:
        raise ValueError(f"No chromosomes found in bigWig while filtering {source} records from {path}")

    filtered: list[_FallbackRegionRecord] = []
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


def _quantile_token(quantile: float) -> str:
    return f"{int(round(float(quantile) * 100.0)):02d}"


def _quantile_string(quantile: float) -> str:
    return f"{float(quantile):.2f}"


def _window_total(bigwig: Any, chrom: str, center: int) -> float:
    chrom_sizes = bigwig.chroms() or {}
    chrom_size = int(chrom_sizes[chrom])
    half_window = WINDOW_BP // 2
    start = max(0, int(center) - half_window)
    end = min(chrom_size, start + WINDOW_BP)
    start = max(0, end - WINDOW_BP)
    values = np.asarray(bigwig.values(chrom, start, end, numpy=True), dtype=np.float64)
    return float(np.nan_to_num(values, nan=0.0).sum())


def _pick_quantile_index(
    sorted_rows: Sequence[dict[str, Any]], quantile: float, used_indices: set[int]
) -> int:
    if not sorted_rows:
        raise ValueError("Cannot choose panel rows from an empty source group")

    target_index = float(quantile) * float(len(sorted_rows) - 1)
    ranked_indices = sorted(
        range(len(sorted_rows)),
        key=lambda idx: (abs(float(idx) - target_index), idx),
    )
    for idx in ranked_indices:
        if idx not in used_indices:
            return idx
    raise ValueError("Ran out of unique rows while selecting quantile anchors")


def choose_panel_rows(
    rows: Sequence[dict[str, Any]], quantiles: Sequence[float]
) -> list[dict[str, Any]]:
    required_sources = ("peak", "nonpeak")
    grouped: dict[str, list[dict[str, Any]]] = {source: [] for source in required_sources}
    for row in rows:
        source = str(row.get("source", "")).strip()
        if source in grouped:
            grouped[source].append(dict(row))

    panel: list[dict[str, Any]] = []
    for source in required_sources:
        source_rows = grouped[source]
        if len(source_rows) < len(quantiles):
            raise ValueError(
                f"Need at least {len(quantiles)} rows for source={source!r}, got {len(source_rows)}"
            )
        ranked_rows = sorted(
            source_rows,
            key=lambda row: (
                float(row["observed_total"]),
                str(row.get("chrom", "")),
                int(row.get("center", 0)),
                str(row.get("label", "")),
            ),
        )
        used_indices: set[int] = set()
        for quantile in quantiles:
            selected_idx = _pick_quantile_index(ranked_rows, quantile, used_indices)
            used_indices.add(selected_idx)
            selected = dict(ranked_rows[selected_idx])
            quantile_text = _quantile_string(quantile)
            selected["quantile"] = quantile_text
            selected["label"] = f"{source}_q{_quantile_token(quantile)}"
            selected["observed_total"] = float(selected["observed_total"])
            selected["notes"] = (
                f"{source} quantile anchor at q={quantile_text} for {WINDOW_BP}bp observed-total window"
            )
            panel.append(selected)

    labels = [row["label"] for row in panel]
    if len(panel) != len(required_sources) * len(quantiles):
        raise ValueError("Panel row count does not match expected peak/nonpeak quantile anchors")
    if len(set(labels)) != len(labels):
        raise ValueError("Panel labels must be unique")
    if {row["source"] for row in panel} != set(required_sources):
        raise ValueError("Panel must include both peak and nonpeak rows")
    return panel


def _load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"Expected mapping in data config: {path}")
    return payload


def _build_observed_rows(data_config: dict[str, Any], split: str) -> list[dict[str, Any]]:
    (
        load_fold_chroms,
        load_regions_from_bed,
        load_bigwig_chrom_sizes,
        filter_records_by_chroms,
    ) = _import_real_data_helpers()

    input_cfg = data_config.get("input")
    if not isinstance(input_cfg, dict):
        raise ValueError("data config is missing `input` mapping")

    folds_json = str(data_config["folds_json"])
    peaks_bed = str(input_cfg["peaks_bed"])
    nonpeaks_bed = str(input_cfg["nonpeaks_bed"])
    bigwig_path = str(input_cfg["bigwig"])

    split_chroms = load_fold_chroms(folds_json, split)
    bigwig_chrom_sizes = load_bigwig_chrom_sizes(bigwig_path)
    peaks = load_regions_from_bed(peaks_bed, split_chroms, source="peak")
    peaks = filter_records_by_chroms(
        peaks,
        bigwig_chrom_sizes.keys(),
        peaks_bed,
        source="peak",
    )
    nonpeaks = load_regions_from_bed(nonpeaks_bed, split_chroms, source="nonpeak")
    nonpeaks = filter_records_by_chroms(
        nonpeaks,
        bigwig_chrom_sizes.keys(),
        nonpeaks_bed,
        source="nonpeak",
    )

    try:
        import pyBigWig
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "pyBigWig is required to build the matched panel; use the TransChromBP runtime environment."
        ) from exc

    rows: list[dict[str, Any]] = []
    bigwig = pyBigWig.open(bigwig_path)
    if bigwig is None:
        raise ValueError(f"Failed to open bigWig: {bigwig_path}")
    try:
        for record in list(peaks) + list(nonpeaks):
            observed_total = _window_total(bigwig, record.chrom, record.center)
            rows.append(
                {
                    "label": f"{record.source}_{record.chrom}_{record.center}",
                    "source": record.source,
                    "chrom": record.chrom,
                    "center": int(record.center),
                    "observed_total": observed_total,
                    "notes": "",
                }
            )
    finally:
        bigwig.close()

    return rows


def _write_rows(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            payload = {column: row[column] for column in OUTPUT_COLUMNS}
            payload["observed_total"] = f"{float(payload['observed_total']):.6f}"
            writer.writerow(payload)


def main() -> None:
    args = parse_args()
    data_config = _load_yaml(args.data_config)
    rows = _build_observed_rows(data_config, split=args.split)
    panel_rows = choose_panel_rows(rows, PANEL_QUANTILES)
    _write_rows(args.output_csv, panel_rows)


if __name__ == "__main__":
    main()
