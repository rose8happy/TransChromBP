#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any, Sequence

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a deterministic matched AlphaGenome tutorial panel."
    )
    parser.add_argument("--data-config", type=Path, required=True)
    parser.add_argument("--split", default="test")
    parser.add_argument("--output-csv", type=Path, required=True)
    return parser.parse_args()


def _import_real_data_helpers() -> tuple[Any, Any, Any, Any]:
    from transchrombp.data.real_data import (
        filter_records_by_chroms,
        load_bigwig_chrom_sizes,
        load_fold_chroms,
        load_regions_from_bed,
    )

    return (
        load_fold_chroms,
        load_regions_from_bed,
        load_bigwig_chrom_sizes,
        filter_records_by_chroms,
    )


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
