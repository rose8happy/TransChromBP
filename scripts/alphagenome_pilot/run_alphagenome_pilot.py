#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np


DEFAULT_API_ENV = Path("~/.config/alphagenome/api.env").expanduser()
DEFAULT_ONTOLOGY = "EFO:0002067"  # K562
DEFAULT_OUTPUT_TYPE = "ATAC"
DEFAULT_TARGET_WIDTH = 1000
API_KEY_ENV_NAMES = ("ALPHAGENOME_API_KEY", "ALPHA_GENOME_API_KEY")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a small AlphaGenome pilot on a list of loci."
    )
    parser.add_argument(
        "--regions-csv",
        type=Path,
        required=True,
        help="CSV with at least columns: label, chrom, center. Optional: source, notes.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory to write summary CSV, JSONL metadata, and per-locus NPZ profiles.",
    )
    parser.add_argument(
        "--api-env",
        type=Path,
        default=DEFAULT_API_ENV,
        help=(
            "Shell env file containing "
            "export ALPHAGENOME_API_KEY=... or export ALPHA_GENOME_API_KEY=..."
        ),
    )
    parser.add_argument(
        "--ontology-term",
        default=DEFAULT_ONTOLOGY,
        help="Ontology term to request/filter tracks to. Default is K562 EFO:0002067.",
    )
    parser.add_argument(
        "--output-type",
        default=DEFAULT_OUTPUT_TYPE,
        choices=["ATAC", "DNASE"],
        help="AlphaGenome output type to request.",
    )
    parser.add_argument(
        "--target-width",
        type=int,
        default=DEFAULT_TARGET_WIDTH,
        help="Central output width in bp to compare against local models.",
    )
    parser.add_argument(
        "--aggregate",
        choices=["mean", "sum"],
        default="mean",
        help="How to merge multiple AlphaGenome tracks after ontology filtering.",
    )
    parser.add_argument(
        "--keep-all-tracks",
        action="store_true",
        help="Also save the raw multi-track matrix into each NPZ output.",
    )
    return parser.parse_args()


def load_api_key(env_path: Path) -> str:
    for env_name in API_KEY_ENV_NAMES:
        direct = os.environ.get(env_name, "").strip()
        if direct:
            return direct
    if not env_path.exists():
        raise FileNotFoundError(f"API env file not found: {env_path}")

    pattern = re.compile(
        r"^\s*(?:export\s+)?(?:ALPHAGENOME_API_KEY|ALPHA_GENOME_API_KEY)=(.*)\s*$"
    )
    for line in env_path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if not match:
            continue
        value = match.group(1).strip()
        if (value.startswith("'") and value.endswith("'")) or (
            value.startswith('"') and value.endswith('"')
        ):
            value = value[1:-1]
        value = value.strip()
        if value:
            return value
    raise RuntimeError(
        f"None of {API_KEY_ENV_NAMES} found in {env_path}"
    )


def import_alphagenome_modules() -> tuple[Any, Any]:
    try:
        from alphagenome.data import genome
        from alphagenome.models import dna_client
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "alphagenome Python package is not installed in this interpreter. "
            "Use Python >=3.10 and install the official SDK first."
        ) from exc
    return genome, dna_client


def read_regions(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"label", "chrom", "center"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns in {path}: {sorted(missing)}")
        for row in reader:
            rows.append(
                {
                    "label": row["label"].strip(),
                    "source": row.get("source", "").strip(),
                    "chrom": row["chrom"].strip(),
                    "center": int(row["center"]),
                    "notes": row.get("notes", "").strip(),
                }
            )
    if not rows:
        raise ValueError(f"No loci found in {path}")
    return rows


def filter_track_data_to_ontology(track_data: Any, ontology_term: str) -> tuple[Any, dict[str, Any]]:
    metadata = track_data.metadata.copy()
    summary: dict[str, Any] = {
        "tracks_before_filter": int(len(metadata)),
        "tracks_after_ontology_filter": int(len(metadata)),
        "tracks_after_strand_filter": int(len(metadata)),
        "used_unstranded_only": False,
        "track_names": metadata.get("name", []).tolist() if "name" in metadata else [],
    }

    if "ontology_curie" in metadata.columns:
        mask = (metadata["ontology_curie"].fillna("") == ontology_term).to_numpy()
        if mask.any():
            track_data = track_data.filter_tracks(mask)
            metadata = track_data.metadata.copy()
            summary["tracks_after_ontology_filter"] = int(len(metadata))

    if "strand" in metadata.columns and (metadata["strand"] == ".").any():
        track_data = track_data.filter_to_unstranded()
        metadata = track_data.metadata.copy()
        summary["tracks_after_strand_filter"] = int(len(metadata))
        summary["used_unstranded_only"] = True

    summary["track_names"] = metadata.get("name", []).tolist() if "name" in metadata else []
    summary["ontology_terms"] = (
        metadata.get("ontology_curie", []).tolist() if "ontology_curie" in metadata else []
    )
    return track_data, summary


def aggregate_track_values(values: np.ndarray, aggregate: str) -> np.ndarray:
    values = np.asarray(values, dtype=np.float32)
    if values.ndim == 1:
        return values
    if values.ndim < 1:
        raise ValueError(f"Unexpected AlphaGenome value shape: {values.shape}")
    if aggregate == "sum":
        return values.sum(axis=-1)
    return values.mean(axis=-1)


def main() -> None:
    args = parse_args()
    genome, dna_client = import_alphagenome_modules()
    api_key = load_api_key(args.api_env)
    loci = read_regions(args.regions_csv)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    profiles_dir = args.output_dir / "profiles"
    profiles_dir.mkdir(parents=True, exist_ok=True)

    model = dna_client.create(api_key)
    output_enum = getattr(dna_client.OutputType, args.output_type)
    summary_rows: list[dict[str, Any]] = []
    metadata_path = args.output_dir / "region_metadata.jsonl"

    with metadata_path.open("w", encoding="utf-8") as metadata_handle:
        for row in loci:
            center = int(row["center"])
            interval = genome.Interval(
                chromosome=row["chrom"],
                start=center,
                end=center + 1,
            ).resize(dna_client.SEQUENCE_LENGTH_1MB)

            outputs = model.predict_interval(
                interval=interval,
                requested_outputs={output_enum},
                ontology_terms=[args.ontology_term],
            )

            track_data = getattr(outputs, args.output_type.lower())
            track_data, track_summary = filter_track_data_to_ontology(
                track_data, args.ontology_term
            )
            target_interval = genome.Interval(
                chromosome=row["chrom"],
                start=center - args.target_width // 2,
                end=center + args.target_width // 2,
            )
            sliced = track_data.slice_by_interval(target_interval, match_resolution=True)
            aggregate_profile = aggregate_track_values(sliced.values, args.aggregate)
            total_count = float(np.asarray(aggregate_profile, dtype=np.float64).sum())
            max_value = float(np.asarray(aggregate_profile).max())
            out_npz = profiles_dir / f"{row['label']}.npz"

            npz_payload: dict[str, Any] = {
                "aggregate_profile": np.asarray(aggregate_profile, dtype=np.float32),
                "track_names": np.asarray(track_summary["track_names"], dtype=object),
                "ontology_terms": np.asarray(track_summary["ontology_terms"], dtype=object),
            }
            if args.keep_all_tracks:
                npz_payload["raw_track_values"] = np.asarray(sliced.values)
            np.savez_compressed(out_npz, **npz_payload)

            metadata_record = {
                **row,
                "ontology_term": args.ontology_term,
                "output_type": args.output_type,
                "interval_start_1mb": int(interval.start),
                "interval_end_1mb": int(interval.end),
                "target_start": int(target_interval.start),
                "target_end": int(target_interval.end),
                "profile_npz": str(out_npz),
                **track_summary,
            }
            metadata_handle.write(json.dumps(metadata_record, ensure_ascii=False) + "\n")

            summary_rows.append(
                {
                    "label": row["label"],
                    "source": row["source"],
                    "chrom": row["chrom"],
                    "center": center,
                    "output_type": args.output_type,
                    "ontology_term": args.ontology_term,
                    "aggregate": args.aggregate,
                    "target_width": args.target_width,
                    "num_tracks_used": track_summary["tracks_after_strand_filter"],
                    "alpha_total_count": total_count,
                    "alpha_max_value": max_value,
                    "profile_npz": str(out_npz),
                }
            )

    summary_path = args.output_dir / "summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    run_meta = {
        "regions_csv": str(args.regions_csv),
        "num_regions": len(loci),
        "output_type": args.output_type,
        "ontology_term": args.ontology_term,
        "target_width": args.target_width,
        "aggregate": args.aggregate,
        "summary_csv": str(summary_path),
        "region_metadata_jsonl": str(metadata_path),
    }
    (args.output_dir / "run_meta.json").write_text(
        json.dumps(run_meta, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(f"Wrote summary: {summary_path}")
    print(f"Wrote metadata: {metadata_path}")


if __name__ == "__main__":
    main()
