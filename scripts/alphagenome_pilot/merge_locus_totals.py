#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge AlphaGenome pilot totals with existing locus totals."
    )
    parser.add_argument("--alpha-summary", type=Path, required=True)
    parser.add_argument("--local-totals", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    alpha = pd.read_csv(args.alpha_summary)
    local = pd.read_csv(args.local_totals)
    merged = local.merge(
        alpha[["label", "source", "alpha_total_count", "alpha_max_value"]],
        on=["label", "source"],
        how="left",
    )
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output_csv, index=False)
    print(f"Wrote merged comparison: {args.output_csv}")


if __name__ == "__main__":
    main()
