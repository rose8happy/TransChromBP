#!/usr/bin/env python3
"""Generate fold_*.json files for paper-aligned ChromBPNet reproduction."""

import argparse
import json
from pathlib import Path


def parse_pair(spec: str):
    """
    Parse a fold spec in the format:
      test_chr[,test_chr...]:valid_chr[,valid_chr...]
    """
    if ":" not in spec:
        raise ValueError(f"Invalid pair spec '{spec}', expected TEST:VALID format")
    test_part, valid_part = spec.split(":", 1)
    test = [x.strip() for x in test_part.split(",") if x.strip()]
    valid = [x.strip() for x in valid_part.split(",") if x.strip()]
    if not test or not valid:
        raise ValueError(f"Invalid pair spec '{spec}', test/valid must be non-empty")
    return test, valid


def read_chroms(chrom_sizes: Path):
    chroms = []
    with chrom_sizes.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            chroms.append(line.split("\t", 1)[0])
    if not chroms:
        raise ValueError(f"No chromosomes found in {chrom_sizes}")
    return chroms


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chrom-sizes", required=True, help="Path to chrom sizes TSV")
    parser.add_argument("--output-dir", required=True, help="Directory for fold_*.json")
    parser.add_argument(
        "--pairs",
        nargs="+",
        default=[
            "chr1:chr2",
            "chr3:chr4",
            "chr5:chr6",
            "chr7:chr8",
            "chr9:chr10",
        ],
        help=(
            "Fold specs in TEST:VALID format. "
            "Each side can be comma-separated (e.g. chr1,chr11:chr2,chr12)."
        ),
    )
    args = parser.parse_args()

    chrom_sizes = Path(args.chrom_sizes).resolve()
    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    chroms = read_chroms(chrom_sizes)
    chrom_set = set(chroms)
    used_test = set()
    used_valid = set()
    manifest = []

    for idx, spec in enumerate(args.pairs):
        test, valid = parse_pair(spec)
        test_set = set(test)
        valid_set = set(valid)

        missing = (test_set | valid_set) - chrom_set
        if missing:
            raise ValueError(f"Fold {idx} has chromosomes not in chrom sizes: {sorted(missing)}")
        if test_set & valid_set:
            raise ValueError(f"Fold {idx} has overlapping test/valid chromosomes: {sorted(test_set & valid_set)}")
        if used_test & test_set:
            raise ValueError(f"Test chromosome reused across folds: {sorted(used_test & test_set)}")
        if used_valid & valid_set:
            raise ValueError(f"Valid chromosome reused across folds: {sorted(used_valid & valid_set)}")

        used_test |= test_set
        used_valid |= valid_set

        train = [c for c in chroms if c not in test_set and c not in valid_set]
        fold = {"test": test, "valid": valid, "train": train}

        fold_path = out_dir / f"fold_{idx}.json"
        with fold_path.open("w", encoding="utf-8") as f:
            json.dump(fold, f, indent=4)

        manifest.append(
            {
                "fold": idx,
                "path": str(fold_path),
                "test": test,
                "valid": valid,
                "n_train_chroms": len(train),
            }
        )

    with (out_dir / "manifest.json").open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=4)

    print(f"Wrote {len(manifest)} fold files to {out_dir}")
    for item in manifest:
        print(f"fold_{item['fold']}: test={item['test']} valid={item['valid']} train_chroms={item['n_train_chroms']}")


if __name__ == "__main__":
    main()

