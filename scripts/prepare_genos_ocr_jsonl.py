#!/usr/bin/env python3
"""Convert a genomic-benchmarks full-seq dataset into Genos benchmark JSONL."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


VALID_BASES = set("ACGTN")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert a full-seq dataset laid out as split/class/*.txt "
            "into Genos benchmark JSONL files."
        )
    )
    parser.add_argument(
        "--source-dir",
        required=True,
        help="Path to the full-seq dataset root, e.g. .../human_ocr_ensembl",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        help="Output directory that will contain train.jsonl/test.jsonl/label_map.json",
    )
    parser.add_argument(
        "--positive-class",
        default=None,
        help=(
            "Optional class directory name to force as label 1 for a binary task. "
            "If omitted, labels follow sorted class-directory order."
        ),
    )
    return parser.parse_args()


def read_sequence(path: Path) -> str:
    seq = "".join(path.read_text(encoding="utf-8").split()).upper()
    invalid = sorted(set(seq) - VALID_BASES)
    if invalid:
        raise ValueError(f"{path} contains invalid bases: {invalid}")
    return seq


def list_class_dirs(split_dir: Path) -> list[Path]:
    class_dirs = sorted(p for p in split_dir.iterdir() if p.is_dir())
    if not class_dirs:
        raise ValueError(f"No class directories found in {split_dir}")
    return class_dirs


def build_label_map(train_dir: Path, test_dir: Path, positive_class: str | None) -> dict[str, int]:
    train_classes = [p.name for p in list_class_dirs(train_dir)]
    test_classes = [p.name for p in list_class_dirs(test_dir)]
    if train_classes != test_classes:
        raise ValueError(
            "Train/test class directories do not match: "
            f"train={train_classes}, test={test_classes}"
        )
    if positive_class is not None:
        if len(train_classes) != 2:
            raise ValueError("--positive-class only supports binary datasets")
        if positive_class not in train_classes:
            raise ValueError(
                f"--positive-class={positive_class!r} not found in classes {train_classes}"
            )
        negative_class = next(name for name in train_classes if name != positive_class)
        return {negative_class: 0, positive_class: 1}
    return {class_name: idx for idx, class_name in enumerate(train_classes)}


def write_split_jsonl(split_dir: Path, out_path: Path, label_map: dict[str, int]) -> dict[str, int]:
    counts = {class_name: 0 for class_name in label_map}
    total = 0

    with out_path.open("w", encoding="utf-8") as out_f:
        for class_name, label in label_map.items():
            class_dir = split_dir / class_name
            files = sorted(p for p in class_dir.iterdir() if p.is_file())
            if not files:
                raise ValueError(f"No sequence files found in {class_dir}")
            for seq_path in files:
                row = {
                    "seq": read_sequence(seq_path),
                    "label": label,
                }
                out_f.write(json.dumps(row, ensure_ascii=True) + "\n")
                counts[class_name] += 1
                total += 1

    counts["total"] = total
    return counts


def main() -> None:
    args = parse_args()
    source_dir = Path(args.source_dir).resolve()
    out_dir = Path(args.out_dir).resolve()

    train_dir = source_dir / "train"
    test_dir = source_dir / "test"
    if not train_dir.is_dir() or not test_dir.is_dir():
        raise FileNotFoundError(
            f"Expected train/ and test/ under {source_dir}, got train={train_dir.exists()} test={test_dir.exists()}"
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    label_map = build_label_map(train_dir, test_dir, args.positive_class)

    train_counts = write_split_jsonl(train_dir, out_dir / "train.jsonl", label_map)
    test_counts = write_split_jsonl(test_dir, out_dir / "test.jsonl", label_map)

    summary = {
        "source_dir": str(source_dir),
        "label_map": label_map,
        "positive_class": args.positive_class,
        "train_counts": train_counts,
        "test_counts": test_counts,
    }

    (out_dir / "label_map.json").write_text(
        json.dumps(label_map, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
