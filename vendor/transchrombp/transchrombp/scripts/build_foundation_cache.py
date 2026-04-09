#!/usr/bin/env python3
"""Build dataset-aligned cached foundation-model features.

Current backend support:
  - nt_v2: Nucleotide Transformer v2 (HF local directory)

The cache layout intentionally mirrors the Genos cached-summary path:
  - manifest_<split>.json
  - <split>_<feature>.f16.npy / .f32.npy
Optional probe metadata:
  - records_<split>.jsonl
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import torch
import yaml


_HERE = Path(__file__).resolve()
_IMPORT_ROOTS = [
    _HERE.parents[1] / "src",
    _HERE.parents[2],
]
for _root in _IMPORT_ROOTS:
    if (_root / "transchrombp").exists():
        if str(_root) not in sys.path:
            sys.path.insert(0, str(_root))
        break


FEATURE_BIN_COUNTS = {
    "global_mean": 0,
    "bins4_mean": 4,
    "bins16_mean": 16,
}


def load_yaml(path: str) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def compute_record_sha1(records) -> str:
    h = hashlib.sha1()
    for record in records:
        h.update(f"{record.chrom}:{record.center}:{record.source}\n".encode("utf-8"))
    return h.hexdigest()


def resolve_data_value(train_data_cfg: dict[str, Any], data_source_cfg: dict[str, Any], key: str, default: Any = None) -> Any:
    if key in train_data_cfg and train_data_cfg[key] not in ("", None):
        return train_data_cfg[key]
    if key in data_source_cfg and data_source_cfg[key] not in ("", None):
        return data_source_cfg[key]
    for section in ("input", "window", "sampling"):
        section_cfg = data_source_cfg.get(section, {})
        if key in section_cfg and section_cfg[key] not in ("", None):
            return section_cfg[key]
    return default


def normalize_split_name(split: str) -> str:
    normalized = str(split).strip().lower()
    return "valid" if normalized == "val" else normalized


def resolve_split_region_source(data_cfg: dict[str, Any], split: str) -> str:
    default_region_source = str(data_cfg.get("region_source", "both"))
    normalized = normalize_split_name(split)
    if normalized == "train":
        return str(data_cfg.get("train_region_source", default_region_source))
    if normalized == "valid":
        return str(
            data_cfg.get(
                "valid_region_source",
                data_cfg.get("val_region_source", default_region_source),
            )
        )
    return str(
        data_cfg.get(
            f"{normalized}_region_source",
            data_cfg.get("val_region_source", default_region_source),
        )
    )


def resolve_split_max_records(data_cfg: dict[str, Any], split: str) -> int:
    normalized = normalize_split_name(split)
    candidate_keys: list[str]
    if normalized == "train":
        candidate_keys = ["max_train_regions"]
    elif normalized == "valid":
        candidate_keys = ["max_valid_regions", "max_val_regions"]
    else:
        candidate_keys = [f"max_{normalized}_regions"]
    candidate_keys.append("max_records")

    for key in candidate_keys:
        value = data_cfg.get(key, None)
        if value not in (None, ""):
            return int(value)
    return 0


def build_dataset_for_cache(split: str, train_config: dict[str, Any], data_source_config: dict[str, Any]):
    from transchrombp.data import ChromBPNetBigWigDataset, resolve_dataset_seed

    data_cfg = train_config.get("data", {})
    seed = int(train_config.get("seed", 1234))
    normalized_split = normalize_split_name(split)

    input_len = int(data_cfg.get("input_len", 2114))
    supervised_bp = int(data_cfg.get("supervised_bp", 1000))
    profile_bin_size = int(data_cfg.get("profile_bin_size", 1))
    max_records = resolve_split_max_records(data_cfg, normalized_split)
    nonpeak_ratio = float(resolve_data_value(data_cfg, data_source_config, "nonpeak_ratio", 1.0))
    region_source = resolve_split_region_source(data_cfg, normalized_split)

    ds = ChromBPNetBigWigDataset(
        genome_fasta=str(resolve_data_value(data_cfg, data_source_config, "genome_fasta", "")),
        bigwig_path=str(resolve_data_value(data_cfg, data_source_config, "bigwig", "")),
        peaks_bed=str(resolve_data_value(data_cfg, data_source_config, "peaks_bed", "")),
        nonpeaks_bed=str(resolve_data_value(data_cfg, data_source_config, "nonpeaks_bed", "")),
        folds_json=str(resolve_data_value(data_cfg, data_source_config, "folds_json", "")),
        split=normalized_split,
        input_len=input_len,
        supervised_bp=supervised_bp,
        profile_bin_size=profile_bin_size,
        max_jitter=0,
        peak_max_jitter=0,
        nonpeak_max_jitter=0,
        seed=resolve_dataset_seed(seed, normalized_split),
        nonpeak_ratio=nonpeak_ratio,
        max_records=max_records,
        region_source=region_source,
        random_revcomp=False,
        revcomp_prob=0.0,
        foundation_cache_dir="",
        foundation_cache_features=(),
        teacher_cache_dir="",
        teacher_target_names=(),
    )
    return ds


def resolve_layers(layer_spec: str, num_hidden_layers: int) -> list[int]:
    if layer_spec == "auto_quartiles":
        raw = [
            max(1, round(num_hidden_layers * 0.25)),
            max(1, round(num_hidden_layers * 0.50)),
            max(1, round(num_hidden_layers * 0.75)),
            num_hidden_layers,
        ]
        ordered: list[int] = []
        for value in raw:
            if value not in ordered:
                ordered.append(value)
        return ordered
    layers = [int(part.strip()) for part in layer_spec.split(",") if part.strip()]
    if not layers:
        raise ValueError(f"Invalid --layers: {layer_spec!r}")
    return layers


def feature_key(layer: int, feature_type: str) -> str:
    return f"layer_{layer:02d}__{feature_type}"


def fetch_sequence_string(ds, record) -> str:
    start = int(record.center) - (ds.input_len // 2)
    end = start + ds.input_len
    clipped_start = max(0, start)
    clipped_end = min(end, ds._chrom_len(record.chrom))
    prefix = "N" * max(0, clipped_start - start)
    suffix = "N" * max(0, end - clipped_end)
    seq = ""
    if clipped_start < clipped_end:
        seq = ds._get_fasta()[record.chrom][clipped_start:clipped_end]
        if not isinstance(seq, str):
            seq = str(seq)
    return (prefix + seq + suffix).upper()


def sum_into_bins(values: np.ndarray, n_bins: int) -> np.ndarray:
    chunks = np.array_split(np.asarray(values, dtype=np.float32), int(n_bins))
    return np.asarray([float(chunk.sum(dtype=np.float32)) for chunk in chunks], dtype=np.float32)


def masked_global_mean(hidden: torch.Tensor, mask: torch.Tensor) -> np.ndarray:
    masked = hidden * mask.unsqueeze(-1)
    denom = mask.sum(dim=1, keepdim=True).clamp_min(1)
    return (masked.sum(dim=1) / denom).detach().cpu().numpy()


def masked_bins_mean(hidden: torch.Tensor, mask: torch.Tensor, bin_count: int) -> np.ndarray:
    outputs: list[np.ndarray] = []
    hidden_np = hidden.detach().cpu().numpy()
    mask_np = mask.detach().cpu().numpy()
    for sample_hidden, sample_mask in zip(hidden_np, mask_np):
        valid_len = int(sample_mask.sum())
        valid_hidden = sample_hidden[:valid_len]
        if valid_len == 0:
            outputs.append(np.zeros((bin_count, sample_hidden.shape[-1]), dtype=np.float32))
            continue
        bins: list[np.ndarray] = []
        for index in range(bin_count):
            start = math.floor(index * valid_len / bin_count)
            end = math.floor((index + 1) * valid_len / bin_count)
            chunk = valid_hidden[start:end]
            if len(chunk) == 0:
                chunk = valid_hidden[max(0, min(start, valid_len - 1)) : max(1, min(start + 1, valid_len))]
            bins.append(chunk.mean(axis=0, dtype=np.float32))
        outputs.append(np.stack(bins, axis=0))
    return np.stack(outputs, axis=0).astype(np.float32)


def build_nt_v2_runtime(model_dir: str, device: str):
    from transformers import AutoModelForMaskedLM, AutoTokenizer

    tokenizer = AutoTokenizer.from_pretrained(model_dir, trust_remote_code=True, local_files_only=True)
    model = AutoModelForMaskedLM.from_pretrained(model_dir, trust_remote_code=True, local_files_only=True)
    model = model.to(device).eval()
    return tokenizer, model


def estimate_feature_shape(n_records: int, hidden_size: int, feature_type: str) -> tuple[int, ...]:
    bins = FEATURE_BIN_COUNTS[feature_type]
    if bins <= 0:
        return (n_records, hidden_size)
    return (n_records, bins, hidden_size)


def shard_bounds(n_records: int, shard_index: int, num_shards: int) -> tuple[int, int]:
    start = (n_records * shard_index) // num_shards
    end = (n_records * (shard_index + 1)) // num_shards
    return start, end


def maybe_init_dist(device_arg: str) -> tuple[int, int, int, str]:
    world_size = int(os.environ.get("WORLD_SIZE", "1"))
    rank = int(os.environ.get("RANK", "0"))
    local_rank = int(os.environ.get("LOCAL_RANK", str(rank)))
    device = device_arg

    if world_size <= 1:
        return rank, world_size, local_rank, device

    import torch.distributed as dist

    if device_arg.startswith("cuda"):
        torch.cuda.set_device(local_rank)
        device = f"cuda:{local_rank}"
        backend = "nccl"
    else:
        backend = "gloo"
    if not dist.is_initialized():
        dist.init_process_group(backend=backend, init_method="env://")
    return rank, world_size, local_rank, device


def maybe_barrier(world_size: int) -> None:
    if world_size <= 1:
        return
    import torch.distributed as dist

    dist.barrier()


def maybe_destroy_dist(world_size: int) -> None:
    if world_size <= 1:
        return
    import torch.distributed as dist

    if dist.is_initialized():
        dist.destroy_process_group()


def cleanup_paths(paths: list[Path]) -> None:
    for path in paths:
        try:
            if path.exists():
                path.unlink()
        except FileNotFoundError:
            pass


def merge_sharded_outputs(
    *,
    output_dir: Path,
    shard_dir: Path,
    split: str,
    rank_count: int,
    n_records: int,
    hidden_size: int,
    layers: list[int],
    feature_types: list[str],
    dtype: np.dtype,
    suffix: str,
    write_records: bool,
    profile_probe_bins: int,
    manifest_base: dict[str, Any],
) -> None:
    shard_meta_paths = [shard_dir / f"{split}.shard{rank:02d}.meta.json" for rank in range(rank_count)]
    shard_metas = [json.loads(path.read_text(encoding="utf-8")) for path in shard_meta_paths]
    shard_metas.sort(key=lambda item: int(item["shard_index"]))

    feature_keys = [feature_key(layer, feature_type) for layer in layers for feature_type in feature_types]
    for key in feature_keys:
        feature_type = key.split("__", 1)[1]
        final_shape = estimate_feature_shape(n_records, hidden_size, feature_type)
        final_path = output_dir / f"{split}_{key}.{suffix}.npy"
        final_arr = np.lib.format.open_memmap(str(final_path), mode="w+", dtype=dtype, shape=final_shape)
        for meta in shard_metas:
            shard_path = Path(meta["feature_files"][key])
            shard_arr = np.load(shard_path, mmap_mode="r")
            final_arr[int(meta["start_idx"]) : int(meta["end_idx"])] = shard_arr
        del final_arr

    records_jsonl_path = output_dir / f"records_{split}.jsonl"
    if write_records:
        with records_jsonl_path.open("w", encoding="utf-8") as out_handle:
            for meta in shard_metas:
                shard_records_path = meta.get("records_jsonl", "")
                if not shard_records_path:
                    continue
                with Path(shard_records_path).open("r", encoding="utf-8") as in_handle:
                    for line in in_handle:
                        out_handle.write(line)

    token_mins = [int(meta["token_length_min"]) for meta in shard_metas if int(meta["shard_rows"]) > 0]
    token_maxs = [int(meta["token_length_max"]) for meta in shard_metas if int(meta["shard_rows"]) > 0]
    weighted_token_sum = sum(float(meta["token_length_sum"]) for meta in shard_metas)
    total_rows = sum(int(meta["shard_rows"]) for meta in shard_metas)

    manifest = dict(manifest_base)
    manifest["features"] = feature_keys
    manifest["feature_types"] = feature_types
    manifest["hidden_size"] = hidden_size
    manifest["dtype"] = "float16" if dtype == np.float16 else "float32"
    manifest["token_length_min"] = int(min(token_mins)) if token_mins else 0
    manifest["token_length_max"] = int(max(token_maxs)) if token_maxs else 0
    manifest["token_length_mean"] = float(weighted_token_sum / max(total_rows, 1))
    manifest["elapsed_sec_total"] = float(max(float(meta["elapsed_sec_total"]) for meta in shard_metas))
    manifest["worker_elapsed_sec_sum"] = float(sum(float(meta["elapsed_sec_total"]) for meta in shard_metas))
    manifest["world_size"] = int(rank_count)
    if write_records:
        manifest["records_jsonl"] = str(records_jsonl_path)
        manifest["profile_probe_bins"] = int(profile_probe_bins)

    manifest_path = output_dir / f"manifest_{split}.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    extra_cleanup: list[Path] = []
    for meta in shard_metas:
        extra_cleanup.extend(Path(path_str) for path_str in meta["feature_files"].values())
        if meta.get("records_jsonl"):
            extra_cleanup.append(Path(meta["records_jsonl"]))
    extra_cleanup.extend(shard_meta_paths)
    cleanup_paths(extra_cleanup)
    try:
        shard_dir.rmdir()
    except OSError:
        pass


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data_config", required=True)
    parser.add_argument("--train_config", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--splits", nargs="+", default=["train", "valid"])
    parser.add_argument("--backend", default="nt_v2", choices=["nt_v2"])
    parser.add_argument("--model_dir", required=True)
    parser.add_argument("--layers", default="auto_quartiles")
    parser.add_argument("--feature_types", nargs="+", default=["global_mean", "bins4_mean"])
    parser.add_argument("--batch_size", type=int, default=8)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--dtype", default="float16", choices=["float16", "float32"])
    parser.add_argument("--record_splits", nargs="*", default=["valid"])
    parser.add_argument("--profile_probe_bins", type=int, default=16)
    parser.add_argument("--dry_run", action="store_true")
    args = parser.parse_args()

    feature_types = [str(item) for item in args.feature_types]
    unknown = [item for item in feature_types if item not in FEATURE_BIN_COUNTS]
    if unknown:
        raise ValueError(f"Unsupported feature_types: {unknown!r}")

    train_config = load_yaml(args.train_config)
    data_source_config = load_yaml(args.data_config)
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    rank, world_size, local_rank, resolved_device = maybe_init_dist(args.device)

    if args.backend != "nt_v2":
        raise ValueError(f"Unsupported backend={args.backend!r}")

    try:
        tokenizer, model = build_nt_v2_runtime(args.model_dir, resolved_device)
        num_hidden_layers = int(getattr(model.config, "num_hidden_layers", 0))
        hidden_size = int(getattr(model.config, "hidden_size", 0))
        if num_hidden_layers <= 0 or hidden_size <= 0:
            raise ValueError(
                f"Unable to infer model depth/width from {args.model_dir}: "
                f"num_hidden_layers={num_hidden_layers}, hidden_size={hidden_size}"
            )
        layers = resolve_layers(args.layers, num_hidden_layers)
        dtype = np.float16 if args.dtype == "float16" else np.float32
        suffix = "f16" if args.dtype == "float16" else "f32"

        for split in args.splits:
            if rank == 0:
                print(f"\n{'=' * 60}")
                print(f"Processing split: {split}")
                print(f"{'=' * 60}")

            ds = build_dataset_for_cache(split, train_config, data_source_config)
            n_records = len(ds)
            if getattr(ds, "supports_epoch_resampling", False):
                epoch_regions = int(getattr(ds, "peak_record_count", 0)) + min(
                    int(getattr(ds, "nonpeak_record_count", 0)),
                    int(getattr(ds, "target_nonpeak_count", 0)),
                )
            else:
                epoch_regions = n_records
            record_sha1 = compute_record_sha1(ds.records)
            write_records = split in set(args.record_splits)

            if rank == 0:
                print(f"  n_records = {n_records}")
                if split == "train":
                    print(f"  train_epoch_regions ≈ {epoch_regions}")
                print(f"  record_sha1 = {record_sha1}")
                print(f"  layers = {layers}")
                print(f"  feature_types = {feature_types}")

                est_total_bytes = 0
                for layer in layers:
                    for feature_type in feature_types:
                        shape = estimate_feature_shape(n_records, hidden_size, feature_type)
                        bytes_per = np.dtype(dtype).itemsize
                        est_bytes = int(np.prod(shape)) * bytes_per
                        est_total_bytes += est_bytes
                        print(
                            f"  estimate {feature_key(layer, feature_type)}: "
                            f"shape={shape} size={est_bytes / (1024 ** 3):.3f} GiB"
                        )
                stat = os.statvfs(str(output_dir))
                free_gib = (stat.f_bavail * stat.f_frsize) / (1024 ** 3)
                print(f"  estimated total size: {est_total_bytes / (1024 ** 3):.3f} GiB")
                print(f"  output_dir free space: {free_gib:.1f} GiB")
            maybe_barrier(world_size)

            if args.dry_run:
                if rank == 0:
                    print("  [dry_run] Skipping extraction.")
                continue

            shard_start, shard_end = shard_bounds(n_records, rank, world_size)
            shard_records = ds.records[shard_start:shard_end]
            shard_rows = shard_end - shard_start
            shard_dir = output_dir / f".cache_shards_{split}"
            if rank == 0:
                shard_dir.mkdir(parents=True, exist_ok=True)
            maybe_barrier(world_size)

            feature_files: dict[str, str] = {}
            memmaps: dict[str, np.memmap] = {}
            for layer in layers:
                for feature_type in feature_types:
                    key = feature_key(layer, feature_type)
                    shape = estimate_feature_shape(shard_rows, hidden_size, feature_type)
                    path = shard_dir / f"{split}.shard{rank:02d}_{key}.{suffix}.npy"
                    feature_files[key] = str(path)
                    memmaps[key] = np.lib.format.open_memmap(
                        str(path),
                        mode="w+",
                        dtype=dtype,
                        shape=shape,
                    )

            records_jsonl_path = shard_dir / f"records_{split}.shard{rank:02d}.jsonl"
            record_handle = None
            if write_records:
                record_handle = records_jsonl_path.open("w", encoding="utf-8")

            token_lengths: list[int] = []
            t0 = time.time()
            batch_size = int(args.batch_size)
            n_batches = (shard_rows + batch_size - 1) // batch_size if shard_rows > 0 else 0

            try:
                for batch_idx in range(n_batches):
                    local_start = batch_idx * batch_size
                    local_end = min(local_start + batch_size, shard_rows)
                    global_start = shard_start + local_start
                    global_end = shard_start + local_end
                    batch_records = shard_records[local_start:local_end]
                    sequences = [fetch_sequence_string(ds, record) for record in batch_records]

                    batch = tokenizer.batch_encode_plus(sequences, return_tensors="pt", padding=True, truncation=True)
                    mask = batch["attention_mask"].bool().to(resolved_device)
                    batch = {name: tensor.to(resolved_device) for name, tensor in batch.items()}

                    with torch.no_grad():
                        outputs = model(
                            batch["input_ids"],
                            attention_mask=mask,
                            encoder_attention_mask=mask,
                            output_hidden_states=True,
                        )
                    token_lengths.extend(mask.sum(dim=1).detach().cpu().tolist())

                    for layer in layers:
                        hidden = outputs.hidden_states[layer]
                        for feature_type in feature_types:
                            key = feature_key(layer, feature_type)
                            if feature_type == "global_mean":
                                values = masked_global_mean(hidden, mask)
                            else:
                                values = masked_bins_mean(hidden, mask, FEATURE_BIN_COUNTS[feature_type])
                            memmaps[key][local_start:local_end] = values.astype(dtype, copy=False)

                    if record_handle is not None:
                        for offset, (record, sequence) in enumerate(zip(batch_records, sequences)):
                            profile_start = int(record.center) - (ds.supervised_bp // 2)
                            profile = ds._fetch_profile(record.chrom, profile_start, ds.supervised_bp)
                            profile16 = sum_into_bins(profile, args.profile_probe_bins)
                            row = {
                                "index": global_start + offset,
                                "chrom": record.chrom,
                                "center": int(record.center),
                                "source": record.source,
                                "sequence": sequence,
                                "true_logcount": float(np.log1p(profile.sum(dtype=np.float32))),
                                "true_profile16": [float(x) for x in profile16.tolist()],
                            }
                            record_handle.write(json.dumps(row, ensure_ascii=True) + "\n")

                    if (batch_idx + 1) % 50 == 0 or batch_idx == n_batches - 1:
                        elapsed = time.time() - t0
                        processed = global_end - shard_start
                        rate = processed / max(elapsed, 1e-6)
                        eta = (shard_rows - processed) / max(rate, 1e-6)
                        print(
                            f"  [{split}][rank {rank}/{world_size}] "
                            f"{global_end}/{n_records} ({elapsed:.0f}s elapsed, ETA {eta:.0f}s)"
                        )
            finally:
                if record_handle is not None:
                    record_handle.close()
                for arr in memmaps.values():
                    del arr

            elapsed = time.time() - t0
            shard_meta = {
                "split": split,
                "shard_index": rank,
                "world_size": world_size,
                "start_idx": shard_start,
                "end_idx": shard_end,
                "shard_rows": shard_rows,
                "record_sha1": record_sha1,
                "feature_files": feature_files,
                "records_jsonl": str(records_jsonl_path) if write_records else "",
                "elapsed_sec_total": float(elapsed),
                "token_length_min": int(min(token_lengths)) if token_lengths else 0,
                "token_length_max": int(max(token_lengths)) if token_lengths else 0,
                "token_length_sum": float(sum(token_lengths)),
            }
            (shard_dir / f"{split}.shard{rank:02d}.meta.json").write_text(
                json.dumps(shard_meta, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            maybe_barrier(world_size)

            if rank == 0:
                manifest_base = {
                    "backend": args.backend,
                    "data_config_path": os.path.abspath(args.data_config),
                    "train_config_path": os.path.abspath(args.train_config),
                    "split": split,
                    "dataset_split": normalize_split_name(split),
                    "input_len": int(ds.input_len),
                    "supervised_bp": int(ds.supervised_bp),
                    "n_records": n_records,
                    "max_records": int(getattr(ds, "max_records", 0)),
                    "nonpeak_ratio": float(getattr(ds, "nonpeak_ratio", 0.0)),
                    "region_source": str(getattr(ds, "region_source", "")),
                    "model_dir": str(Path(args.model_dir).resolve()),
                    "layers": layers,
                    "record_sha1": record_sha1,
                }
                if split == "train":
                    manifest_base["train_epoch_regions"] = epoch_regions
                    manifest_base["supports_epoch_resampling"] = bool(getattr(ds, "supports_epoch_resampling", False))

                merge_sharded_outputs(
                    output_dir=output_dir,
                    shard_dir=shard_dir,
                    split=split,
                    rank_count=world_size,
                    n_records=n_records,
                    hidden_size=hidden_size,
                    layers=layers,
                    feature_types=feature_types,
                    dtype=dtype,
                    suffix=suffix,
                    write_records=write_records,
                    profile_probe_bins=args.profile_probe_bins,
                    manifest_base=manifest_base,
                )
                print(f"  [{split}] done with world_size={world_size}")
                print(f"  Manifest written: {output_dir / f'manifest_{split}.json'}")
            maybe_barrier(world_size)

        if rank == 0:
            print("\nAll done.")
    finally:
        maybe_destroy_dist(world_size)


if __name__ == "__main__":
    main()
