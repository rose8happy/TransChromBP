"""Export record-aligned teacher caches directly from a model checkpoint."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import torch
from torch.utils.data import DataLoader

from transchrombp.data import ChromBPNetBigWigDataset, resolve_dataset_seed
from transchrombp.data.real_data import compute_record_sha1
from transchrombp.evaluation.evaluate_checkpoint import load_model_state_with_count_pool_compat
from transchrombp.training.train_ddp import (
    _pool_profile_probs_to_bins as training_pool_profile_probs_to_bins,
    _resolve_real_data_value,
    apply_external_data_config_defaults,
    apply_training_mode_defaults,
    build_model,
    build_online_feature_runtime,
    check_profile_resolution_compatibility,
    extract_foundation_feature_kwargs,
    extract_genos_summary_kwargs,
    extract_online_feature_kwargs,
    load_yaml,
    resolve_data_config_path,
    validate_model_config,
    validate_semantics_profile,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SPLITS: tuple[str, ...] = ("train", "valid")
DEFAULT_TARGETS: tuple[str, ...] = ("profile16", "logcount")


def _resolve_path(path_value: str, base_dir: Path) -> str:
    path = Path(path_value)
    if path.is_absolute():
        return str(path)
    return str((base_dir / path).resolve())


def _get_device(device_name: str) -> torch.device:
    if device_name == "auto":
        return torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    return torch.device(device_name)


def _to_numpy_f32(array: Any) -> np.ndarray:
    if torch.is_tensor(array):
        return array.detach().cpu().numpy().astype(np.float32, copy=False)
    return np.asarray(array, dtype=np.float32)


def _resolve_region_source(data_cfg: Mapping[str, Any], split: str) -> str:
    default_region_source = str(data_cfg.get("region_source", "both"))
    if split == "train":
        return str(data_cfg.get("train_region_source", default_region_source))
    return str(data_cfg.get(f"{split}_region_source", data_cfg.get("val_region_source", default_region_source)))


def _build_export_dataset(
    train_cfg: Mapping[str, Any],
    data_source_cfg: Mapping[str, Any],
    split: str,
    max_regions_override: int,
) -> ChromBPNetBigWigDataset:
    data_cfg = dict(train_cfg.get("data", {}))
    input_len = int(data_cfg.get("input_len", 2114))
    profile_bin_size = int(data_cfg.get("profile_bin_size", 1))
    supervised_bp = int(data_cfg.get("supervised_bp", data_cfg.get("output_len", 1000)))
    nonpeak_ratio = float(data_cfg.get("nonpeak_ratio", 1.0))
    track_total_count_target = float(data_cfg.get("track_total_count_target", 0.0))
    region_source = _resolve_region_source(data_cfg, split)
    max_records = int(max_regions_override) if max_regions_override > 0 else int(data_cfg.get(f"max_{split}_regions", 0))

    genome_fasta = str(_resolve_real_data_value(data_cfg, data_source_cfg, "genome_fasta", ""))
    folds_json = str(_resolve_real_data_value(data_cfg, data_source_cfg, "folds_json", ""))
    peaks_bed = str(_resolve_real_data_value(data_cfg, data_source_cfg, "peaks_bed", ""))
    nonpeaks_bed = str(_resolve_real_data_value(data_cfg, data_source_cfg, "nonpeaks_bed", ""))
    bigwig_path = str(_resolve_real_data_value(data_cfg, data_source_cfg, "bigwig", ""))

    required = {
        "genome_fasta": genome_fasta,
        "folds_json": folds_json,
        "peaks_bed": peaks_bed,
        "bigwig": bigwig_path,
    }
    missing = [name for name, value in required.items() if not value]
    if missing:
        raise ValueError("Real-data backend requires the following paths: " + ", ".join(missing))

    return ChromBPNetBigWigDataset(
        genome_fasta=genome_fasta,
        bigwig_path=bigwig_path,
        peaks_bed=peaks_bed,
        nonpeaks_bed=nonpeaks_bed,
        folds_json=folds_json,
        split=split,
        input_len=input_len,
        supervised_bp=supervised_bp,
        profile_bin_size=profile_bin_size,
        max_jitter=0,
        peak_max_jitter=0,
        nonpeak_max_jitter=0,
        seed=resolve_dataset_seed(int(train_cfg.get("seed", 1234)), split),
        nonpeak_ratio=nonpeak_ratio,
        max_records=max_records,
        region_source=region_source,
        random_revcomp=False,
        revcomp_prob=float(data_cfg.get("revcomp_prob", 0.5)),
        track_total_count_target=track_total_count_target,
        genos_cache_dir=str(data_cfg.get("genos_cache_dir", "")),
        genos_cache_features=list(data_cfg.get("genos_cache_features", [])),
        foundation_cache_dir=str(data_cfg.get("foundation_cache_dir", "")),
        foundation_cache_features=list(data_cfg.get("foundation_cache_features", [])),
        teacher_cache_dir="",
        teacher_target_names=[],
    )


def pool_profile_probs_to_bins(profile_probs: Any, n_bins: int = 16) -> np.ndarray:
    """Pool per-base profile probabilities using the training distill contract."""

    probs = _to_numpy_f32(profile_probs)
    if probs.ndim != 2:
        raise ValueError(f"Expected profile_probs [N, L], got shape={tuple(probs.shape)!r}")
    if probs.shape[1] <= 0:
        raise ValueError("profile_probs must have a non-empty sequence axis")
    n_bins = int(n_bins)
    if n_bins <= 0:
        raise ValueError(f"n_bins must be positive, got {n_bins}")
    pooled = training_pool_profile_probs_to_bins(torch.from_numpy(probs), n_bins)
    return pooled.detach().cpu().numpy().astype(np.float32, copy=False)


def export_arrays(
    *,
    output_dir: Path | str,
    split: str,
    profile16: Any,
    logcount: Any,
) -> None:
    output_root = Path(output_dir).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    np.save(output_root / f"{split}_profile16.f32.npy", _to_numpy_f32(profile16))
    np.save(output_root / f"{split}_logcount.f32.npy", _to_numpy_f32(logcount).reshape(-1))


def write_teacher_manifest(
    *,
    output_dir: Path | str,
    split: str,
    n_records: int,
    record_sha1: str,
    model_checkpoint: Path | str,
    targets: Sequence[str] = DEFAULT_TARGETS,
    extra_metadata: Mapping[str, Any] | None = None,
) -> Path:
    manifest = {
        "split": str(split),
        "n_records": int(n_records),
        "record_sha1": str(record_sha1),
        "targets": [str(name) for name in targets],
        "model_checkpoint": str(Path(model_checkpoint).resolve()),
    }
    if extra_metadata:
        manifest.update(dict(extra_metadata))

    output_root = Path(output_dir).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    manifest_path = output_root / f"teacher_manifest_{split}.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest_path


def export_model_teacher_cache(
    *,
    checkpoint: Path | str,
    output_dir: Path | str,
    data_config: Path | str = "",
    splits: Sequence[str] = DEFAULT_SPLITS,
    batch_size: int = 0,
    num_workers: int = -1,
    max_regions: int = 0,
    device: str = "auto",
    log_every: int = 100,
) -> None:
    checkpoint_path = Path(checkpoint).resolve()
    payload = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    if not isinstance(payload, dict):
        raise ValueError(f"Unsupported checkpoint payload at {checkpoint_path}")

    model_cfg = copy.deepcopy(payload.get("model_config", {}))
    train_cfg = copy.deepcopy(payload.get("train_config", {}))
    model_state = payload.get("model_state")
    if not model_cfg or not train_cfg or not isinstance(model_state, dict):
        raise ValueError(f"Checkpoint {checkpoint_path} is missing model/train config payloads or model_state")

    cli_data_config = _resolve_path(str(data_config), REPO_ROOT) if str(data_config).strip() else ""
    data_config_path = resolve_data_config_path(train_cfg, cli_data_config)
    if data_config_path:
        data_config_path = _resolve_path(str(data_config_path), REPO_ROOT)
    data_source_cfg = load_yaml(data_config_path) if data_config_path else {}

    if data_source_cfg:
        apply_external_data_config_defaults(train_cfg, data_source_cfg)
    apply_training_mode_defaults(train_cfg)
    validate_semantics_profile(train_cfg, data_source_cfg, rank=0)
    check_profile_resolution_compatibility(model_cfg, train_cfg, rank=0)
    validate_model_config(model_cfg)

    bias_cfg = model_cfg.setdefault("bias_branch", {})
    if bool(bias_cfg.get("freeze_bias_core", False)) and not str(bias_cfg.get("pretrained_path", "")):
        bias_cfg["freeze_bias_core"] = False

    resolved_device = _get_device(device)
    model, _ = build_model(model_cfg, train_cfg, rank=0)
    load_model_state_with_count_pool_compat(model, model_state)
    model.to(resolved_device)
    model.eval()

    online_runtime = build_online_feature_runtime(model_cfg, resolved_device, rank=0)
    uses_cached_genos = bool(model_cfg.get("genos_cached", {}).get("enabled", False))

    data_cfg = train_cfg.get("data", {})
    loader_batch_size = int(batch_size) if batch_size > 0 else int(data_cfg.get("batch_size_per_gpu", 16))
    loader_num_workers = int(num_workers) if num_workers >= 0 else int(data_cfg.get("num_workers", 2))
    loader_kwargs: dict[str, Any] = {
        "batch_size": loader_batch_size,
        "shuffle": False,
        "num_workers": loader_num_workers,
        "pin_memory": bool(data_cfg.get("pin_memory", True)) and (resolved_device.type == "cuda"),
        "drop_last": False,
    }
    if loader_num_workers > 0:
        loader_kwargs["prefetch_factor"] = int(data_cfg.get("prefetch_factor", 2))
        loader_kwargs["persistent_workers"] = bool(data_cfg.get("persistent_workers", True))

    output_root = Path(output_dir).resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    with torch.no_grad():
        for split in splits:
            dataset = _build_export_dataset(train_cfg, data_source_cfg, str(split), int(max_regions))
            loader = DataLoader(dataset, **loader_kwargs)
            profile_batches: list[np.ndarray] = []
            logcount_batches: list[np.ndarray] = []

            for batch_idx, batch in enumerate(loader, start=1):
                seq = batch["seq"].to(resolved_device, non_blocking=True)
                online_feature_kwargs = extract_online_feature_kwargs(seq, online_runtime)
                foundation_feature_kwargs = extract_foundation_feature_kwargs(batch, model_cfg, resolved_device)
                genos_summary_kwargs = extract_genos_summary_kwargs(
                    batch,
                    model_cfg,
                    resolved_device,
                    require_genos_summary=uses_cached_genos,
                )

                outputs = model(
                    seq,
                    **genos_summary_kwargs,
                    **online_feature_kwargs,
                    **foundation_feature_kwargs,
                )
                profile_probs = torch.softmax(outputs.profile_logits_debiased.detach().float(), dim=-1)
                profile_batches.append(pool_profile_probs_to_bins(profile_probs, n_bins=16))
                logcount_batches.append(_to_numpy_f32(outputs.logcount_debiased).reshape(-1))

                if log_every > 0 and (batch_idx % log_every == 0):
                    print(f"[export] split={split} batch={batch_idx} records={sum(x.shape[0] for x in logcount_batches)}")

            if profile_batches:
                profile16 = np.concatenate(profile_batches, axis=0).astype(np.float32, copy=False)
                logcount_array = np.concatenate(logcount_batches, axis=0).astype(np.float32, copy=False)
            else:
                profile16 = np.zeros((0, 16), dtype=np.float32)
                logcount_array = np.zeros((0,), dtype=np.float32)

            if profile16.shape[0] != len(dataset) or logcount_array.shape[0] != len(dataset):
                raise ValueError(
                    f"Export row mismatch for split={split!r}: profile16={profile16.shape[0]} "
                    f"logcount={logcount_array.shape[0]} dataset={len(dataset)}"
                )

            export_arrays(output_dir=output_root, split=str(split), profile16=profile16, logcount=logcount_array)
            write_teacher_manifest(
                output_dir=output_root,
                split=str(split),
                n_records=len(dataset),
                record_sha1=compute_record_sha1(dataset.records),
                model_checkpoint=checkpoint_path,
                extra_metadata={
                    "data_config": data_config_path,
                    "profile_bins": 16,
                },
            )


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", required=True, help="Teacher checkpoint .pt path")
    parser.add_argument("--output-dir", required=True, help="Output directory for cache arrays and manifests")
    parser.add_argument("--data-config", default="", help="Optional override for external data config YAML")
    parser.add_argument("--splits", nargs="+", default=list(DEFAULT_SPLITS), help="Splits to export")
    parser.add_argument("--batch-size", type=int, default=0, help="Override batch size when > 0")
    parser.add_argument("--num-workers", type=int, default=-1, help="Override num_workers when >= 0")
    parser.add_argument("--max-regions", type=int, default=0, help="Override max records when > 0")
    parser.add_argument("--device", default="auto", help="cuda, cpu, or auto")
    parser.add_argument("--log-every", type=int, default=100, help="Progress print interval in steps")
    return parser.parse_args(list(argv) if argv is not None else None)


def main(argv: Iterable[str] | None = None) -> None:
    args = parse_args(argv)
    export_model_teacher_cache(
        checkpoint=args.checkpoint,
        output_dir=args.output_dir,
        data_config=args.data_config,
        splits=args.splits,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        max_regions=args.max_regions,
        device=args.device,
        log_every=args.log_every,
    )


if __name__ == "__main__":
    main()
