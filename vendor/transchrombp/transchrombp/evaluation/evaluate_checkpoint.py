"""Evaluate a TransChromBP checkpoint on a held-out split."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Sequence

import numpy as np
import torch
from torch.utils.data import DataLoader

from transchrombp.data import ChromBPNetBigWigDataset
from transchrombp.training.train_ddp import (
    BIAS_DIAGNOSTIC_KEYS,
    LOSS_METRIC_KEYS,
    VALIDATION_REGION_KEYS,
    apply_external_data_config_defaults,
    apply_training_mode_defaults,
    build_online_feature_runtime,
    build_model,
    compute_bias_diagnostics,
    compute_losses,
    compute_profile_jsd,
    extract_foundation_feature_kwargs,
    extract_online_feature_kwargs,
    extract_teacher_targets,
    load_yaml,
    resolve_amp_dtype,
    resolve_data_config_path,
    slice_outputs,
    validate_model_config,
)


SCOPE_KEYS: tuple[str, ...] = ("overall",) + VALIDATION_REGION_KEYS
REPO_ROOT = Path(__file__).resolve().parents[3]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate a TransChromBP checkpoint on a held-out split")
    parser.add_argument("--checkpoint", required=True, help="Checkpoint .pt path")
    parser.add_argument("--split", default="test", help="Fold split to evaluate, default: test")
    parser.add_argument("--data-config", default="", help="Optional override for data config YAML")
    parser.add_argument("--batch-size", type=int, default=0, help="Override batch size when > 0")
    parser.add_argument("--num-workers", type=int, default=-1, help="Override num_workers when >= 0")
    parser.add_argument("--max-regions", type=int, default=0, help="Override max records when > 0")
    parser.add_argument("--region-source", default="", help="Optional override: both, peak, or nonpeak")
    parser.add_argument(
        "--nonpeak-ratio",
        type=float,
        default=-1.0,
        help="Optional override for nonpeak_ratio; use 1.0 for full held-out evaluation",
    )
    parser.add_argument("--device", default="auto", help="cuda, cpu, or auto")
    parser.add_argument("--output", default="", help="Where to write JSON results")
    parser.add_argument("--log-every", type=int, default=100, help="Progress print interval in steps")
    return parser.parse_args()


def resolve_path(path_value: str, base_dir: Path) -> str:
    path = Path(path_value)
    if path.is_absolute():
        return str(path)
    return str((base_dir / path).resolve())


def get_device(device_name: str) -> torch.device:
    if device_name == "auto":
        return torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    return torch.device(device_name)


def rankdata_average(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    ranks = np.empty(values.shape[0], dtype=np.float64)

    start = 0
    while start < sorted_values.shape[0]:
        end = start + 1
        while end < sorted_values.shape[0] and sorted_values[end] == sorted_values[start]:
            end += 1
        average_rank = 0.5 * (start + end - 1) + 1.0
        ranks[order[start:end]] = average_rank
        start = end
    return ranks


def pearsonr_safe(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2 or y.size < 2:
        return float("nan")
    x = x.astype(np.float64, copy=False)
    y = y.astype(np.float64, copy=False)
    x = x - x.mean()
    y = y - y.mean()
    denom = np.sqrt(np.sum(x * x) * np.sum(y * y))
    if denom <= 0.0:
        return float("nan")
    return float(np.sum(x * y) / denom)


def spearmanr_safe(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2 or y.size < 2:
        return float("nan")
    return pearsonr_safe(rankdata_average(x), rankdata_average(y))


def auroc_numpy(labels: np.ndarray, scores: np.ndarray) -> float:
    """Compute AUROC using numpy only (no sklearn dependency)."""
    if labels.size < 2:
        return float("nan")
    n_pos = int(labels.sum())
    n_neg = labels.size - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(-scores, kind="mergesort")
    sorted_labels = labels[order].astype(np.float64)
    tps = np.cumsum(sorted_labels)
    fps = np.cumsum(1.0 - sorted_labels)
    tpr = np.concatenate([[0.0], tps / n_pos])
    fpr = np.concatenate([[0.0], fps / n_neg])
    return float(np.trapz(tpr, fpr))


def auprc_numpy(labels: np.ndarray, scores: np.ndarray) -> float:
    """Compute Average Precision (area under precision-recall curve) using numpy only."""
    if labels.size < 2:
        return float("nan")
    n_pos = int(labels.sum())
    if n_pos == 0 or n_pos == labels.size:
        return float("nan")
    order = np.argsort(-scores, kind="mergesort")
    sorted_labels = labels[order].astype(np.float64)
    tps = np.cumsum(sorted_labels)
    precision = tps / np.arange(1, labels.size + 1, dtype=np.float64)
    recall = tps / n_pos
    recall_change = np.diff(recall, prepend=0.0)
    return float(np.sum(precision * recall_change))


def best_f1_metrics(labels: np.ndarray, scores: np.ndarray) -> Dict[str, float]:
    """Find threshold that maximizes F1, return precision/recall/f1/threshold."""
    n_pos = int(labels.sum())
    if labels.size < 2 or n_pos == 0 or n_pos == labels.size:
        return {"best_f1": float("nan"), "best_f1_precision": float("nan"),
                "best_f1_recall": float("nan"), "best_f1_threshold": float("nan")}
    order = np.argsort(-scores, kind="mergesort")
    sorted_labels = labels[order].astype(np.float64)
    sorted_scores = scores[order]
    tps = np.cumsum(sorted_labels)
    fps = np.cumsum(1.0 - sorted_labels)
    precision = tps / (tps + fps)
    recall = tps / n_pos
    f1 = 2 * precision * recall / np.maximum(precision + recall, 1e-8)
    best_idx = int(np.argmax(f1))
    return {
        "best_f1": float(f1[best_idx]),
        "best_f1_precision": float(precision[best_idx]),
        "best_f1_recall": float(recall[best_idx]),
        "best_f1_threshold": float(sorted_scores[best_idx]),
    }


def profile_pearson_batch(target: torch.Tensor, pred: torch.Tensor) -> torch.Tensor:
    """Per-region Pearson correlation between target and predicted profiles. Returns (batch,)."""
    t = target.float()
    p = pred.float()
    t = t - t.mean(dim=-1, keepdim=True)
    p = p - p.mean(dim=-1, keepdim=True)
    num = (t * p).sum(dim=-1)
    denom = torch.sqrt((t * t).sum(dim=-1) * (p * p).sum(dim=-1))
    return torch.where(denom > 0, num / denom, torch.zeros_like(num))


def summarize_values(values: list[np.ndarray], prefix: str) -> Dict[str, float]:
    if not values:
        return {
            f"{prefix}_mean": float("nan"),
            f"{prefix}_median": float("nan"),
        }
    merged = np.concatenate(values, axis=0).astype(np.float64, copy=False)
    return {
        f"{prefix}_mean": float(np.mean(merged)),
        f"{prefix}_median": float(np.median(merged)),
    }


def finalize_scope(metric_sums: Dict[str, float], total_weight: float) -> Dict[str, float]:
    if total_weight <= 0.0:
        return {}
    return {key: float(metric_sums[key] / total_weight) for key in metric_sums}


def append_array(store: Dict[str, list[np.ndarray]], key: str, value: torch.Tensor) -> None:
    store[key].append(value.detach().cpu().numpy().reshape(-1))


def make_output_path(args: argparse.Namespace, checkpoint_path: Path) -> Path:
    if args.output:
        return Path(args.output)
    label = f"{checkpoint_path.parent.name}__{checkpoint_path.stem}__{args.split}.json"
    return REPO_ROOT / "outputs" / "metrics" / label


def load_model_state_with_count_pool_compat(model: torch.nn.Module, state_dict: Dict[str, Any]) -> None:
    try:
        model.load_state_dict(state_dict, strict=True)
        return
    except RuntimeError as exc:
        missing_allowed = {
            "count_pool_proj.weight",
            "count_pool_proj.bias",
        }
        model_keys = set(model.state_dict().keys())
        state_keys = set(state_dict.keys())
        missing = sorted(model_keys - state_keys)
        unexpected = sorted(state_keys - model_keys)
        if any(key not in missing_allowed for key in missing):
            raise exc
        if any(not key.startswith("count_pool_proj.") for key in unexpected):
            raise exc
        model.load_state_dict(state_dict, strict=False)
        print(
            "[warn] Loaded checkpoint with count_pool_proj backward-compat mode; "
            "old attention-pooling checkpoints were trained before the pooling layer was implemented."
        )


def build_eval_dataset(
    train_cfg: Dict[str, Any],
    data_source_cfg: Dict[str, Any],
    split: str,
    max_regions_override: int = 0,
    region_source_override: str = "",
    nonpeak_ratio_override: float = -1.0,
) -> ChromBPNetBigWigDataset:
    data_cfg = train_cfg.get("data", {})
    input_len = int(data_cfg.get("input_len", 2114))
    output_len = int(data_cfg.get("output_len", 1000))
    supervised_bp = int(data_cfg.get("supervised_bp", output_len))
    profile_bin_size = int(data_cfg.get("profile_bin_size", 1))
    nonpeak_ratio = (
        float(nonpeak_ratio_override)
        if nonpeak_ratio_override >= 0.0
        else float(data_cfg.get("nonpeak_ratio", 1.0))
    )
    region_source = str(
        region_source_override
        or data_cfg.get(f"{split}_region_source", data_cfg.get("val_region_source", data_cfg.get("region_source", "both")))
    )
    track_total_count_target = float(data_cfg.get("track_total_count_target", 0.0))
    genos_cache_dir = str(data_cfg.get("genos_cache_dir", ""))
    genos_cache_features = list(data_cfg.get("genos_cache_features", []))
    foundation_cache_dir = str(data_cfg.get("foundation_cache_dir", ""))
    foundation_cache_features = list(data_cfg.get("foundation_cache_features", []))
    teacher_cache_dir = str(data_cfg.get("teacher_cache_dir", ""))
    teacher_target_names = list(data_cfg.get("teacher_target_names", []))

    max_records = int(max_regions_override) if max_regions_override > 0 else int(data_cfg.get(f"max_{split}_regions", 0))

    return ChromBPNetBigWigDataset(
        genome_fasta=str(data_source_cfg["genome_fasta"]),
        bigwig_path=str(data_source_cfg["input"]["bigwig"]),
        peaks_bed=str(data_source_cfg["input"]["peaks_bed"]),
        nonpeaks_bed=str(data_source_cfg["input"]["nonpeaks_bed"]),
        folds_json=str(data_source_cfg["folds_json"]),
        split=split,
        input_len=input_len,
        supervised_bp=supervised_bp,
        profile_bin_size=profile_bin_size,
        max_jitter=0,
        peak_max_jitter=0,
        nonpeak_max_jitter=0,
        seed=int(train_cfg.get("seed", 1234)) + 20_000,
        nonpeak_ratio=nonpeak_ratio,
        max_records=max_records,
        region_source=region_source,
        random_revcomp=False,
        revcomp_prob=float(data_cfg.get("revcomp_prob", 0.5)),
        track_total_count_target=track_total_count_target,
        genos_cache_dir=genos_cache_dir,
        genos_cache_features=genos_cache_features,
        foundation_cache_dir=foundation_cache_dir,
        foundation_cache_features=foundation_cache_features,
        teacher_cache_dir=teacher_cache_dir,
        teacher_target_names=teacher_target_names,
    )


def main() -> None:
    args = parse_args()
    checkpoint_path = Path(args.checkpoint).resolve()
    payload = torch.load(checkpoint_path, map_location="cpu", weights_only=False)

    model_cfg = payload.get("model_config", {})
    train_cfg = payload.get("train_config", {})
    if not model_cfg or not train_cfg:
        raise ValueError(f"Checkpoint {checkpoint_path} is missing model/train config payloads")

    repo_root = REPO_ROOT
    cli_data_config = resolve_path(args.data_config, repo_root) if args.data_config else ""
    data_config_path = resolve_data_config_path(train_cfg, cli_data_config)
    if data_config_path:
        data_config_path = resolve_path(data_config_path, repo_root)
    data_source_cfg = load_yaml(data_config_path) if data_config_path else {}
    if data_source_cfg:
        apply_external_data_config_defaults(train_cfg, data_source_cfg)
    apply_training_mode_defaults(train_cfg)

    bias_cfg = model_cfg.setdefault("bias_branch", {})
    if bool(bias_cfg.get("freeze_bias_core", False)) and not str(bias_cfg.get("pretrained_path", "")):
        bias_cfg["freeze_bias_core"] = False
        print("[warn] checkpoint has freeze_bias_core=true without pretrained_path; overriding to false for evaluation")

    validate_model_config(model_cfg)

    device = get_device(args.device)
    model, training_mode = build_model(model_cfg, train_cfg, rank=0)
    load_model_state_with_count_pool_compat(model, payload["model_state"])
    model.to(device)
    model.eval()
    online_runtime = build_online_feature_runtime(model_cfg, device, rank=0)
    uses_cached_genos = bool(model_cfg.get("genos_cached", {}).get("enabled", False))

    data_cfg = train_cfg.get("data", {})
    batch_size = int(args.batch_size) if args.batch_size > 0 else int(data_cfg.get("batch_size_per_gpu", 16))
    num_workers = int(args.num_workers) if args.num_workers >= 0 else int(data_cfg.get("num_workers", 2))
    pin_memory = bool(data_cfg.get("pin_memory", True)) and (device.type == "cuda")
    prefetch_factor = int(data_cfg.get("prefetch_factor", 2))
    persistent_workers = bool(data_cfg.get("persistent_workers", True)) if num_workers > 0 else False

    eval_ds = build_eval_dataset(
        train_cfg,
        data_source_cfg,
        split=args.split,
        max_regions_override=args.max_regions,
        region_source_override=args.region_source,
        nonpeak_ratio_override=args.nonpeak_ratio,
    )
    loader_kwargs: Dict[str, Any] = {
        "batch_size": batch_size,
        "shuffle": False,
        "num_workers": num_workers,
        "pin_memory": pin_memory,
        "drop_last": False,
    }
    if num_workers > 0:
        loader_kwargs["prefetch_factor"] = prefetch_factor
        loader_kwargs["persistent_workers"] = persistent_workers
    eval_loader = DataLoader(eval_ds, **loader_kwargs)

    trainer_cfg = train_cfg.get("trainer", {})
    precision = str(trainer_cfg.get("precision", "bf16"))
    amp_dtype = resolve_amp_dtype(precision)
    use_amp = (device.type == "cuda") and (amp_dtype != torch.float32)

    metric_template = {key: 0.0 for key in LOSS_METRIC_KEYS + BIAS_DIAGNOSTIC_KEYS}
    metric_sums = {scope: dict(metric_template) for scope in SCOPE_KEYS}
    total_examples = {scope: 0.0 for scope in SCOPE_KEYS}
    arrays = {
        scope: {
            "target_total": [],
            "target_log": [],
            "pred_total_full": [],
            "pred_total_debiased": [],
            "pred_log_full": [],
            "pred_log_debiased": [],
            "profile_jsd_full": [],
            "profile_jsd_debiased": [],
            "profile_pearson_full": [],
            "profile_pearson_debiased": [],
        }
        for scope in SCOPE_KEYS
    }
    region_labels: list[np.ndarray] = []  # 1=peak, 0=nonpeak, for overall scope classification metrics

    def update_scope(
        scope: str,
        outputs: Any,
        counts: torch.Tensor,
        batch_weight: float,
        teacher_targets: Dict[str, torch.Tensor] | None = None,
    ) -> None:
        _, loss_metrics = compute_losses(
            outputs,
            counts,
            train_cfg.get("loss", {}),
            teacher_targets=teacher_targets,
        )
        diag_metrics = compute_bias_diagnostics(model, outputs)
        combined = {**loss_metrics, **diag_metrics}
        total_examples[scope] += float(batch_weight)
        for key in metric_sums[scope]:
            metric_sums[scope][key] += float(combined[key].detach().item()) * float(batch_weight)

        counts = counts.detach().float()
        target_total = counts.sum(dim=-1)
        target_log = torch.log1p(target_total)
        pred_log_full = outputs.logcount_full.detach().float().reshape(-1)
        pred_log_debiased = outputs.logcount_debiased.detach().float().reshape(-1)
        pred_total_full = torch.clamp(torch.expm1(pred_log_full), min=0.0)
        pred_total_debiased = torch.clamp(torch.expm1(pred_log_debiased), min=0.0)

        target_probs = counts + 1e-6
        target_probs = target_probs / target_probs.sum(dim=-1, keepdim=True)
        full_probs = torch.softmax(outputs.profile_logits_full.detach().float(), dim=-1)
        debiased_probs = torch.softmax(outputs.profile_logits_debiased.detach().float(), dim=-1)
        jsd_full = compute_profile_jsd(full_probs, target_probs)
        jsd_debiased = compute_profile_jsd(debiased_probs, target_probs)

        # Per-region profile Pearson correlation
        pr_full = profile_pearson_batch(target_probs, full_probs)
        pr_debiased = profile_pearson_batch(target_probs, debiased_probs)

        append_array(arrays[scope], "target_total", target_total)
        append_array(arrays[scope], "target_log", target_log)
        append_array(arrays[scope], "pred_total_full", pred_total_full)
        append_array(arrays[scope], "pred_total_debiased", pred_total_debiased)
        append_array(arrays[scope], "pred_log_full", pred_log_full)
        append_array(arrays[scope], "pred_log_debiased", pred_log_debiased)
        append_array(arrays[scope], "profile_jsd_full", jsd_full)
        append_array(arrays[scope], "profile_jsd_debiased", jsd_debiased)
        append_array(arrays[scope], "profile_pearson_full", pr_full)
        append_array(arrays[scope], "profile_pearson_debiased", pr_debiased)

    with torch.no_grad():
        for step_idx, batch in enumerate(eval_loader, start=1):
            seq = batch["seq"].to(device, non_blocking=True)
            profile_counts = batch["profile_counts"].to(device, non_blocking=True)
            region_sources = batch.get("region_source")
            online_feature_kwargs = extract_online_feature_kwargs(seq, online_runtime)
            foundation_feature_kwargs = extract_foundation_feature_kwargs(batch, model_cfg, device)
            teacher_targets = extract_teacher_targets(batch, device)
            genos_summary = None
            if "genos_global_mean" in batch:
                genos_summary = batch["genos_global_mean"].to(device, non_blocking=True)
            elif uses_cached_genos:
                raise RuntimeError(
                    "Cached Genos mode is enabled, but evaluation batch is missing "
                    "'genos_global_mean'. Check train_config.data.genos_cache_dir / genos_cache_features."
                )

            with torch.autocast(device_type=device.type, dtype=amp_dtype, enabled=use_amp):
                outputs = model(
                    seq,
                    genos_summary=genos_summary,
                    **online_feature_kwargs,
                    **foundation_feature_kwargs,
                )

            batch_weight = float(seq.shape[0])
            update_scope("overall", outputs, profile_counts, batch_weight, teacher_targets=teacher_targets)

            # Collect region labels for classification metrics (1=peak, 0=nonpeak)
            if isinstance(region_sources, Sequence) and not isinstance(region_sources, (str, bytes)):
                batch_labels = np.array([1.0 if rs == "peak" else 0.0 for rs in region_sources], dtype=np.float64)
                region_labels.append(batch_labels)

            if isinstance(region_sources, Sequence) and not isinstance(region_sources, (str, bytes)):
                for source in VALIDATION_REGION_KEYS:
                    mask_values = [region_source == source for region_source in region_sources]
                    local_count = sum(mask_values)
                    if local_count <= 0:
                        continue
                    mask = torch.tensor(mask_values, device=device, dtype=torch.bool)
                    source_outputs = slice_outputs(outputs, mask)
                    source_counts = profile_counts[mask]
                    source_teacher_targets = {
                        name: value[mask]
                        for name, value in teacher_targets.items()
                    }
                    update_scope(
                        source,
                        source_outputs,
                        source_counts,
                        float(local_count),
                        teacher_targets=source_teacher_targets,
                    )

            if args.log_every > 0 and (step_idx % args.log_every == 0):
                print(f"[eval] step={step_idx}/{len(eval_loader)} examples={int(total_examples['overall'])}")

    results: Dict[str, Any] = {
        "checkpoint": str(checkpoint_path),
        "checkpoint_epoch": int(payload.get("epoch", 0)),
        "checkpoint_global_step": int(payload.get("global_step", 0)),
        "split": args.split,
        "device": str(device),
        "training_mode": training_mode,
        "batch_size": batch_size,
        "num_workers": num_workers,
        "dataset_size": int(len(eval_ds)),
        "nonpeak_ratio": float(args.nonpeak_ratio) if args.nonpeak_ratio >= 0.0 else float(data_cfg.get("nonpeak_ratio", 1.0)),
        "data_config_path": data_config_path,
        "model_run_name": checkpoint_path.parent.name,
        "results": {},
    }

    for scope in SCOPE_KEYS:
        n_examples = total_examples[scope]
        if n_examples <= 0:
            continue
        target_total = np.concatenate(arrays[scope]["target_total"], axis=0)
        target_log = np.concatenate(arrays[scope]["target_log"], axis=0)
        pred_total_full = np.concatenate(arrays[scope]["pred_total_full"], axis=0)
        pred_total_debiased = np.concatenate(arrays[scope]["pred_total_debiased"], axis=0)
        pred_log_full = np.concatenate(arrays[scope]["pred_log_full"], axis=0)
        pred_log_debiased = np.concatenate(arrays[scope]["pred_log_debiased"], axis=0)

        scope_result = {
            "n_examples": int(round(n_examples)),
            **finalize_scope(metric_sums[scope], n_examples),
            "count_pearson_full": pearsonr_safe(target_total, pred_total_full),
            "count_spearman_full": spearmanr_safe(target_total, pred_total_full),
            "count_pearson_debiased": pearsonr_safe(target_total, pred_total_debiased),
            "count_spearman_debiased": spearmanr_safe(target_total, pred_total_debiased),
            "logcount_pearson_full": pearsonr_safe(target_log, pred_log_full),
            "logcount_spearman_full": spearmanr_safe(target_log, pred_log_full),
            "logcount_pearson_debiased": pearsonr_safe(target_log, pred_log_debiased),
            "logcount_spearman_debiased": spearmanr_safe(target_log, pred_log_debiased),
            "count_mae_full": float(np.mean(np.abs(pred_total_full - target_total))),
            "count_mae_debiased": float(np.mean(np.abs(pred_total_debiased - target_total))),
            **summarize_values(arrays[scope]["profile_jsd_full"], "profile_target_jsd_full"),
            **summarize_values(arrays[scope]["profile_jsd_debiased"], "profile_target_jsd_debiased"),
            **summarize_values(arrays[scope]["profile_pearson_full"], "profile_pearson_full"),
            **summarize_values(arrays[scope]["profile_pearson_debiased"], "profile_pearson_debiased"),
        }

        # Classification metrics: peak vs non-peak using predicted count as score
        if scope == "overall" and region_labels:
            all_labels = np.concatenate(region_labels, axis=0)
            if all_labels.sum() > 0 and all_labels.sum() < all_labels.size:
                scope_result["peak_auroc_count_full"] = auroc_numpy(all_labels, pred_total_full)
                scope_result["peak_auroc_count_debiased"] = auroc_numpy(all_labels, pred_total_debiased)
                scope_result["peak_auroc_logcount_full"] = auroc_numpy(all_labels, pred_log_full)
                scope_result["peak_auroc_logcount_debiased"] = auroc_numpy(all_labels, pred_log_debiased)
                scope_result["peak_auprc_count_full"] = auprc_numpy(all_labels, pred_total_full)
                scope_result["peak_auprc_count_debiased"] = auprc_numpy(all_labels, pred_total_debiased)
                scope_result["peak_auprc_logcount_full"] = auprc_numpy(all_labels, pred_log_full)
                scope_result["peak_auprc_logcount_debiased"] = auprc_numpy(all_labels, pred_log_debiased)
                f1_full = best_f1_metrics(all_labels, pred_total_full)
                f1_deb = best_f1_metrics(all_labels, pred_total_debiased)
                f1_log_full = best_f1_metrics(all_labels, pred_log_full)
                f1_log_deb = best_f1_metrics(all_labels, pred_log_debiased)
                for k, v in f1_full.items():
                    scope_result[f"peak_{k}_count_full"] = v
                for k, v in f1_deb.items():
                    scope_result[f"peak_{k}_count_debiased"] = v
                for k, v in f1_log_full.items():
                    scope_result[f"peak_{k}_logcount_full"] = v
                for k, v in f1_log_deb.items():
                    scope_result[f"peak_{k}_logcount_debiased"] = v
                scope_result["peak_classification_n_pos"] = int(all_labels.sum())
                scope_result["peak_classification_n_neg"] = int(all_labels.size - all_labels.sum())

        results["results"][scope] = scope_result

    output_path = make_output_path(args, checkpoint_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"[done] wrote {output_path}")
    for scope in SCOPE_KEYS:
        scope_result = results["results"].get(scope)
        if not scope_result:
            continue
        summary_parts = [
            f"[summary:{scope}]",
            f"n={scope_result['n_examples']}",
            f"loss={scope_result['loss_total']:.5f}",
            f"count_r={scope_result['count_pearson_full']:.4f}",
            f"jsd={scope_result['profile_target_jsd_full_mean']:.4f}",
            f"prof_r={scope_result.get('profile_pearson_full_mean', float('nan')):.4f}",
        ]
        if "peak_auroc_count_full" in scope_result:
            summary_parts.extend([
                f"auroc={scope_result['peak_auroc_count_full']:.4f}",
                f"auprc={scope_result['peak_auprc_count_full']:.4f}",
                f"f1={scope_result['peak_best_f1_count_full']:.4f}",
            ])
        print(" ".join(summary_parts))


if __name__ == "__main__":
    main()
