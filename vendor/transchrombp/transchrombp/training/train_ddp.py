"""DDP training entrypoint for TransChromBP.

Single-node multi-GPU training features:
- torch.distributed + DistributedDataParallel (DDP)
- AMP (bf16/fp16/fp32)
- DistributedSampler
- rank0-only logging/checkpointing

Current data backend supports:
- synthetic data for infrastructure validation
- FASTA + bigWig + BED regions for tutorial/real-data smoke runs
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import shutil
import time
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Sequence, Tuple

import numpy as np
import torch
import torch.distributed as dist
import torch.nn.functional as F
import yaml
from torch import Tensor, nn
from torch.nn.parallel import DistributedDataParallel as DDP
from torch.utils.data import DataLoader, Dataset, DistributedSampler, Sampler

from transchrombp.data import ChromBPNetBigWigDataset
from transchrombp.models import (
    TransChromBPOutput,
    build_bias_branch_from_config,
    build_transchrombp_from_config,
)
from transchrombp.utils.foundation_contract import (
    foundation_batch_key,
    resolve_foundation_contract,
    validate_foundation_cache_config,
)


NORM_LAYER_TYPES: Tuple[type, ...] = (
    nn.LayerNorm,
    nn.BatchNorm1d,
    nn.BatchNorm2d,
    nn.BatchNorm3d,
    nn.SyncBatchNorm,
    nn.GroupNorm,
    nn.InstanceNorm1d,
    nn.InstanceNorm2d,
    nn.InstanceNorm3d,
)

LOSS_METRIC_KEYS: Tuple[str, ...] = (
    "loss_total",
    "loss_profile",
    "loss_count",
    "loss_debiased_profile",
    "loss_debiased_count",
    "loss_distill_profile",
    "loss_distill_count",
    "loss_distill_rank",
)

BIAS_DIAGNOSTIC_KEYS: Tuple[str, ...] = (
    "effective_profile_scale",
    "effective_count_scale",
    "profile_bias_rms_over_signal_rms",
    "profile_full_debiased_jsd",
    "count_full_debiased_abs",
)

VALIDATION_REGION_KEYS: Tuple[str, ...] = ("peak", "nonpeak")

SELECTION_SUM_KEYS: Tuple[str, ...] = (
    "n_examples",
    "target_total_sum",
    "target_total_sq_sum",
    "pred_total_full_sum",
    "pred_total_full_sq_sum",
    "target_pred_total_full_sum",
    "pred_total_debiased_sum",
    "pred_total_debiased_sq_sum",
    "target_pred_total_debiased_sum",
    "count_mae_full_sum",
    "count_mae_debiased_sum",
    "profile_jsd_full_sum",
    "profile_jsd_debiased_sum",
)


def load_yaml(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def set_seed(seed: int, rank: int) -> None:
    s = int(seed) + int(rank)
    random.seed(s)
    np.random.seed(s)
    torch.manual_seed(s)
    torch.cuda.manual_seed_all(s)


@dataclass
class DistEnv:
    rank: int
    local_rank: int
    world_size: int
    distributed: bool
    device: torch.device


@dataclass
class OnlineFeatureRuntime:
    kind: str
    extractor: Any


class PeakNonpeakResamplingSampler(Sampler[int]):
    """Keep all peaks and draw a fresh nonpeak subset each epoch."""

    def __init__(
        self,
        dataset: Dataset,
        num_replicas: int = 1,
        rank: int = 0,
        shuffle: bool = True,
        seed: int = 0,
        drop_last: bool = False,
    ) -> None:
        self.dataset = dataset
        self.num_replicas = int(num_replicas)
        self.rank = int(rank)
        self.shuffle = bool(shuffle)
        self.seed = int(seed)
        self.drop_last = bool(drop_last)
        self.epoch = 0

        self.peak_count = int(getattr(dataset, "peak_record_count", 0))
        self.nonpeak_count = int(getattr(dataset, "nonpeak_record_count", 0))
        self.target_nonpeak_count = int(getattr(dataset, "target_nonpeak_count", 0))
        if self.num_replicas <= 0:
            raise ValueError(f"num_replicas must be positive, got {self.num_replicas}")
        if not (0 <= self.rank < self.num_replicas):
            raise ValueError(f"rank must lie in [0, {self.num_replicas}), got {self.rank}")
        if self.peak_count <= 0 or self.nonpeak_count <= 0 or self.target_nonpeak_count <= 0:
            raise ValueError("PeakNonpeakResamplingSampler requires both peak and nonpeak pools to be non-empty")

        self.epoch_size = self.peak_count + min(self.nonpeak_count, self.target_nonpeak_count)
        if self.drop_last and (self.epoch_size % self.num_replicas != 0):
            self.num_samples = self.epoch_size // self.num_replicas
        else:
            self.num_samples = int(math.ceil(self.epoch_size / self.num_replicas))
        self.total_size = self.num_samples * self.num_replicas

    def __iter__(self) -> Iterable[int]:
        generator = torch.Generator()
        generator.manual_seed(self.seed + self.epoch)

        peak_indices = torch.arange(self.peak_count, dtype=torch.int64)
        if self.target_nonpeak_count >= self.nonpeak_count:
            nonpeak_indices = torch.arange(self.nonpeak_count, dtype=torch.int64)
        else:
            nonpeak_indices = torch.randperm(self.nonpeak_count, generator=generator)[: self.target_nonpeak_count]
        nonpeak_indices = nonpeak_indices + self.peak_count

        indices = torch.cat([peak_indices, nonpeak_indices], dim=0)
        if self.shuffle:
            perm = torch.randperm(indices.numel(), generator=generator)
            indices = indices[perm]

        if not self.drop_last and indices.numel() < self.total_size:
            padding_size = self.total_size - indices.numel()
            repeats = int(math.ceil(padding_size / max(1, indices.numel())))
            padding = indices.repeat(repeats)[:padding_size]
            indices = torch.cat([indices, padding], dim=0)
        else:
            indices = indices[: self.total_size]

        indices = indices[self.rank : self.total_size : self.num_replicas]
        return iter(indices.tolist())

    def __len__(self) -> int:
        return self.num_samples

    def set_epoch(self, epoch: int) -> None:
        self.epoch = int(epoch)


def init_distributed(backend: str = "nccl") -> DistEnv:
    rank = int(os.environ.get("RANK", "0"))
    local_rank = int(os.environ.get("LOCAL_RANK", "0"))
    world_size = int(os.environ.get("WORLD_SIZE", "1"))
    distributed = world_size > 1

    if distributed:
        if torch.cuda.is_available():
            torch.cuda.set_device(local_rank)
        dist.init_process_group(backend=backend)

    if torch.cuda.is_available():
        device = torch.device(f"cuda:{local_rank}")
    else:
        device = torch.device("cpu")

    return DistEnv(
        rank=rank,
        local_rank=local_rank,
        world_size=world_size,
        distributed=distributed,
        device=device,
    )


def cleanup_distributed(distributed: bool) -> None:
    if distributed and dist.is_initialized():
        dist.destroy_process_group()


def is_main_process(rank: int) -> bool:
    return rank == 0


def to_scalar(x: Any) -> float:
    if isinstance(x, (float, int)):
        return float(x)
    if torch.is_tensor(x):
        return float(x.detach().item())
    return float(x)


def unwrap_model(model: nn.Module) -> nn.Module:
    return model.module if isinstance(model, DDP) else model


def reduce_mean(value: Tensor, distributed: bool) -> Tensor:
    if not distributed:
        return value
    dist.all_reduce(value, op=dist.ReduceOp.SUM)
    value = value / dist.get_world_size()
    return value


def compute_profile_jsd(probs_a: Tensor, probs_b: Tensor, eps: float = 1e-8) -> Tensor:
    probs_a = torch.clamp(probs_a, min=eps)
    probs_b = torch.clamp(probs_b, min=eps)
    probs_a = probs_a / probs_a.sum(dim=-1, keepdim=True)
    probs_b = probs_b / probs_b.sum(dim=-1, keepdim=True)
    mean_probs = 0.5 * (probs_a + probs_b)
    kl_a = (probs_a * (torch.log(probs_a) - torch.log(mean_probs))).sum(dim=-1)
    kl_b = (probs_b * (torch.log(probs_b) - torch.log(mean_probs))).sum(dim=-1)
    return torch.sqrt(torch.clamp(0.5 * (kl_a + kl_b), min=0.0))


def compute_bias_diagnostics(model: nn.Module, outputs: TransChromBPOutput) -> Dict[str, Tensor]:
    raw_model = unwrap_model(model)
    zero = outputs.profile_logits_full.detach().new_tensor(0.0)
    diagnostics = {key: zero for key in BIAS_DIAGNOSTIC_KEYS}

    if not (hasattr(raw_model, "profile_scale") and hasattr(raw_model, "_resolve_scale")):
        return diagnostics

    profile_signal = outputs.profile_logits_debiased.detach().float()
    profile_full = outputs.profile_logits_full.detach().float()
    profile_bias = outputs.profile_bias.detach().float()
    count_full = outputs.logcount_full.detach().float()
    count_debiased = outputs.logcount_debiased.detach().float()

    profile_scale = raw_model._resolve_scale(raw_model.profile_scale, profile_signal).detach().float()
    count_scale = raw_model._resolve_scale(raw_model.count_scale, count_debiased).detach().float()
    effective_profile_bias = profile_scale * profile_bias

    centered_signal = profile_signal - profile_signal.mean(dim=-1, keepdim=True)
    centered_bias = effective_profile_bias - effective_profile_bias.mean(dim=-1, keepdim=True)
    signal_rms = torch.sqrt(torch.mean(centered_signal**2, dim=-1))
    bias_rms = torch.sqrt(torch.mean(centered_bias**2, dim=-1))

    diagnostics["effective_profile_scale"] = profile_scale.mean()
    diagnostics["effective_count_scale"] = count_scale.mean()
    diagnostics["profile_bias_rms_over_signal_rms"] = (
        bias_rms / torch.clamp(signal_rms, min=1e-8)
    ).mean()

    full_probs = torch.softmax(profile_full, dim=-1)
    debiased_probs = torch.softmax(profile_signal, dim=-1)
    diagnostics["profile_full_debiased_jsd"] = compute_profile_jsd(full_probs, debiased_probs).mean()
    diagnostics["count_full_debiased_abs"] = torch.mean(torch.abs(count_full - count_debiased))
    return diagnostics


GENOS_FILM_KEYS: Tuple[str, ...] = (
    "film_gamma_mean",
    "film_gamma_std",
    "film_beta_rms",
    "film_delta_rms",
    "genos_count_delta_rms",
    "caduceus_gate_mean",
    "caduceus_gate_std",
    "caduceus_delta_rms",
    "foundation_cross_gate_mean",
    "foundation_cross_gate_std",
    "foundation_cross_delta_rms",
    "foundation_residual_profile_rms",
    "foundation_residual_count_rms",
)

TRAIN_RUNTIME_KEYS: Tuple[str, ...] = (
    "step_time_total",
    "step_time_data_to_device",
    "step_time_genos_fetch",
    "step_time_forward_backward",
    "step_time_optimizer",
)


def compute_genos_film_diagnostics(model: nn.Module) -> Dict[str, float]:
    """Collect current-batch online-feature usage statistics.

    The helper keeps its legacy name because existing Genos cached logs already
    reference these fields.
    """
    raw_model = unwrap_model(model)
    diag: Dict[str, float] = {k: 0.0 for k in GENOS_FILM_KEYS}

    caduceus_adapter = getattr(raw_model, "caduceus_adapter", None)
    if caduceus_adapter is not None and hasattr(caduceus_adapter, "_last_forward_stats"):
        for key, value in getattr(caduceus_adapter, "_last_forward_stats", {}).items():
            if key in diag:
                diag[key] = to_scalar(value)

    film = getattr(raw_model, "genos_film", None)
    if film is not None and hasattr(film, "_last_forward_stats"):
        for key, value in getattr(film, "_last_forward_stats", {}).items():
            if key in diag:
                diag[key] = to_scalar(value)

    count_proj = getattr(raw_model, "genos_count_proj", None)
    if count_proj is not None and hasattr(count_proj, "_last_forward_stats"):
        for key, value in getattr(count_proj, "_last_forward_stats", {}).items():
            if key in diag:
                diag[key] = to_scalar(value)

    foundation_cross = getattr(raw_model, "foundation_cross_adapter", None)
    if foundation_cross is not None and hasattr(foundation_cross, "_last_forward_stats"):
        for key, value in getattr(foundation_cross, "_last_forward_stats", {}).items():
            if key in diag:
                diag[key] = to_scalar(value)

    foundation_residual = getattr(raw_model, "foundation_residual_head", None)
    if foundation_residual is not None and hasattr(foundation_residual, "_last_forward_stats"):
        for key, value in getattr(foundation_residual, "_last_forward_stats", {}).items():
            if key in diag:
                diag[key] = to_scalar(value)

    return diag


def format_bias_diag_suffix(metrics: Dict[str, float]) -> str:
    keys = [
        ("effective_profile_scale", "p_scale"),
        ("effective_count_scale", "c_scale"),
        ("profile_bias_rms_over_signal_rms", "p_ratio"),
        ("profile_full_debiased_jsd", "p_jsd"),
        ("count_full_debiased_abs", "c_shift"),
    ]
    parts = []
    for key, short_name in keys:
        if key in metrics:
            parts.append(f"{short_name}={metrics[key]:.4g}")
    return " ".join(parts)


def append_epoch_metrics(log_dir: Path, payload: Dict[str, Any]) -> Path:
    log_dir.mkdir(parents=True, exist_ok=True)
    path = log_dir / "epoch_metrics.jsonl"
    with open(path, "a", encoding="utf-8") as f:
        f.write(json.dumps(payload, sort_keys=True) + "\n")
    return path


def slice_outputs(outputs: TransChromBPOutput, mask: Tensor) -> TransChromBPOutput:
    return TransChromBPOutput(
        profile_logits_full=outputs.profile_logits_full[mask],
        logcount_full=outputs.logcount_full[mask],
        profile_logits_debiased=outputs.profile_logits_debiased[mask],
        logcount_debiased=outputs.logcount_debiased[mask],
        profile_bias=outputs.profile_bias[mask],
        count_bias=outputs.count_bias[mask],
    )


def reduce_weighted_metrics(
    metric_sums: Dict[str, float],
    total_weight: float,
    dist_env: DistEnv,
) -> Dict[str, float]:
    keys = list(metric_sums.keys())
    sum_tensor = torch.tensor(
        [metric_sums[key] for key in keys],
        device=dist_env.device,
        dtype=torch.float32,
    )
    weight_tensor = torch.tensor(float(total_weight), device=dist_env.device, dtype=torch.float32)
    if dist_env.distributed:
        dist.all_reduce(sum_tensor, op=dist.ReduceOp.SUM)
        dist.all_reduce(weight_tensor, op=dist.ReduceOp.SUM)
    denom = max(to_scalar(weight_tensor), 1.0)
    return {
        key: to_scalar(sum_tensor[idx]) / denom
        for idx, key in enumerate(keys)
    }


def reduce_sum_metrics(metric_sums: Dict[str, float], dist_env: DistEnv) -> Dict[str, float]:
    keys = list(metric_sums.keys())
    sum_tensor = torch.tensor(
        [metric_sums[key] for key in keys],
        device=dist_env.device,
        dtype=torch.float64,
    )
    if dist_env.distributed:
        dist.all_reduce(sum_tensor, op=dist.ReduceOp.SUM)
    return {
        key: float(sum_tensor[idx].item())
        for idx, key in enumerate(keys)
    }


def safe_pearson_from_sums(
    n: float,
    sum_x: float,
    sum_y: float,
    sum_x2: float,
    sum_y2: float,
    sum_xy: float,
) -> float:
    if n < 2.0:
        return float("nan")
    numerator = n * sum_xy - sum_x * sum_y
    denom_x = n * sum_x2 - sum_x * sum_x
    denom_y = n * sum_y2 - sum_y * sum_y
    denom = math.sqrt(max(denom_x, 0.0) * max(denom_y, 0.0))
    if denom <= 0.0:
        return float("nan")
    return float(numerator / denom)


def compute_selection_metric_sums(outputs: TransChromBPOutput, profile_counts: Tensor) -> Dict[str, Tensor]:
    counts = profile_counts.detach().double()
    target_total = counts.sum(dim=-1).reshape(-1)

    pred_log_full = outputs.logcount_full.detach().double().reshape(-1)
    pred_log_debiased = outputs.logcount_debiased.detach().double().reshape(-1)
    pred_total_full = torch.clamp(torch.expm1(pred_log_full), min=0.0)
    pred_total_debiased = torch.clamp(torch.expm1(pred_log_debiased), min=0.0)

    target_probs = counts + 1e-6
    target_probs = target_probs / target_probs.sum(dim=-1, keepdim=True)
    full_probs = torch.softmax(outputs.profile_logits_full.detach().double(), dim=-1)
    debiased_probs = torch.softmax(outputs.profile_logits_debiased.detach().double(), dim=-1)
    jsd_full = compute_profile_jsd(full_probs, target_probs)
    jsd_debiased = compute_profile_jsd(debiased_probs, target_probs)

    return {
        "n_examples": target_total.new_tensor(float(target_total.numel())),
        "target_total_sum": target_total.sum(),
        "target_total_sq_sum": torch.square(target_total).sum(),
        "pred_total_full_sum": pred_total_full.sum(),
        "pred_total_full_sq_sum": torch.square(pred_total_full).sum(),
        "target_pred_total_full_sum": (target_total * pred_total_full).sum(),
        "pred_total_debiased_sum": pred_total_debiased.sum(),
        "pred_total_debiased_sq_sum": torch.square(pred_total_debiased).sum(),
        "target_pred_total_debiased_sum": (target_total * pred_total_debiased).sum(),
        "count_mae_full_sum": torch.abs(pred_total_full - target_total).sum(),
        "count_mae_debiased_sum": torch.abs(pred_total_debiased - target_total).sum(),
        "profile_jsd_full_sum": jsd_full.double().sum(),
        "profile_jsd_debiased_sum": jsd_debiased.double().sum(),
    }


def finalize_selection_metrics(metric_sums: Dict[str, float]) -> Dict[str, float]:
    n = float(metric_sums.get("n_examples", 0.0))
    if n <= 0.0:
        return {}
    return {
        "count_pearson_full": safe_pearson_from_sums(
            n=n,
            sum_x=metric_sums["target_total_sum"],
            sum_y=metric_sums["pred_total_full_sum"],
            sum_x2=metric_sums["target_total_sq_sum"],
            sum_y2=metric_sums["pred_total_full_sq_sum"],
            sum_xy=metric_sums["target_pred_total_full_sum"],
        ),
        "count_pearson_debiased": safe_pearson_from_sums(
            n=n,
            sum_x=metric_sums["target_total_sum"],
            sum_y=metric_sums["pred_total_debiased_sum"],
            sum_x2=metric_sums["target_total_sq_sum"],
            sum_y2=metric_sums["pred_total_debiased_sq_sum"],
            sum_xy=metric_sums["target_pred_total_debiased_sum"],
        ),
        "count_mae_full": float(metric_sums["count_mae_full_sum"] / n),
        "count_mae_debiased": float(metric_sums["count_mae_debiased_sum"] / n),
        "profile_target_jsd_full_mean": float(metric_sums["profile_jsd_full_sum"] / n),
        "profile_target_jsd_debiased_mean": float(metric_sums["profile_jsd_debiased_sum"] / n),
    }


def format_metric_summary(metrics: Dict[str, float]) -> str:
    return " ".join(
        f"{k}={v:.5f}" if k in LOSS_METRIC_KEYS else f"{k}={v:.4g}"
        for k, v in metrics.items()
    )


def _metric_path_variants(path: str) -> list[str]:
    variants = [path]
    if "." in path:
        prefix, suffix = path.split(".", 1)
        if prefix == "overall":
            variants.append(suffix)
    else:
        variants.append(f"overall.{path}")
    return variants


def extract_metric_value(metrics: Dict[str, Any], path: str) -> Optional[float]:
    for candidate in _metric_path_variants(path):
        current: Any = metrics
        found = True
        for part in candidate.split("."):
            if not isinstance(current, dict) or part not in current:
                found = False
                break
            current = current[part]
        if found:
            return float(current)
    return None


def resolve_best_metric_name(
    trainer_cfg: Dict[str, Any],
    val_metrics: Dict[str, Any],
) -> str:
    configured = str(trainer_cfg.get("best_metric", "")).strip()
    if configured:
        return configured
    if isinstance(val_metrics, dict) and "peak" in val_metrics:
        return "peak.loss_total"
    if isinstance(val_metrics, dict) and "overall" in val_metrics:
        return "overall.loss_total"
    return "loss_total"


def is_better_metric(candidate: float, current_best: Optional[float], mode: str, min_delta: float = 0.0) -> bool:
    if current_best is None:
        return True
    if mode == "max":
        return candidate > (current_best + min_delta)
    return candidate < (current_best - min_delta)


def check_sequence_length_compatibility(model_cfg: Dict[str, Any], train_cfg: Dict[str, Any], rank: int) -> None:
    seq_cfg = model_cfg.get("sequence_encoder", {})
    data_cfg = train_cfg.get("data", {})
    trainer_cfg = train_cfg.get("trainer", {})

    model_max_len = int(seq_cfg.get("max_len", 0))
    input_len = int(data_cfg.get("input_len", 0))
    preprocessing_max_len = int(data_cfg.get("max_seq_len", input_len if input_len > 0 else 0))
    strict_check = bool(trainer_cfg.get("strict_max_len_check", False))

    if input_len > 0 and preprocessing_max_len > 0 and input_len > preprocessing_max_len:
        raise ValueError(
            f"data.input_len ({input_len}) cannot exceed data.max_seq_len ({preprocessing_max_len})"
        )

    if model_max_len > 0 and input_len > 0 and input_len > model_max_len:
        msg = (
            f"data.input_len ({input_len}) > sequence_encoder.max_len ({model_max_len}). "
            "RoPE cache can grow dynamically, but this may increase memory usage unexpectedly."
        )
        if strict_check:
            raise ValueError(msg)
        if is_main_process(rank):
            print(f"[warn] {msg}")

    if model_max_len > 0 and preprocessing_max_len > 0 and model_max_len != preprocessing_max_len:
        msg = (
            "sequence length mismatch: "
            f"sequence_encoder.max_len={model_max_len}, "
            f"data.max_seq_len={preprocessing_max_len}."
        )
        if strict_check:
            raise ValueError(msg)
        if is_main_process(rank):
            print(f"[warn] {msg}")



def check_profile_resolution_compatibility(model_cfg: Dict[str, Any], train_cfg: Dict[str, Any], rank: int) -> None:
    data_cfg = train_cfg.setdefault("data", {})
    heads_cfg = model_cfg.setdefault("heads", {})

    profile_bin_size = int(data_cfg.get("profile_bin_size", 1))
    if profile_bin_size <= 0:
        raise ValueError(f"data.profile_bin_size must be positive, got {profile_bin_size}")

    data_output_len = int(data_cfg.get("output_len", 0))
    model_output_len = int(heads_cfg.get("profile_output_len", 0))
    supervised_bp = int(data_cfg.get("supervised_bp", 0))

    if supervised_bp <= 0:
        base_output_len = data_output_len if data_output_len > 0 else model_output_len
        if base_output_len <= 0:
            raise ValueError(
                "Unable to infer data.supervised_bp: require either data.supervised_bp > 0 "
                "or one of data.output_len / heads.profile_output_len > 0"
            )
        supervised_bp = base_output_len * profile_bin_size
        data_cfg["supervised_bp"] = supervised_bp

    if supervised_bp % profile_bin_size != 0:
        raise ValueError(
            "data.supervised_bp must be divisible by data.profile_bin_size, got "
            f"{supervised_bp} and {profile_bin_size}"
        )

    derived_output_len = supervised_bp // profile_bin_size

    if data_output_len <= 0:
        data_cfg["output_len"] = derived_output_len
        data_output_len = derived_output_len
    if model_output_len <= 0:
        heads_cfg["profile_output_len"] = derived_output_len
        model_output_len = derived_output_len

    if data_output_len != derived_output_len:
        raise ValueError(
            "data.output_len mismatch with supervised_bp/profile_bin_size: "
            f"output_len={data_output_len}, expected={derived_output_len}"
        )

    if model_output_len != derived_output_len:
        raise ValueError(
            "heads.profile_output_len mismatch with supervised_bp/profile_bin_size: "
            f"profile_output_len={model_output_len}, expected={derived_output_len}"
        )

    if is_main_process(rank):
        input_len = int(data_cfg.get("input_len", 0))
        print(
            "[data] "
            f"input_len={input_len} "
            f"profile_bin_size={profile_bin_size} supervised_bp={supervised_bp} "
            f"output_len={derived_output_len}"
        )


def apply_external_data_config_defaults(train_cfg: Dict[str, Any], data_source_cfg: Dict[str, Any]) -> None:
    data_cfg = train_cfg.setdefault("data", {})
    window_cfg = data_source_cfg.get("window", {})
    sampling_cfg = data_source_cfg.get("sampling", {})

    if "input_len" not in data_cfg and "input_len" in window_cfg:
        data_cfg["input_len"] = int(window_cfg["input_len"])
    if "output_len" not in data_cfg and "output_len" in window_cfg:
        data_cfg["output_len"] = int(window_cfg["output_len"])
    if "max_seq_len" not in data_cfg and "input_len" in data_cfg:
        data_cfg["max_seq_len"] = int(data_cfg["input_len"])
    if "max_jitter" not in data_cfg and "max_jitter" in window_cfg:
        data_cfg["max_jitter"] = int(window_cfg["max_jitter"])
    if "nonpeak_ratio" not in data_cfg and "nonpeak_ratio" in sampling_cfg:
        data_cfg["nonpeak_ratio"] = float(sampling_cfg["nonpeak_ratio"])


def get_training_mode(train_cfg: Dict[str, Any]) -> str:
    return str(train_cfg.get("training_mode", "full")).lower()


def apply_training_mode_defaults(train_cfg: Dict[str, Any]) -> None:
    data_cfg = train_cfg.setdefault("data", {})
    training_mode = get_training_mode(train_cfg)
    if training_mode == "bias_pretrain":
        data_cfg.setdefault("region_source", "nonpeak")
        data_cfg.setdefault("train_region_source", "nonpeak")
        data_cfg.setdefault("val_region_source", "nonpeak")
        data_cfg.setdefault("nonpeak_ratio", 1.0)


def _semantics_value_matches(actual: Any, expected: Any) -> bool:
    if isinstance(expected, float):
        try:
            return math.isclose(float(actual), float(expected), rel_tol=0.0, abs_tol=1e-8)
        except (TypeError, ValueError):
            return False
    return actual == expected


def validate_semantics_profile(
    train_cfg: Dict[str, Any],
    data_source_cfg: Dict[str, Any],
    rank: int,
) -> None:
    semantics_profile = str(train_cfg.get("semantics_profile", "")).strip()
    if not semantics_profile:
        return

    normalized_profile = semantics_profile.lower()
    if normalized_profile != "chrombpnet_atac_canonical_v1":
        raise ValueError(
            "Unsupported semantics_profile="
            f"{semantics_profile!r}; expected 'chrombpnet_atac_canonical_v1'"
        )

    training_mode = get_training_mode(train_cfg)
    data_cfg = train_cfg.setdefault("data", {})
    loss_cfg = train_cfg.setdefault("loss", {})
    input_cfg = data_source_cfg.get("input", {})
    mismatches: list[str] = []

    def expect(path: str, actual: Any, expected: Any) -> None:
        if not _semantics_value_matches(actual, expected):
            mismatches.append(f"{path}={actual!r} (expected {expected!r})")

    expect("data.source", str(data_cfg.get("source", "")).strip().lower(), "chrombpnet_bigwig")
    expect(
        "loss.count_weight_strategy",
        str(loss_cfg.get("count_weight_strategy", "fixed")).strip().lower(),
        "chrombpnet_auto",
    )
    expect("data.train_revcomp", bool(data_cfg.get("train_revcomp", False)), True)
    expect("data.val_revcomp", bool(data_cfg.get("val_revcomp", False)), False)
    expect("data.revcomp_prob", float(data_cfg.get("revcomp_prob", 0.5)), 0.5)
    expect("data.track_total_count_target", float(data_cfg.get("track_total_count_target", 0.0)), 0.0)

    bigwig_path = str(input_cfg.get("bigwig", "")).strip()
    peaks_bed = str(input_cfg.get("peaks_bed", "")).strip()
    nonpeaks_bed = str(input_cfg.get("nonpeaks_bed", "")).strip()
    if not bigwig_path.endswith("_unstranded.bw"):
        mismatches.append(
            f"input.bigwig={bigwig_path!r} (expected shifted unstranded bigWig ending with '_unstranded.bw')"
        )
    if not peaks_bed:
        mismatches.append("input.peaks_bed is empty")
    if not nonpeaks_bed:
        mismatches.append("input.nonpeaks_bed is empty")

    if training_mode == "bias_pretrain":
        expect("data.region_source", str(data_cfg.get("region_source", "")).strip().lower(), "nonpeak")
        expect("data.train_region_source", str(data_cfg.get("train_region_source", "")).strip().lower(), "nonpeak")
        expect("data.val_region_source", str(data_cfg.get("val_region_source", "")).strip().lower(), "nonpeak")
        expect("data.nonpeak_ratio", float(data_cfg.get("nonpeak_ratio", 0.0)), 1.0)
        expect("data.max_jitter", int(data_cfg.get("max_jitter", -1)), 0)
        expect("data.peak_max_jitter", int(data_cfg.get("peak_max_jitter", -1)), 0)
        expect("data.nonpeak_max_jitter", int(data_cfg.get("nonpeak_max_jitter", -1)), 0)
    elif training_mode == "full":
        expect("data.region_source", str(data_cfg.get("region_source", "")).strip().lower(), "both")
        expect("data.train_region_source", str(data_cfg.get("train_region_source", "")).strip().lower(), "both")
        expect("data.val_region_source", str(data_cfg.get("val_region_source", "")).strip().lower(), "both")
        expect("data.nonpeak_ratio", float(data_cfg.get("nonpeak_ratio", 0.0)), 0.1)
        expect("data.max_jitter", int(data_cfg.get("max_jitter", -1)), 500)
        expect("data.peak_max_jitter", int(data_cfg.get("peak_max_jitter", -1)), 500)
        expect("data.nonpeak_max_jitter", int(data_cfg.get("nonpeak_max_jitter", -1)), 0)
    else:
        mismatches.append(
            f"training_mode={training_mode!r} (expected 'full' or 'bias_pretrain' for canonical semantics)"
        )

    if mismatches:
        joined = "; ".join(mismatches)
        raise ValueError(
            "semantics_profile='chrombpnet_atac_canonical_v1' validation failed: "
            f"{joined}"
        )

    if is_main_process(rank):
        print(
            "[config] semantics_profile=chrombpnet_atac_canonical_v1 "
            f"validated for training_mode={training_mode}"
        )


def resolve_data_config_path(train_cfg: Dict[str, Any], cli_data_config: str) -> str:
    if cli_data_config:
        return cli_data_config
    data_cfg = train_cfg.get("data", {})
    return str(data_cfg.get("config_path", ""))


def _resolve_real_data_value(
    train_data_cfg: Dict[str, Any],
    data_source_cfg: Dict[str, Any],
    key: str,
    default: Any = None,
) -> Any:
    if key in train_data_cfg and train_data_cfg[key] not in ("", None):
        return train_data_cfg[key]
    if key in data_source_cfg and data_source_cfg[key] not in ("", None):
        return data_source_cfg[key]

    for section in ("input", "window", "sampling"):
        section_cfg = data_source_cfg.get(section, {})
        if key in section_cfg and section_cfg[key] not in ("", None):
            return section_cfg[key]
    return default


def auto_configure_loss_weights(
    train_cfg: Dict[str, Any],
    train_dataset: Dataset,
    dist_env: DistEnv,
) -> Dict[str, float]:
    loss_cfg = train_cfg.setdefault("loss", {})
    strategy = str(loss_cfg.get("count_weight_strategy", "fixed")).strip().lower()
    if strategy in {"", "fixed", "manual", "none"}:
        return {}
    if strategy not in {"chrombpnet_auto", "median_over_10"}:
        raise ValueError(
            "Unsupported loss.count_weight_strategy="
            f"{strategy!r}; expected one of fixed/manual/none/chrombpnet_auto"
        )
    if not hasattr(train_dataset, "estimate_total_count_statistics"):
        raise ValueError(
            "Automatic count-weight calibration requires the dataset to expose "
            "`estimate_total_count_statistics(...)`"
        )

    sample_size = int(loss_cfg.get("count_weight_sample_size", 32768))
    quantile_low = float(loss_cfg.get("count_weight_quantile_low", 1e-4))
    quantile_high = float(loss_cfg.get("count_weight_quantile_high", 0.9999))
    divisor = float(loss_cfg.get("count_weight_divisor", 10.0))
    floor = float(loss_cfg.get("count_weight_floor", 1.0))
    if divisor <= 0.0:
        raise ValueError(f"loss.count_weight_divisor must be positive, got {divisor}")
    if floor <= 0.0:
        raise ValueError(f"loss.count_weight_floor must be positive, got {floor}")

    payload = torch.zeros(8, device=dist_env.device, dtype=torch.float32)
    if is_main_process(dist_env.rank):
        stats = train_dataset.estimate_total_count_statistics(
            sample_size=sample_size,
            seed=int(train_cfg.get("seed", 1234)) + 17_291,
            quantile_low=quantile_low,
            quantile_high=quantile_high,
        )
        recommended = max(floor, float(stats["median_total_count"]) / divisor)
        payload = torch.tensor(
            [
                recommended,
                float(stats["median_total_count"]),
                float(stats["mean_total_count"]),
                float(stats["lower_threshold"]),
                float(stats["upper_threshold"]),
                float(stats["sample_size"]),
                float(stats["retained_size"]),
                float(stats.get("track_scale_factor", 1.0)),
            ],
            device=dist_env.device,
            dtype=torch.float32,
        )
    if dist_env.distributed:
        dist.broadcast(payload, src=0)

    recommended, median_total, mean_total, lower_threshold, upper_threshold, used_n, retained_n, track_scale = (
        float(x) for x in payload.tolist()
    )
    loss_cfg["count_weight"] = recommended
    train_cfg["loss"] = loss_cfg

    calibration = {
        "strategy": strategy,
        "count_weight": recommended,
        "median_total_count": median_total,
        "mean_total_count": mean_total,
        "lower_threshold": lower_threshold,
        "upper_threshold": upper_threshold,
        "sample_size": used_n,
        "retained_size": retained_n,
        "count_weight_divisor": divisor,
        "count_weight_floor": floor,
        "track_scale_factor": track_scale,
    }
    if is_main_process(dist_env.rank):
        print(
            "[loss] "
            f"count_weight_strategy={strategy} count_weight={recommended:.4f} "
            f"median_total={median_total:.4f} mean_total={mean_total:.4f} "
            f"trim=[{lower_threshold:.4f}, {upper_threshold:.4f}] "
            f"used={int(round(used_n))} retained={int(round(retained_n))} "
            f"track_scale={track_scale:.6g}"
        )
    return calibration


class SyntheticChromDataset(Dataset):
    """Deterministic synthetic dataset for training pipeline validation."""

    def __init__(
        self,
        n_samples: int,
        input_len: int,
        output_len: int,
        seed: int = 1234,
        count_scale: float = 20.0,
    ) -> None:
        self.n_samples = int(n_samples)
        self.input_len = int(input_len)
        self.output_len = int(output_len)
        self.seed = int(seed)
        self.count_scale = float(count_scale)

    def __len__(self) -> int:
        return self.n_samples

    def __getitem__(self, idx: int) -> Dict[str, Tensor]:
        g = torch.Generator().manual_seed(self.seed + idx)

        bases = torch.randint(0, 4, (self.input_len,), generator=g)
        seq = F.one_hot(bases, num_classes=4).to(dtype=torch.float32)

        lam = torch.rand(self.output_len, generator=g) * self.count_scale + 1.0
        profile_counts = torch.poisson(lam)

        return {
            "seq": seq,
            "profile_counts": profile_counts.to(dtype=torch.float32),
        }


class BiasOnlyWrapper(nn.Module):
    """Train only the bias branch while keeping the loss interface unchanged."""

    def __init__(self, model_cfg: Dict[str, Any]) -> None:
        super().__init__()
        self.bias_branch = build_bias_branch_from_config(model_cfg)

    def forward(self, seq_onehot: Tensor) -> TransChromBPOutput:
        profile_bias, count_bias = self.bias_branch(seq_onehot)
        return TransChromBPOutput(
            profile_logits_full=profile_bias,
            logcount_full=count_bias,
            profile_logits_debiased=profile_bias,
            logcount_debiased=count_bias,
            profile_bias=profile_bias,
            count_bias=count_bias,
        )


def build_dataloaders(
    train_cfg: Dict[str, Any],
    dist_env: DistEnv,
    data_source_cfg: Optional[Dict[str, Any]] = None,
) -> tuple[DataLoader, Optional[DataLoader], Optional[DistributedSampler]]:
    data_cfg = train_cfg.get("data", {})
    source = str(data_cfg.get("source", "synthetic")).lower()

    input_len = int(data_cfg.get("input_len", 2114))
    output_len = int(data_cfg.get("output_len", 1000))
    batch_size = int(data_cfg.get("batch_size_per_gpu", 8))
    num_workers = int(data_cfg.get("num_workers", 2))
    pin_memory = bool(data_cfg.get("pin_memory", True))
    prefetch_factor = int(data_cfg.get("prefetch_factor", 2))
    persistent_workers = bool(data_cfg.get("persistent_workers", True))
    seed = int(train_cfg.get("seed", 1234))
    if source == "synthetic":
        train_samples = int(data_cfg.get("train_samples", 8192))
        val_samples = int(data_cfg.get("val_samples", 1024))
        count_scale = float(data_cfg.get("synthetic_count_scale", 20.0))

        train_ds = SyntheticChromDataset(
            n_samples=train_samples,
            input_len=input_len,
            output_len=output_len,
            seed=seed,
            count_scale=count_scale,
        )
        val_ds = SyntheticChromDataset(
            n_samples=val_samples,
            input_len=input_len,
            output_len=output_len,
            seed=seed + 10_000,
            count_scale=count_scale,
        )
    elif source in {"chrombpnet_bigwig", "bigwig"}:
        data_source_cfg = data_source_cfg or {}
        supervised_bp = int(data_cfg.get("supervised_bp", output_len))
        profile_bin_size = int(data_cfg.get("profile_bin_size", 1))
        max_jitter = int(_resolve_real_data_value(data_cfg, data_source_cfg, "max_jitter", 0))
        peak_max_jitter = int(_resolve_real_data_value(data_cfg, data_source_cfg, "peak_max_jitter", max_jitter))
        nonpeak_max_jitter = int(
            _resolve_real_data_value(data_cfg, data_source_cfg, "nonpeak_max_jitter", max_jitter)
        )
        nonpeak_ratio = float(_resolve_real_data_value(data_cfg, data_source_cfg, "nonpeak_ratio", 1.0))
        train_revcomp = bool(_resolve_real_data_value(data_cfg, data_source_cfg, "train_revcomp", False))
        val_revcomp = bool(_resolve_real_data_value(data_cfg, data_source_cfg, "val_revcomp", False))
        revcomp_prob = float(_resolve_real_data_value(data_cfg, data_source_cfg, "revcomp_prob", 0.5))
        track_total_count_target = float(
            _resolve_real_data_value(data_cfg, data_source_cfg, "track_total_count_target", 0.0)
        )
        max_train_regions = int(data_cfg.get("max_train_regions", 0))
        max_val_regions = int(data_cfg.get("max_val_regions", 0))
        default_region_source = str(data_cfg.get("region_source", "both"))
        train_region_source = str(data_cfg.get("train_region_source", default_region_source))
        val_region_source = str(data_cfg.get("val_region_source", default_region_source))

        genos_cache_dir = str(data_cfg.get("genos_cache_dir", ""))
        genos_cache_features = list(data_cfg.get("genos_cache_features", []))
        foundation_cache_dir = str(data_cfg.get("foundation_cache_dir", ""))
        foundation_cache_features = list(data_cfg.get("foundation_cache_features", []))
        teacher_cache_dir = str(data_cfg.get("teacher_cache_dir", ""))
        teacher_target_names = list(data_cfg.get("teacher_target_names", []))

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
        missing = [key for key, value in required.items() if not value]
        if missing:
            raise ValueError(
                "Real-data backend requires the following paths: "
                + ", ".join(missing)
            )

        train_ds = ChromBPNetBigWigDataset(
            genome_fasta=genome_fasta,
            bigwig_path=bigwig_path,
            peaks_bed=peaks_bed,
            nonpeaks_bed=nonpeaks_bed,
            folds_json=folds_json,
            split="train",
            input_len=input_len,
            supervised_bp=supervised_bp,
            profile_bin_size=profile_bin_size,
            max_jitter=max_jitter,
            peak_max_jitter=peak_max_jitter,
            nonpeak_max_jitter=nonpeak_max_jitter,
            seed=seed,
            nonpeak_ratio=nonpeak_ratio,
            max_records=max_train_regions,
            region_source=train_region_source,
            random_revcomp=train_revcomp,
            revcomp_prob=revcomp_prob,
            track_total_count_target=track_total_count_target,
            genos_cache_dir=genos_cache_dir,
            genos_cache_features=genos_cache_features,
            foundation_cache_dir=foundation_cache_dir,
            foundation_cache_features=foundation_cache_features,
            teacher_cache_dir=teacher_cache_dir,
            teacher_target_names=teacher_target_names,
        )
        val_ds = ChromBPNetBigWigDataset(
            genome_fasta=genome_fasta,
            bigwig_path=bigwig_path,
            peaks_bed=peaks_bed,
            nonpeaks_bed=nonpeaks_bed,
            folds_json=folds_json,
            split="valid",
            input_len=input_len,
            supervised_bp=supervised_bp,
            profile_bin_size=profile_bin_size,
            max_jitter=0,
            peak_max_jitter=0,
            nonpeak_max_jitter=0,
            seed=seed + 10_000,
            nonpeak_ratio=nonpeak_ratio,
            max_records=max_val_regions,
            region_source=val_region_source,
            random_revcomp=val_revcomp,
            revcomp_prob=revcomp_prob,
            track_total_count_target=track_total_count_target,
            genos_cache_dir=genos_cache_dir,
            genos_cache_features=genos_cache_features,
            foundation_cache_dir=foundation_cache_dir,
            foundation_cache_features=foundation_cache_features,
            teacher_cache_dir=teacher_cache_dir,
            teacher_target_names=teacher_target_names,
        )
    else:
        raise ValueError(f"Unsupported data.source: {source}")

    if getattr(train_ds, "supports_epoch_resampling", False):
        train_sampler = PeakNonpeakResamplingSampler(
            train_ds,
            num_replicas=dist_env.world_size,
            rank=dist_env.rank,
            shuffle=True,
            seed=seed,
            drop_last=False,
        )
        if dist_env.distributed:
            val_sampler = DistributedSampler(val_ds, num_replicas=dist_env.world_size, rank=dist_env.rank, shuffle=False)
        else:
            val_sampler = None
    elif dist_env.distributed:
        train_sampler = DistributedSampler(train_ds, num_replicas=dist_env.world_size, rank=dist_env.rank, shuffle=True)
        val_sampler = DistributedSampler(val_ds, num_replicas=dist_env.world_size, rank=dist_env.rank, shuffle=False)
    else:
        train_sampler = None
        val_sampler = None

    if source in {"chrombpnet_bigwig", "bigwig"} and is_main_process(dist_env.rank):
        effective_train_regions = int(getattr(train_sampler, "epoch_size", len(train_ds)))
        print(
            f"[data] source={source} train_regions={len(train_ds)} train_epoch_regions={effective_train_regions} "
            f"val_regions={len(val_ds)} "
            f"nonpeak_ratio={nonpeak_ratio} train_region_source={train_region_source} "
            f"val_region_source={val_region_source} peak_jitter={peak_max_jitter} "
            f"nonpeak_jitter={nonpeak_max_jitter} train_revcomp={train_revcomp} "
            f"track_total_count_target={track_total_count_target}"
        )

    train_loader_kwargs: Dict[str, Any] = dict(
        batch_size=batch_size,
        num_workers=num_workers,
        pin_memory=pin_memory,
        drop_last=True,
    )
    if num_workers > 0:
        train_loader_kwargs["prefetch_factor"] = prefetch_factor
        train_loader_kwargs["persistent_workers"] = persistent_workers

    val_loader_kwargs = dict(train_loader_kwargs)
    val_loader_kwargs["drop_last"] = False

    train_loader = DataLoader(
        train_ds,
        sampler=train_sampler,
        shuffle=(train_sampler is None),
        **train_loader_kwargs,
    )

    val_loader = DataLoader(
        val_ds,
        sampler=val_sampler,
        shuffle=False,
        **val_loader_kwargs,
    )

    return train_loader, val_loader, train_sampler


def profile_multinomial_nll(logits: Tensor, counts: Tensor, eps: float = 1e-6) -> Tensor:
    if logits.shape != counts.shape:
        raise ValueError(f"profile shape mismatch: logits {tuple(logits.shape)} vs counts {tuple(counts.shape)}")

    total_counts = counts.sum(dim=-1, keepdim=True)
    safe_totals = torch.where(total_counts > 0, total_counts, torch.ones_like(total_counts))
    target = counts / safe_totals
    log_probs = torch.log_softmax(logits, dim=-1)
    per_example = -(target * log_probs).sum(dim=-1)
    nonzero_mask = (total_counts.squeeze(-1) > 0).to(dtype=per_example.dtype)
    per_example = per_example * nonzero_mask
    return per_example.mean()


def _pool_profile_probs_to_bins(probs: Tensor, n_bins: int) -> Tensor:
    if probs.dim() != 2:
        raise ValueError(f"Expected probs [B, L], got {tuple(probs.shape)}")
    pooled = F.adaptive_avg_pool1d(probs.unsqueeze(1), output_size=int(n_bins)).squeeze(1)
    pooled = torch.clamp(pooled, min=1e-8)
    return pooled / pooled.sum(dim=-1, keepdim=True)


def standardized_mse_loss(pred: Tensor, target: Tensor) -> Tensor:
    pred = pred.reshape(pred.size(0), -1)
    target = target.reshape(target.size(0), -1)
    pred = (pred - pred.mean(dim=0, keepdim=True)) / torch.clamp(pred.std(dim=0, keepdim=True, unbiased=False), min=1e-6)
    target = (target - target.mean(dim=0, keepdim=True)) / torch.clamp(target.std(dim=0, keepdim=True, unbiased=False), min=1e-6)
    return F.mse_loss(pred, target)


def compute_losses(
    outputs: TransChromBPOutput,
    profile_counts: Tensor,
    loss_cfg: Dict[str, Any],
    teacher_targets: Optional[Dict[str, Tensor]] = None,
) -> tuple[Tensor, Dict[str, Tensor]]:
    profile_w = float(loss_cfg.get("profile_weight", 1.0))
    count_w = float(loss_cfg.get("count_weight", 0.1))
    debiased_profile_w = float(loss_cfg.get("debiased_profile_weight", 0.0))
    debiased_count_w = float(loss_cfg.get("debiased_count_weight", 0.0))
    distill_profile_w = float(loss_cfg.get("distill_profile_weight", 0.0))
    distill_count_w = float(loss_cfg.get("distill_count_weight", 0.0))
    distill_rank_w = float(loss_cfg.get("distill_rank_weight", 0.0))

    target_logcount = torch.log1p(profile_counts.sum(dim=-1, keepdim=True))

    profile_loss = profile_multinomial_nll(outputs.profile_logits_full, profile_counts)
    count_loss = F.mse_loss(outputs.logcount_full, target_logcount)

    debiased_profile_loss = profile_multinomial_nll(outputs.profile_logits_debiased, profile_counts)
    debiased_count_loss = F.mse_loss(outputs.logcount_debiased, target_logcount)

    total = profile_w * profile_loss + count_w * count_loss
    if debiased_profile_w > 0:
        total = total + debiased_profile_w * debiased_profile_loss
    if debiased_count_w > 0:
        total = total + debiased_count_w * debiased_count_loss

    zero = profile_loss.detach().new_tensor(0.0)
    distill_profile_loss = zero
    distill_count_loss = zero
    distill_rank_loss = zero

    teacher_targets = teacher_targets or {}
    if distill_profile_w > 0 and "profile16" in teacher_targets:
        teacher_profile = teacher_targets["profile16"].to(device=outputs.profile_logits_debiased.device, dtype=torch.float32)
        teacher_profile = torch.clamp(teacher_profile, min=1e-8)
        teacher_profile = teacher_profile / teacher_profile.sum(dim=-1, keepdim=True)
        pred_profile = torch.softmax(outputs.profile_logits_debiased, dim=-1)
        pred_profile = _pool_profile_probs_to_bins(pred_profile, teacher_profile.size(-1))
        distill_profile_loss = F.mse_loss(pred_profile, teacher_profile)
        total = total + distill_profile_w * distill_profile_loss
    if distill_count_w > 0 and "logcount" in teacher_targets:
        teacher_logcount = teacher_targets["logcount"].to(device=outputs.logcount_debiased.device, dtype=torch.float32)
        distill_count_loss = F.mse_loss(outputs.logcount_debiased, teacher_logcount.reshape_as(outputs.logcount_debiased))
        total = total + distill_count_w * distill_count_loss
    if distill_rank_w > 0 and "rank" in teacher_targets:
        teacher_rank = teacher_targets["rank"].to(device=outputs.logcount_debiased.device, dtype=torch.float32)
        distill_rank_loss = standardized_mse_loss(
            outputs.logcount_debiased,
            teacher_rank.reshape_as(outputs.logcount_debiased),
        )
        total = total + distill_rank_w * distill_rank_loss

    metrics = {
        "loss_total": total.detach(),
        "loss_profile": profile_loss.detach(),
        "loss_count": count_loss.detach(),
        "loss_debiased_profile": debiased_profile_loss.detach(),
        "loss_debiased_count": debiased_count_loss.detach(),
        "loss_distill_profile": distill_profile_loss.detach(),
        "loss_distill_count": distill_count_loss.detach(),
        "loss_distill_rank": distill_rank_loss.detach(),
    }
    return total, metrics


def build_model(model_cfg: Dict[str, Any], train_cfg: Dict[str, Any], rank: int) -> tuple[nn.Module, str]:
    training_mode = get_training_mode(train_cfg)
    if training_mode == "full":
        model = build_transchrombp_from_config(model_cfg)
    elif training_mode == "bias_pretrain":
        model = BiasOnlyWrapper(model_cfg)
    else:
        raise ValueError(f"Unsupported training_mode={training_mode!r}; expected 'full' or 'bias_pretrain'")

    if is_main_process(rank):
        print(f"[model] training_mode={training_mode}")
    return model, training_mode


def build_online_feature_runtime(
    model_cfg: Dict[str, Any],
    device: torch.device,
    rank: int,
) -> Optional[OnlineFeatureRuntime]:
    """Build a frozen online feature extractor if an online foundation branch is enabled."""
    genos_cfg = model_cfg.get("genos_branch", {})
    caduceus_cfg = model_cfg.get("caduceus_branch", {})
    genos_enabled = bool(genos_cfg.get("enabled", False))
    caduceus_enabled = bool(caduceus_cfg.get("enabled", False))

    if genos_enabled and caduceus_enabled:
        raise ValueError("genos_branch.enabled and caduceus_branch.enabled cannot both be true")

    if genos_enabled:
        from transchrombp.models.genos_adapter import GenosFeatureExtractor

        extractor = GenosFeatureExtractor(
            model_path=genos_cfg["model_path"],
            layer=int(genos_cfg.get("layer", 6)),
            device=device,
            attn_implementation=str(genos_cfg.get("attn_implementation", "flash_attention_2")),
            bidirectional=bool(genos_cfg.get("bidirectional", True)),
        )
        if is_main_process(rank):
            print(
                f"[genos] Loaded Genos extractor: layer={extractor.layer} "
                f"hidden_size={extractor.hidden_size} bidir={extractor.bidirectional}"
            )
        return OnlineFeatureRuntime(kind="genos", extractor=extractor)

    if caduceus_enabled:
        from transchrombp.models.caduceus_adapter import CaduceusFeatureExtractor

        extractor = CaduceusFeatureExtractor(
            model_path=str(caduceus_cfg["model_path"]),
            layer=int(caduceus_cfg.get("layer", -1)),
            device=device,
            trust_remote_code=bool(caduceus_cfg.get("trust_remote_code", True)),
            local_files_only=bool(caduceus_cfg.get("local_files_only", True)),
            dtype=str(caduceus_cfg.get("dtype", "bfloat16")),
        )
        cfg_hidden_size = caduceus_cfg.get("hidden_size")
        if cfg_hidden_size is not None and int(cfg_hidden_size) != extractor.hidden_size:
            raise ValueError(
                "caduceus_branch.hidden_size mismatch: "
                f"config={int(cfg_hidden_size)} extractor={extractor.hidden_size}. "
                "For RCPS Caduceus models the effective token width is typically 2x d_model."
            )
        if is_main_process(rank):
            print(
                f"[caduceus] Loaded Caduceus extractor: layer={extractor.layer} "
                f"hidden_size={extractor.hidden_size} rcps={extractor.rcps} rc_equivariant=true"
            )
        return OnlineFeatureRuntime(kind="caduceus", extractor=extractor)

    return None


def extract_online_feature_kwargs(
    seq: Tensor,
    online_runtime: Optional[OnlineFeatureRuntime],
) -> Dict[str, Tensor]:
    if online_runtime is None:
        return {}

    feat = online_runtime.extractor.extract(seq)
    if online_runtime.kind == "genos":
        return {"genos_feat": feat}
    if online_runtime.kind == "caduceus":
        return {"caduceus_feat": feat}
    raise ValueError(f"Unsupported online runtime kind={online_runtime.kind!r}")


def extract_foundation_feature_kwargs(
    batch: Dict[str, Any],
    model_cfg: Dict[str, Any],
    device: torch.device,
) -> Dict[str, Tensor]:
    contract = resolve_foundation_contract(model_cfg)
    if contract is None:
        return {}

    key = foundation_batch_key(contract.feature_name)
    if key not in batch:
        raise RuntimeError(
            f"Foundation mode requires batch field {key!r}; check data.foundation_cache_dir / "
            "data.foundation_cache_features."
        )

    feat = batch[key].to(device, non_blocking=True)
    kwargs: Dict[str, Tensor] = {}

    if contract.feature_layout == "summary":
        kwargs["foundation_summary"] = feat
        return kwargs

    if feat.dim() == 3:
        foundation_tokens = feat
    else:
        expected_dim = contract.feature_tokens * contract.hidden_size
        if feat.dim() != 2 or int(feat.shape[1]) != expected_dim:
            raise ValueError(
                f"Foundation token feature {contract.feature_name!r} has shape {tuple(feat.shape)}, expected "
                f"[B, {expected_dim}] for {contract.feature_tokens} tokens x hidden_size {contract.hidden_size}"
            )
        foundation_tokens = feat.view(feat.size(0), contract.feature_tokens, contract.hidden_size)

    kwargs["foundation_tokens"] = foundation_tokens
    if contract.mode == "residual_head":
        kwargs["foundation_summary"] = foundation_tokens.mean(dim=1)
    return kwargs


def extract_teacher_targets(batch: Dict[str, Any], device: torch.device) -> Dict[str, Tensor]:
    out: Dict[str, Tensor] = {}
    for key, value in batch.items():
        if not key.startswith("teacher_"):
            continue
        out[key.split("teacher_", 1)[1]] = value.to(device, non_blocking=True)
    return out


def validate_model_config(model_cfg: Dict[str, Any]) -> None:
    heads_cfg = model_cfg.get("heads", {})
    count_head_type = str(heads_cfg.get("count_head", "linear")).lower()
    if count_head_type != "linear":
        raise ValueError(f"Unsupported heads.count_head={count_head_type!r}; expected 'linear'")
    count_pool_mode = str(heads_cfg.get("count_pool_mode", "full")).lower()
    if count_pool_mode not in {"full", "center", "attention"}:
        raise ValueError(
            f"Unsupported heads.count_pool_mode={count_pool_mode!r}; "
            "expected 'full', 'center', or 'attention'"
        )
    if bool(model_cfg.get("genos_branch", {}).get("enabled", False)) and bool(
        model_cfg.get("caduceus_branch", {}).get("enabled", False)
    ):
        raise ValueError("genos_branch.enabled and caduceus_branch.enabled cannot both be true")
    enabled_count = sum(
        int(bool(model_cfg.get(name, {}).get("enabled", False)))
        for name in ("genos_branch", "caduceus_branch", "genos_cached", "foundation_model")
    )
    if enabled_count > 1:
        raise ValueError(
            "At most one of genos_branch, caduceus_branch, genos_cached, foundation_model may be enabled"
        )


def freeze_stopgrad_bias_profile_head(model: nn.Module, model_cfg: Dict[str, Any], rank: int) -> None:
    """Freeze bias profile-head params when profile bias is explicitly detached.

    With `fusion.profile_bias_stop_gradient=true`, `bias_branch.profile_head` does not
    contribute gradients to the loss. Single-rank training tolerates this silently,
    but multi-rank DDP treats these permanently-unused trainable params as an error.
    Freezing them here keeps semantics aligned with the configured stop-gradient path
    and lets optimizer/DDP parameter groups reflect reality.
    """
    fusion_cfg = model_cfg.get("fusion", {})
    if not bool(fusion_cfg.get("profile_bias_stop_gradient", False)):
        return

    bias_branch = getattr(model, "bias_branch", None)
    profile_head = getattr(bias_branch, "profile_head", None) if bias_branch is not None else None
    if profile_head is None:
        return

    frozen: list[str] = []
    for name, param in profile_head.named_parameters():
        if param.requires_grad:
            param.requires_grad_(False)
            frozen.append(f"bias_branch.profile_head.{name}")

    if frozen and is_main_process(rank):
        print(
            "[model] profile_bias_stop_gradient=true -> froze bias profile head params: "
            + ", ".join(frozen)
        )


def _set_module_requires_grad(module: Optional[nn.Module], requires_grad: bool) -> None:
    if module is None:
        return
    for param in module.parameters():
        param.requires_grad_(requires_grad)


def apply_foundation_freeze_policy(model: nn.Module, model_cfg: Dict[str, Any], epoch: int, rank: int) -> None:
    foundation_cfg = model_cfg.get("foundation_model", {})
    if not bool(foundation_cfg.get("enabled", False)):
        return
    if str(foundation_cfg.get("mode", "")).strip().lower() != "residual_head":
        return

    freeze_until = int(foundation_cfg.get("freeze_backbone_until_epoch", 0))
    unfreeze_last_n = int(foundation_cfg.get("unfreeze_transformer_last_n_blocks", 0))
    raw_model = unwrap_model(model)

    foundation_head = getattr(raw_model, "foundation_residual_head", None)
    _set_module_requires_grad(foundation_head, True)

    if freeze_until > 0 and epoch <= freeze_until:
        _set_module_requires_grad(getattr(raw_model, "conv_stem", None), False)
        _set_module_requires_grad(getattr(raw_model, "local_tower", None), False)
        _set_module_requires_grad(getattr(raw_model, "transformer", None), False)
        _set_module_requires_grad(getattr(raw_model, "profile_signal_head", None), False)
        _set_module_requires_grad(getattr(raw_model, "count_signal_head", None), False)
        if is_main_process(rank):
            print(f"[foundation] epoch={epoch}: residual head only")
        return

    _set_module_requires_grad(getattr(raw_model, "profile_signal_head", None), True)
    _set_module_requires_grad(getattr(raw_model, "count_signal_head", None), True)
    _set_module_requires_grad(getattr(raw_model, "conv_stem", None), False)
    _set_module_requires_grad(getattr(raw_model, "local_tower", None), False)
    transformer = getattr(raw_model, "transformer", None)
    _set_module_requires_grad(transformer, False)
    if transformer is not None and hasattr(transformer, "layers"):
        layers = list(getattr(transformer, "layers"))
        if unfreeze_last_n <= 0 or unfreeze_last_n >= len(layers):
            _set_module_requires_grad(transformer, True)
        else:
            for layer in layers[-unfreeze_last_n:]:
                _set_module_requires_grad(layer, True)
            if hasattr(transformer, "final_norm"):
                _set_module_requires_grad(getattr(transformer, "final_norm"), True)
    if is_main_process(rank):
        print(
            f"[foundation] epoch={epoch}: residual head + last {max(unfreeze_last_n, 0)} transformer blocks"
        )


def resolve_amp_dtype(dtype_name: str) -> torch.dtype:
    name = dtype_name.lower()
    if name == "bf16":
        return torch.bfloat16
    if name == "fp16":
        return torch.float16
    if name in {"fp32", "float32"}:
        return torch.float32
    raise ValueError(f"Unsupported precision dtype: {dtype_name}")


def build_grad_scaler(device_type: str, enabled: bool) -> Any:
    # torch>=2.3 exposes torch.amp.GradScaler(device=...), while torch 2.2 still
    # uses torch.cuda.amp.GradScaler without the device argument.
    amp_grad_scaler = getattr(torch.amp, "GradScaler", None)
    if amp_grad_scaler is not None:
        return amp_grad_scaler(device=device_type, enabled=enabled)
    if device_type == "cuda":
        return torch.cuda.amp.GradScaler(enabled=enabled)
    return torch.cuda.amp.GradScaler(enabled=False)


def iter_named_params_without_duplicates(model: nn.Module) -> Iterable[tuple[str, nn.Module, str, nn.Parameter]]:
    seen: set[int] = set()
    for module_name, module in model.named_modules():
        for param_name, param in module.named_parameters(recurse=False):
            if not param.requires_grad:
                continue
            pid = id(param)
            if pid in seen:
                continue
            seen.add(pid)
            full_name = f"{module_name}.{param_name}" if module_name else param_name
            yield full_name, module, param_name, param


def split_decay_param_groups(model: nn.Module) -> tuple[list[nn.Parameter], list[nn.Parameter], list[str], list[str]]:
    decay_params: list[nn.Parameter] = []
    no_decay_params: list[nn.Parameter] = []
    decay_names: list[str] = []
    no_decay_names: list[str] = []

    for full_name, module, param_name, param in iter_named_params_without_duplicates(model):
        is_bias_param = param_name.endswith("bias")
        is_norm_param = isinstance(module, NORM_LAYER_TYPES)
        is_scalar_or_vector = param.ndim <= 1

        if is_bias_param or is_norm_param or is_scalar_or_vector:
            no_decay_params.append(param)
            no_decay_names.append(full_name)
        else:
            decay_params.append(param)
            decay_names.append(full_name)

    return decay_params, no_decay_params, decay_names, no_decay_names


def build_optimizer(model: nn.Module, train_cfg: Dict[str, Any], rank: int) -> torch.optim.Optimizer:
    opt_cfg = train_cfg.get("optimizer", {})
    lr = float(opt_cfg.get("learning_rate", 5e-4))
    weight_decay = float(opt_cfg.get("weight_decay", 0.01))
    betas = opt_cfg.get("betas", [0.9, 0.95])

    name = str(opt_cfg.get("name", "adamw")).lower()
    if name != "adamw":
        raise ValueError(f"Only adamw is supported currently, got: {name}")

    decay_params, no_decay_params, decay_names, no_decay_names = split_decay_param_groups(model)

    rope_buffers = [
        buf_name
        for buf_name, _ in model.named_buffers()
        if ("inv_freq" in buf_name) or ("cos_cached" in buf_name) or ("sin_cached" in buf_name)
    ]

    param_groups = []
    if decay_params:
        param_groups.append({"params": decay_params, "weight_decay": weight_decay})
    if no_decay_params:
        param_groups.append({"params": no_decay_params, "weight_decay": 0.0})

    if is_main_process(rank):
        print(
            "[opt] "
            f"decay_tensors={len(decay_params)} no_decay_tensors={len(no_decay_params)} "
            f"rope_buffers={len(rope_buffers)}"
        )
        if rope_buffers:
            print("[opt] rope buffers (not optimized): " + ", ".join(rope_buffers))
        # Optional short preview
        print("[opt] decay sample: " + ", ".join(decay_names[:3]))
        print("[opt] no_decay sample: " + ", ".join(no_decay_names[:3]))

    return torch.optim.AdamW(
        param_groups,
        lr=lr,
        betas=(float(betas[0]), float(betas[1])),
    )


def build_scheduler(
    optimizer: torch.optim.Optimizer,
    train_cfg: Dict[str, Any],
    steps_per_epoch: int,
) -> Optional[torch.optim.lr_scheduler.LambdaLR]:
    sched_cfg = train_cfg.get("schedule", {})
    name = str(sched_cfg.get("name", "cosine")).lower()

    if name == "none":
        return None
    if name != "cosine":
        raise ValueError(f"Unsupported scheduler: {name}")

    max_epochs = int(train_cfg.get("max_epochs", 10))
    total_steps = max(1, max_epochs * max(1, steps_per_epoch))

    warmup_steps_cfg = int(sched_cfg.get("warmup_steps", 0))
    warmup_ratio = float(sched_cfg.get("warmup_ratio", 0.06))
    warmup_steps = warmup_steps_cfg if warmup_steps_cfg > 0 else int(total_steps * warmup_ratio)
    warmup_steps = max(0, min(warmup_steps, total_steps - 1 if total_steps > 1 else 0))

    min_lr_ratio = float(sched_cfg.get("min_lr_ratio", 0.1))

    def lr_lambda(step: int) -> float:
        if warmup_steps > 0 and step < warmup_steps:
            return float(step + 1) / float(warmup_steps)

        progress = float(step - warmup_steps) / float(max(1, total_steps - warmup_steps))
        progress = min(max(progress, 0.0), 1.0)
        cosine = 0.5 * (1.0 + math.cos(math.pi * progress))
        return min_lr_ratio + (1.0 - min_lr_ratio) * cosine

    return torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda=lr_lambda)


def maybe_compile(model: nn.Module, train_cfg: Dict[str, Any], rank: int) -> nn.Module:
    trainer_cfg = train_cfg.get("trainer", {})
    use_compile = bool(trainer_cfg.get("compile", False))
    if not use_compile:
        return model

    if not hasattr(torch, "compile"):
        if is_main_process(rank):
            print("[warn] torch.compile requested but unavailable in this torch build")
        return model

    return torch.compile(model)


def save_checkpoint(
    output_dir: Path,
    run_name: str,
    epoch: int,
    global_step: int,
    model: nn.Module,
    optimizer: torch.optim.Optimizer,
    scheduler: Optional[torch.optim.lr_scheduler.LambdaLR],
    model_cfg: Dict[str, Any],
    train_cfg: Dict[str, Any],
) -> Path:
    ckpt_dir = output_dir / "checkpoints" / run_name
    ckpt_dir.mkdir(parents=True, exist_ok=True)

    raw_model = model.module if isinstance(model, DDP) else model
    ckpt_path = ckpt_dir / f"epoch_{epoch:03d}.pt"
    payload = {
        "epoch": epoch,
        "global_step": global_step,
        "model_state": raw_model.state_dict(),
        "optimizer_state": optimizer.state_dict(),
        "scheduler_state": scheduler.state_dict() if scheduler is not None else None,
        "model_config": model_cfg,
        "train_config": train_cfg,
    }
    torch.save(payload, ckpt_path)
    return ckpt_path


def run_validation(
    model: nn.Module,
    val_loader: DataLoader,
    dist_env: DistEnv,
    loss_cfg: Dict[str, Any],
    use_amp: bool,
    amp_dtype: torch.dtype,
    online_runtime: Optional[OnlineFeatureRuntime] = None,
    require_genos_summary: bool = False,
    max_batches: int = 0,
) -> Dict[str, Any]:
    model.eval()

    totals: Dict[str, float] = {key: 0.0 for key in LOSS_METRIC_KEYS + BIAS_DIAGNOSTIC_KEYS}
    selection_totals: Dict[str, float] = {key: 0.0 for key in SELECTION_SUM_KEYS}
    n_examples = 0.0
    source_totals: Dict[str, Dict[str, float]] = {
        source: {key: 0.0 for key in LOSS_METRIC_KEYS + BIAS_DIAGNOSTIC_KEYS}
        for source in VALIDATION_REGION_KEYS
    }
    source_selection_totals: Dict[str, Dict[str, float]] = {
        source: {key: 0.0 for key in SELECTION_SUM_KEYS}
        for source in VALIDATION_REGION_KEYS
    }
    source_examples = {source: 0.0 for source in VALIDATION_REGION_KEYS}

    with torch.no_grad():
        for batch_idx, batch in enumerate(val_loader, start=1):
            if max_batches > 0 and batch_idx > max_batches:
                break
            seq = batch["seq"].to(dist_env.device, non_blocking=True)
            profile_counts = batch["profile_counts"].to(dist_env.device, non_blocking=True)
            region_sources = batch.get("region_source")

            online_feature_kwargs = extract_online_feature_kwargs(seq, online_runtime)
            foundation_feature_kwargs = extract_foundation_feature_kwargs(batch, model_cfg, dist_env.device)
            teacher_targets = extract_teacher_targets(batch, dist_env.device)

            genos_summary = None
            if "genos_global_mean" in batch:
                genos_summary = batch["genos_global_mean"].to(dist_env.device, non_blocking=True)
            elif require_genos_summary:
                raise RuntimeError(
                    "Cached Genos mode is enabled, but validation batch is missing "
                    "'genos_global_mean'. Check genos_cache_dir and genos_cache_features."
                )

            with torch.autocast(
                device_type=dist_env.device.type,
                dtype=amp_dtype,
                enabled=use_amp,
            ):
                outputs = model(
                    seq,
                    genos_summary=genos_summary,
                    **online_feature_kwargs,
                    **foundation_feature_kwargs,
                )
                _, metrics = compute_losses(outputs, profile_counts, loss_cfg, teacher_targets=teacher_targets)
                diagnostics = compute_bias_diagnostics(model, outputs)

            batch_weight = float(seq.shape[0])
            n_examples += batch_weight
            combined = {**metrics, **diagnostics}
            selection_metrics = compute_selection_metric_sums(outputs, profile_counts)
            for k in totals:
                totals[k] += to_scalar(combined[k]) * batch_weight
            for key in selection_totals:
                selection_totals[key] += to_scalar(selection_metrics[key])

            if isinstance(region_sources, Sequence) and not isinstance(region_sources, (str, bytes)):
                for source in VALIDATION_REGION_KEYS:
                    mask_values = [region_source == source for region_source in region_sources]
                    local_count = sum(mask_values)
                    if local_count <= 0:
                        continue
                    mask = torch.tensor(mask_values, device=seq.device, dtype=torch.bool)
                    source_outputs = slice_outputs(outputs, mask)
                    source_counts = profile_counts[mask]
                    source_teacher_targets = {
                        name: value[mask]
                        for name, value in teacher_targets.items()
                    }
                    _, source_loss_metrics = compute_losses(
                        source_outputs,
                        source_counts,
                        loss_cfg,
                        teacher_targets=source_teacher_targets,
                    )
                    source_diag_metrics = compute_bias_diagnostics(model, source_outputs)
                    source_selection_metrics = compute_selection_metric_sums(source_outputs, source_counts)
                    source_combined = {**source_loss_metrics, **source_diag_metrics}
                    source_examples[source] += float(local_count)
                    for key in source_totals[source]:
                        source_totals[source][key] += to_scalar(source_combined[key]) * float(local_count)
                    for key in source_selection_totals[source]:
                        source_selection_totals[source][key] += to_scalar(source_selection_metrics[key])

    overall_metrics = reduce_weighted_metrics(totals, n_examples, dist_env)
    overall_metrics.update(finalize_selection_metrics(reduce_sum_metrics(selection_totals, dist_env)))
    out: Dict[str, Any] = {
        "overall": overall_metrics,
    }
    for source in VALIDATION_REGION_KEYS:
        metrics = reduce_weighted_metrics(source_totals[source], source_examples[source], dist_env)
        metrics.update(finalize_selection_metrics(reduce_sum_metrics(source_selection_totals[source], dist_env)))
        if any(abs(value) > 0.0 for value in metrics.values()):
            out[source] = metrics
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="DDP trainer for TransChromBP")
    parser.add_argument("--model-config", type=str, default="configs/model/transchrombp_base.yaml")
    parser.add_argument("--train-config", type=str, default="configs/train/train_base.yaml")
    parser.add_argument("--data-config", type=str, default="")
    parser.add_argument("--output-dir", type=str, default="")
    parser.add_argument("--run-name", type=str, default="")
    parser.add_argument("--max-epochs", type=int, default=0, help="Override config when > 0")
    parser.add_argument("--dry-run-steps", type=int, default=0, help="Stop after N optimizer steps")
    parser.add_argument("--batch-size-per-gpu", type=int, default=0, help="Override config when > 0")
    parser.add_argument("--num-workers", type=int, default=-1, help="Override config when >= 0")
    parser.add_argument("--genos-cache-dir", type=str, default="",
                        help="Override data.genos_cache_dir in train config")
    parser.add_argument("--foundation-cache-dir", type=str, default="",
                        help="Override data.foundation_cache_dir in train config")
    parser.add_argument("--teacher-cache-dir", type=str, default="",
                        help="Override data.teacher_cache_dir in train config")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    model_cfg = load_yaml(args.model_config)
    train_cfg = load_yaml(args.train_config)
    data_config_path = resolve_data_config_path(train_cfg, args.data_config)
    data_source_cfg = load_yaml(data_config_path) if data_config_path else {}
    if data_source_cfg:
        apply_external_data_config_defaults(train_cfg, data_source_cfg)
    apply_training_mode_defaults(train_cfg)

    if args.max_epochs > 0:
        train_cfg["max_epochs"] = args.max_epochs
    if args.batch_size_per_gpu > 0:
        train_cfg.setdefault("data", {})["batch_size_per_gpu"] = args.batch_size_per_gpu
    if args.num_workers >= 0:
        train_cfg.setdefault("data", {})["num_workers"] = args.num_workers
    if args.genos_cache_dir:
        train_cfg.setdefault("data", {})["genos_cache_dir"] = args.genos_cache_dir
    if args.foundation_cache_dir:
        train_cfg.setdefault("data", {})["foundation_cache_dir"] = args.foundation_cache_dir
    if args.teacher_cache_dir:
        train_cfg.setdefault("data", {})["teacher_cache_dir"] = args.teacher_cache_dir

    trainer_cfg = train_cfg.get("trainer", {})
    backend = str(trainer_cfg.get("backend", "nccl"))

    dist_env = init_distributed(backend=backend)

    try:
        set_seed(int(train_cfg.get("seed", 1234)), dist_env.rank)

        validate_semantics_profile(train_cfg, data_source_cfg, dist_env.rank)
        validate_model_config(model_cfg)
        check_sequence_length_compatibility(model_cfg, train_cfg, dist_env.rank)
        check_profile_resolution_compatibility(model_cfg, train_cfg, dist_env.rank)

        output_dir = Path(args.output_dir or train_cfg.get("logging", {}).get("output_dir", "outputs"))
        run_name = args.run_name or train_cfg.get("logging", {}).get("run_name", "transchrombp_ddp")

        model, training_mode = build_model(model_cfg, train_cfg, dist_env.rank)
        freeze_stopgrad_bias_profile_head(model, model_cfg, dist_env.rank)
        model = model.to(dist_env.device)

        # Build frozen online foundation extractor (outside DDP / optimizer).
        genos_cached_cfg = model_cfg.get("genos_cached", {})
        uses_cached_genos = bool(genos_cached_cfg.get("enabled", False))
        cached_fusion_mode = str(genos_cached_cfg.get("fusion_mode", "")).strip().lower()
        foundation_cfg = model_cfg.get("foundation_model", {})
        uses_foundation_model = bool(foundation_cfg.get("enabled", False))
        online_runtime = build_online_feature_runtime(model_cfg, dist_env.device, dist_env.rank)

        # Prevent online fallback when cached mode is active.
        if uses_cached_genos and online_runtime is not None:
            raise RuntimeError(
                "genos_cached.enabled cannot be combined with an online foundation runtime. "
                "Disable genos_branch.enabled / caduceus_branch.enabled when using cached summaries."
            )
        required_cache_features: list[str] = []
        required_teacher_targets: list[str] = []
        if uses_cached_genos:
            data_cfg = train_cfg.setdefault("data", {})
            cache_dir = str(data_cfg.get("genos_cache_dir", "")).strip()
            configured_features = [str(x) for x in data_cfg.get("genos_cache_features", [])]
            if not cache_dir:
                raise ValueError(
                    "genos_cached.enabled=true requires data.genos_cache_dir or --genos-cache-dir"
                )
            if cached_fusion_mode in {"late_film", "count_only", "late_film_and_count"}:
                required_cache_features.append("global_mean")
            else:
                raise ValueError(
                    "Unsupported genos_cached.fusion_mode="
                    f"{cached_fusion_mode!r}; expected one of "
                    "'late_film', 'count_only', 'late_film_and_count'"
                )
            missing_features = [feat for feat in required_cache_features if feat not in configured_features]
            if missing_features:
                raise ValueError(
                    "genos_cached.enabled=true but data.genos_cache_features is missing required "
                    f"features: {missing_features!r}"
                )
        if uses_foundation_model:
            data_cfg = train_cfg.setdefault("data", {})
            validate_foundation_cache_config(model_cfg, data_cfg)
            if float(train_cfg.get("loss", {}).get("distill_profile_weight", 0.0)) > 0:
                required_teacher_targets.append("profile16")
            if float(train_cfg.get("loss", {}).get("distill_count_weight", 0.0)) > 0:
                required_teacher_targets.append("logcount")
            if float(train_cfg.get("loss", {}).get("distill_rank_weight", 0.0)) > 0:
                required_teacher_targets.append("rank")
            if required_teacher_targets:
                teacher_cache_dir = str(data_cfg.get("teacher_cache_dir", "")).strip()
                configured_targets = [str(x) for x in data_cfg.get("teacher_target_names", [])]
                if not teacher_cache_dir:
                    raise ValueError(
                        "Distillation weights are enabled but data.teacher_cache_dir / "
                        "--teacher-cache-dir is missing"
                    )
                missing_targets = [
                    name for name in required_teacher_targets if name not in configured_targets
                ]
                if missing_targets:
                    raise ValueError(
                        "Distillation weights are enabled but data.teacher_target_names is missing "
                        f"{missing_targets!r}"
                    )

        # Build optimizer before DDP wrapping so parameter grouping sees raw module types.
        optimizer = build_optimizer(model, train_cfg, dist_env.rank)

        model = maybe_compile(model, train_cfg, dist_env.rank)

        if dist_env.distributed:
            if dist_env.device.type == "cuda":
                model = DDP(
                    model,
                    device_ids=[dist_env.local_rank],
                    output_device=dist_env.local_rank,
                    find_unused_parameters=bool(trainer_cfg.get("find_unused_parameters", False)),
                )
            else:
                model = DDP(model, find_unused_parameters=bool(trainer_cfg.get("find_unused_parameters", False)))

        train_loader, val_loader, train_sampler = build_dataloaders(train_cfg, dist_env, data_source_cfg)
        scheduler = build_scheduler(optimizer, train_cfg, steps_per_epoch=len(train_loader))

        precision = str(trainer_cfg.get("precision", "bf16"))
        amp_dtype = resolve_amp_dtype(precision)
        use_amp = (dist_env.device.type == "cuda") and (amp_dtype != torch.float32)
        use_scaler = use_amp and (amp_dtype == torch.float16)
        scaler = build_grad_scaler(dist_env.device.type, enabled=use_scaler)

        grad_accum_steps = int(trainer_cfg.get("grad_accum_steps", 1))
        grad_clip = float(trainer_cfg.get("clip_grad_norm", 0.0))
        log_every = int(trainer_cfg.get("log_every_steps", 20))
        max_epochs = int(train_cfg.get("max_epochs", 10))
        validate_every = int(trainer_cfg.get("validate_every_epochs", 1))
        save_every = int(trainer_cfg.get("checkpoint_every_epochs", 1))
        best_metric_mode = str(trainer_cfg.get("best_metric_mode", "min")).lower()
        early_stop_patience = int(trainer_cfg.get("early_stop_patience", 0))
        early_stop_min_delta = float(trainer_cfg.get("early_stop_min_delta", 0.0))
        if best_metric_mode not in {"min", "max"}:
            raise ValueError(f"trainer.best_metric_mode must be 'min' or 'max', got {best_metric_mode!r}")
        if early_stop_patience < 0:
            raise ValueError(f"trainer.early_stop_patience must be >= 0, got {early_stop_patience}")
        if early_stop_min_delta < 0.0:
            raise ValueError(f"trainer.early_stop_min_delta must be >= 0, got {early_stop_min_delta}")

        if is_main_process(dist_env.rank):
            global_batch = (
                int(train_cfg.get("data", {}).get("batch_size_per_gpu", 8))
                * max(1, dist_env.world_size)
                * max(1, grad_accum_steps)
            )
            current_lr = optimizer.param_groups[0]["lr"] if optimizer.param_groups else float("nan")
            print(
                f"[setup] device={dist_env.device} world_size={dist_env.world_size} "
                f"precision={precision} use_amp={use_amp} global_batch={global_batch} lr={current_lr:.3e} "
                f"early_stop_patience={early_stop_patience} early_stop_min_delta={early_stop_min_delta:.3g}"
            )

        global_step = 0
        dry_limit = int(args.dry_run_steps)
        latest_val_metrics: Dict[str, Any] = {}
        latest_train_metrics: Dict[str, float] = {}
        epoch_metrics_path: Optional[Path] = None
        best_metric_name = ""
        best_metric_value: Optional[float] = None
        best_epoch = 0
        best_checkpoint_path = ""
        epochs_without_improvement = 0
        stopped_early = False
        stop_reason = ""
        loss_calibration = auto_configure_loss_weights(train_cfg, train_loader.dataset, dist_env)

        for epoch in range(1, max_epochs + 1):
            model.train()
            apply_foundation_freeze_policy(model, model_cfg, epoch, dist_env.rank)
            if train_sampler is not None:
                train_sampler.set_epoch(epoch)

            t0 = time.time()
            running = {key: 0.0 for key in LOSS_METRIC_KEYS + BIAS_DIAGNOSTIC_KEYS}
            running_genos_diag = {key: 0.0 for key in GENOS_FILM_KEYS}
            running_runtime = {key: 0.0 for key in TRAIN_RUNTIME_KEYS}
            running_batches = 0
            running_examples = 0.0

            optimizer.zero_grad(set_to_none=True)

            for step_idx, batch in enumerate(train_loader, start=1):
                batch_t0 = time.perf_counter()
                seq = batch["seq"].to(dist_env.device, non_blocking=True)
                profile_counts = batch["profile_counts"].to(dist_env.device, non_blocking=True)
                t_after_data = time.perf_counter()

                genos_fetch_t0 = time.perf_counter()
                online_feature_kwargs = extract_online_feature_kwargs(seq, online_runtime)
                foundation_feature_kwargs = extract_foundation_feature_kwargs(batch, model_cfg, dist_env.device)
                teacher_targets = extract_teacher_targets(batch, dist_env.device)

                genos_summary = None
                if "genos_global_mean" in batch:
                    genos_summary = batch["genos_global_mean"].to(dist_env.device, non_blocking=True)
                elif uses_cached_genos:
                    raise RuntimeError(
                        "Cached Genos mode is enabled, but training batch is missing "
                        "'genos_global_mean'. Check genos_cache_dir and genos_cache_features."
                    )
                t_after_genos = time.perf_counter()

                should_step = (step_idx % grad_accum_steps == 0)
                no_sync_ctx = model.no_sync() if (isinstance(model, DDP) and not should_step) else nullcontext()

                with no_sync_ctx:
                    with torch.autocast(
                        device_type=dist_env.device.type,
                        dtype=amp_dtype,
                        enabled=use_amp,
                    ):
                        outputs = model(
                            seq,
                            genos_summary=genos_summary,
                            **online_feature_kwargs,
                            **foundation_feature_kwargs,
                        )
                        loss, metrics = compute_losses(
                            outputs,
                            profile_counts,
                            train_cfg.get("loss", {}),
                            teacher_targets=teacher_targets,
                        )
                        diagnostics = compute_bias_diagnostics(model, outputs)
                        loss = loss / grad_accum_steps

                    if use_scaler:
                        scaler.scale(loss).backward()
                    else:
                        loss.backward()
                t_after_forward_backward = time.perf_counter()

                optimizer_time = 0.0
                if should_step:
                    opt_t0 = time.perf_counter()
                    if grad_clip > 0:
                        if use_scaler:
                            scaler.unscale_(optimizer)
                        nn.utils.clip_grad_norm_(model.parameters(), grad_clip)

                    if use_scaler:
                        scaler.step(optimizer)
                        scaler.update()
                    else:
                        optimizer.step()

                    optimizer.zero_grad(set_to_none=True)

                    if scheduler is not None:
                        scheduler.step()

                    global_step += 1
                    optimizer_time = time.perf_counter() - opt_t0

                batch_total_time = time.perf_counter() - batch_t0
                batch_runtime_metrics = {
                    "step_time_total": batch_total_time,
                    "step_time_data_to_device": t_after_data - batch_t0,
                    "step_time_genos_fetch": t_after_genos - genos_fetch_t0,
                    "step_time_forward_backward": t_after_forward_backward - t_after_genos,
                    "step_time_optimizer": optimizer_time,
                }

                batch_weight = float(seq.shape[0])
                running_batches += 1
                running_examples += batch_weight
                combined_metrics = {**metrics, **diagnostics}
                genos_diag_metrics = compute_genos_film_diagnostics(model)
                for k in running:
                    running[k] += to_scalar(combined_metrics[k]) * batch_weight
                for k in running_genos_diag:
                    running_genos_diag[k] += to_scalar(genos_diag_metrics[k]) * batch_weight
                for k in running_runtime:
                    running_runtime[k] += float(batch_runtime_metrics[k]) * batch_weight

                if (step_idx % log_every == 0) and is_main_process(dist_env.rank):
                    elapsed = time.time() - t0
                    avg_loss = running["loss_total"] / max(1.0, running_examples)
                    lr = optimizer.param_groups[0]["lr"]
                    avg_metrics = {key: running[key] / max(1.0, running_examples) for key in running}
                    avg_genos_diag = {
                        key: running_genos_diag[key] / max(1.0, running_examples)
                        for key in running_genos_diag
                    }
                    avg_runtime = {
                        key: running_runtime[key] / max(1.0, running_examples)
                        for key in running_runtime
                    }
                    diag_suffix = format_bias_diag_suffix(avg_metrics)
                    runtime_suffix = (
                        f"t_total={avg_runtime['step_time_total']:.4f}s "
                        f"t_data={avg_runtime['step_time_data_to_device']:.4f}s "
                        f"t_fm={avg_runtime['step_time_genos_fetch']:.4f}s "
                        f"t_fb={avg_runtime['step_time_forward_backward']:.4f}s"
                    )
                    genos_suffix = " ".join(
                        f"{k}={v:.4g}" for k, v in avg_genos_diag.items() if abs(v) > 0.0
                    )
                    print(
                        f"[train] epoch={epoch} step={step_idx}/{len(train_loader)} "
                        f"global_step={global_step} loss={avg_loss:.5f} lr={lr:.3e} time={elapsed:.1f}s "
                        f"{diag_suffix} {runtime_suffix} {genos_suffix}".rstrip()
                    )

                if dry_limit > 0 and global_step >= dry_limit:
                    break

            # Flush remainder grads when steps not divisible by grad_accum_steps.
            if (running_batches > 0) and (running_batches % grad_accum_steps != 0):
                if grad_clip > 0:
                    if use_scaler:
                        scaler.unscale_(optimizer)
                    nn.utils.clip_grad_norm_(model.parameters(), grad_clip)

                if use_scaler:
                    scaler.step(optimizer)
                    scaler.update()
                else:
                    optimizer.step()

                optimizer.zero_grad(set_to_none=True)
                if scheduler is not None:
                    scheduler.step()
                global_step += 1

            train_epoch_metrics = reduce_weighted_metrics(running, running_examples, dist_env)
            train_genos_diag = reduce_weighted_metrics(running_genos_diag, running_examples, dist_env)
            train_runtime_metrics = reduce_weighted_metrics(running_runtime, running_examples, dist_env)
            train_epoch_metrics = {
                **train_epoch_metrics,
                **train_runtime_metrics,
                **train_genos_diag,
            }
            latest_train_metrics = train_epoch_metrics

            if is_main_process(dist_env.rank):
                elapsed = time.time() - t0
                diag_suffix = format_bias_diag_suffix(train_epoch_metrics)
                runtime_suffix = (
                    f"t_total={train_epoch_metrics['step_time_total']:.4f}s "
                    f"t_data={train_epoch_metrics['step_time_data_to_device']:.4f}s "
                    f"t_fm={train_epoch_metrics['step_time_genos_fetch']:.4f}s "
                    f"t_fb={train_epoch_metrics['step_time_forward_backward']:.4f}s "
                    f"t_opt={train_epoch_metrics['step_time_optimizer']:.4f}s"
                )
                film_suffix = " ".join(
                    f"{k}={train_epoch_metrics[k]:.4g}" for k in GENOS_FILM_KEYS if abs(train_epoch_metrics[k]) > 0.0
                )
                print(
                    f"[epoch] {epoch} "
                    f"train_loss={train_epoch_metrics['loss_total']:.5f} "
                    f"profile={train_epoch_metrics['loss_profile']:.5f} "
                    f"count={train_epoch_metrics['loss_count']:.5f} elapsed={elapsed:.1f}s "
                    f"{diag_suffix} {runtime_suffix} {film_suffix}".rstrip()
                )

            val_metrics: Dict[str, Any] = {}
            if val_loader is not None and (epoch % validate_every == 0):
                val_metrics = run_validation(
                    model=model,
                    val_loader=val_loader,
                    dist_env=dist_env,
                    loss_cfg=train_cfg.get("loss", {}),
                    use_amp=use_amp,
                    amp_dtype=amp_dtype,
                    online_runtime=online_runtime,
                    require_genos_summary=uses_cached_genos,
                    max_batches=2 if dry_limit > 0 else 0,
                )
                latest_val_metrics = val_metrics
                if is_main_process(dist_env.rank):
                    overall_val_metrics = val_metrics.get("overall", {})
                    if overall_val_metrics:
                        print("[val] " + format_metric_summary(overall_val_metrics))
                    for source in VALIDATION_REGION_KEYS:
                        source_metrics = val_metrics.get(source, {})
                        if source_metrics:
                            print(f"[val:{source}] " + format_metric_summary(source_metrics))
            else:
                latest_val_metrics = {}

            current_ckpt: Optional[Path] = None
            candidate_metric_name = resolve_best_metric_name(trainer_cfg, val_metrics) if val_metrics else ""
            candidate_metric_value = extract_metric_value(val_metrics, candidate_metric_name) if candidate_metric_name else None
            is_best = candidate_metric_value is not None and is_better_metric(
                candidate=candidate_metric_value,
                current_best=best_metric_value,
                mode=best_metric_mode,
                min_delta=early_stop_min_delta,
            )
            if is_main_process(dist_env.rank):
                save_due = (epoch % save_every == 0)

                if save_due or is_best:
                    current_ckpt = save_checkpoint(
                        output_dir=output_dir,
                        run_name=run_name,
                        epoch=epoch,
                        global_step=global_step,
                        model=model,
                        optimizer=optimizer,
                        scheduler=scheduler,
                        model_cfg=model_cfg,
                        train_cfg=train_cfg,
                    )
                    print(f"[ckpt] saved {current_ckpt}")

                if is_best and current_ckpt is not None:
                    ckpt_dir = output_dir / "checkpoints" / run_name
                    best_alias_path = ckpt_dir / "best.pt"
                    shutil.copy2(current_ckpt, best_alias_path)
                    best_metric_name = candidate_metric_name
                    best_metric_value = candidate_metric_value
                    best_epoch = epoch
                    best_checkpoint_path = str(best_alias_path)
                    print(
                        f"[best] epoch={best_epoch} {best_metric_name}={best_metric_value:.5f} "
                        f"ckpt={best_alias_path}"
                    )

            if is_best and candidate_metric_value is not None:
                best_metric_name = candidate_metric_name
                best_metric_value = candidate_metric_value
                best_epoch = epoch
            if candidate_metric_value is not None and best_metric_value is not None:
                if is_best:
                    epochs_without_improvement = 0
                else:
                    epochs_without_improvement += 1

            if is_main_process(dist_env.rank):
                epoch_payload = {
                    "epoch": epoch,
                    "global_step": global_step,
                    "training_mode": training_mode,
                    "train": train_epoch_metrics,
                    "val": val_metrics,
                    "best_metric_name": best_metric_name,
                    "best_metric_value": best_metric_value,
                    "best_epoch": best_epoch,
                    "best_checkpoint_path": best_checkpoint_path,
                    "epochs_without_improvement": epochs_without_improvement,
                }
                epoch_metrics_path = append_epoch_metrics(output_dir / "logs" / run_name, epoch_payload)
                print(f"[metrics] wrote {epoch_metrics_path}")

            should_early_stop = (
                early_stop_patience > 0
                and candidate_metric_value is not None
                and best_metric_value is not None
                and not is_best
                and epochs_without_improvement >= early_stop_patience
            )
            if should_early_stop:
                stopped_early = True
                stop_reason = (
                    f"no improvement in {candidate_metric_name} for {epochs_without_improvement} validations "
                    f"(patience={early_stop_patience}, min_delta={early_stop_min_delta:.3g})"
                )
                if is_main_process(dist_env.rank):
                    print(f"[early-stop] epoch={epoch} {stop_reason}")
                break

            if dry_limit > 0 and global_step >= dry_limit:
                if is_main_process(dist_env.rank):
                    print(f"[done] dry run limit reached: {dry_limit} steps")
                break

        if is_main_process(dist_env.rank):
            meta_dir = output_dir / "logs" / run_name
            meta_dir.mkdir(parents=True, exist_ok=True)
            meta_path = meta_dir / "run_meta.json"
            meta = {
                "world_size": dist_env.world_size,
                "device": str(dist_env.device),
                "precision": precision,
                "training_mode": training_mode,
                "model_config": args.model_config,
                "train_config": args.train_config,
                "data_config": data_config_path,
                "global_step": global_step,
                "latest_train_metrics": latest_train_metrics,
                "latest_val_metrics": latest_val_metrics,
                "epoch_metrics_path": str(epoch_metrics_path) if epoch_metrics_path is not None else "",
                "bias_pretrained_path": str(model_cfg.get("bias_branch", {}).get("pretrained_path", "")),
                "best_metric_name": best_metric_name,
                "best_metric_value": best_metric_value,
                "best_metric_mode": best_metric_mode,
                "best_epoch": best_epoch,
                "best_checkpoint_path": best_checkpoint_path,
                "early_stop_patience": early_stop_patience,
                "early_stop_min_delta": early_stop_min_delta,
                "stopped_early": stopped_early,
                "early_stop_reason": stop_reason,
                "loss_calibration": loss_calibration,
            }
            with open(meta_path, "w", encoding="utf-8") as f:
                json.dump(meta, f, indent=2)
            print(f"[meta] wrote {meta_path}")

    finally:
        cleanup_distributed(dist_env.distributed)


if __name__ == "__main__":
    main()
