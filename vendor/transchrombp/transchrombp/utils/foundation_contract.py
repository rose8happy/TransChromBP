"""Shared contract helpers for cached foundation-model features."""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Any, Dict, Optional, Sequence, Tuple


FEATURE_NAME_RE = re.compile(r"^layer_(\d+)__(.+)$")
SUPPORTED_FOUNDATION_MODES = {"distill_only", "cross_attention", "residual_head"}


@dataclass(frozen=True)
class FoundationContract:
    mode: str
    feature_name: str
    feature_layout: str
    feature_tokens: int
    hidden_size: int


def foundation_batch_key(feature_name: str) -> str:
    return f"foundation_{feature_name}"


def parse_foundation_feature_name(feature_name: str) -> Tuple[int, str]:
    normalized = str(feature_name).strip()
    match = FEATURE_NAME_RE.fullmatch(normalized)
    if match is None:
        raise ValueError(
            "foundation_model.feature_name must look like 'layer_07__bins4_mean', "
            f"got {feature_name!r}"
        )
    return int(match.group(1)), match.group(2)


def resolve_foundation_contract(model_cfg: Dict[str, Any]) -> Optional[FoundationContract]:
    foundation_cfg = model_cfg.get("foundation_model", {})
    if not bool(foundation_cfg.get("enabled", False)):
        return None

    mode = str(foundation_cfg.get("mode", "distill_only")).strip().lower()
    if mode not in SUPPORTED_FOUNDATION_MODES:
        raise ValueError(
            "Unsupported foundation_model.mode="
            f"{mode!r}; expected one of {sorted(SUPPORTED_FOUNDATION_MODES)!r}"
        )
    if mode == "distill_only":
        return None

    feature_name = str(foundation_cfg.get("feature_name", "")).strip()
    if not feature_name:
        raise ValueError(
            "foundation_model.enabled=true requires foundation_model.feature_name for "
            f"mode={mode!r}"
        )
    parse_foundation_feature_name(feature_name)

    feature_layout = str(foundation_cfg.get("feature_layout", "summary")).strip().lower()
    if feature_layout not in {"summary", "token"}:
        raise ValueError(
            "Unsupported foundation_model.feature_layout="
            f"{feature_layout!r}; expected 'summary' or 'token'"
        )

    hidden_size = int(foundation_cfg.get("hidden_size", foundation_cfg.get("feature_hidden_size", 0)))
    feature_tokens = int(foundation_cfg.get("feature_tokens", 0))
    if feature_layout == "token":
        if hidden_size <= 0:
            raise ValueError(
                "Token foundation features require foundation_model.hidden_size > 0 "
                "(or feature_hidden_size > 0)"
            )
        if feature_tokens <= 0:
            raise ValueError(
                "Token foundation features require foundation_model.feature_tokens > 0"
            )

    return FoundationContract(
        mode=mode,
        feature_name=feature_name,
        feature_layout=feature_layout,
        feature_tokens=feature_tokens,
        hidden_size=hidden_size,
    )


def infer_required_foundation_cache_features(model_cfg: Dict[str, Any]) -> Tuple[str, ...]:
    contract = resolve_foundation_contract(model_cfg)
    if contract is None:
        return ()
    return (contract.feature_name,)


def infer_foundation_cache_build_request(
    model_cfg: Dict[str, Any],
) -> Tuple[Tuple[int, ...], Tuple[str, ...]]:
    required_features = infer_required_foundation_cache_features(model_cfg)
    ordered_layers: list[int] = []
    ordered_feature_types: list[str] = []
    for feature_name in required_features:
        layer, feature_type = parse_foundation_feature_name(feature_name)
        if layer not in ordered_layers:
            ordered_layers.append(layer)
        if feature_type not in ordered_feature_types:
            ordered_feature_types.append(feature_type)
    return tuple(ordered_layers), tuple(ordered_feature_types)


def validate_foundation_cache_build_request(
    model_cfg: Dict[str, Any],
    requested_layers: Sequence[int],
    requested_feature_types: Sequence[str],
) -> Optional[FoundationContract]:
    contract = resolve_foundation_contract(model_cfg)
    if contract is None:
        return None

    required_layer, required_feature_type = parse_foundation_feature_name(contract.feature_name)
    if int(required_layer) not in {int(x) for x in requested_layers}:
        raise ValueError(
            "foundation cache build request must include required layer "
            f"{required_layer} for feature {contract.feature_name!r}"
        )
    if str(required_feature_type) not in {str(x) for x in requested_feature_types}:
        raise ValueError(
            "foundation cache build request must include required feature type "
            f"{required_feature_type!r} for feature {contract.feature_name!r}"
        )
    return contract


def normalize_cache_split_name(split: str) -> str:
    normalized = str(split).strip().lower()
    return "valid" if normalized == "val" else normalized


def resolve_cache_split_region_source(data_cfg: Dict[str, Any], split: str) -> str:
    default_region_source = str(data_cfg.get("region_source", "both"))
    normalized = normalize_cache_split_name(split)
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


def resolve_cache_split_max_records(data_cfg: Dict[str, Any], split: str) -> int:
    normalized = normalize_cache_split_name(split)
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


def resolve_cache_split_seed(base_seed: int, split: str) -> int:
    normalized = normalize_cache_split_name(split)
    if normalized == "train":
        return int(base_seed)
    return int(base_seed) + 10_000


def validate_foundation_cache_manifest(
    manifest: Dict[str, Any],
    *,
    model_cfg: Dict[str, Any],
    data_cfg: Dict[str, Any],
    split: str,
    base_seed: int,
) -> Optional[FoundationContract]:
    normalized_split = normalize_cache_split_name(split)
    manifest_split = normalize_cache_split_name(str(manifest.get("split", "")))
    if manifest_split != normalized_split:
        raise ValueError(
            f"manifest split mismatch: expected {normalized_split!r}, got {manifest_split!r}"
        )

    contract = validate_foundation_cache_build_request(
        model_cfg,
        requested_layers=tuple(int(x) for x in manifest.get("layers", ())),
        requested_feature_types=tuple(str(x) for x in manifest.get("feature_types", ())),
    )
    if contract is None:
        return None

    if "region_source" not in manifest:
        raise ValueError("manifest is missing required field 'region_source'")
    expected_region_source = resolve_cache_split_region_source(data_cfg, normalized_split)
    if str(manifest.get("region_source")) != expected_region_source:
        raise ValueError(
            f"manifest region_source mismatch: expected {expected_region_source!r}, "
            f"got {manifest.get('region_source')!r}"
        )

    if "max_records" not in manifest:
        raise ValueError("manifest is missing required field 'max_records'")
    expected_max_records = resolve_cache_split_max_records(data_cfg, normalized_split)
    if int(manifest.get("max_records")) != expected_max_records:
        raise ValueError(
            f"manifest max_records mismatch: expected {expected_max_records}, "
            f"got {manifest.get('max_records')!r}"
        )

    if "split_seed" not in manifest:
        raise ValueError("manifest is missing required field 'split_seed'")
    expected_seed = resolve_cache_split_seed(base_seed, normalized_split)
    if int(manifest.get("split_seed")) != expected_seed:
        raise ValueError(
            f"manifest split_seed mismatch: expected {expected_seed}, got {manifest.get('split_seed')!r}"
        )
    return contract


def validate_foundation_cache_config(
    model_cfg: Dict[str, Any],
    data_cfg: Dict[str, Any],
) -> Optional[FoundationContract]:
    contract = resolve_foundation_contract(model_cfg)
    if contract is None:
        return None

    cache_dir = str(data_cfg.get("foundation_cache_dir", "")).strip()
    if not cache_dir:
        raise ValueError(
            "foundation_model.enabled=true requires data.foundation_cache_dir or "
            "--foundation-cache-dir"
        )

    configured_features = tuple(str(x) for x in data_cfg.get("foundation_cache_features", ()))
    required_features = infer_required_foundation_cache_features(model_cfg)
    missing_features = [feat for feat in required_features if feat not in configured_features]
    if missing_features:
        raise ValueError(
            "foundation_model requires data.foundation_cache_features to include "
            f"{missing_features!r}"
        )
    return contract
