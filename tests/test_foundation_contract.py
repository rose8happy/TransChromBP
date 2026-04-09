from pathlib import Path
import unittest

import yaml

from transchrombp.utils.foundation_contract import (
    foundation_batch_key,
    infer_foundation_cache_build_request,
    infer_required_foundation_cache_features,
    parse_foundation_feature_name,
    normalize_cache_split_name,
    resolve_cache_split_max_records,
    resolve_cache_split_region_source,
    resolve_cache_split_seed,
    resolve_foundation_contract,
    validate_foundation_cache_manifest,
    validate_foundation_cache_build_request,
    validate_foundation_cache_config,
)


CONFIG_ROOT = Path(__file__).resolve().parents[1] / "vendor" / "transchrombp" / "transchrombp" / "configs"


def load_yaml_config(relative_path: str) -> dict:
    path = CONFIG_ROOT / relative_path
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def make_bins16_residual_cfg() -> dict:
    model_cfg = load_yaml_config("model/transchrombp_teacher_v2_center_pool_ntv2_residual.yaml")
    foundation_cfg = model_cfg["foundation_model"]
    foundation_cfg["feature_name"] = "layer_07__bins16_mean"
    foundation_cfg["feature_tokens"] = 16
    foundation_cfg["profile_bin_count"] = 128
    foundation_cfg["aligned_token_count"] = 16
    foundation_cfg["alignment_mode"] = "center_token_bins"
    return model_cfg


class FoundationContractTests(unittest.TestCase):
    def test_resolve_foundation_contract_from_model_configs(self) -> None:
        residual_cfg = load_yaml_config("model/transchrombp_teacher_v2_center_pool_ntv2_residual.yaml")
        residual_contract = resolve_foundation_contract(residual_cfg)
        self.assertIsNotNone(residual_contract)
        self.assertEqual(residual_contract.mode, "residual_head")
        self.assertEqual(residual_contract.feature_name, "layer_07__bins4_mean")
        self.assertEqual(residual_contract.feature_layout, "token")
        self.assertEqual(residual_contract.feature_tokens, 4)
        self.assertEqual(residual_contract.hidden_size, 1024)

        cross_attn_cfg = load_yaml_config("model/transchrombp_teacher_v2_center_pool_ntv2_cross_attention.yaml")
        cross_attn_contract = resolve_foundation_contract(cross_attn_cfg)
        self.assertIsNotNone(cross_attn_contract)
        self.assertEqual(cross_attn_contract.mode, "cross_attention")
        self.assertEqual(cross_attn_contract.feature_name, "layer_07__bins4_mean")
        self.assertEqual(cross_attn_contract.feature_tokens, 4)

        bins16_contract = resolve_foundation_contract(make_bins16_residual_cfg())
        self.assertIsNotNone(bins16_contract)
        self.assertEqual(bins16_contract.feature_name, "layer_07__bins16_mean")
        self.assertEqual(bins16_contract.feature_tokens, 16)

    def test_parse_feature_name_and_batch_key(self) -> None:
        self.assertEqual(parse_foundation_feature_name("layer_07__bins4_mean"), (7, "bins4_mean"))
        self.assertEqual(parse_foundation_feature_name("layer_14__global_mean"), (14, "global_mean"))
        self.assertEqual(
            foundation_batch_key("layer_07__bins16_mean"),
            "foundation_layer_07__bins16_mean",
        )

    def test_infer_required_cache_features(self) -> None:
        residual_cfg = load_yaml_config("model/transchrombp_teacher_v2_center_pool_ntv2_residual.yaml")
        self.assertEqual(
            infer_required_foundation_cache_features(residual_cfg),
            ("layer_07__bins4_mean",),
        )

        bins16_cfg = make_bins16_residual_cfg()
        self.assertEqual(
            infer_required_foundation_cache_features(bins16_cfg),
            ("layer_07__bins16_mean",),
        )

    def test_infer_foundation_cache_build_request(self) -> None:
        residual_cfg = load_yaml_config("model/transchrombp_teacher_v2_center_pool_ntv2_residual.yaml")
        self.assertEqual(
            infer_foundation_cache_build_request(residual_cfg),
            ((7,), ("bins4_mean",)),
        )

        bins16_cfg = make_bins16_residual_cfg()
        self.assertEqual(
            infer_foundation_cache_build_request(bins16_cfg),
            ((7,), ("bins16_mean",)),
        )

    def test_validate_foundation_cache_build_request(self) -> None:
        model_cfg = make_bins16_residual_cfg()
        validate_foundation_cache_build_request(
            model_cfg,
            requested_layers=(7, 14),
            requested_feature_types=("global_mean", "bins16_mean"),
        )
        with self.assertRaisesRegex(ValueError, "must include"):
            validate_foundation_cache_build_request(
                model_cfg,
                requested_layers=(7, 14),
                requested_feature_types=("global_mean", "bins4_mean"),
            )

    def test_cache_split_helpers(self) -> None:
        data_cfg = {
            "region_source": "both",
            "train_region_source": "peaks_only",
            "val_region_source": "nonpeaks_only",
            "max_records": 17,
            "max_train_regions": 101,
            "max_valid_regions": 23,
            "max_test_regions": 31,
        }
        self.assertEqual(normalize_cache_split_name("val"), "valid")
        self.assertEqual(normalize_cache_split_name("test"), "test")
        self.assertEqual(resolve_cache_split_region_source(data_cfg, "train"), "peaks_only")
        self.assertEqual(resolve_cache_split_region_source(data_cfg, "val"), "nonpeaks_only")
        self.assertEqual(resolve_cache_split_region_source(data_cfg, "test"), "nonpeaks_only")
        self.assertEqual(resolve_cache_split_max_records(data_cfg, "train"), 101)
        self.assertEqual(resolve_cache_split_max_records(data_cfg, "valid"), 23)
        self.assertEqual(resolve_cache_split_max_records(data_cfg, "test"), 31)
        self.assertEqual(resolve_cache_split_seed(42, "train"), 42)
        self.assertEqual(resolve_cache_split_seed(42, "valid"), 10042)
        self.assertEqual(resolve_cache_split_seed(42, "test"), 10042)

    def test_validate_foundation_cache_manifest(self) -> None:
        model_cfg = make_bins16_residual_cfg()
        data_cfg = {
            "foundation_cache_dir": "/tmp/foundation-cache",
            "foundation_cache_features": ["layer_07__bins16_mean"],
            "region_source": "both",
            "train_region_source": "both",
            "val_region_source": "both",
            "max_valid_regions": 64,
        }
        manifest = {
            "split": "valid",
            "layers": [7],
            "feature_types": ["bins16_mean"],
            "region_source": "both",
            "max_records": 64,
            "split_seed": 10042,
        }
        validate_foundation_cache_manifest(
            manifest,
            model_cfg=model_cfg,
            data_cfg=data_cfg,
            split="valid",
            base_seed=42,
        )
        bad_manifest = dict(manifest)
        bad_manifest["feature_types"] = ["bins4_mean"]
        with self.assertRaisesRegex(ValueError, "required feature type"):
            validate_foundation_cache_manifest(
                bad_manifest,
                model_cfg=model_cfg,
                data_cfg=data_cfg,
                split="valid",
                base_seed=42,
            )

    def test_validate_cache_config_rejects_missing_feature(self) -> None:
        model_cfg = make_bins16_residual_cfg()
        train_cfg = load_yaml_config("train/train_tutorial_foundation_short10.yaml")
        data_cfg = train_cfg["data"]
        data_cfg["foundation_cache_dir"] = "/tmp/foundation-cache"
        data_cfg["foundation_cache_features"] = ["layer_07__bins4_mean"]

        with self.assertRaisesRegex(ValueError, "foundation_cache_features"):
            validate_foundation_cache_config(model_cfg, data_cfg)

    def test_validate_cache_config_accepts_matching_feature(self) -> None:
        model_cfg = make_bins16_residual_cfg()
        train_cfg = load_yaml_config("train/train_tutorial_foundation_short10.yaml")
        data_cfg = train_cfg["data"]
        data_cfg["foundation_cache_dir"] = "/tmp/foundation-cache"
        data_cfg["foundation_cache_features"] = ["layer_07__bins16_mean"]

        contract = validate_foundation_cache_config(model_cfg, data_cfg)
        self.assertEqual(contract.feature_name, "layer_07__bins16_mean")


if __name__ == "__main__":
    unittest.main()
