from pathlib import Path
import sys
import types

import pytest
import torch
import yaml

if "pyBigWig" not in sys.modules:
    sys.modules["pyBigWig"] = types.SimpleNamespace(open=lambda *args, **kwargs: None)

import transchrombp.models as model_lib
from transchrombp.models import TransChromBPOutput, build_transchrombp_from_config
from transchrombp.training.train_ddp import extract_genos_summary_kwargs


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ARCHIVE_ROOT = PROJECT_ROOT / "vendor/transchrombp"
MODEL_CONFIG = ARCHIVE_ROOT / "configs/model/transchrombp_teacher_v2_hierdec4096.yaml"


class CustomCheckpointMetadata:
    def __init__(self, note: str) -> None:
        self.note = note


def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def test_hierarchical_config_dispatch_builds_model() -> None:
    cfg = load_yaml(MODEL_CONFIG)

    assert cfg["model_name"] == "transchrombp_teacher_v2_hierdec4096"
    assert cfg["architecture"]["variant"] == "hierarchical_encoder_decoder_v1"
    assert cfg["sequence_encoder"]["max_len"] == 4096
    assert hasattr(model_lib, "HierarchicalTransChromBP")
    assert hasattr(model_lib, "build_hierarchical_transchrombp_from_config")

    model = build_transchrombp_from_config(cfg)

    assert isinstance(model, model_lib.HierarchicalTransChromBP)
    assert hasattr(model, "bias_branch")
    assert model.bias_branch is not None


def test_hierarchical_config_avoids_unused_baseline_only_fields() -> None:
    cfg = load_yaml(MODEL_CONFIG)

    assert "local_tower" not in cfg
    assert "n_conv_layers" not in cfg["conv_stem"]
    assert "type" not in cfg["sequence_encoder"]
    assert "type" not in cfg["bias_branch"]


def test_hierarchical_forward_matches_transchrombp_output_contract() -> None:
    cfg = load_yaml(MODEL_CONFIG)
    model = build_transchrombp_from_config(cfg)

    outputs = model(torch.randn(2, 4096, 4))

    assert isinstance(outputs, TransChromBPOutput)
    assert outputs.profile_logits_full.shape == (2, 1000)
    assert outputs.profile_logits_debiased.shape == (2, 1000)
    assert outputs.logcount_full.shape == (2, 1)
    assert outputs.logcount_debiased.shape == (2, 1)
    assert outputs.profile_bias.shape == (2, 1000)
    assert outputs.count_bias.shape == (2, 1)


def test_hierarchical_build_rejects_unsupported_foundation_config() -> None:
    cfg = load_yaml(MODEL_CONFIG)
    cfg["foundation_model"] = {"enabled": True, "mode": "distill_only"}

    with pytest.raises(ValueError, match="foundation_model"):
        build_transchrombp_from_config(cfg)


def test_hierarchical_forward_rejects_unexpected_adapter_kwargs() -> None:
    cfg = load_yaml(MODEL_CONFIG)
    model = build_transchrombp_from_config(cfg)
    model.eval()

    with pytest.raises(ValueError, match="foundation_summary|foundation_tokens"):
        model(
            torch.randn(2, 4096, 4),
            foundation_tokens=torch.randn(2, 16, 32),
            foundation_summary=torch.randn(2, 32),
        )


def test_hierarchical_runtime_skips_genos_summary_when_cached_genos_disabled() -> None:
    batch = {"genos_global_mean": torch.randn(2, 4)}

    kwargs = extract_genos_summary_kwargs(
        batch=batch,
        model_cfg={"genos_cached": {"enabled": False}},
        device=torch.device("cpu"),
        require_genos_summary=False,
    )

    assert kwargs == {}


def test_hierarchical_bias_checkpoint_loading_matches_baseline_contract(tmp_path: Path) -> None:
    cfg = load_yaml(MODEL_CONFIG)
    model = build_transchrombp_from_config(cfg)
    state_dict = model.bias_branch.state_dict()
    first_key = next(iter(state_dict))
    patched_state = {name: tensor.clone() for name, tensor in state_dict.items()}
    patched_state[first_key] = torch.full_like(patched_state[first_key], 0.125)

    checkpoint_path = tmp_path / "hier_bias_with_metadata.pt"
    torch.save(
        {
            "model_state": {f"bias_branch.{name}": tensor for name, tensor in patched_state.items()},
            "metadata": CustomCheckpointMetadata("baseline-compatible"),
        },
        checkpoint_path,
    )

    reloaded_model = build_transchrombp_from_config(cfg)
    reloaded_model.load_bias_branch_weights(str(checkpoint_path))

    loaded_state = reloaded_model.bias_branch.state_dict()
    assert torch.allclose(loaded_state[first_key], patched_state[first_key])
