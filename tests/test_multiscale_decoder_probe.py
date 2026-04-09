"""Tests for the multiscale profile decoder probe."""

from __future__ import annotations

import sys
from pathlib import Path

import torch
import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]
IMPORT_ROOTS = [
    REPO_ROOT / "vendor" / "transchrombp",
    REPO_ROOT / "src",
    REPO_ROOT,
]
for import_root in IMPORT_ROOTS:
    if (import_root / "transchrombp").exists() and str(import_root) not in sys.path:
        sys.path.insert(0, str(import_root))
        break

from transchrombp.models.transchrombp import TransChromBP, build_transchrombp_from_config


def test_multiscale_profile_decoder_emits_profile_and_count_shapes() -> None:
    model = TransChromBP(
        d_model=16,
        n_heads=4,
        n_layers=2,
        ff_mult=2,
        output_len=128,
        signal_n_dil_layers=2,
        use_bias_decomposition=False,
        count_pool_mode="center",
        profile_readout_mode="multiscale_decoder_v1",
        profile_decoder_hidden_channels=24,
        profile_decoder_num_scales=3,
    )

    seq = torch.randn(2, 256, 4)
    outputs = model(seq)

    assert outputs.profile_logits_full.shape == (2, 128)
    assert outputs.profile_logits_debiased.shape == (2, 128)
    assert outputs.logcount_full.shape == (2, 1)
    assert outputs.logcount_debiased.shape == (2, 1)


def test_multiscale_profile_decoder_preserves_full_debiased_additive_contract() -> None:
    model = TransChromBP(
        d_model=16,
        n_heads=4,
        n_layers=2,
        ff_mult=2,
        output_len=128,
        signal_n_dil_layers=2,
        use_bias_decomposition=False,
        count_pool_mode="center",
        count_fusion="add",
        learnable_scales=False,
        positive_scales=False,
        profile_scale_init=1.0,
        count_scale_init=1.0,
        profile_readout_mode="multiscale_decoder_v1",
        profile_decoder_hidden_channels=24,
        profile_decoder_num_scales=3,
    )

    seq = torch.randn(2, 256, 4)
    bias_profile = torch.randn(2, 128)
    bias_count = torch.randn(2, 1)
    outputs = model(seq, bias_profile=bias_profile, bias_count=bias_count)

    torch.testing.assert_close(
        outputs.profile_logits_full,
        outputs.profile_logits_debiased + bias_profile,
    )
    torch.testing.assert_close(
        outputs.logcount_full,
        outputs.logcount_debiased + bias_count,
    )


def test_skip_probe_profile_readout_emits_profile_and_count_shapes() -> None:
    model = TransChromBP(
        d_model=16,
        n_heads=4,
        n_layers=2,
        ff_mult=2,
        output_len=128,
        signal_n_dil_layers=2,
        use_bias_decomposition=False,
        count_pool_mode="center",
        profile_readout_mode="skip_probe_v1",
        profile_decoder_hidden_channels=24,
    )

    seq = torch.randn(2, 256, 4)
    outputs = model(seq)

    assert outputs.profile_logits_full.shape == (2, 128)
    assert outputs.profile_logits_debiased.shape == (2, 128)
    assert outputs.logcount_full.shape == (2, 1)
    assert outputs.logcount_debiased.shape == (2, 1)


def test_build_from_config_attaches_multiscale_profile_decoder() -> None:
    config = {
        "sequence_encoder": {
            "enabled": True,
            "d_model": 16,
            "n_heads": 4,
            "n_layers": 2,
            "ff_mult": 2,
            "dropout": 0.0,
            "use_rope": True,
            "rope_theta": 10000.0,
            "max_len": 256,
            "use_sdpa": True,
        },
        "conv_stem": {
            "stem_kernel_size": 15,
            "conv_kernel_size": 5,
            "n_conv_layers": 2,
        },
        "local_tower": {
            "kernel_size": 3,
            "n_dil_layers": 2,
            "dilation_cycle_length": 2,
        },
        "bias_branch": {
            "enabled": False,
            "profile_pool_factor": 32,
        },
        "fusion": {
            "profile_fusion": "add",
            "count_fusion": "logsumexp",
            "profile_scale_init": 1.0,
            "count_scale_init": 1.0,
            "learnable_scales": True,
            "positive_scales": True,
            "profile_bias_stop_gradient": True,
        },
        "heads": {
            "profile_output_len": 128,
            "count_head": "linear",
            "count_pool_mode": "center",
        },
        "profile_decoder": {
            "enabled": True,
            "mode": "multiscale_decoder_v1",
            "hidden_channels": 24,
            "num_scales": 3,
        },
    }

    model = build_transchrombp_from_config(config)

    assert model.profile_readout_mode == "multiscale_decoder_v1"
    assert model.profile_decoder is not None


def test_model_configs_build_narrow_and_skip_probe_variants() -> None:
    config_dir = REPO_ROOT / "vendor" / "transchrombp" / "transchrombp" / "configs" / "model"

    narrow_cfg = yaml.safe_load(
        (config_dir / "transchrombp_teacher_v2_center_pool_msdec_v1_narrow.yaml").read_text(
            encoding="utf-8"
        )
    )
    narrow_model = build_transchrombp_from_config(narrow_cfg)
    assert narrow_model.profile_readout_mode == "multiscale_decoder_v1"
    assert narrow_model.profile_decoder is not None
    assert narrow_model.profile_decoder.num_scales == 3
    assert narrow_model.profile_decoder.input_proj.out_channels == 128

    skip_cfg = yaml.safe_load(
        (config_dir / "transchrombp_teacher_v2_center_pool_skipprobe_v1.yaml").read_text(
            encoding="utf-8"
        )
    )
    skip_model = build_transchrombp_from_config(skip_cfg)
    assert skip_model.profile_readout_mode == "skip_probe_v1"
    assert skip_model.profile_decoder is not None

    skip_wide_cfg = yaml.safe_load(
        (config_dir / "transchrombp_teacher_v2_center_pool_skipprobe_v1_wide.yaml").read_text(
            encoding="utf-8"
        )
    )
    skip_wide_model = build_transchrombp_from_config(skip_wide_cfg)
    assert skip_wide_model.profile_readout_mode == "skip_probe_v1"
    assert skip_wide_model.profile_decoder is not None
    assert skip_wide_model.profile_decoder.readout[0].out_channels == 256
