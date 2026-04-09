from __future__ import annotations

from copy import deepcopy

import pytest
import torch
from torch import nn

from transchrombp.models.bias_branch import center_crop_1d
from transchrombp.models.profile_decoder import MultiScaleLocalSkipDecoderV2
from transchrombp.models.transchrombp import build_transchrombp_from_config


def _base_config() -> dict:
    return {
        "sequence_encoder": {
            "enabled": True,
            "d_model": 64,
            "n_heads": 4,
            "n_layers": 2,
            "ff_mult": 2,
            "dropout": 0.0,
            "use_rope": True,
            "rope_theta": 10000.0,
            "max_len": 2114,
            "use_sdpa": True,
        },
        "conv_stem": {
            "stem_kernel_size": 21,
            "conv_kernel_size": 7,
            "n_conv_layers": 2,
        },
        "local_tower": {
            "kernel_size": 3,
            "n_dil_layers": 4,
            "dilation_cycle_length": 4,
        },
        "bias_branch": {
            "enabled": True,
            "hidden_channels": 32,
            "kernel_size": 21,
            "n_dil_layers": 2,
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
            "profile_output_len": 1000,
            "count_head": "linear",
            "count_pool_mode": "center",
        },
        "profile_decoder": {
            "mode": "multiscale_local_skip_v2",
            "hidden_channels": 64,
            "decoder_channels": [64, 48, 32],
            "dropout": 0.0,
            "upsample_mode": "linear",
        },
    }


def test_multiscale_local_skip_v2_builds_and_preserves_output_shapes():
    model = build_transchrombp_from_config(_base_config())
    seq = torch.randn(2, 2114, 4)
    out = model(seq)
    assert out.profile_logits_full.shape == (2, 1000)
    assert out.profile_logits_debiased.shape == (2, 1000)
    assert out.logcount_full.shape == (2, 1)
    assert out.logcount_debiased.shape == (2, 1)


def test_multiscale_local_skip_v2_keeps_count_pool_mode_center():
    model = build_transchrombp_from_config(_base_config())
    assert model.count_pool_mode == "center"
    assert model.profile_bias_stop_gradient is True


def test_unknown_profile_decoder_mode_raises():
    cfg = deepcopy(_base_config())
    cfg["profile_decoder"]["mode"] = "does_not_exist"
    with pytest.raises(ValueError, match="profile_decoder.mode"):
        build_transchrombp_from_config(cfg)


class _IdentityRefine(nn.Module):
    def forward(self, x: torch.Tensor, skip: torch.Tensor) -> torch.Tensor:
        return x


class _FirstChannelHead(nn.Module):
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x[:, :1, :]


def test_multiscale_local_skip_v2_decodes_on_full_sequence_grid_before_crop():
    seq_len = 2114
    output_len = 1000
    decoder = MultiScaleLocalSkipDecoderV2(
        encoded_channels=4,
        local_channels=4,
        output_len=output_len,
        hidden_channels=4,
        decoder_channels=[4],
        dropout=0.0,
        upsample_mode="linear",
    )
    decoder.input_proj = nn.Identity()
    decoder.refine_blocks = nn.ModuleList([_IdentityRefine()])
    decoder.profile_head = _FirstChannelHead()

    encoded = torch.zeros(1, seq_len, 4)
    encoded[0, :, 0] = torch.arange(seq_len, dtype=torch.float32)
    local_feat = torch.zeros(1, 4, seq_len)

    profile_logits = decoder(encoded, local_feat)

    expected = center_crop_1d(torch.arange(seq_len, dtype=torch.float32).unsqueeze(0), output_len)
    assert torch.equal(profile_logits, expected)
