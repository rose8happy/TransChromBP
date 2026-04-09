"""Unit tests for the center-aligned NT v2 residual head."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
import torch


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

from transchrombp.models.foundation_adapter import FoundationResidualHead


def test_center_token_bins_residual_head_emits_profile_and_count_shapes() -> None:
    head = FoundationResidualHead(
        d_model=8,
        foundation_hidden_size=12,
        hidden_size=16,
        profile_bin_count=128,
        output_len=1000,
        alignment_mode="center_token_bins",
        aligned_token_count=16,
    )

    encoded = torch.randn(2, 2114, 8)
    foundation_tokens = torch.randn(2, 16, 12)

    profile_delta, count_delta = head(encoded, foundation_tokens=foundation_tokens)

    assert profile_delta.shape == (2, 1000)
    assert count_delta.shape == (2, 1)
    assert torch.allclose(profile_delta, torch.zeros_like(profile_delta))
    assert torch.allclose(count_delta, torch.zeros_like(count_delta))


def test_center_token_bins_residual_head_requires_token_inputs() -> None:
    head = FoundationResidualHead(
        d_model=8,
        foundation_hidden_size=12,
        hidden_size=16,
        profile_bin_count=128,
        output_len=1000,
        alignment_mode="center_token_bins",
        aligned_token_count=16,
    )

    encoded = torch.randn(1, 2114, 8)
    foundation_summary = torch.randn(1, 12)

    with pytest.raises(ValueError, match="foundation_tokens"):
        head(encoded, foundation_summary=foundation_summary)
