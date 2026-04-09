"""Profile readout modules for TransChromBP."""

from __future__ import annotations

from typing import List

import torch
import torch.nn.functional as F
from torch import Tensor, nn

from .bias_branch import center_crop_1d


class _FusionBlock(nn.Module):
    """Lightweight conv block used inside the multiscale decoder."""

    def __init__(self, in_channels: int, hidden_channels: int, dropout: float) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv1d(in_channels, hidden_channels, kernel_size=3, padding=1),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Conv1d(hidden_channels, hidden_channels, kernel_size=3, padding=1),
            nn.GELU(),
        )

    def forward(self, x: Tensor) -> Tensor:
        return self.net(x)


class MultiScaleProfileDecoder(nn.Module):
    """Coarse-to-fine profile decoder with local skip features."""

    def __init__(
        self,
        d_model: int,
        output_len: int,
        hidden_channels: int,
        num_scales: int = 3,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        if num_scales < 2:
            raise ValueError(f"num_scales must be >= 2, got {num_scales}")
        if output_len <= 0:
            raise ValueError(f"output_len must be positive, got {output_len}")

        self.output_len = output_len
        self.num_scales = num_scales
        self.input_proj = nn.Conv1d(d_model, hidden_channels, kernel_size=1)
        self.stage_blocks = nn.ModuleList(
            [
                _FusionBlock(
                    in_channels=hidden_channels + (2 * d_model),
                    hidden_channels=hidden_channels,
                    dropout=dropout,
                )
                for _ in range(num_scales)
            ]
        )
        self.output_proj = nn.Conv1d(hidden_channels, 1, kernel_size=1)

    def _scale_lengths(self) -> List[int]:
        lengths: List[int] = []
        max_divisor = 2 ** (self.num_scales - 1)
        for scale_idx in range(self.num_scales):
            divisor = max_divisor // (2**scale_idx)
            length = max(1, self.output_len // max(1, divisor))
            lengths.append(length)
        lengths[-1] = self.output_len
        return lengths

    @staticmethod
    def _resize(x: Tensor, target_len: int) -> Tensor:
        if x.size(-1) == target_len:
            return x
        if x.size(-1) < target_len:
            return F.interpolate(x, size=target_len, mode="linear", align_corners=False)
        return F.adaptive_avg_pool1d(x, target_len)

    def forward(self, encoded: Tensor, local_feat: Tensor) -> Tensor:
        if encoded.dim() != 3:
            raise ValueError(f"Expected encoded shape [B, L, D], got {tuple(encoded.shape)}")
        if local_feat.dim() != 3:
            raise ValueError(f"Expected local_feat shape [B, D, L], got {tuple(local_feat.shape)}")

        encoded_cf = center_crop_1d(encoded.transpose(1, 2), self.output_len)
        local_cf = center_crop_1d(local_feat, self.output_len)

        scale_lengths = self._scale_lengths()
        x = self.input_proj(self._resize(encoded_cf, scale_lengths[0]))

        for stage_idx, target_len in enumerate(scale_lengths):
            if stage_idx > 0:
                x = F.interpolate(x, size=target_len, mode="linear", align_corners=False)
            encoded_scaled = self._resize(encoded_cf, target_len)
            local_scaled = self._resize(local_cf, target_len)
            x = self.stage_blocks[stage_idx](torch.cat([x, encoded_scaled, local_scaled], dim=1))

        return self.output_proj(x).squeeze(1)


class LocalSkipProfileProbe(nn.Module):
    """Minimal skip probe that fuses local high-frequency features at readout."""

    def __init__(
        self,
        d_model: int,
        output_len: int,
        hidden_channels: int,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        if output_len <= 0:
            raise ValueError(f"output_len must be positive, got {output_len}")
        if hidden_channels <= 0:
            raise ValueError(f"hidden_channels must be positive, got {hidden_channels}")

        self.output_len = output_len
        self.readout = nn.Sequential(
            nn.Conv1d(2 * d_model, hidden_channels, kernel_size=1),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Conv1d(hidden_channels, 1, kernel_size=1),
        )

    def forward(self, encoded: Tensor, local_feat: Tensor) -> Tensor:
        if encoded.dim() != 3:
            raise ValueError(f"Expected encoded shape [B, L, D], got {tuple(encoded.shape)}")
        if local_feat.dim() != 3:
            raise ValueError(f"Expected local_feat shape [B, D, L], got {tuple(local_feat.shape)}")

        encoded_cf = center_crop_1d(encoded.transpose(1, 2), self.output_len)
        local_cf = center_crop_1d(local_feat, self.output_len)
        fused = torch.cat([encoded_cf, local_cf], dim=1)
        return self.readout(fused).squeeze(1)
