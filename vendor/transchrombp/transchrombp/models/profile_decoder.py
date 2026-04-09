"""Profile decoder variants for TransChromBP."""

from __future__ import annotations

import math
from typing import Sequence

import torch.nn.functional as F
from torch import Tensor, nn

from .bias_branch import center_crop_1d


class ConvRefineBlock(nn.Module):
    """Refine a decoder state with a same-resolution local skip tensor."""

    def __init__(
        self,
        in_channels: int,
        skip_channels: int,
        out_channels: int,
        kernel_size: int = 3,
        dropout: float = 0.0,
    ) -> None:
        super().__init__()
        if kernel_size % 2 == 0:
            raise ValueError("kernel_size must be odd for same-length padding")

        padding = kernel_size // 2
        self.input_proj = (
            nn.Identity()
            if in_channels == out_channels
            else nn.Conv1d(in_channels, out_channels, kernel_size=1)
        )
        self.skip_proj = nn.Conv1d(skip_channels, out_channels, kernel_size=1)
        self.block = nn.Sequential(
            nn.GELU(),
            nn.Conv1d(out_channels, out_channels, kernel_size=kernel_size, padding=padding),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Conv1d(out_channels, out_channels, kernel_size=kernel_size, padding=padding),
        )
        self.norm = nn.BatchNorm1d(out_channels)
        self.act = nn.GELU()

    def forward(self, x: Tensor, skip: Tensor) -> Tensor:
        if x.size(-1) != skip.size(-1):
            raise ValueError(
                f"ConvRefineBlock requires matched sequence lengths, got {x.size(-1)} and {skip.size(-1)}"
            )

        fused = self.input_proj(x) + self.skip_proj(skip)
        refined = self.block(fused)
        return self.act(self.norm(refined + fused))


class MultiScaleLocalSkipDecoderV2(nn.Module):
    """Decode transformer tokens with multiscale local skip refinement."""

    def __init__(
        self,
        encoded_channels: int,
        local_channels: int,
        output_len: int,
        hidden_channels: int,
        decoder_channels: Sequence[int],
        dropout: float = 0.0,
        upsample_mode: str = "linear",
    ) -> None:
        super().__init__()
        if output_len <= 0:
            raise ValueError(f"output_len must be positive, got {output_len}")
        if hidden_channels <= 0:
            raise ValueError(f"hidden_channels must be positive, got {hidden_channels}")

        stage_channels = [int(channels) for channels in decoder_channels]
        if not stage_channels:
            raise ValueError("decoder_channels must contain at least one stage")
        if any(channels <= 0 for channels in stage_channels):
            raise ValueError(f"decoder_channels must be positive, got {stage_channels}")

        mode = str(upsample_mode).strip().lower()
        if mode not in {"linear", "nearest"}:
            raise ValueError(f"Unsupported upsample_mode={upsample_mode!r}; expected 'linear' or 'nearest'")

        self.output_len = output_len
        self.upsample_mode = mode
        self.input_proj = nn.Conv1d(encoded_channels, hidden_channels, kernel_size=1)
        self.num_stages = len(stage_channels)
        self.refine_blocks = nn.ModuleList()

        in_channels = hidden_channels
        for out_channels in stage_channels:
            self.refine_blocks.append(
                ConvRefineBlock(
                    in_channels=in_channels,
                    skip_channels=local_channels,
                    out_channels=out_channels,
                    dropout=dropout,
                )
            )
            in_channels = out_channels

        self.profile_head = nn.Sequential(
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Conv1d(in_channels, 1, kernel_size=1),
        )

    def _resize(self, x: Tensor, target_len: int) -> Tensor:
        if x.size(-1) == target_len:
            return x
        if target_len < x.size(-1):
            return F.adaptive_avg_pool1d(x, target_len)
        if self.upsample_mode == "linear":
            return F.interpolate(x, size=target_len, mode="linear", align_corners=False)
        return F.interpolate(x, size=target_len, mode=self.upsample_mode)

    def _stage_lengths(self, full_len: int) -> list[int]:
        return [
            int(math.ceil(full_len / (2 ** (self.num_stages - index - 1))))
            for index in range(self.num_stages)
        ]

    def forward(self, encoded: Tensor, local_feat: Tensor) -> Tensor:
        if encoded.dim() != 3:
            raise ValueError(f"Expected encoded shape [B, L, D], got {tuple(encoded.shape)}")
        if local_feat.dim() != 3:
            raise ValueError(f"Expected local_feat shape [B, D, L], got {tuple(local_feat.shape)}")
        if encoded.size(1) != local_feat.size(-1):
            raise ValueError(
                "MultiScaleLocalSkipDecoderV2 requires encoded/local_feat to share sequence length, "
                f"got {encoded.size(1)} and {local_feat.size(-1)}"
            )

        full_len = encoded.size(1)
        current = self.input_proj(encoded.transpose(1, 2))
        for target_len, block in zip(self._stage_lengths(full_len), self.refine_blocks):
            current = self._resize(current, target_len)
            local_stage = self._resize(local_feat, target_len)
            current = block(current, local_stage)

        current = self._resize(current, full_len)
        profile_logits = self.profile_head(current).squeeze(1)
        return center_crop_1d(profile_logits, self.output_len)
