"""ChromBPNet-style bias branch for TransChromBP."""

from __future__ import annotations

from typing import Tuple

from torch import Tensor, nn


def center_crop_1d(x: Tensor, target_len: int) -> Tensor:
    """Center-crop last dimension to target length."""
    if x.dim() < 1:
        raise ValueError("center_crop_1d expects tensor with at least 1 dimension")
    if target_len <= 0:
        raise ValueError(f"target_len must be positive, got {target_len}")

    seq_len = x.size(-1)
    if seq_len < target_len:
        raise ValueError(f"Cannot crop from length {seq_len} to larger target {target_len}")
    if seq_len == target_len:
        return x

    start = (seq_len - target_len) // 2
    end = start + target_len
    return x[..., start:end]


def _to_channel_first(seq_onehot: Tensor) -> Tensor:
    """Convert sequence tensor to [B, C, L] where C=4."""
    if seq_onehot.dim() != 3:
        raise ValueError(
            f"Expected sequence tensor with 3 dims [B, L, 4] or [B, 4, L], got {tuple(seq_onehot.shape)}"
        )

    if seq_onehot.size(1) == 4:
        return seq_onehot
    if seq_onehot.size(2) == 4:
        return seq_onehot.transpose(1, 2)

    raise ValueError(
        f"Unable to infer channel axis in sequence tensor with shape {tuple(seq_onehot.shape)}; expected a 4-channel axis"
    )


class ResidualDilatedConvBlock(nn.Module):
    """Simple residual dilated conv block."""

    def __init__(self, channels: int, kernel_size: int, dilation: int, dropout: float) -> None:
        super().__init__()
        if kernel_size % 2 == 0:
            raise ValueError("kernel_size must be odd for same-length padding")

        padding = (kernel_size // 2) * dilation
        self.block = nn.Sequential(
            nn.Conv1d(channels, channels, kernel_size=kernel_size, padding=padding, dilation=dilation),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Conv1d(channels, channels, kernel_size=1),
        )
        self.norm = nn.BatchNorm1d(channels)
        self.act = nn.GELU()

    def forward(self, x: Tensor) -> Tensor:
        residual = x
        x = self.block(x)
        x = self.norm(x + residual)
        return self.act(x)


class ChromBPNetBiasBranch(nn.Module):
    """Bias branch that predicts profile/count bias logits from sequence."""

    def __init__(
        self,
        input_channels: int = 4,
        hidden_channels: int = 128,
        kernel_size: int = 21,
        n_dil_layers: int = 4,
        output_len: int = 1000,
        dropout: float = 0.1,
        profile_pool_factor: int = 0,
    ) -> None:
        super().__init__()
        if kernel_size % 2 == 0:
            raise ValueError("kernel_size must be odd for same-length padding")
        if n_dil_layers <= 0:
            raise ValueError(f"n_dil_layers must be positive, got {n_dil_layers}")

        self.output_len = output_len
        self.stem = nn.Sequential(
            nn.Conv1d(
                in_channels=input_channels,
                out_channels=hidden_channels,
                kernel_size=kernel_size,
                padding=kernel_size // 2,
            ),
            nn.GELU(),
            nn.Dropout(dropout),
        )

        blocks = []
        for i in range(n_dil_layers):
            dilation = 2**i
            blocks.append(
                ResidualDilatedConvBlock(
                    channels=hidden_channels,
                    kernel_size=3,
                    dilation=dilation,
                    dropout=dropout,
                )
            )
        self.dilated_tower = nn.Sequential(*blocks)

        self.profile_head = nn.Conv1d(hidden_channels, 1, kernel_size=1)
        self.count_pool = nn.AdaptiveAvgPool1d(1)
        self.count_head = nn.Linear(hidden_channels, 1)

        # Resolution bottleneck: limits profile output to low-frequency patterns.
        # When profile_pool_factor > 0, the profile is downsampled then upsampled,
        # preventing the bias branch from representing fine-grained motif shapes.
        self.profile_pool_factor = int(profile_pool_factor)
        if self.profile_pool_factor > 0:
            bottleneck_len = max(1, output_len // self.profile_pool_factor)
            self.profile_bottleneck = nn.Sequential(
                nn.AdaptiveAvgPool1d(bottleneck_len),
                nn.Upsample(size=output_len, mode="linear", align_corners=False),
            )
        else:
            self.profile_bottleneck = None

    def forward(self, seq_onehot: Tensor) -> Tuple[Tensor, Tensor]:
        """Forward pass.

        Args:
            seq_onehot: [B, L, 4] or [B, 4, L]

        Returns:
            profile_bias_logits: [B, output_len]
            count_bias: [B, 1]
        """
        x = _to_channel_first(seq_onehot)
        x = self.stem(x)
        x = self.dilated_tower(x)

        profile_bias_logits = self.profile_head(x).squeeze(1)
        profile_bias_logits = center_crop_1d(profile_bias_logits, self.output_len)

        if self.profile_bottleneck is not None:
            profile_bias_logits = self.profile_bottleneck(
                profile_bias_logits.unsqueeze(1)
            ).squeeze(1)

        pooled = self.count_pool(x).squeeze(-1)
        count_bias = self.count_head(pooled)

        return profile_bias_logits, count_bias

    def freeze_core(self) -> None:
        """Freeze all bias branch parameters."""
        for p in self.parameters():
            p.requires_grad = False

    def unfreeze_core(self) -> None:
        """Unfreeze all bias branch parameters."""
        for p in self.parameters():
            p.requires_grad = True
