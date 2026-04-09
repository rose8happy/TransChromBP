"""Generic foundation-model adapters for TransChromBP."""

from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import Tensor, nn

from .bias_branch import center_crop_1d


class FoundationCrossAttentionAdapter(nn.Module):
    """Residual gated cross-attention from local tokens to external foundation tokens."""

    def __init__(
        self,
        d_model: int = 256,
        foundation_hidden_size: int = 1024,
        n_heads: int = 4,
        dropout: float = 0.1,
        gate_bias_init: float = -4.0,
    ) -> None:
        super().__init__()
        if d_model <= 0 or foundation_hidden_size <= 0:
            raise ValueError("d_model and foundation_hidden_size must be positive")
        if n_heads <= 0 or d_model % n_heads != 0:
            raise ValueError(f"d_model={d_model} must be divisible by n_heads={n_heads}")

        self.foundation_proj = nn.Linear(foundation_hidden_size, d_model)
        self.query_norm = nn.LayerNorm(d_model)
        self.kv_norm = nn.LayerNorm(d_model)
        self.attn = nn.MultiheadAttention(
            embed_dim=d_model,
            num_heads=n_heads,
            dropout=dropout,
            batch_first=True,
        )
        self.out_proj = nn.Linear(d_model, d_model)
        self.gate_proj = nn.Linear(d_model, d_model)

        nn.init.zeros_(self.out_proj.weight)
        nn.init.zeros_(self.out_proj.bias)
        nn.init.constant_(self.gate_proj.bias, gate_bias_init)
        self._last_forward_stats: dict[str, Tensor] = {}

    def forward(self, tokens: Tensor, foundation_tokens: Tensor) -> Tensor:
        if tokens.dim() != 3 or foundation_tokens.dim() != 3:
            raise ValueError(
                "Expected [B, L, D] tokens and [B, M, H] foundation_tokens, got "
                f"{tuple(tokens.shape)} and {tuple(foundation_tokens.shape)}"
            )
        if tokens.size(0) != foundation_tokens.size(0):
            raise ValueError(
                "Batch size mismatch between local and foundation tokens: "
                f"{tokens.size(0)} vs {foundation_tokens.size(0)}"
            )

        q = self.query_norm(tokens)
        kv = self.kv_norm(self.foundation_proj(foundation_tokens))
        attn_out, _ = self.attn(q, kv, kv, need_weights=False)
        delta = self.out_proj(attn_out)
        gate = torch.sigmoid(self.gate_proj(tokens))
        fused = tokens + gate * delta

        with torch.no_grad():
            self._last_forward_stats = {
                "foundation_cross_gate_mean": gate.detach().float().mean(),
                "foundation_cross_gate_std": gate.detach().float().std(unbiased=False),
                "foundation_cross_delta_rms": delta.detach().float().pow(2).mean().sqrt(),
            }
        return fused


class FoundationResidualHead(nn.Module):
    """Predict coarse residual corrections from encoded tokens and external foundation features."""

    def __init__(
        self,
        d_model: int = 256,
        foundation_hidden_size: int = 1024,
        hidden_size: int = 256,
        profile_bin_count: int = 16,
        output_len: int = 1000,
        alignment_mode: str = "summary",
        aligned_token_count: int = 0,
    ) -> None:
        super().__init__()
        if d_model <= 0 or foundation_hidden_size <= 0 or hidden_size <= 0:
            raise ValueError("All hidden sizes must be positive")
        if profile_bin_count <= 0 or output_len <= 0:
            raise ValueError("profile_bin_count and output_len must be positive")
        alignment_mode = str(alignment_mode).strip().lower()
        if alignment_mode not in {"summary", "center_token_bins"}:
            raise ValueError(
                "alignment_mode must be one of {'summary', 'center_token_bins'}, got "
                f"{alignment_mode!r}"
            )
        if alignment_mode == "center_token_bins" and aligned_token_count <= 0:
            raise ValueError("center_token_bins alignment requires aligned_token_count > 0")

        self.output_len = int(output_len)
        self.profile_bin_count = int(profile_bin_count)
        self.alignment_mode = alignment_mode
        self.aligned_token_count = int(aligned_token_count)
        self.profile_bins_per_token = 0
        if self.alignment_mode == "center_token_bins":
            self.profile_bins_per_token = max(
                1,
                (self.profile_bin_count + self.aligned_token_count - 1) // self.aligned_token_count,
            )

        self.encoded_norm = nn.LayerNorm(d_model)
        self.foundation_norm = nn.LayerNorm(foundation_hidden_size)
        self.encoded_proj = nn.Linear(d_model, hidden_size)
        self.foundation_proj = nn.Linear(foundation_hidden_size, hidden_size)
        self.fuse = nn.Linear(hidden_size * 2, hidden_size)
        profile_out_dim = self.profile_bins_per_token if self.profile_bins_per_token > 0 else profile_bin_count
        self.profile_out = nn.Linear(hidden_size, profile_out_dim)
        self.count_out = nn.Linear(hidden_size, 1)

        nn.init.zeros_(self.profile_out.weight)
        nn.init.zeros_(self.profile_out.bias)
        nn.init.zeros_(self.count_out.weight)
        nn.init.zeros_(self.count_out.bias)
        self._last_forward_stats: dict[str, Tensor] = {}

    def forward(
        self,
        encoded: Tensor,
        foundation_summary: Tensor | None = None,
        foundation_tokens: Tensor | None = None,
    ) -> tuple[Tensor, Tensor]:
        if encoded.dim() != 3:
            raise ValueError(f"Expected encoded shape [B, L, D], got {tuple(encoded.shape)}")

        if self.alignment_mode == "center_token_bins":
            if foundation_tokens is None:
                raise ValueError("center_token_bins residual requires foundation_tokens")
            if foundation_tokens.dim() != 3:
                raise ValueError(f"Expected foundation_tokens [B, M, H], got {tuple(foundation_tokens.shape)}")
            if foundation_tokens.size(1) != self.aligned_token_count:
                raise ValueError(
                    "center_token_bins residual expects foundation_tokens with "
                    f"{self.aligned_token_count} bins, got {foundation_tokens.size(1)}"
                )

            encoded_bins = self.encoded_norm(encoded).transpose(1, 2)
            encoded_bins = center_crop_1d(encoded_bins, self.output_len)
            encoded_bins = F.adaptive_avg_pool1d(encoded_bins, output_size=self.aligned_token_count).transpose(1, 2)
            foundation_tokens = self.foundation_norm(foundation_tokens)

            fused_tokens = torch.cat(
                [
                    self.encoded_proj(encoded_bins),
                    self.foundation_proj(foundation_tokens),
                ],
                dim=-1,
            )
            fused_tokens = F.gelu(self.fuse(fused_tokens))
            profile_bins = self.profile_out(fused_tokens).reshape(encoded.size(0), -1)
            profile_bins = profile_bins[:, : self.profile_bin_count]
            count_delta = self.count_out(fused_tokens.mean(dim=1))
        else:
            if foundation_summary is None:
                if foundation_tokens is None:
                    raise ValueError("FoundationResidualHead requires foundation_summary or foundation_tokens")
                if foundation_tokens.dim() != 3:
                    raise ValueError(f"Expected foundation_tokens [B, M, H], got {tuple(foundation_tokens.shape)}")
                foundation_summary = foundation_tokens.mean(dim=1)
            elif foundation_summary.dim() != 2:
                raise ValueError(f"Expected foundation_summary [B, H], got {tuple(foundation_summary.shape)}")

            encoded_summary = self.encoded_norm(encoded).mean(dim=1)
            foundation_summary = self.foundation_norm(foundation_summary)
            fused = torch.cat(
                [
                    self.encoded_proj(encoded_summary),
                    self.foundation_proj(foundation_summary),
                ],
                dim=-1,
            )
            fused = F.gelu(self.fuse(fused))
            profile_bins = self.profile_out(fused)
            count_delta = self.count_out(fused)

        profile_delta = F.interpolate(
            profile_bins.unsqueeze(1),
            size=self.output_len,
            mode="linear",
            align_corners=False,
        ).squeeze(1)

        with torch.no_grad():
            self._last_forward_stats = {
                "foundation_residual_profile_rms": profile_delta.detach().float().pow(2).mean().sqrt(),
                "foundation_residual_count_rms": count_delta.detach().float().pow(2).mean().sqrt(),
            }
        return profile_delta, count_delta
