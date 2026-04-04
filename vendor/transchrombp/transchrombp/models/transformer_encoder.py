"""Transformer encoder modules with Rotary Positional Embeddings (RoPE).

Design notes:
- Inputs are batch-first: [B, L, D].
- Attention is non-causal by default (genomic context is bidirectional).
- RoPE cache dynamically expands when sequence length exceeds cached size.
- Optional SDPA path uses torch.scaled_dot_product_attention for speed.
"""

from __future__ import annotations

from typing import Optional

import torch
import torch.nn.functional as F
from torch import Tensor, nn


class RotaryEmbedding(nn.Module):
    """Rotary positional embedding with dynamic cache growth."""

    def __init__(self, dim: int, max_seq_len: int = 32768, theta: float = 10000.0) -> None:
        super().__init__()
        if dim <= 0:
            raise ValueError(f"dim must be positive, got {dim}")
        if dim % 2 != 0:
            raise ValueError(f"RoPE dim must be even, got {dim}")
        if max_seq_len <= 0:
            raise ValueError(f"max_seq_len must be positive, got {max_seq_len}")
        if theta <= 0:
            raise ValueError(f"theta must be positive, got {theta}")

        self.dim = dim
        self.theta = float(theta)
        inv_freq = 1.0 / (self.theta ** (torch.arange(0, dim, 2, dtype=torch.float32) / dim))
        self.register_buffer("inv_freq", inv_freq, persistent=False)

        self.max_seq_len_cached = 0
        self.register_buffer("cos_cached", torch.empty(0), persistent=False)
        self.register_buffer("sin_cached", torch.empty(0), persistent=False)
        self._build_cache(max_seq_len)

    def _build_cache(self, seq_len: int) -> None:
        t = torch.arange(seq_len, device=self.inv_freq.device, dtype=self.inv_freq.dtype)
        freqs = torch.einsum("i,j->ij", t, self.inv_freq)
        emb = torch.cat((freqs, freqs), dim=-1)
        self.cos_cached = emb.cos()[None, None, :, :]
        self.sin_cached = emb.sin()[None, None, :, :]
        self.max_seq_len_cached = seq_len

    def forward(self, seq_len: int, device: torch.device, dtype: torch.dtype) -> tuple[Tensor, Tensor]:
        """Return cos/sin with shape [1, 1, seq_len, dim]."""
        if seq_len <= 0:
            raise ValueError(f"seq_len must be positive, got {seq_len}")
        if seq_len > self.max_seq_len_cached:
            self._build_cache(seq_len)

        cos = self.cos_cached[:, :, :seq_len, :].to(device=device, dtype=dtype)
        sin = self.sin_cached[:, :, :seq_len, :].to(device=device, dtype=dtype)
        return cos, sin


def rotate_half(x: Tensor) -> Tensor:
    """Rotate last-dim halves: [x1, x2] -> [-x2, x1]."""
    x1, x2 = x.chunk(2, dim=-1)
    return torch.cat((-x2, x1), dim=-1)


def apply_rotary_pos_emb(q: Tensor, k: Tensor, cos: Tensor, sin: Tensor) -> tuple[Tensor, Tensor]:
    """Apply RoPE to query and key tensors.

    Args:
        q, k: [B, H, L, Dh]
        cos, sin: [1, 1, L, Dh]
    """
    q_out = (q * cos) + (rotate_half(q) * sin)
    k_out = (k * cos) + (rotate_half(k) * sin)
    return q_out, k_out


class RoPEMultiheadAttention(nn.Module):
    """Batch-first multi-head self-attention with RoPE.

    Notes:
    - Non-causal by default.
    - key_padding_mask uses PyTorch convention: True means "ignore this token".
    """

    def __init__(
        self,
        d_model: int,
        n_heads: int,
        dropout: float = 0.1,
        rotary_emb: Optional[RotaryEmbedding] = None,
        use_sdpa: bool = True,
    ) -> None:
        super().__init__()
        if d_model <= 0:
            raise ValueError(f"d_model must be positive, got {d_model}")
        if n_heads <= 0:
            raise ValueError(f"n_heads must be positive, got {n_heads}")
        if d_model % n_heads != 0:
            raise ValueError(f"d_model ({d_model}) must be divisible by n_heads ({n_heads})")

        self.d_model = d_model
        self.n_heads = n_heads
        self.head_dim = d_model // n_heads
        self.scale = self.head_dim**-0.5

        self.q_proj = nn.Linear(d_model, d_model)
        self.k_proj = nn.Linear(d_model, d_model)
        self.v_proj = nn.Linear(d_model, d_model)
        self.out_proj = nn.Linear(d_model, d_model)

        self.dropout_p = float(dropout)
        self.dropout = nn.Dropout(dropout)
        self.rotary_emb = rotary_emb
        self.use_sdpa = bool(use_sdpa) and hasattr(F, "scaled_dot_product_attention")

    def _shape_qkv(self, x: Tensor) -> tuple[Tensor, Tensor, Tensor]:
        bsz, seq_len, _ = x.shape
        q = self.q_proj(x).view(bsz, seq_len, self.n_heads, self.head_dim).transpose(1, 2)
        k = self.k_proj(x).view(bsz, seq_len, self.n_heads, self.head_dim).transpose(1, 2)
        v = self.v_proj(x).view(bsz, seq_len, self.n_heads, self.head_dim).transpose(1, 2)
        return q, k, v

    def forward(self, x: Tensor, key_padding_mask: Optional[Tensor] = None) -> Tensor:
        """Forward self-attention.

        Args:
            x: [B, L, D]
            key_padding_mask: optional [B, L], True means token is padding/ignored.
        """
        if x.dim() != 3:
            raise ValueError(f"Expected x with shape [B, L, D], got {tuple(x.shape)}")

        bsz, seq_len, _ = x.shape
        if key_padding_mask is not None:
            if key_padding_mask.shape != (bsz, seq_len):
                raise ValueError(
                    "key_padding_mask shape mismatch: "
                    f"expected {(bsz, seq_len)}, got {tuple(key_padding_mask.shape)}"
                )
            key_padding_mask = key_padding_mask.to(dtype=torch.bool)
            if torch.any(torch.all(key_padding_mask, dim=1)):
                raise ValueError("Found sequence with all tokens masked in key_padding_mask")

        q, k, v = self._shape_qkv(x)

        if self.rotary_emb is not None:
            cos, sin = self.rotary_emb(seq_len=seq_len, device=q.device, dtype=q.dtype)
            q, k = apply_rotary_pos_emb(q, k, cos, sin)

        if self.use_sdpa:
            attn_mask = None
            if key_padding_mask is not None:
                # SDPA bool mask uses True = keep, False = mask.
                keep_mask = (~key_padding_mask).unsqueeze(1).unsqueeze(2)  # [B,1,1,L]
                attn_mask = keep_mask

            out = F.scaled_dot_product_attention(
                q,
                k,
                v,
                attn_mask=attn_mask,
                dropout_p=self.dropout_p if self.training else 0.0,
                is_causal=False,
            )
        else:
            attn_weights = torch.matmul(q, k.transpose(-2, -1)) * self.scale
            if key_padding_mask is not None:
                mask = key_padding_mask.unsqueeze(1).unsqueeze(2)  # [B,1,1,L]
                attn_weights = attn_weights.masked_fill(mask, float("-inf"))
            attn_weights = F.softmax(attn_weights, dim=-1)
            attn_weights = self.dropout(attn_weights)
            out = torch.matmul(attn_weights, v)

        out = out.transpose(1, 2).contiguous().view(bsz, seq_len, self.d_model)
        return self.out_proj(out)


class RoPETransformerLayer(nn.Module):
    """Pre-Norm transformer encoder layer with RoPE attention."""

    def __init__(
        self,
        d_model: int,
        n_heads: int,
        dim_feedforward: int,
        dropout: float,
        activation: str,
        rotary_emb: Optional[RotaryEmbedding] = None,
        use_sdpa: bool = True,
    ) -> None:
        super().__init__()
        self.norm1 = nn.LayerNorm(d_model)
        self.attn = RoPEMultiheadAttention(
            d_model=d_model,
            n_heads=n_heads,
            dropout=dropout,
            rotary_emb=rotary_emb,
            use_sdpa=use_sdpa,
        )
        self.dropout1 = nn.Dropout(dropout)

        self.norm2 = nn.LayerNorm(d_model)
        self.linear1 = nn.Linear(d_model, dim_feedforward)
        self.dropout = nn.Dropout(dropout)
        self.linear2 = nn.Linear(dim_feedforward, d_model)
        self.dropout2 = nn.Dropout(dropout)

        if activation == "relu":
            self.act = nn.ReLU()
        elif activation == "gelu":
            self.act = nn.GELU()
        else:
            raise ValueError(f"Unsupported activation: {activation}")

    def forward(self, x: Tensor, key_padding_mask: Optional[Tensor] = None) -> Tensor:
        x_norm = self.norm1(x)
        attn_out = self.attn(x_norm, key_padding_mask=key_padding_mask)
        x = x + self.dropout1(attn_out)

        x_norm = self.norm2(x)
        ff_out = self.linear2(self.dropout(self.act(self.linear1(x_norm))))
        x = x + self.dropout2(ff_out)
        return x


class SequenceTransformerEncoder(nn.Module):
    """Transformer encoder with optional RoPE and optional SDPA acceleration."""

    def __init__(
        self,
        d_model: int,
        n_heads: int,
        n_layers: int,
        ff_mult: int = 4,
        dropout: float = 0.1,
        activation: str = "gelu",
        use_positional_encoding: bool = True,
        max_len: int = 32768,
        rope_theta: float = 10000.0,
        use_sdpa: bool = True,
    ) -> None:
        super().__init__()
        if d_model <= 0:
            raise ValueError(f"d_model must be positive, got {d_model}")
        if n_heads <= 0:
            raise ValueError(f"n_heads must be positive, got {n_heads}")
        if d_model % n_heads != 0:
            raise ValueError(f"d_model ({d_model}) must be divisible by n_heads ({n_heads})")
        if n_layers <= 0:
            raise ValueError(f"n_layers must be positive, got {n_layers}")
        if ff_mult <= 0:
            raise ValueError(f"ff_mult must be positive, got {ff_mult}")

        head_dim = d_model // n_heads
        self.rotary_emb = (
            RotaryEmbedding(dim=head_dim, max_seq_len=max_len, theta=rope_theta)
            if use_positional_encoding
            else None
        )

        self.layers = nn.ModuleList(
            [
                RoPETransformerLayer(
                    d_model=d_model,
                    n_heads=n_heads,
                    dim_feedforward=d_model * ff_mult,
                    dropout=dropout,
                    activation=activation,
                    rotary_emb=self.rotary_emb,
                    use_sdpa=use_sdpa,
                )
                for _ in range(n_layers)
            ]
        )
        self.final_norm = nn.LayerNorm(d_model)

    def forward(self, x: Tensor, key_padding_mask: Optional[Tensor] = None) -> Tensor:
        if x.dim() != 3:
            raise ValueError(f"Expected x with shape [B, L, D], got {tuple(x.shape)}")

        for layer in self.layers:
            x = layer(x, key_padding_mask=key_padding_mask)
        return self.final_norm(x)
