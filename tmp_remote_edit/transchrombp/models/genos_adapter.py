"""Genos-1.2B feature extraction and gated fusion adapter for TransChromBP.

Provides:
- GenosFeatureExtractor: frozen Genos wrapper that extracts hidden-state
  features from DNA sequences (runs outside the trainable model, under
  torch.no_grad).  NOT an nn.Module — lives outside DDP and optimizer.
- GenosGatedAdapter: small trainable module that fuses Genos features into
  the TransChromBP signal branch via learned gating.
"""

from __future__ import annotations

import torch
from torch import Tensor, nn


# ─── One-hot → DNA string helpers ───────────────────────────

_ONEHOT_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T"}
_RC_MAP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def onehot_to_dna(seq_onehot: Tensor) -> list[str]:
    """Convert a batch of one-hot sequences to DNA strings.

    Args:
        seq_onehot: [B, L, 4] one-hot encoded DNA.
                    All-zero positions map to 'N'.
    """
    indices = seq_onehot.argmax(dim=-1).cpu().numpy()       # [B, L]
    valid = (seq_onehot.sum(dim=-1) > 0).cpu().numpy()      # [B, L]
    result: list[str] = []
    for b in range(indices.shape[0]):
        chars = []
        for i in range(indices.shape[1]):
            chars.append(_ONEHOT_TO_BASE[indices[b, i]] if valid[b, i] else "N")
        result.append("".join(chars))
    return result


def reverse_complement(seq: str) -> str:
    return "".join(_RC_MAP.get(c, "N") for c in reversed(seq))


# ─── Frozen Genos Feature Extractor ─────────────────────────

class GenosFeatureExtractor:
    """Frozen Genos-1.2B wrapper for extracting intermediate hidden states.

    NOT an nn.Module.  Lives outside the trainable model and DDP.
    Each DDP rank holds its own copy in eval() + no_grad mode.
    """

    def __init__(
        self,
        model_path: str,
        layer: int = 6,
        device: torch.device | str = "cuda",
        attn_implementation: str = "flash_attention_2",
        bidirectional: bool = True,
    ) -> None:
        from transformers import AutoModelForCausalLM, AutoTokenizer

        self.layer = layer
        self.device = torch.device(device)
        self.bidirectional = bidirectional

        self.tokenizer = AutoTokenizer.from_pretrained(
            model_path, local_files_only=True,
        )
        self.model = AutoModelForCausalLM.from_pretrained(
            model_path,
            local_files_only=True,
            dtype=torch.bfloat16,
            attn_implementation=attn_implementation,
        ).to(self.device).eval()

        for p in self.model.parameters():
            p.requires_grad_(False)

        self.hidden_size: int = self.model.config.hidden_size

    # ── internal helpers ──

    def _tokenize(self, sequences: list[str]) -> dict[str, Tensor]:
        tokens = self.tokenizer(
            sequences,
            return_tensors="pt",
            padding=True,
            add_special_tokens=False,
        )
        return {k: v.to(self.device) for k, v in tokens.items()}

    def _forward_hidden(self, token_inputs: dict[str, Tensor]) -> Tensor:
        outputs = self.model(**token_inputs, output_hidden_states=True)
        # hidden_states: tuple of (num_layers + 1) tensors
        # index 0 = embedding output, index i = after transformer layer i
        hidden = outputs.hidden_states[self.layer]  # [B, L, H]
        return hidden.float()

    # ── public API ──

    @torch.no_grad()
    def extract(self, seq_onehot: Tensor) -> Tensor:
        """Extract Genos features from one-hot DNA sequences.

        Args:
            seq_onehot: [B, L, 4] or [B, 4, L] one-hot sequences on any device.
        Returns:
            [B, L, hidden_size] float32 features on the same device as self.device.
        """
        # ensure [B, L, 4]
        if seq_onehot.dim() == 3 and seq_onehot.size(1) == 4 and seq_onehot.size(2) != 4:
            seq_onehot = seq_onehot.transpose(1, 2)

        dna_strings = onehot_to_dna(seq_onehot)

        fwd_inputs = self._tokenize(dna_strings)
        fwd_hidden = self._forward_hidden(fwd_inputs)

        if not self.bidirectional:
            return fwd_hidden

        rc_strings = [reverse_complement(s) for s in dna_strings]
        rc_inputs = self._tokenize(rc_strings)
        rc_hidden = self._forward_hidden(rc_inputs)
        rc_hidden = rc_hidden.flip(dims=[1])

        return (fwd_hidden + rc_hidden) * 0.5


# ─── Trainable Gated Fusion Adapter ─────────────────────────

class GenosGatedAdapter(nn.Module):
    """Gated additive fusion of Genos features into the signal branch.

    ``fused = tokens + gate * proj(genos_feat)``

    The gate is computed from *local tokens* and initialized near zero
    (``gate_bias_init = -2.0`` ⇒ sigmoid ≈ 0.12), so the model starts
    close to the Genos-free baseline.
    """

    def __init__(
        self,
        d_model: int = 256,
        genos_hidden_size: int = 1024,
        gate_bias_init: float = -2.0,
        pool_mode: str = "none",
    ) -> None:
        super().__init__()
        if pool_mode not in {"none", "mean"}:
            raise ValueError(f"Unsupported genos pool_mode={pool_mode!r}")
        self.pool_mode = pool_mode

        self.proj = nn.Linear(genos_hidden_size, d_model)
        self.ln = nn.LayerNorm(d_model)
        self.gate_proj = nn.Linear(d_model, d_model)
        nn.init.constant_(self.gate_proj.bias, gate_bias_init)

    def forward(self, tokens: Tensor, genos_feat: Tensor) -> Tensor:
        """
        Args:
            tokens:     [B, L, d_model]  from conv_stem / local_tower.
            genos_feat: [B, L, genos_hidden_size]  from GenosFeatureExtractor.
        Returns:
            [B, L, d_model] fused tokens.
        """
        if self.pool_mode == "mean":
            genos_feat = genos_feat.mean(dim=1, keepdim=True).expand_as(genos_feat)

        projected = self.ln(self.proj(genos_feat))          # [B, L, d_model]
        gate = torch.sigmoid(self.gate_proj(tokens))        # [B, L, d_model]
        return tokens + gate * projected


# ─── Cached Summary FiLM Adapter ──────────────────────────

class GenosSummaryFiLM(nn.Module):
    """Feature-wise Linear Modulation (FiLM) using a Genos global summary.

    Applies channel-wise affine modulation to the Transformer output:
        ``output = (1 + gamma) * encoded + beta``

    Zero-initialized so the model starts identical to baseline (gamma=0, beta=0).
    """

    def __init__(self, d_model: int = 256, genos_dim: int = 1024) -> None:
        super().__init__()
        self.proj = nn.Linear(genos_dim, d_model * 2)
        nn.init.zeros_(self.proj.weight)
        nn.init.zeros_(self.proj.bias)
        self._last_forward_stats: dict[str, Tensor] = {}

    def forward(self, encoded: Tensor, genos_summary: Tensor) -> Tensor:
        """
        Args:
            encoded:       [B, L, d_model]  Transformer output.
            genos_summary: [B, genos_dim]   Global-mean Genos feature.
        Returns:
            [B, L, d_model] modulated encoded.
        """
        params = self.proj(genos_summary)               # [B, d_model * 2]
        gamma_delta, beta = params.chunk(2, dim=-1)      # [B, d_model] each
        gamma = 1.0 + gamma_delta.unsqueeze(1)           # [B, 1, d_model]
        beta = beta.unsqueeze(1)                         # [B, 1, d_model]
        modulated = gamma * encoded + beta

        with torch.no_grad():
            delta = modulated - encoded
            self._last_forward_stats = {
                "film_gamma_mean": gamma.detach().float().mean(),
                "film_gamma_std": gamma.detach().float().std(unbiased=False),
                "film_beta_rms": beta.detach().float().pow(2).mean().sqrt(),
                "film_delta_rms": delta.detach().float().pow(2).mean().sqrt(),
            }

        return modulated


class GenosCountProj(nn.Module):
    """Additive projection of Genos global summary into the count head hidden layer.

    Zero-initialized so the model starts identical to baseline.
    """

    def __init__(self, d_count_hidden: int = 128, genos_dim: int = 1024) -> None:
        super().__init__()
        self.proj = nn.Linear(genos_dim, d_count_hidden)
        nn.init.zeros_(self.proj.weight)
        nn.init.zeros_(self.proj.bias)
        self._last_forward_stats: dict[str, Tensor] = {}

    def forward(self, count_hidden: Tensor, genos_summary: Tensor) -> Tensor:
        """
        Args:
            count_hidden:  [B, d_count_hidden]  Intermediate count head activation.
            genos_summary: [B, genos_dim]       Global-mean Genos feature.
        Returns:
            [B, d_count_hidden] with additive Genos contribution.
        """
        delta = self.proj(genos_summary)
        with torch.no_grad():
            self._last_forward_stats = {
                "genos_count_delta_rms": delta.detach().float().pow(2).mean().sqrt(),
            }
        return count_hidden + delta
