"""Caduceus-PS feature extraction and token-level fusion for TransChromBP."""

from __future__ import annotations

from typing import Any

import torch
from torch import Tensor, nn

from .genos_adapter import onehot_to_dna


def _resolve_hidden_size(config: Any) -> int:
    base_hidden_size = None
    for attr in ("hidden_size", "d_model", "n_embd"):
        value = getattr(config, attr, None)
        if value is not None:
            base_hidden_size = int(value)
            break
    if base_hidden_size is None:
        raise AttributeError("Could not infer hidden size from Caduceus config")

    # RCPS models concatenate forward/reverse streams in the exposed token states.
    if bool(getattr(config, "rcps", False)):
        return base_hidden_size * 2
    return base_hidden_size


class CaduceusFeatureExtractor:
    """Frozen Caduceus wrapper for extracting token hidden states.

    Caduceus-PS is RC-equivariant, so we intentionally do not perform
    forward + reverse-complement averaging here.
    """

    def __init__(
        self,
        model_path: str,
        layer: int = -1,
        device: torch.device | str = "cuda",
        trust_remote_code: bool = True,
        local_files_only: bool = True,
        dtype: str = "bfloat16",
    ) -> None:
        from transformers import AutoModel, AutoModelForMaskedLM, AutoTokenizer

        self.layer = int(layer)
        self.device = torch.device(device)
        self.trust_remote_code = bool(trust_remote_code)
        self.local_files_only = bool(local_files_only)

        dtype_name = str(dtype).lower()
        if dtype_name in {"bf16", "bfloat16"}:
            torch_dtype = torch.bfloat16
        elif dtype_name in {"fp16", "float16", "half"}:
            torch_dtype = torch.float16
        elif dtype_name in {"fp32", "float32"}:
            torch_dtype = torch.float32
        else:
            raise ValueError(f"Unsupported Caduceus dtype={dtype!r}")

        self.tokenizer = AutoTokenizer.from_pretrained(
            model_path,
            local_files_only=self.local_files_only,
            trust_remote_code=self.trust_remote_code,
        )
        if self.tokenizer.pad_token_id is None:
            if self.tokenizer.eos_token_id is not None:
                self.tokenizer.pad_token = self.tokenizer.eos_token
            elif self.tokenizer.unk_token_id is not None:
                self.tokenizer.pad_token = self.tokenizer.unk_token
            else:
                raise ValueError("Caduceus tokenizer is missing pad/eos/unk tokens; cannot batch-pad inputs")

        model_kwargs = dict(
            local_files_only=self.local_files_only,
            trust_remote_code=self.trust_remote_code,
            torch_dtype=torch_dtype,
        )
        try:
            self.model = AutoModel.from_pretrained(model_path, **model_kwargs)
        except Exception:
            self.model = AutoModelForMaskedLM.from_pretrained(model_path, **model_kwargs)
        self.model = self.model.to(self.device).eval()

        for p in self.model.parameters():
            p.requires_grad_(False)

        self.rcps = bool(getattr(self.model.config, "rcps", False))
        self.hidden_size = _resolve_hidden_size(self.model.config)

    def _tokenize(self, sequences: list[str]) -> dict[str, Tensor]:
        tokens = self.tokenizer(
            sequences,
            return_tensors="pt",
            padding=True,
            add_special_tokens=False,
        )
        return {k: v.to(self.device) for k, v in tokens.items()}

    @torch.no_grad()
    def extract(self, seq_onehot: Tensor) -> Tensor:
        """Extract token-level hidden states with shape [B, L, H]."""
        if seq_onehot.dim() != 3:
            raise ValueError(f"Expected [B, L, 4] or [B, 4, L], got {tuple(seq_onehot.shape)}")
        if seq_onehot.size(1) == 4 and seq_onehot.size(2) != 4:
            seq_onehot = seq_onehot.transpose(1, 2)
        if seq_onehot.size(-1) != 4:
            raise ValueError(f"Expected one-hot channel dim=4, got {tuple(seq_onehot.shape)}")

        seq_len = int(seq_onehot.size(1))
        token_inputs = self._tokenize(onehot_to_dna(seq_onehot))
        outputs = self.model(**token_inputs, output_hidden_states=True)

        hidden_states = getattr(outputs, "hidden_states", None)
        if hidden_states is None:
            last_hidden = getattr(outputs, "last_hidden_state", None)
            if last_hidden is None:
                raise RuntimeError("Caduceus model output does not expose hidden_states or last_hidden_state")
            hidden = last_hidden
        else:
            hidden = hidden_states[self.layer]

        if hidden.size(1) != seq_len:
            raise ValueError(
                "Caduceus token length mismatch: "
                f"expected {seq_len}, got {int(hidden.size(1))}. "
                "Check tokenizer special-token behavior."
            )
        return hidden.float()


class CaduceusTokenAdapter(nn.Module):
    """Residual gated fusion of Caduceus token features before the transformer."""

    def __init__(
        self,
        d_model: int = 256,
        caduceus_hidden_size: int = 256,
        gate_bias_init: float = -2.0,
    ) -> None:
        super().__init__()
        self.proj = nn.Linear(caduceus_hidden_size, d_model)
        self.norm = nn.LayerNorm(d_model)
        self.gate_proj = nn.Linear(d_model, d_model)

        nn.init.zeros_(self.proj.weight)
        nn.init.zeros_(self.proj.bias)
        nn.init.constant_(self.gate_proj.bias, gate_bias_init)

        self._last_forward_stats: dict[str, Tensor] = {}

    def forward(self, tokens: Tensor, caduceus_feat: Tensor) -> Tensor:
        if tokens.shape[:2] != caduceus_feat.shape[:2]:
            raise ValueError(
                "Caduceus feature shape mismatch: "
                f"tokens={tuple(tokens.shape)} caduceus_feat={tuple(caduceus_feat.shape)}"
            )

        projected = self.norm(self.proj(caduceus_feat))
        gate = torch.sigmoid(self.gate_proj(tokens))
        delta = gate * projected
        fused = tokens + delta

        with torch.no_grad():
            self._last_forward_stats = {
                "caduceus_gate_mean": gate.detach().float().mean(),
                "caduceus_gate_std": gate.detach().float().std(unbiased=False),
                "caduceus_delta_rms": delta.detach().float().pow(2).mean().sqrt(),
            }
        return fused
