#!/usr/bin/env python3
"""Minimal local smoke test for NT v2 loading and hidden-state extraction."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import torch


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-dir", required=True, help="Local NT v2 model directory.")
    parser.add_argument("--genome-fasta", default="", help="Reference FASTA for a real-window smoke.")
    parser.add_argument("--chrom", default="", help="Chromosome for a real-window smoke.")
    parser.add_argument("--center", type=int, default=0, help="Center coordinate for a real-window smoke.")
    parser.add_argument("--input-len", type=int, default=2114, help="Length of the real genomic window.")
    parser.add_argument("--device", default="cuda", help="Torch device.")
    parser.add_argument("--output-json", default="", help="Optional path to write the smoke summary JSON.")
    return parser.parse_args()


def read_real_window(genome_fasta: str, chrom: str, center: int, input_len: int) -> str:
    if not genome_fasta or not chrom or center <= 0:
        return ("ACGT" * ((input_len // 4) + 1))[:input_len]

    from pyfaidx import Fasta

    fasta = Fasta(genome_fasta)
    start = max(0, center - input_len // 2)
    end = start + input_len
    seq = str(fasta[chrom][start:end]).upper()
    if len(seq) < input_len:
        seq = seq + ("N" * (input_len - len(seq)))
    return seq


def inspect_intermediate_shapes(model_dir: Path) -> dict[str, Any]:
    summary: dict[str, Any] = {}

    safe_path = model_dir / "model.safetensors"
    if safe_path.exists():
        try:
            from safetensors import safe_open

            with safe_open(str(safe_path), framework="pt", device="cpu") as f:
                for key in f.keys():
                    if "intermediate" in key and key.endswith("weight"):
                        summary["safetensors_intermediate_key"] = key
                        summary["safetensors_intermediate_shape"] = list(f.get_tensor(key).shape)
                        break
        except Exception as exc:  # pragma: no cover - best effort diagnostics
            summary["safetensors_error"] = repr(exc)

    bin_path = model_dir / "pytorch_model.bin"
    if bin_path.exists():
        try:
            state = torch.load(bin_path, map_location="cpu")
            state_dict = state.get("state_dict", state)
            if isinstance(state_dict, dict):
                for key, value in state_dict.items():
                    if "intermediate" in key and key.endswith("weight"):
                        summary["bin_intermediate_key"] = key
                        summary["bin_intermediate_shape"] = list(value.shape)
                        break
        except Exception as exc:  # pragma: no cover - best effort diagnostics
            summary["bin_error"] = repr(exc)

    return summary


def run_forward(tokenizer, model, sequences: list[str], label: str) -> dict[str, Any]:
    batch = tokenizer.batch_encode_plus(
        sequences,
        return_tensors="pt",
        padding=True,
        truncation=False,
    )
    attention_mask = batch["attention_mask"].bool()

    device = next(model.parameters()).device
    batch = {k: v.to(device) for k, v in batch.items()}
    attention_mask = attention_mask.to(device)

    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    with torch.no_grad():
        outputs = model(
            batch["input_ids"],
            attention_mask=attention_mask,
            encoder_attention_mask=attention_mask,
            output_hidden_states=True,
        )
    elapsed = time.perf_counter() - started
    hidden = outputs["hidden_states"][-1]
    mask = attention_mask.unsqueeze(-1)
    pooled = (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp_min(1)

    result = {
        "label": label,
        "input_ids_shape": list(batch["input_ids"].shape),
        "attention_mask_shape": list(attention_mask.shape),
        "last_hidden_state_shape": list(hidden.shape),
        "mean_sequence_embedding_shape": list(pooled.shape),
        "elapsed_sec": elapsed,
    }
    if device.type == "cuda":
        result["peak_memory_mb"] = round(torch.cuda.max_memory_allocated(device) / (1024 ** 2), 3)
    return result


def main() -> None:
    args = parse_args()
    model_dir = Path(args.model_dir).resolve()
    summary: dict[str, Any] = {
        "model_dir": str(model_dir),
        "device": args.device,
        "status": "started",
    }

    try:
        from transformers import AutoModelForMaskedLM, AutoTokenizer

        tokenizer = AutoTokenizer.from_pretrained(
            model_dir,
            trust_remote_code=True,
            local_files_only=True,
        )
        model = AutoModelForMaskedLM.from_pretrained(
            model_dir,
            trust_remote_code=True,
            local_files_only=True,
        )
        model = model.to(args.device).eval()

        short_sequences = [
            "ATTCCGATTCCGATTCCG",
            "ATTTCTCTCTCTCTCTGAGATCGATCGATCGAT",
        ]
        real_window = read_real_window(args.genome_fasta, args.chrom, args.center, args.input_len)

        summary["short_smoke"] = run_forward(tokenizer, model, short_sequences, "short_sequences")
        summary["real_window_smoke"] = run_forward(tokenizer, model, [real_window], "real_window")
        summary["model_hidden_size"] = int(getattr(model.config, "hidden_size", -1))
        summary["model_layers"] = int(getattr(model.config, "num_hidden_layers", -1))
        summary["status"] = "ok"
    except Exception as exc:
        summary["status"] = "error"
        summary["error"] = repr(exc)
        summary["config_intermediate_size"] = json.loads((model_dir / "config.json").read_text())["intermediate_size"]
        summary["weight_inspection"] = inspect_intermediate_shapes(model_dir)

    if args.output_json:
        output_path = Path(args.output_json)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
