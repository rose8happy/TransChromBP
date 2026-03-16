from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import torch

from transchrombp.data import ChromBPNetBigWigDataset
from transchrombp.training.train_ddp import (
    apply_external_data_config_defaults,
    apply_training_mode_defaults,
    build_model,
    load_yaml,
    resolve_data_config_path,
    validate_model_config,
)


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
FIG_DIR = ROOT / "figures"
REPO_ROOT = Path("/data1/zhoujiazhen/bylw_atac/TransChromBP")
DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

MODEL_SPECS = {
    "bias2main_best": {
        "checkpoint": Path("/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/checkpoints/tutorial_bias2main_20260315_023539_main_learnable/best.pt"),
        "color": "#bf5b17",
        "label": "bias→main best full",
        "debiased_label": "bias→main best debiased",
        "linestyle": "-",
    },
    "bias2main_epoch20": {
        "checkpoint": Path("/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/checkpoints/tutorial_bias2main_20260315_023539_main_learnable/epoch_020.pt"),
        "color": "#1b9e77",
        "label": "bias→main epoch20 full",
        "debiased_label": "bias→main epoch20 debiased",
        "linestyle": "-",
    },
    "baseline_epoch20": {
        "checkpoint": Path("/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/checkpoints/tutorial_real_baseline_20260313_1740/epoch_020.pt"),
        "color": "#4c78a8",
        "label": "baseline epoch20 full",
        "debiased_label": "baseline epoch20 debiased",
        "linestyle": "-",
    },
    "bias_pretrain_best": {
        "checkpoint": Path("/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/checkpoints/tutorial_bias2main_20260315_023539_bias_pretrain/best.pt"),
        "color": "#7f7f7f",
        "label": "bias-pretrain only",
        "debiased_label": "bias-pretrain only",
        "linestyle": "--",
    },
}

SELECTORS = [
    ("peak_high", "peak", 0.90),
    ("peak_mid", "peak", 0.50),
    ("nonpeak_high", "nonpeak", 0.90),
    ("nonpeak_mid", "nonpeak", 0.50),
]


def resolve_path(path_value: str, base_dir: Path) -> str:
    path = Path(path_value)
    return str(path if path.is_absolute() else (base_dir / path).resolve())


def load_checkpoint_model(checkpoint_path: Path) -> tuple[torch.nn.Module, dict[str, Any], dict[str, Any], dict[str, Any]]:
    payload = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    model_cfg = payload["model_config"]
    train_cfg = payload["train_config"]

    data_config_path = resolve_path(resolve_data_config_path(train_cfg, ""), REPO_ROOT)
    data_source_cfg = load_yaml(data_config_path)
    apply_external_data_config_defaults(train_cfg, data_source_cfg)
    apply_training_mode_defaults(train_cfg)

    bias_cfg = model_cfg.setdefault("bias_branch", {})
    if bool(bias_cfg.get("freeze_bias_core", False)) and not str(bias_cfg.get("pretrained_path", "")):
        bias_cfg["freeze_bias_core"] = False

    validate_model_config(model_cfg)
    model, _ = build_model(model_cfg, train_cfg, rank=0)
    model.load_state_dict(payload["model_state"], strict=True)
    model.to(DEVICE)
    model.eval()
    return model, payload, train_cfg, data_source_cfg


def build_test_dataset(train_cfg: dict[str, Any], data_source_cfg: dict[str, Any]) -> ChromBPNetBigWigDataset:
    data_cfg = train_cfg["data"]
    return ChromBPNetBigWigDataset(
        genome_fasta=str(data_source_cfg["genome_fasta"]),
        bigwig_path=str(data_source_cfg["input"]["bigwig"]),
        peaks_bed=str(data_source_cfg["input"]["peaks_bed"]),
        nonpeaks_bed=str(data_source_cfg["input"]["nonpeaks_bed"]),
        folds_json=str(data_source_cfg["folds_json"]),
        split="test",
        input_len=int(data_cfg["input_len"]),
        supervised_bp=int(data_cfg.get("supervised_bp", data_cfg["output_len"])),
        profile_bin_size=int(data_cfg.get("profile_bin_size", 1)),
        max_jitter=0,
        peak_max_jitter=0,
        nonpeak_max_jitter=0,
        seed=int(train_cfg.get("seed", 1234)) + 20000,
        nonpeak_ratio=float(data_cfg.get("nonpeak_ratio", 1.0)),
        max_records=0,
        region_source="both",
        random_revcomp=False,
        revcomp_prob=float(data_cfg.get("revcomp_prob", 0.5)),
        track_total_count_target=float(data_cfg.get("track_total_count_target", 0.0)),
    )


def profile_from_outputs(outputs: Any, use_debiased: bool) -> np.ndarray:
    if use_debiased:
        logits = outputs.profile_logits_debiased.detach().float()
        logcount = outputs.logcount_debiased.detach().float().reshape(-1, 1)
    else:
        logits = outputs.profile_logits_full.detach().float()
        logcount = outputs.logcount_full.detach().float().reshape(-1, 1)
    probs = torch.softmax(logits, dim=-1)
    totals = torch.clamp(torch.expm1(logcount), min=0.0)
    return (probs * totals).cpu().numpy()[0]


def compute_region_totals(dataset: ChromBPNetBigWigDataset) -> np.ndarray:
    totals = np.zeros(len(dataset.records), dtype=np.float64)
    profile_half_width = dataset.supervised_bp // 2
    for idx, record in enumerate(dataset.records):
        profile_start = int(record.center) - profile_half_width
        totals[idx] = float(dataset._fetch_profile(record.chrom, profile_start, dataset.supervised_bp).sum())
    return totals


def select_indices(dataset: ChromBPNetBigWigDataset, totals: np.ndarray) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    used: set[int] = set()
    source_to_indices: dict[str, np.ndarray] = {}
    for source in ("peak", "nonpeak"):
        source_to_indices[source] = np.array([idx for idx, record in enumerate(dataset.records) if record.source == source], dtype=np.int64)

    for label, source, quantile in SELECTORS:
        source_indices = source_to_indices[source]
        source_totals = totals[source_indices]
        target = float(np.quantile(source_totals, quantile))
        order = np.argsort(np.abs(source_totals - target))
        chosen = None
        for pos in order:
            candidate = int(source_indices[int(pos)])
            if candidate not in used:
                chosen = candidate
                used.add(candidate)
                break
        if chosen is None:
            chosen = int(source_indices[int(order[0])])
        record = dataset.records[chosen]
        selected.append(
            {
                "label": label,
                "source": source,
                "quantile": quantile,
                "index": chosen,
                "chrom": record.chrom,
                "center": int(record.center),
                "observed_total_count": float(totals[chosen]),
            }
        )
    return selected


def fetch_sample(dataset: ChromBPNetBigWigDataset, idx: int) -> tuple[np.ndarray, np.ndarray]:
    sample = dataset[idx]
    seq = sample["seq"].unsqueeze(0).to(DEVICE)
    profile_counts = sample["profile_counts"].numpy()
    return seq, profile_counts


def make_predictions(seq: torch.Tensor, models: dict[str, torch.nn.Module]) -> dict[str, dict[str, np.ndarray]]:
    outputs: dict[str, dict[str, np.ndarray]] = {}
    with torch.no_grad():
        for key, model in models.items():
            model_out = model(seq)
            outputs[key] = {
                "full": profile_from_outputs(model_out, use_debiased=False),
                "debiased": profile_from_outputs(model_out, use_debiased=True),
            }
    return outputs


def plot_selected_loci(dataset: ChromBPNetBigWigDataset, selected: list[dict[str, Any]], models: dict[str, torch.nn.Module]) -> None:
    fig, axes = plt.subplots(len(selected), 2, figsize=(15, 3.4 * len(selected)), constrained_layout=True)
    x = np.arange(dataset.supervised_bp)

    for row_idx, info in enumerate(selected):
        seq, truth = fetch_sample(dataset, info["index"])
        preds = make_predictions(seq, models)
        center_pos = dataset.supervised_bp // 2

        ax_full = axes[row_idx, 0]
        ax_deb = axes[row_idx, 1]

        ax_full.plot(x, truth, color="black", linewidth=1.7, label="observed")
        ax_deb.plot(x, truth, color="black", linewidth=1.7, label="observed")

        for model_key, spec in MODEL_SPECS.items():
            ax_full.plot(
                x,
                preds[model_key]["full"],
                color=spec["color"],
                linestyle=spec["linestyle"],
                linewidth=1.5,
                label=spec["label"],
                alpha=0.95 if model_key != "bias_pretrain_best" else 0.9,
            )
            if model_key != "bias_pretrain_best":
                ax_deb.plot(
                    x,
                    preds[model_key]["debiased"],
                    color=spec["color"],
                    linestyle=spec["linestyle"],
                    linewidth=1.5,
                    label=spec["debiased_label"],
                )

        title = (
            f"{info['label']} | {info['source']} | {info['chrom']}:{info['center'] - dataset.supervised_bp//2}-"
            f"{info['center'] + dataset.supervised_bp//2} | observed total={info['observed_total_count']:.1f}"
        )
        ax_full.set_title(title, fontsize=10.5)
        ax_deb.set_title(f"{info['label']} debiased comparison", fontsize=10.5)
        for ax in (ax_full, ax_deb):
            ax.axvline(center_pos, color="#999999", linestyle=":", linewidth=1.0, alpha=0.7)
            ax.grid(alpha=0.18)
            ax.set_xlim(0, dataset.supervised_bp - 1)
            ax.set_ylabel("Predicted counts")

        if row_idx == len(selected) - 1:
            ax_full.set_xlabel("Position within 1000 bp output window")
            ax_deb.set_xlabel("Position within 1000 bp output window")

    axes[0, 0].legend(frameon=False, fontsize=9, ncol=3, loc="upper right")
    axes[0, 1].legend(frameon=False, fontsize=9, ncol=2, loc="upper right")
    out = FIG_DIR / "locus_track_comparison.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_selected_csv(selected: list[dict[str, Any]]) -> None:
    path = ROOT / "selected_loci.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["label", "source", "quantile", "index", "chrom", "center", "observed_total_count"],
        )
        writer.writeheader()
        writer.writerows(selected)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    models: dict[str, torch.nn.Module] = {}

    primary_model, _, primary_train_cfg, primary_data_cfg = load_checkpoint_model(MODEL_SPECS["bias2main_best"]["checkpoint"])
    models["bias2main_best"] = primary_model
    dataset = build_test_dataset(primary_train_cfg, primary_data_cfg)

    for key in ("bias2main_epoch20", "baseline_epoch20", "bias_pretrain_best"):
        model, _, _, _ = load_checkpoint_model(MODEL_SPECS[key]["checkpoint"])
        models[key] = model

    totals = compute_region_totals(dataset)
    selected = select_indices(dataset, totals)
    write_selected_csv(selected)
    plot_selected_loci(dataset, selected, models)

    summary = {
        "device": str(DEVICE),
        "n_test_regions": int(len(dataset.records)),
        "selected": selected,
    }
    (ROOT / "locus_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
