from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
METRICS_DIR = ROOT / "metrics"
FIGURES_DIR = ROOT / "figures"


RUNS = {
    "orig_best_full": {
        "label": "Orig best",
        "path": "tutorial_bias2main_best_test.json",
        "family": "orig",
        "view": "full",
    },
    "orig_epoch20_full": {
        "label": "Orig epoch20",
        "path": "tutorial_bias2main_epoch020_test.json",
        "family": "orig",
        "view": "full",
    },
    "baseline_epoch20_full": {
        "label": "Baseline epoch20",
        "path": "tutorial_baseline_epoch020_test.json",
        "family": "baseline",
        "view": "full",
    },
    "hybrid_learnable_best_trimmed": {
        "label": "Hybrid learnable best",
        "path": "tutorial_hybrid_learnable_best_test.json",
        "family": "hybrid_learnable",
        "view": "trimmed",
    },
    "hybrid_learnable_epoch20_trimmed": {
        "label": "Hybrid learnable epoch20",
        "path": "tutorial_hybrid_learnable_epoch020_test.json",
        "family": "hybrid_learnable",
        "view": "trimmed",
    },
    "hybrid_fixed_best_trimmed": {
        "label": "Hybrid fixed best",
        "path": "tutorial_hybrid_fixed_best_test.json",
        "family": "hybrid_fixed",
        "view": "trimmed",
    },
    "hybrid_fixed_epoch20_trimmed": {
        "label": "Hybrid fixed epoch20",
        "path": "tutorial_hybrid_fixed_epoch020_test.json",
        "family": "hybrid_fixed",
        "view": "trimmed",
    },
    "hybrid_learnable_best_full": {
        "label": "Hybrid learnable best",
        "path": "tutorial_hybrid_learnable_best_test_full.json",
        "family": "hybrid_learnable",
        "view": "full",
    },
    "hybrid_learnable_epoch20_full": {
        "label": "Hybrid learnable epoch20",
        "path": "tutorial_hybrid_learnable_epoch020_test_full.json",
        "family": "hybrid_learnable",
        "view": "full",
    },
    "hybrid_fixed_best_full": {
        "label": "Hybrid fixed best",
        "path": "tutorial_hybrid_fixed_best_test_full.json",
        "family": "hybrid_fixed",
        "view": "full",
    },
    "hybrid_fixed_epoch20_full": {
        "label": "Hybrid fixed epoch20",
        "path": "tutorial_hybrid_fixed_epoch020_test_full.json",
        "family": "hybrid_fixed",
        "view": "full",
    },
    "orig_best_trimmed": {
        "label": "Orig best",
        "path": "tutorial_orig_best_test_npr025.json",
        "family": "orig",
        "view": "trimmed",
    },
    "orig_epoch20_trimmed": {
        "label": "Orig epoch20",
        "path": "tutorial_orig_epoch020_test_npr025.json",
        "family": "orig",
        "view": "trimmed",
    },
    "baseline_epoch20_trimmed": {
        "label": "Baseline epoch20",
        "path": "tutorial_baseline_epoch020_test_npr025.json",
        "family": "baseline",
        "view": "trimmed",
    },
}


BIAS_RUNS = {
    "orig_bias_peak": {
        "label": "Orig bias_pretrain",
        "path": "tutorial_bias_pretrain_best_test_peak.json",
        "scope": "peak",
    },
    "orig_bias_nonpeak": {
        "label": "Orig bias_pretrain",
        "path": "tutorial_bias_pretrain_best_test_nonpeak.json",
        "scope": "nonpeak",
    },
    "hybrid_bias_peak": {
        "label": "Hybrid bias_pretrain",
        "path": "tutorial_hybrid_bias_pretrain_best_test_peak_full.json",
        "scope": "peak",
    },
    "hybrid_bias_nonpeak": {
        "label": "Hybrid bias_pretrain",
        "path": "tutorial_hybrid_bias_pretrain_best_test_nonpeak_full.json",
        "scope": "nonpeak",
    },
}


SCOPES = ("overall", "peak", "nonpeak")
COUNT_KEY = "count_pearson_full"
JSD_KEY = "profile_target_jsd_full_mean"
LOSS_KEY = "loss_total"


def load_json(name: str) -> dict:
    return json.loads((METRICS_DIR / name).read_text(encoding="utf-8"))


def build_main_rows() -> list[dict]:
    rows: list[dict] = []
    for run_id, meta in RUNS.items():
        payload = load_json(meta["path"])
        for scope in SCOPES:
            if scope not in payload["results"]:
                continue
            metrics = payload["results"][scope]
            rows.append(
                {
                    "run_id": run_id,
                    "label": meta["label"],
                    "family": meta["family"],
                    "view": meta["view"],
                    "dataset_size": payload.get("dataset_size"),
                    "nonpeak_ratio": payload.get("nonpeak_ratio"),
                    "scope": scope,
                    "n_examples": metrics["n_examples"],
                    "count_r": metrics[COUNT_KEY],
                    "profile_jsd": metrics[JSD_KEY],
                    "loss_total": metrics[LOSS_KEY],
                }
            )
    return rows


def build_bias_rows() -> list[dict]:
    rows: list[dict] = []
    for run_id, meta in BIAS_RUNS.items():
        payload = load_json(meta["path"])
        metrics = payload["results"][meta["scope"]]
        rows.append(
            {
                "run_id": run_id,
                "label": meta["label"],
                "scope": meta["scope"],
                "dataset_size": payload.get("dataset_size"),
                "count_r": metrics[COUNT_KEY],
                "profile_jsd": metrics[JSD_KEY],
                "loss_total": metrics[LOSS_KEY],
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot_full_heatmaps(rows: list[dict]) -> None:
    model_order = [
        "Orig best",
        "Orig epoch20",
        "Baseline epoch20",
        "Hybrid learnable best",
        "Hybrid learnable epoch20",
        "Hybrid fixed best",
        "Hybrid fixed epoch20",
    ]
    scope_order = ["overall", "peak", "nonpeak"]
    value_maps = {
        "Count Pearson r": {},
        "Profile JSD": {},
    }
    for row in rows:
        if row["view"] != "full":
            continue
        key = (row["label"], row["scope"])
        value_maps["Count Pearson r"][key] = row["count_r"]
        value_maps["Profile JSD"][key] = row["profile_jsd"]

    matrices = []
    for title in ("Count Pearson r", "Profile JSD"):
        matrix = np.array(
            [[value_maps[title][(label, scope)] for scope in scope_order] for label in model_order],
            dtype=float,
        )
        matrices.append((title, matrix))

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    cmap_specs = {
        "Count Pearson r": ("viridis", None, None),
        "Profile JSD": ("magma_r", None, None),
    }
    for ax, (title, matrix) in zip(axes, matrices):
        cmap, vmin, vmax = cmap_specs[title]
        image = ax.imshow(matrix, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)
        ax.set_title(title)
        ax.set_xticks(np.arange(len(scope_order)))
        ax.set_xticklabels(scope_order)
        ax.set_yticks(np.arange(len(model_order)))
        ax.set_yticklabels(model_order)
        for row_idx in range(matrix.shape[0]):
            for col_idx in range(matrix.shape[1]):
                ax.text(col_idx, row_idx, f"{matrix[row_idx, col_idx]:.3f}", ha="center", va="center", color="white")
        fig.colorbar(image, ax=ax, shrink=0.85)

    fig.suptitle("Full test(chr1) scale-invariant metrics", fontsize=14)
    fig.savefig(FIGURES_DIR / "full_test_scale_invariant_heatmaps.png", dpi=180)
    plt.close(fig)


def plot_trimmed_vs_full(rows: list[dict]) -> None:
    model_order = [
        ("Orig best", "orig"),
        ("Orig epoch20", "orig"),
        ("Baseline epoch20", "baseline"),
        ("Hybrid learnable best", "hybrid_learnable"),
        ("Hybrid learnable epoch20", "hybrid_learnable"),
        ("Hybrid fixed best", "hybrid_fixed"),
        ("Hybrid fixed epoch20", "hybrid_fixed"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(14, 7), constrained_layout=True)
    metric_keys = [("count_r", "Overall count Pearson r"), ("profile_jsd", "Overall profile JSD")]
    colors = {"trimmed": "#d95f02", "full": "#1b9e77"}

    for ax, (metric_key, title) in zip(axes, metric_keys):
        for idx, (label, _) in enumerate(model_order):
            subset = [row for row in rows if row["label"] == label and row["scope"] == "overall"]
            by_view = {row["view"]: row for row in subset}
            if "trimmed" not in by_view or "full" not in by_view:
                continue
            y = len(model_order) - idx - 1
            x1 = by_view["trimmed"][metric_key]
            x2 = by_view["full"][metric_key]
            ax.plot([x1, x2], [y, y], color="#999999", linewidth=1.5, zorder=1)
            ax.scatter(x1, y, color=colors["trimmed"], s=50, zorder=2)
            ax.scatter(x2, y, color=colors["full"], s=50, zorder=2)
        ax.set_title(title)
        ax.set_yticks(np.arange(len(model_order)))
        ax.set_yticklabels([label for label, _ in reversed(model_order)])
        ax.grid(axis="x", alpha=0.25)

    axes[0].legend(["trimmed (nonpeak_ratio=0.25)", "full (nonpeak_ratio=1.0)"], loc="lower right")
    fig.suptitle("Overall metric shift from trimmed to full held-out test", fontsize=14)
    fig.savefig(FIGURES_DIR / "trimmed_vs_full_overall.png", dpi=180)
    plt.close(fig)


def plot_bias_sanity(rows: list[dict]) -> None:
    labels = ["Orig bias_pretrain", "Hybrid bias_pretrain"]
    scopes = ["peak", "nonpeak"]
    count_matrix = np.array(
        [
            [next(row["count_r"] for row in rows if row["label"] == label and row["scope"] == scope) for scope in scopes]
            for label in labels
        ],
        dtype=float,
    )
    jsd_matrix = np.array(
        [
            [next(row["profile_jsd"] for row in rows if row["label"] == label and row["scope"] == scope) for scope in scopes]
            for label in labels
        ],
        dtype=float,
    )

    x = np.arange(len(scopes))
    width = 0.35
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    colors = ["#4c72b0", "#dd8452"]

    for idx, label in enumerate(labels):
        axes[0].bar(x + (idx - 0.5) * width, count_matrix[idx], width=width, label=label, color=colors[idx])
        axes[1].bar(x + (idx - 0.5) * width, jsd_matrix[idx], width=width, label=label, color=colors[idx])

    axes[0].set_title("Bias-only count Pearson r")
    axes[1].set_title("Bias-only profile JSD")
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(scopes)
        ax.grid(axis="y", alpha=0.25)
        ax.legend()

    fig.suptitle("Bias pretrain sanity check", fontsize=14)
    fig.savefig(FIGURES_DIR / "bias_pretrain_sanity.png", dpi=180)
    plt.close(fig)


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    main_rows = build_main_rows()
    bias_rows = build_bias_rows()

    write_csv(
        ROOT / "main_metric_summary.csv",
        main_rows,
        ["run_id", "label", "family", "view", "dataset_size", "nonpeak_ratio", "scope", "n_examples", "count_r", "profile_jsd", "loss_total"],
    )
    write_csv(
        ROOT / "bias_pretrain_summary.csv",
        bias_rows,
        ["run_id", "label", "scope", "dataset_size", "count_r", "profile_jsd", "loss_total"],
    )

    plot_full_heatmaps(main_rows)
    plot_trimmed_vs_full(main_rows)
    plot_bias_sanity(bias_rows)


if __name__ == "__main__":
    main()
