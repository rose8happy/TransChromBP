from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
FIG_DIR = ROOT / "figures"

COLORS = {
    "bias2main_best": "#bf5b17",
    "bias2main_epoch020": "#1b9e77",
    "baseline_epoch020": "#4c78a8",
}


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_jsonl(path: Path) -> list[dict]:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line:
            rows.append(json.loads(line))
    return rows


def parse_baseline_log(path: Path) -> list[dict[str, float]]:
    epoch_pat = re.compile(
        r"^\[epoch\]\s+(?P<epoch>\d+)\s+train_loss=(?P<train_loss>[0-9.]+)\s+profile=(?P<profile>[0-9.]+)\s+count=(?P<count>[0-9.]+)\s+elapsed=(?P<elapsed>[0-9.]+)s$"
    )
    val_pat = re.compile(
        r"^\[val\]\s+loss_total=(?P<loss_total>[0-9.]+)\s+loss_profile=(?P<loss_profile>[0-9.]+)\s+loss_count=(?P<loss_count>[0-9.]+)\s+loss_debiased_profile=(?P<loss_debiased_profile>[0-9.]+)\s+loss_debiased_count=(?P<loss_debiased_count>[0-9.]+)$"
    )

    rows: list[dict[str, float]] = []
    current: dict[str, float] | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        m_epoch = epoch_pat.match(line.strip())
        if m_epoch:
            current = {k: float(v) for k, v in m_epoch.groupdict().items()}
            current["epoch"] = int(current["epoch"])
            rows.append(current)
            continue

        m_val = val_pat.match(line.strip())
        if m_val and current is not None:
            for k, v in m_val.groupdict().items():
                current[f"val_{k}"] = float(v)
    return rows


def build_test_metric_figure() -> None:
    checkpoints = [
        ("bias2main_best", "bias→main best"),
        ("bias2main_epoch020", "bias→main epoch20"),
        ("baseline_epoch020", "baseline epoch20"),
    ]
    metrics = {
        "loss_total": ("Loss Total", False),
        "count_pearson_full": ("Count Pearson r", True),
        "profile_target_jsd_full_mean": ("Profile JSD", False),
    }
    scopes = [("overall", "Overall"), ("peak", "Peak"), ("nonpeak", "Nonpeak")]

    payloads = {
        key: load_json(DATA_DIR / f"tutorial_{key}_test.json" if key != "baseline_epoch020" else DATA_DIR / "tutorial_baseline_epoch020_test.json")
        for key, _ in checkpoints
    }
    # Normalize naming above for the two bias2main files.
    payloads["bias2main_best"] = load_json(DATA_DIR / "tutorial_bias2main_best_test.json")
    payloads["bias2main_epoch020"] = load_json(DATA_DIR / "tutorial_bias2main_epoch020_test.json")
    payloads["baseline_epoch020"] = load_json(DATA_DIR / "tutorial_baseline_epoch020_test.json")

    fig, axes = plt.subplots(3, 1, figsize=(8.8, 10), constrained_layout=True)
    width = 0.22
    x = np.arange(len(scopes))

    for row_idx, (metric_key, (title, higher_better)) in enumerate(metrics.items()):
        ax = axes[row_idx]
        for idx, (checkpoint_key, checkpoint_label) in enumerate(checkpoints):
            values = [payloads[checkpoint_key]["results"][scope_key][metric_key] for scope_key, _ in scopes]
            ax.bar(
                x + (idx - 1) * width,
                values,
                width=width,
                color=COLORS[checkpoint_key],
                label=checkpoint_label,
            )

        ax.set_xticks(x)
        ax.set_xticklabels([label for _, label in scopes])
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.25, linewidth=0.8)
        if not higher_better:
            best_values = [
                min(payloads[checkpoint_key]["results"][scope_key][metric_key] for checkpoint_key, _ in checkpoints)
                for scope_key, _ in scopes
            ]
        else:
            best_values = [
                max(payloads[checkpoint_key]["results"][scope_key][metric_key] for checkpoint_key, _ in checkpoints)
                for scope_key, _ in scopes
            ]
        for pos, best_value in enumerate(best_values):
            ax.axhline(best_value, xmin=(pos - 0.35) / 3, xmax=(pos + 0.35) / 3, color="#777777", alpha=0.15)

    axes[0].legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.22))
    fig.suptitle("Held-out test (chr1) comparison across checkpoints", fontsize=15, y=1.02)
    out = FIG_DIR / "test_metric_comparison.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)


def build_validation_curve_figure() -> None:
    baseline_rows = parse_baseline_log(DATA_DIR / "transchrombp_tutorial_real_baseline_20260313_1740.log")
    bias_rows = load_jsonl(DATA_DIR / "epoch_metrics.jsonl")

    baseline_epochs = [row["epoch"] for row in baseline_rows]
    baseline_val = [row["val_loss_total"] for row in baseline_rows]

    bias_epochs = [row["epoch"] for row in bias_rows]
    bias_overall = [row["val"]["overall"]["loss_total"] for row in bias_rows]
    bias_peak = [row["val"]["peak"]["loss_total"] for row in bias_rows]
    bias_nonpeak = [row["val"]["nonpeak"]["loss_total"] for row in bias_rows]

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), constrained_layout=True)

    ax = axes[0]
    ax.plot(baseline_epochs, baseline_val, marker="o", linewidth=2.2, color=COLORS["baseline_epoch020"], label="baseline overall val")
    ax.plot(bias_epochs, bias_overall, marker="o", linewidth=2.2, color=COLORS["bias2main_epoch020"], label="bias→main overall val")
    ax.axvline(10, color=COLORS["bias2main_best"], linestyle="--", linewidth=1.6, alpha=0.9, label="bias→main best epoch")
    ax.set_title("Overall validation loss by epoch")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Loss Total")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)

    ax = axes[1]
    ax.plot(bias_epochs, bias_overall, marker="o", linewidth=2.0, color="#444444", label="overall")
    ax.plot(bias_epochs, bias_peak, marker="o", linewidth=2.0, color=COLORS["bias2main_best"], label="peak")
    ax.plot(bias_epochs, bias_nonpeak, marker="o", linewidth=2.0, color=COLORS["bias2main_epoch020"], label="nonpeak")
    ax.axvline(10, color=COLORS["bias2main_best"], linestyle="--", linewidth=1.4, alpha=0.8)
    ax.set_title("Bias→main split validation loss by epoch")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Loss Total")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)

    out = FIG_DIR / "validation_curves.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)


def build_bias_diagnostic_figure() -> None:
    rows = load_jsonl(DATA_DIR / "epoch_metrics.jsonl")
    epochs = [row["epoch"] for row in rows]

    metrics = [
        ("effective_profile_scale", "Val overall effective_profile_scale"),
        ("effective_count_scale", "Val overall effective_count_scale"),
        ("profile_full_debiased_jsd", "Val overall full-vs-debiased JSD"),
        ("count_full_debiased_abs", "Val overall |count_full-debiased|"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.2), constrained_layout=True)
    axes = axes.ravel()

    for ax, (metric_key, title) in zip(axes, metrics):
        values = [row["val"]["overall"][metric_key] for row in rows]
        ax.plot(epochs, values, marker="o", linewidth=2.1, color="#5b6cb8")
        ax.axvline(10, color=COLORS["bias2main_best"], linestyle="--", linewidth=1.4, alpha=0.85)
        ax.set_title(title)
        ax.set_xlabel("Epoch")
        ax.grid(alpha=0.25)

    out = FIG_DIR / "bias_diagnostics.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)


def build_summary_json() -> None:
    bias_rows = load_jsonl(DATA_DIR / "epoch_metrics.jsonl")
    baseline_rows = parse_baseline_log(DATA_DIR / "transchrombp_tutorial_real_baseline_20260313_1740.log")
    payload_best = load_json(DATA_DIR / "tutorial_bias2main_best_test.json")
    payload_epoch20 = load_json(DATA_DIR / "tutorial_bias2main_epoch020_test.json")
    payload_base = load_json(DATA_DIR / "tutorial_baseline_epoch020_test.json")

    summary = {
        "bias2main_best_epoch": min(bias_rows, key=lambda row: row["val"]["peak"]["loss_total"])["epoch"],
        "bias2main_final_epoch": bias_rows[-1]["epoch"],
        "baseline_final_epoch": baseline_rows[-1]["epoch"],
        "test_dataset_size": payload_best["dataset_size"],
        "deltas_vs_baseline": {
            "best_peak_loss_total": payload_best["results"]["peak"]["loss_total"] - payload_base["results"]["peak"]["loss_total"],
            "epoch20_overall_loss_total": payload_epoch20["results"]["overall"]["loss_total"] - payload_base["results"]["overall"]["loss_total"],
            "epoch20_nonpeak_loss_total": payload_epoch20["results"]["nonpeak"]["loss_total"] - payload_base["results"]["nonpeak"]["loss_total"],
            "best_peak_count_r": payload_best["results"]["peak"]["count_pearson_full"] - payload_base["results"]["peak"]["count_pearson_full"],
            "epoch20_nonpeak_count_r": payload_epoch20["results"]["nonpeak"]["count_pearson_full"] - payload_base["results"]["nonpeak"]["count_pearson_full"],
        },
    }
    (ROOT / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    build_test_metric_figure()
    build_validation_curve_figure()
    build_bias_diagnostic_figure()
    build_summary_json()


if __name__ == "__main__":
    main()
