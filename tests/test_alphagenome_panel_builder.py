from __future__ import annotations

from scripts.alphagenome_pilot.build_matched_panel_v2 import choose_panel_rows


def test_choose_panel_rows_returns_unique_quantile_anchors():
    rows = [
        {"label": f"peak_{i}", "source": "peak", "observed_total": float(i)}
        for i in range(1, 21)
    ] + [
        {"label": f"nonpeak_{i}", "source": "nonpeak", "observed_total": float(i)}
        for i in range(1, 21)
    ]

    panel = choose_panel_rows(
        rows,
        quantiles=[0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9, 0.98],
    )

    assert len(panel) == 16
    assert len({row["label"] for row in panel}) == 16
    assert {row["source"] for row in panel} == {"peak", "nonpeak"}
