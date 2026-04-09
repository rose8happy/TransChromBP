from __future__ import annotations

from scripts.alphagenome_pilot.build_matched_panel_v2 import (
    _import_real_data_helpers,
    choose_panel_rows,
)


def test_import_real_data_helpers_resolves_vendor_layout():
    helpers = _import_real_data_helpers()

    assert [helper.__name__ for helper in helpers] == [
        "load_fold_chroms",
        "load_regions_from_bed",
        "load_bigwig_chrom_sizes",
        "filter_records_by_chroms",
    ]


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
