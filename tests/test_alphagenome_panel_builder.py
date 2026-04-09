from __future__ import annotations

import builtins
import json
import sys
import types

import pyBigWig
import pytest

from scripts.alphagenome_pilot.build_matched_panel_v2 import (
    _import_real_data_helpers,
    choose_panel_rows,
)


def test_import_real_data_helpers_falls_back_without_transchrombp_module(tmp_path, monkeypatch):
    original_import = builtins.__import__

    def controlled_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "transchrombp.data.real_data":
            raise ModuleNotFoundError("No module named 'transchrombp'", name="transchrombp")
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", controlled_import)

    load_fold_chroms, load_regions_from_bed, load_bigwig_chrom_sizes, filter_records_by_chroms = (
        _import_real_data_helpers()
    )

    folds_path = tmp_path / "folds.json"
    folds_path.write_text(json.dumps({"test": ["chr1"]}), encoding="utf-8")
    assert load_fold_chroms(str(folds_path), "test") == ["chr1"]

    bed_path = tmp_path / "regions.bed"
    bed_path.write_text("chr1\t10\t30\t.\t.\t.\t.\t.\t.\t7\nchr2\t50\t80\n", encoding="utf-8")
    records = load_regions_from_bed(str(bed_path), ["chr1", "chr2"], source="peak")
    assert [(record.chrom, record.center, record.source) for record in records] == [
        ("chr1", 17, "peak"),
        ("chr2", 65, "peak"),
    ]

    with pytest.warns(RuntimeWarning, match="Filtered 1 peak records"):
        filtered = filter_records_by_chroms(records, ["chr1"], str(bed_path), source="peak")
    assert [(record.chrom, record.center, record.source) for record in filtered] == [
        ("chr1", 17, "peak")
    ]

    bigwig_path = tmp_path / "toy.bw"
    writer = pyBigWig.open(str(bigwig_path), "w")
    try:
        writer.addHeader([("chr1", 100)])
        writer.addEntries(["chr1"], [0], ends=[10], values=[1.0])
    finally:
        writer.close()

    assert load_bigwig_chrom_sizes(str(bigwig_path)) == {"chr1": 100}


def test_import_real_data_helpers_prefers_real_module_when_available(monkeypatch):
    package = types.ModuleType("transchrombp")
    data_package = types.ModuleType("transchrombp.data")
    real_module = types.ModuleType("transchrombp.data.real_data")

    def real_load_fold_chroms(*args, **kwargs):
        return ["real-folds"]

    def real_load_regions_from_bed(*args, **kwargs):
        return ["real-regions"]

    def real_load_bigwig_chrom_sizes(*args, **kwargs):
        return {"chrReal": 123}

    def real_filter_records_by_chroms(*args, **kwargs):
        return ["real-filtered"]

    real_module.load_fold_chroms = real_load_fold_chroms
    real_module.load_regions_from_bed = real_load_regions_from_bed
    real_module.load_bigwig_chrom_sizes = real_load_bigwig_chrom_sizes
    real_module.filter_records_by_chroms = real_filter_records_by_chroms
    data_package.real_data = real_module
    package.data = data_package

    monkeypatch.setitem(sys.modules, "transchrombp", package)
    monkeypatch.setitem(sys.modules, "transchrombp.data", data_package)
    monkeypatch.setitem(sys.modules, "transchrombp.data.real_data", real_module)

    helpers = _import_real_data_helpers()

    assert helpers == (
        real_load_fold_chroms,
        real_load_regions_from_bed,
        real_load_bigwig_chrom_sizes,
        real_filter_records_by_chroms,
    )


def test_import_real_data_helpers_reraises_unrelated_missing_dependency(monkeypatch):
    original_import = builtins.__import__

    def controlled_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "transchrombp.data.real_data":
            raise ModuleNotFoundError("No module named 'numpy'", name="numpy")
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", controlled_import)

    with pytest.raises(ModuleNotFoundError, match="No module named 'numpy'"):
        _import_real_data_helpers()


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
