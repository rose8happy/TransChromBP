"""Tests for model teacher cache export."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import types

import numpy as np
import torch

if "pyBigWig" not in sys.modules:
    sys.modules["pyBigWig"] = types.SimpleNamespace(open=lambda *args, **kwargs: None)

from transchrombp.data.real_data import ChromBPNetBigWigDataset, RegionRecord, compute_record_sha1
from transchrombp.evaluation.model_teacher_cache_export import (
    export_arrays,
    export_model_teacher_cache,
    parse_args,
    pool_profile_probs_to_bins,
    write_teacher_manifest,
)
from transchrombp.training.train_ddp import _pool_profile_probs_to_bins as training_pool_profile_probs_to_bins


def test_pool_profile_probs_to_bins_matches_training_distill_pooling() -> None:
    profile_probs = np.full((1, 1000), 1.0 / 1000.0, dtype=np.float32)

    pooled = pool_profile_probs_to_bins(profile_probs)
    expected = training_pool_profile_probs_to_bins(torch.from_numpy(profile_probs), 16).numpy()
    assert pooled.shape == (1, 16)
    np.testing.assert_allclose(pooled.sum(axis=1), np.ones(1, dtype=np.float32), atol=1e-6)
    np.testing.assert_allclose(pooled, expected, atol=1e-6)
    np.testing.assert_allclose(pooled, np.full((1, 16), 1.0 / 16.0, dtype=np.float32), atol=1e-6)


def test_export_arrays_and_manifest_are_accepted_by_teacher_cache_loader(tmp_path: Path) -> None:
    profile16 = np.arange(32, dtype=np.float64).reshape(2, 16)
    profile16 = profile16 / profile16.sum(axis=1, keepdims=True)
    logcount = np.asarray([1.5, 2.5], dtype=np.float64)
    records = [
        RegionRecord(chrom="chr1", center=100, source="peak"),
        RegionRecord(chrom="chr1", center=200, source="nonpeak"),
    ]

    export_arrays(output_dir=tmp_path, split="valid", profile16=profile16, logcount=logcount)
    write_teacher_manifest(
        output_dir=tmp_path,
        split="valid",
        n_records=2,
        record_sha1=compute_record_sha1(records),
        model_checkpoint="/tmp/fake.ckpt",
    )

    dataset = ChromBPNetBigWigDataset.__new__(ChromBPNetBigWigDataset)
    dataset.split = "valid"
    dataset.records = records
    dataset._teacher_targets = {}
    dataset._bigwig = None

    dataset._load_teacher_cache(str(tmp_path), ["profile16", "logcount"])

    assert set(dataset._teacher_targets) == {"profile16", "logcount"}
    assert dataset._teacher_targets["profile16"].dtype == np.float32
    assert dataset._teacher_targets["logcount"].dtype == np.float32
    assert dataset._teacher_targets["profile16"].shape == (2, 16)
    assert dataset._teacher_targets["logcount"].shape == (2,)


def test_parse_args_supports_export_invocation() -> None:
    args = parse_args(
        [
            "--checkpoint",
            "teacher.pt",
            "--output-dir",
            "/tmp/out",
            "--splits",
            "train",
            "test",
            "--data-config",
            "data.yaml",
            "--batch-size",
            "4",
            "--num-workers",
            "1",
            "--max-regions",
            "9",
            "--device",
            "cpu",
        ]
    )

    assert args.checkpoint == "teacher.pt"
    assert args.output_dir == "/tmp/out"
    assert args.splits == ["train", "test"]
    assert args.data_config == "data.yaml"
    assert args.batch_size == 4
    assert args.num_workers == 1
    assert args.max_regions == 9
    assert args.device == "cpu"


def test_export_model_teacher_cache_writes_split_outputs(tmp_path: Path, monkeypatch) -> None:
    import transchrombp.evaluation.model_teacher_cache_export as module

    calls: list[tuple[str, object]] = []

    class FakeDataset(torch.utils.data.Dataset):
        def __init__(self, split: str) -> None:
            self.split = split
            self.records = [
                RegionRecord(chrom="chr1", center=100 + idx * 50, source="peak")
                for idx in range(2 if split == "train" else 3)
            ]

        def __len__(self) -> int:
            return len(self.records)

        def __getitem__(self, idx: int) -> dict[str, torch.Tensor]:
            seq = torch.zeros(32, 4, dtype=torch.float32)
            seq[:, 0] = float(idx + 1)
            profile_counts = torch.full((1000,), float(idx + 1), dtype=torch.float32)
            return {
                "seq": seq,
                "profile_counts": profile_counts,
            }

    class FakeModel(torch.nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.forward_calls = 0

        def forward(self, seq: torch.Tensor, **_: torch.Tensor):
            self.forward_calls += 1
            batch = seq.shape[0]
            sample_ids = seq[:, 0, 0].to(dtype=torch.float32)
            logits = torch.full((batch, 1000), -8.0, device=seq.device, dtype=torch.float32)
            logits[:, :500] = 0.0
            logits[:, 500:] = sample_ids.unsqueeze(1)
            logcount = sample_ids.unsqueeze(1) + 0.25
            return type(
                "FakeOutput",
                (),
                {
                    "profile_logits_debiased": logits,
                    "logcount_debiased": logcount,
                },
            )()

    fake_model = FakeModel()
    checkpoint_path = tmp_path / "teacher.pt"
    checkpoint_path.write_bytes(b"fake")
    data_cfg_path = tmp_path / "data.yaml"
    data_cfg_path.write_text("input: {}\n", encoding="utf-8")

    payload = {
        "model_config": {"heads": {"profile_output_len": 1000}},
        "train_config": {
            "seed": 7,
            "data": {
                "batch_size_per_gpu": 2,
                "num_workers": 0,
                "teacher_cache_dir": "/tmp/preexisting_teacher_cache",
                "teacher_target_names": ["profile16", "logcount"],
            },
        },
        "model_state": {"weight": torch.tensor([1.0])},
    }

    def _record(name: str):
        def _inner(*args, **kwargs):
            calls.append((name, {"args": args, "kwargs": kwargs}))
        return _inner

    def fake_load_yaml(path: str) -> dict:
        calls.append(("load_yaml", path))
        return {
            "genome_fasta": "/tmp/fake.fa",
            "folds_json": "/tmp/folds.json",
            "input": {
                "bigwig": "/tmp/fake.bw",
                "peaks_bed": "/tmp/peaks.bed",
                "nonpeaks_bed": "/tmp/nonpeaks.bed",
            },
            "window": {},
            "sampling": {},
        }

    def fake_resolve_data_config_path(train_cfg: dict, cli_data_config: str) -> str:
        calls.append(("resolve_data_config_path", cli_data_config))
        return cli_data_config

    def fake_dataset_factory(**kwargs):
        calls.append(("dataset_kwargs", dict(kwargs)))
        calls.append(("dataset", kwargs["split"]))
        return FakeDataset(kwargs["split"])

    monkeypatch.setattr(module.torch, "load", lambda *args, **kwargs: payload)
    monkeypatch.setattr(module, "load_yaml", fake_load_yaml)
    monkeypatch.setattr(module, "resolve_data_config_path", fake_resolve_data_config_path)
    monkeypatch.setattr(module, "apply_external_data_config_defaults", _record("apply_external_data_config_defaults"))
    monkeypatch.setattr(module, "apply_training_mode_defaults", _record("apply_training_mode_defaults"))
    monkeypatch.setattr(module, "validate_semantics_profile", _record("validate_semantics_profile"))
    monkeypatch.setattr(module, "check_profile_resolution_compatibility", _record("check_profile_resolution_compatibility"))
    monkeypatch.setattr(module, "validate_model_config", _record("validate_model_config"))
    monkeypatch.setattr(module, "build_model", lambda model_cfg, train_cfg, rank: (fake_model, "full"))
    monkeypatch.setattr(module, "load_model_state_with_count_pool_compat", _record("load_model_state_with_count_pool_compat"))
    monkeypatch.setattr(module, "build_online_feature_runtime", lambda model_cfg, device, rank: None)
    monkeypatch.setattr(module, "extract_online_feature_kwargs", lambda seq, runtime: {})
    monkeypatch.setattr(module, "extract_foundation_feature_kwargs", lambda batch, model_cfg, device: {})
    monkeypatch.setattr(module, "ChromBPNetBigWigDataset", fake_dataset_factory)

    export_model_teacher_cache(
        checkpoint=checkpoint_path,
        output_dir=tmp_path,
        data_config=data_cfg_path,
        splits=("train", "valid"),
        batch_size=2,
        num_workers=0,
        max_regions=0,
        device="cpu",
    )

    train_profile = np.load(tmp_path / "train_profile16.f32.npy")
    train_logcount = np.load(tmp_path / "train_logcount.f32.npy")
    valid_manifest = json.loads((tmp_path / "teacher_manifest_valid.json").read_text(encoding="utf-8"))
    expected_valid_sha1 = compute_record_sha1(FakeDataset("valid").records)

    assert train_profile.shape == (2, 16)
    assert train_logcount.shape == (2,)
    np.testing.assert_allclose(train_profile.sum(axis=1), np.ones(2, dtype=np.float32), atol=1e-6)
    assert valid_manifest["n_records"] == 3
    assert valid_manifest["record_sha1"] == expected_valid_sha1
    assert valid_manifest["targets"] == ["profile16", "logcount"]
    assert valid_manifest["model_checkpoint"] == str(checkpoint_path.resolve())
    assert fake_model.forward_calls == 3
    assert ("resolve_data_config_path", str(data_cfg_path)) in calls
    assert ("dataset", "train") in calls
    assert ("dataset", "valid") in calls
    dataset_kwargs = [payload for name, payload in calls if name == "dataset_kwargs"]
    assert dataset_kwargs
    for kwargs in dataset_kwargs:
        assert kwargs["teacher_cache_dir"] == ""
        assert kwargs["teacher_target_names"] == []


def test_export_model_teacher_cache_skips_genos_summary_for_non_cached_models(
    tmp_path: Path,
    monkeypatch,
) -> None:
    import transchrombp.evaluation.model_teacher_cache_export as module

    class FakeDataset(torch.utils.data.Dataset):
        def __init__(self, split: str) -> None:
            self.split = split
            self.records = [
                RegionRecord(chrom="chr1", center=100 + idx * 50, source="peak")
                for idx in range(2)
            ]

        def __len__(self) -> int:
            return len(self.records)

        def __getitem__(self, idx: int) -> dict[str, torch.Tensor]:
            seq = torch.zeros(32, 4, dtype=torch.float32)
            seq[:, 0] = float(idx + 1)
            profile_counts = torch.full((1000,), float(idx + 1), dtype=torch.float32)
            return {
                "seq": seq,
                "profile_counts": profile_counts,
                "genos_global_mean": torch.full((4,), float(idx + 1), dtype=torch.float32),
            }

    class StrictNoGenosModel(torch.nn.Module):
        def forward(self, seq: torch.Tensor, **kwargs: torch.Tensor):
            assert "genos_summary" not in kwargs
            batch = seq.shape[0]
            logits = torch.zeros((batch, 1000), device=seq.device, dtype=torch.float32)
            logcount = torch.zeros((batch, 1), device=seq.device, dtype=torch.float32)
            return type(
                "FakeOutput",
                (),
                {
                    "profile_logits_debiased": logits,
                    "logcount_debiased": logcount,
                },
            )()

    checkpoint_path = tmp_path / "teacher.pt"
    checkpoint_path.write_bytes(b"fake")
    data_cfg_path = tmp_path / "data.yaml"
    data_cfg_path.write_text("input: {}\n", encoding="utf-8")

    payload = {
        "model_config": {
            "heads": {"profile_output_len": 1000},
            "genos_cached": {"enabled": False},
        },
        "train_config": {
            "seed": 7,
            "data": {
                "batch_size_per_gpu": 2,
                "num_workers": 0,
            },
        },
        "model_state": {"weight": torch.tensor([1.0])},
    }

    def fake_load_yaml(path: str) -> dict:
        return {
            "genome_fasta": "/tmp/fake.fa",
            "folds_json": "/tmp/folds.json",
            "input": {
                "bigwig": "/tmp/fake.bw",
                "peaks_bed": "/tmp/peaks.bed",
                "nonpeaks_bed": "/tmp/nonpeaks.bed",
            },
            "window": {},
            "sampling": {},
        }

    monkeypatch.setattr(module.torch, "load", lambda *args, **kwargs: payload)
    monkeypatch.setattr(module, "load_yaml", fake_load_yaml)
    monkeypatch.setattr(module, "resolve_data_config_path", lambda train_cfg, cli_data_config: cli_data_config)
    monkeypatch.setattr(module, "apply_external_data_config_defaults", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "apply_training_mode_defaults", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "validate_semantics_profile", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "check_profile_resolution_compatibility", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "validate_model_config", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "build_model", lambda model_cfg, train_cfg, rank: (StrictNoGenosModel(), "full"))
    monkeypatch.setattr(module, "load_model_state_with_count_pool_compat", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "build_online_feature_runtime", lambda model_cfg, device, rank: None)
    monkeypatch.setattr(module, "extract_online_feature_kwargs", lambda seq, runtime: {})
    monkeypatch.setattr(module, "extract_foundation_feature_kwargs", lambda batch, model_cfg, device: {})
    monkeypatch.setattr(module, "ChromBPNetBigWigDataset", lambda **kwargs: FakeDataset(kwargs["split"]))

    export_model_teacher_cache(
        checkpoint=checkpoint_path,
        output_dir=tmp_path,
        data_config=data_cfg_path,
        splits=("train",),
        batch_size=2,
        num_workers=0,
        max_regions=0,
        device="cpu",
    )

    train_profile = np.load(tmp_path / "train_profile16.f32.npy")
    assert train_profile.shape == (2, 16)
