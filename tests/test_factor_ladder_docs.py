from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_tracking_mentions_factor_ladder_live_item() -> None:
    tracking_text = (REPO_ROOT / "TRACKING.md").read_text(encoding="utf-8")

    assert "AlphaGenome-like factor ladder" in tracking_text
    assert "AlphaGenome-like factor ladder（E1/E2/E3）" in tracking_text


def test_registry_registers_factor_ladder_family() -> None:
    registry_text = (REPO_ROOT / "docs/experiments/registry.md").read_text(
        encoding="utf-8"
    )

    assert "`alphagenome_factor_ladder`" in registry_text
