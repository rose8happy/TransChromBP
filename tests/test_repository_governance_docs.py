from __future__ import annotations

from pathlib import Path
import subprocess


REPO_ROOT = Path(__file__).resolve().parents[1]


def read_text(path: str) -> str:
    return (REPO_ROOT / path).read_text(encoding="utf-8")


def test_governance_doc_declares_single_source_of_truth() -> None:
    governance_path = REPO_ROOT / "docs/env/repository_governance.md"
    genos_env_path = REPO_ROOT / "docs/env/transchrombp_genos_env.md"
    assert governance_path.exists()
    assert genos_env_path.exists()

    text = governance_path.read_text(encoding="utf-8")
    assert "single source of truth" in text.lower()
    assert "/home/zhengwei/project/python/chromBPNet" in text
    assert "/data1/zhoujiazhen/bylw_atac/chromBPNet" in text
    assert "/data1/zhoujiazhen/bylw_atac/TransChromBP" in text
    assert "/home/zhengwei/bylw_atac/TransChromBP" in text
    assert "单向发布" in text


def test_core_docs_reference_governance_doc_and_explicit_publish_commands() -> None:
    for path in ("README.md", "DEVELOPMENT.md", "AGENTS.md"):
        text = read_text(path)
        assert "docs/env/repository_governance.md" in text

    readme = read_text("README.md")
    development = read_text("DEVELOPMENT.md")
    assert "publish-runtime-6000" in readme
    assert "publish-runtime-6002" in readme
    assert "publish-runtime-6000" in development
    assert "pull-results-6002" in development


def test_legacy_root_files_are_rehomed() -> None:
    assert not (REPO_ROOT / "TRAINING_ANALYSIS.md").exists()
    assert not (REPO_ROOT / "REPRODUCTION_NOTES.md").exists()
    assert not (REPO_ROOT / "bpnet_model_annotated.py").exists()
    assert not (REPO_ROOT / "visualize_bpnet.py").exists()

    assert (REPO_ROOT / "reports/chrombpnet_tutorial_bias_training_analysis_20260123.md").exists()
    assert (REPO_ROOT / "docs/research/chrombpnet_reproduction_notes.md").exists()
    assert (REPO_ROOT / "docs/learning/assets/bpnet_model_annotated.py").exists()
    assert (REPO_ROOT / "scripts/visualize_bpnet.py").exists()


def test_supporting_readmes_exist_for_reports_references_and_staging() -> None:
    reports_readme = REPO_ROOT / "reports/README.md"
    references_readme = REPO_ROOT / "references/README.md"
    local_only_readme = REPO_ROOT / "references/local-only/README.md"
    staging_readme = REPO_ROOT / "tmp_remote_edit/README.md"

    assert reports_readme.exists()
    assert references_readme.exists()
    assert local_only_readme.exists()
    assert staging_readme.exists()

    assert "md/tex/assets" in reports_readme.read_text(encoding="utf-8")
    assert "避免和 LaTeX build 日志混淆" in reports_readme.read_text(encoding="utf-8")
    assert "local-only" in references_readme.read_text(encoding="utf-8")
    assert "staging" in staging_readme.read_text(encoding="utf-8")


def test_curated_artifacts_and_ignore_rules_exist() -> None:
    gitignore = read_text(".gitignore")
    assert ".pytest_cache/" in gitignore
    assert "references/local-only/*" in gitignore
    assert "!references/local-only/README.md" in gitignore
    assert "reports/*.aux" in gitignore
    assert "reports/*.fdb_latexmk" in gitignore
    assert "/bias_training_curve.svg" in gitignore
    assert "/bias_training_curve_detailed.svg" in gitignore
    assert "/bias_training_log.csv" in gitignore

    asset_dir = REPO_ROOT / "reports/assets/chrombpnet_bias_training_20260123"
    assert (asset_dir / "bias_training_curve.svg").exists()
    assert (asset_dir / "bias_training_curve_detailed.svg").exists()
    assert (asset_dir / "bias_training_log.csv").exists()

    asset_readme = (asset_dir / "README.md").read_text(encoding="utf-8")
    assert "reports/chrombpnet_tutorial_bias_training_analysis_20260123.md" in asset_readme

    bundle_dir = REPO_ROOT / "reports/chatgpt_bundle_loss_balance_20260416"
    assert (bundle_dir / "00_readme_upload_order_and_prompt.md").exists()
    assert (bundle_dir / "03_questions_for_external_models_dynamic_loss.md").exists()
    assert (bundle_dir / "10_ready_to_send_chatgpt_pro_prompt.md").exists()

    assert (REPO_ROOT / "reports/repository_hygiene_cleanup_20260417.md").exists()
    assert (REPO_ROOT / "vendor/transchrombp/scripts/summarize_loss_balance_selectors.py").exists()


def test_docs_env_is_not_masked_by_root_env_ignore_rule() -> None:
    for path in (
        "docs/env/repository_governance.md",
        "docs/env/transchrombp_genos_env.md",
    ):
        result = subprocess.run(
            ["git", "check-ignore", path],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        assert result.returncode == 1, (
            f"{path} should stay visible to git, but check-ignore reported: "
            f"{result.stdout}{result.stderr}"
        )

    root_env = subprocess.run(
        ["git", "check-ignore", "env/example.txt"],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert root_env.returncode == 0
