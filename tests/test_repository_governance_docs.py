from __future__ import annotations

from pathlib import Path
import subprocess


REPO_ROOT = Path(__file__).resolve().parents[1]


def read_text(path: str) -> str:
    return (REPO_ROOT / path).read_text(encoding="utf-8")


def parse_markdown_table(path: str) -> list[dict[str, str]]:
    text = read_text(path)
    rows = [line.strip() for line in text.splitlines() if line.strip().startswith("|")]
    header = [cell.strip() for cell in rows[0].strip("|").split("|")]
    data_rows = rows[2:]

    parsed_rows = []
    for row in data_rows:
        cells = [cell.strip() for cell in row.strip("|").split("|")]
        if len(cells) != len(header):
            continue
        parsed_rows.append(dict(zip(header, cells)))

    return parsed_rows


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


def test_governance_doc_uses_transchrombp_local_root() -> None:
    text = read_text("docs/env/repository_governance.md")
    assert "/home/zhengwei/project/python/TransChromBP" in text
    assert "| `/home/zhengwei/project/python/chromBPNet` |" not in text


def test_no_fake_transchrombp_worktree_paths_exist() -> None:
    tree = subprocess.run(
        [
            "rg",
            "-n",
            "/home/zhengwei/\\.config/superpowers/worktrees/TransChromBP",
            "AGENTS.md",
            "README.md",
            "TRACKING.md",
            "docs",
            "reports",
            "scripts",
            "tests",
            "workflows",
        ],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert tree.stdout.strip() == ""


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


def test_docs_plan_is_live_home_for_factor_ladder_and_superpowers_is_archive_only() -> None:
    tracking = read_text("TRACKING.md")
    agents = read_text("AGENTS.md")
    superpowers_readme = read_text("docs/superpowers/README.md")

    assert "docs/plan/alphagenome_like_factor_ladder_design_20260411.md" in tracking
    assert "docs/plan/alphagenome_like_factor_ladder_implementation_20260411.md" in tracking
    assert "legacy archive" in agents
    assert "not a second live planning root" in superpowers_readme
    assert (REPO_ROOT / "docs/plan/alphagenome_like_factor_ladder_design_20260411.md").exists()
    assert (REPO_ROOT / "docs/plan/alphagenome_like_factor_ladder_implementation_20260411.md").exists()


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
    superpowers_readme = REPO_ROOT / "docs/superpowers/README.md"

    assert reports_readme.exists()
    assert references_readme.exists()
    assert local_only_readme.exists()
    assert staging_readme.exists()
    assert superpowers_readme.exists()

    assert "md/tex/assets" in reports_readme.read_text(encoding="utf-8")
    assert "避免和 LaTeX build 日志混淆" in reports_readme.read_text(encoding="utf-8")
    assert "冻结代码/配置快照" in reports_readme.read_text(encoding="utf-8")
    assert "local-only" in references_readme.read_text(encoding="utf-8")
    assert "staging" in staging_readme.read_text(encoding="utf-8")
    assert "archive namespace" in superpowers_readme.read_text(encoding="utf-8")


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
    assert (REPO_ROOT / "reports/chrombpnet_official_externalization_closeout_20260408.md").exists()
    assert (REPO_ROOT / "reports/chrombpnet_official_patch_ledger_20260406.md").exists()
    assert (REPO_ROOT / "reports/foundation_cache_contract_snapshot_20260406.md").exists()
    assert (REPO_ROOT / "reports/dual_track_pivot_snapshot_20260409.md").exists()
    assert (REPO_ROOT / "reports/foundation_cache_contract_closeout_20260419.md").exists()
    assert not (REPO_ROOT / "scripts/deploy_strict_compare_staging_to_6000.sh").exists()


def test_registry_allows_archive_branches_without_mounted_worktrees() -> None:
    rows = parse_markdown_table("docs/experiments/registry.md")
    indexed = {row["family_id"]: row for row in rows}

    assert indexed["`unet_lite_v1`"]["mounted_worktree"] == "`n/a`"
    assert indexed["`ntv2_teacher_distill_tutorial`"]["mounted_worktree"] == "`n/a`"
    assert indexed["`foundation_cache_contract`"]["mounted_worktree"] == "`n/a`"
    assert indexed["`chrombpnet_externalization`"]["mounted_worktree"] == "`n/a`"

    assert "dual_track_pivot_snapshot_20260409.md" in indexed["`unet_lite_v1`"]["notes"]
    assert "closeout/foundation_cache_contract/20260419" in indexed["`foundation_cache_contract`"]["closeout_tags"]
    assert "closeout/chrombpnet_externalization/20260408" in indexed["`chrombpnet_externalization`"]["closeout_tags"]


def test_archive_closeout_tags_exist() -> None:
    tag_list = subprocess.run(
        [
            "git",
            "tag",
            "--list",
            "closeout/chrombpnet_externalization/20260408",
            "closeout/foundation_cache_contract/20260419",
        ],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=True,
    ).stdout

    assert "closeout/chrombpnet_externalization/20260408" in tag_list
    assert "closeout/foundation_cache_contract/20260419" in tag_list


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
