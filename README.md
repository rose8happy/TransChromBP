# TransChromBP Main Repository And Archive

This repository is the local canonical trunk for TransChromBP plus the project archive. Official ChromBPNet source lookup, reproduction, and compare runs no longer live here.

Canonical topology, source-of-truth rules, and sync commands live in `docs/env/repository_governance.md`.

The active local canonical root is `/home/zhengwei/project/python/TransChromBP`. If `/home/zhengwei/project/python/chromBPNet` still exists locally, treat it only as a compatibility symlink / rollback path.

## Canonical Entrypoints

- Live status: [TRACKING.md](TRACKING.md)
- Experiment family / run index: [docs/experiments/registry.md](docs/experiments/registry.md), [docs/experiments/runs.csv](docs/experiments/runs.csv)
- Learning / onboarding navigation: [docs/learning/学习材料索引.md](docs/learning/学习材料索引.md)
- Reusable reports and closeouts: [reports/README.md](reports/README.md)

## What Lives Here

- `vendor/transchrombp/`: versioned local snapshot used by this repo's launchers and experiments
- `docs/`: governance, plans, research notes, learning materials, and archived runbooks
- `reports/`: reusable analysis, closeout, governance, handoff, external, and paper materials
- `scripts/`, `workflows/`: project launchers and publish / pull helpers
- `tests/`: shell and doc-contract checks
- `references/`: local-only reference indices; heavy materials live under `references/local-only/`
- `tmp_remote_edit/`: remote-write staging only, not a final archive

## Official ChromBPNet Source

- Official ChromBPNet source lookup, reproduction, and compare runs use the 6000 external repo:
  `/data1/zhoujiazhen/bylw_atac/chromBPNet`
- The canonical official source is not kept in this repository.
- If you need upstream package layout or packaging files, inspect the 6000 external repo first.

## Practical Commands

- `bash workflows/tutorial/step1_download_bams_and_peaks.sh /path/to/data`
- `bash tests/full_workflow.sh 0`
- `./scripts/sync_project.sh publish-runtime-6000 --dry-run`
- `./scripts/sync_project.sh publish-runtime-6002 --dry-run`
- `./scripts/sync_project.sh pull-results-6000 --dry-run`
- `./scripts/sync_project.sh status-all`

## Notes

- This repo is for project execution, documentation, and evidence retention.
- For official ChromBPNet behavior, compare against the external repo on 6000 before updating local conclusions.
- Remote hotfixes must be recovered into the local canonical trunk before they are published back to 6000/6002.
- Directory-specific hygiene rules live in `reports/README.md`, `references/README.md`, and `tmp_remote_edit/README.md`.
