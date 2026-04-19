# TransChromBP main repository and archive

This repository is the main TransChromBP repository plus the project archive. It is not the home of the official ChromBPNet source anymore.

Canonical repository topology, single-source-of-truth rules, and sync commands live in `docs/env/repository_governance.md`.

The active local canonical root is `/home/zhengwei/project/python/TransChromBP`. If `/home/zhengwei/project/python/chromBPNet` still exists locally, treat it only as a compatibility symlink / rollback path.

## What lives here

- `vendor/transchrombp/`: versioned local snapshot of TransChromBP used by this repo's launchers and experiments.
- `docs/`, `reports/`, `TRACKING.md`: documentation, evidence, and live-status archive.
- `scripts/`, `workflows/`: project launchers and explicit publish/pull helpers.
- `tests/`: shell checks for workflow and data-path regressions.
- `images/`: documentation assets.
- `references/`: local-only reference indices; heavy materials live under `references/local-only/`.
- `tmp_remote_edit/`: remote-write staging only, not a final archive.

## Where to use the official ChromBPNet source

- Official ChromBPNet source lookup, reproduction, and compare runs use the 6000 external repo:
  `/data1/zhoujiazhen/bylw_atac/chromBPNet`
- The canonical official source is not kept in this repository.
- If you need upstream package layout or packaging files, inspect the 6000 external repo first.

## Practical entry points

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
