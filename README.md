# TransChromBP main repository and archive

This repository is the main TransChromBP repository plus the project archive. It is not the home of the official ChromBPNet source anymore.

## What lives here

- `vendor/transchrombp/`: versioned local snapshot of TransChromBP used by this repo's launchers and experiments.
- `docs/`, `reports/`, `TRACKING.md`: documentation, evidence, and live-status archive.
- `scripts/`, `workflows/`: project launchers and official bridge/sync scripts.
- `tests/`: shell checks for workflow and data-path regressions.
- `images/`: documentation assets.

## Where to use the official ChromBPNet source

- Official ChromBPNet source lookup, reproduction, and compare runs use the 6000 external repo:
  `/data1/zhoujiazhen/bylw_atac/chrombpnet_official`
- The canonical official source is not kept in this repository.
- If you need upstream package layout or packaging files, inspect the 6000 external repo first.

## Practical entry points

- `bash workflows/tutorial/step1_download_bams_and_peaks.sh /path/to/data`
- `bash tests/full_workflow.sh 0`
- `./scripts/sync_project.sh deploy`
- `./scripts/sync_project.sh download_results`

## Notes

- This repo is for project execution, documentation, and evidence retention.
- For official ChromBPNet behavior, compare against the external repo on 6000 before updating local conclusions.
