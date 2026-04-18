# ChromBPNet Reproduction Notes

## Remote Environment
- Host: `ssh zhoujiazhen@127.0.0.1 -p 6000`
- Repo: `/data1/zhoujiazhen/bylw_atac/chromBPNet`
- Working dir: `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial`
- Python env: `/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet`
- CLI binary: `/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin/chrombpnet` (use full path if not on `PATH`)

## Existing Reference Genome
- `/data1/zhoujiazhen/bylw_atac/hg38.fa` and `/data1/zhoujiazhen/bylw_atac/hg38.fa.fai` already exist.
- If skipping the tutorial download, you can symlink these into `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/` as `hg38.fa` and `hg38.fa.fai`.

## Step 1: Download tutorial inputs without reference genome (completed)
- Previous full step1 download was stopped to avoid re-downloading GRCh38.
- Command:
  `nohup bash -c "set -euo pipefail; wget -c https://www.encodeproject.org/files/ENCFF077FBI/@@download/ENCFF077FBI.bam -O rep1.bam; wget -c https://www.encodeproject.org/files/ENCFF128WZG/@@download/ENCFF128WZG.bam -O rep2.bam; wget -c https://www.encodeproject.org/files/ENCFF534DCE/@@download/ENCFF534DCE.bam -O rep3.bam; wget -c https://www.encodeproject.org/files/ENCFF333TAT/@@download/ENCFF333TAT.bed.gz -O overlap.bed.gz; wget -c https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz -O blacklist.bed.gz; wget -c https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv -O hg38.chrom.sizes; ln -sf /data1/zhoujiazhen/bylw_atac/hg38.fa hg38.fa; ln -sf /data1/zhoujiazhen/bylw_atac/hg38.fa.fai hg38.fa.fai" > /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step1_download_no_genome.log 2>&1 &`
- Started: `2026-01-21 17:19:53 CST`
- Completed: `2026-01-22 01:02:41 CST` (last log entry)
- Log: `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step1_download_no_genome.log`
- Data downloaded (ENCODE tutorial dataset):
  - `ENCFF077FBI.bam`, `ENCFF128WZG.bam`, `ENCFF534DCE.bam`: ATAC-seq alignments (K562), experiment `ENCSR868FGK` (merged into `merged.bam`).
  - `ENCFF333TAT.bed.gz`: pseudoreplicated narrowPeak peaks for `ENCSR868FGK`.
  - `ENCFF356LFX.bed.gz`: ENCODE blacklist/exclusion regions.
  - `GRCh38_EBV.chrom.sizes.tsv`: chromosome sizes for GRCh38 (saved as `hg38.chrom.sizes`).
  - `hg38.fa`/`hg38.fa.fai`: symlinked from `/data1/zhoujiazhen/bylw_atac/`.

## Download Validation
- Files present in `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data`:
  - `rep1.bam`: 7,496,509,558 bytes (ENCODE metadata: 7,492,959,554)
  - `rep2.bam`: 7,056,946,200 bytes
  - `rep3.bam`: 10,012,174,956 bytes
  - `overlap.bed.gz`: 7,067,871 bytes
  - `blacklist.bed.gz`: 8,211 bytes
  - `hg38.chrom.sizes`: 11,686 bytes
  - `hg38.fa`/`hg38.fa.fai`: symlinks to existing reference
- Integrity check: `samtools quickcheck` on `rep1/rep2/rep3` exited 0.

## Rep1 Re-download for Verification (running)
- Purpose: download `rep1.bam` to a separate location and compare md5 with the existing file.
- Command:
  `nohup bash -c "set -euo pipefail; cd /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/rep1_check; EXPECTED=135e5a65e5eafac97a8a2242b4c13963; echo expected_md5=$EXPECTED; wget -c https://www.encodeproject.org/files/ENCFF077FBI/@@download/ENCFF077FBI.bam -O rep1.bam; md5sum rep1.bam /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/rep1.bam; ls -l rep1.bam /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/rep1.bam" > /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/rep1_redownload.log 2>&1 &`
- PID: `97445`
- Started: `2026-01-22 09:38:27 CST`
- Log: `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/rep1_redownload.log`
- Check progress:
  - `tail -f /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/rep1_redownload.log`
  - `ps -p 97445`

## Step 2: Splits + Nonpeaks (running)
- Command:
  `nohup bash -c "set -euo pipefail; /data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin/chrombpnet prep splits -op /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/folds -c /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes -tcr chr1 -vcr chr2; /data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin/chrombpnet prep nonpeaks -g /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.fa -o /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/nonpeaks -p /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/overlap.bed.gz -c /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes -fl /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/folds.json -br /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/blacklist.bed.gz" > /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step2_splits_nonpeaks.log 2>&1 &`
- PID: `97526`
- Started: `2026-01-22 09:38:35 CST`
- Log: `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step2_splits_nonpeaks.log`
- Check progress:
  - `tail -f /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step2_splits_nonpeaks.log`
  - `ps -p 97526`
  - Outputs: `data/folds.json`, `data/nonpeaks_negatives.bed`

## Step 2 Status (completed)
- `folds.json` and `nonpeaks_negatives.bed` created.
- Log shows no filtered regions due to inputlen constraints.

## Step 3: Merge BAMs + BigWig (failed)
- Command:
  `nohup bash -c "set -euo pipefail; cd /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data; samtools merge -f merged_unsorted.bam rep1.bam rep2.bam rep3.bam; samtools sort -@4 merged_unsorted.bam -o merged.bam; samtools index merged.bam; PATH=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin:$PATH bash /data1/zhoujiazhen/bylw_atac/chromBPNet/workflows/tutorial/step2_make_bigwigs_from_bams.sh /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged.bam /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged ATAC /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.fa /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes" > /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step3_merge_bigwig.log 2>&1 &`
- Started: `2026-01-22 09:57:56 CST`
- Log: `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step3_merge_bigwig.log`
- Failure:
  - `samtools merge: "rep1.bam" is truncated`
  - `Invalid BGZF header at offset 7492959554`
- Action:
  - Re-download `rep1.bam` to a separate location and verify md5.
  - Replace `data/rep1.bam` after verification, then re-run Step 3.

## Rep1 Verification Result (completed)
- `rep1.bam` re-download md5: `135e5a65e5eafac97a8a2242b4c13963` (matches ENCODE).
- Replaced `data/rep1.bam` with the verified file (size 7,492,959,554).
- `samtools quickcheck` passes for the new `rep1.bam`.

## Step 3 Retry: Merge BAMs + BigWig (running)
- Command:
  `nohup bash -c "set -euo pipefail; cd /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data; rm -f merged_unsorted.bam merged.bam merged.bam.bai merged_unstranded.bw merged_bias_pwm.*; samtools merge -f merged_unsorted.bam rep1.bam rep2.bam rep3.bam; samtools sort -@4 merged_unsorted.bam -o merged.bam; samtools index merged.bam; PATH=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin:$PATH bash /data1/zhoujiazhen/bylw_atac/chromBPNet/workflows/tutorial/step2_make_bigwigs_from_bams.sh /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged.bam /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged ATAC /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.fa /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes" > /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step3_merge_bigwig_retry.log 2>&1 &`
- PID: `178657`
- Started: `2026-01-22 14:47:04 CST`
- Log: `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step3_merge_bigwig_retry.log`
- Check progress:
  - `tail -f /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step3_merge_bigwig_retry.log`
  - `ps -p 178657`

## Step 3 Retry Status (failed)
- Reason: `chrombpnet_makebigwig` CLI missing in this environment.
- `merged.bam` and `merged.bam.bai` were created successfully.

## Step 3b: BigWig + PWM via Python modules (restart with PATH)
- Previous run used a PATH without `bedGraphToBigWig`; restarted with env PATH.
- Command:
  `nohup bash -c "set -euo pipefail; export PATH=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin:$PATH; rm -f /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged_unstranded.bw /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged_bias_pwm.png; /data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin/python -m chrombpnet.helpers.preprocessing.reads_to_bigwig -g /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.fa -ibam /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged.bam -c /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes -op /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged -d ATAC; /data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin/python -m chrombpnet.helpers.preprocessing.analysis.build_pwm_from_bigwig -i /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged_unstranded.bw -g /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.fa -op /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged_bias_pwm -cr chr20 -c /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes -pw 24" > /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step3_bigwig_only_retry.log 2>&1 &`
- PID: `208867`
- Started: `2026-01-22 16:09:23 CST`
- Log: `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step3_bigwig_only_retry.log`
- Check progress:
  - `tail -f /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step3_bigwig_only_retry.log`
  - `ps -p 208867`

## Step 3b Status (completed)
- Outputs created:
  - `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged_unstranded.bw`
  - `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged_bias_pwm.png`

## Step 4: Bias Model Training (running)
- Command:
  `nohup bash -c "set -euo pipefail; /data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin/chrombpnet bias pipeline -g /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.fa -c /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes -ibam /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged.bam -o /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/bias -d ATAC -p /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/overlap.bed.gz -n /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/nonpeaks_negatives.bed -fl /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/folds.json -b 0.5" > /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step4_bias_train.log 2>&1 &`
- PID: `214619`
- Started: `2026-01-22 17:13:34 CST`
- Log: `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step4_bias_train.log`
- Check progress:
  - `tail -f /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step4_bias_train.log`
  - `ps -p 214619`

## Step 4: Bias Model Training (restart with GPU)
- Stopped previous CPU run and restarted with `LD_LIBRARY_PATH` pointing to env CUDA libs.
- Command:
  `nohup bash -c "set -euo pipefail; export LD_LIBRARY_PATH=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/lib:$LD_LIBRARY_PATH; export PATH=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin:$PATH; /data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin/python - <<\"PY\"\\nimport tensorflow as tf\\nprint(\"gpus\", tf.config.list_physical_devices(\"GPU\"))\\nPY\\n/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin/chrombpnet bias pipeline -g /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.fa -c /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes -ibam /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged.bam -o /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/bias -d ATAC -p /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/overlap.bed.gz -n /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/nonpeaks_negatives.bed -fl /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/folds.json -b 0.5" > /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step4_bias_train_gpu.log 2>&1 &`
- PID: `217298`
- Started: `2026-01-22 17:43:02 CST`
- Log: `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step4_bias_train_gpu.log`
- Check progress:
  - `tail -f /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step4_bias_train_gpu.log`
  - `ps -p 217298`
- Status: failed due to bash heredoc syntax; training never started.

## Step 4: Bias Model Training (retry with GPU, fixed heredoc)
- Command:
  `nohup bash -lc "set -euo pipefail; export PATH=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin:$PATH; export LD_LIBRARY_PATH=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/lib:${LD_LIBRARY_PATH:-}; export CUDA_VISIBLE_DEVICES=0; /data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin/chrombpnet bias pipeline -g /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.fa -c /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes -ibam /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged.bam -o /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/bias -d ATAC -p /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/overlap.bed.gz -n /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/nonpeaks_negatives.bed -fl /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/folds.json -b 0.5" > /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step4_bias_train_gpu_retry.log 2>&1 &`
- PID: `217701`
- Started: `2026-01-22 18:05:13 CST`
- Log: `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step4_bias_train_gpu_retry.log`
- Check progress:
  - `tail -f /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step4_bias_train_gpu_retry.log`
  - `ps -p 217701`
- Status: failed because `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/bias/logs` already existed.
- Action: archived previous outputs to `outputs/bias_prev_20260122_174740`.

## Step 4: Bias Model Training (retry after archiving outputs)
- Command:
  `nohup bash -lc "set -euo pipefail; export PATH=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin:$PATH; export LD_LIBRARY_PATH=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/lib:${LD_LIBRARY_PATH:-}; export CUDA_VISIBLE_DEVICES=0; /data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin/chrombpnet bias pipeline -g /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.fa -c /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes -ibam /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged.bam -o /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/bias -d ATAC -p /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/overlap.bed.gz -n /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/nonpeaks_negatives.bed -fl /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/folds.json -b 0.5" > /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step4_bias_train_gpu_retry2.log 2>&1 &`
- PID: `218082`
- Started: `2026-01-22 18:07:26 CST`
- Log: `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step4_bias_train_gpu_retry2.log`
- Check progress:
  - `tail -f /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/step4_bias_train_gpu_retry2.log`
  - `ps -p 218082`

## Step 4 Status (completed)
- Log tail ended with: `Saving 'profile' scores`.
- PID `218082` no longer running.
- Outputs created under `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/bias`, including:
  - Model: `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/bias/models/bias.h5`
  - Reports: `evaluation/overall_report.html`, `evaluation/overall_report.pdf`
  - Metrics: `evaluation/bias_metrics.json`
## Step 3b Warning Note
- Warning observed: BAM contains contigs not present in the FASTA.
- Detected extra contigs in BAM: `chrEBV` and `*` (unmapped). This is expected and can be ignored.
- FASTA does not include `chrEBV`, so those reads are filtered out during preprocessing.

## Environment Fixes (completed)
- Added wrapper CLIs in `/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin`:
  - `chrombpnet_makebigwig` -> `python -m chrombpnet.helpers.preprocessing.reads_to_bigwig`
  - `chrombpnet_pwm_from_bigwig` -> `python -m chrombpnet.helpers.preprocessing.analysis.build_pwm_from_bigwig`
- `bedGraphToBigWig` exists in the same env bin; ensure commands run with `PATH=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin:$PATH` when calling tutorial scripts.

## GPU Status (needs fix)
- TensorFlow version in env: `2.8.0` (built with CUDA, but expects CUDA 11 libraries).
- CUDA 11 + cuDNN 8 are installed inside the micromamba env, but not visible unless `LD_LIBRARY_PATH` includes the env `lib`.
- Fix: `LD_LIBRARY_PATH=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/lib:$LD_LIBRARY_PATH`
- With this env var set, `tf.config.list_physical_devices('GPU')` returns both GPUs.

## Additional ATAC Downloads (running)
## Dataset Roles
- chrombpnet tutorial dataset (K562 ENCSR868FGK) is used only to reproduce the official ChromBPNet tutorial.
  - Location: `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data`
- Additional datasets (GM12878, K562 ENCSR859USB, HeLa-S3 SRA) are for downstream experiments after the tutorial.
  - Location: `/data1/zhoujiazhen/bylw_atac/chrombpnet_datasets`

### GM12878 (ENCSR637XSC)
- BAMs: `ENCFF962FMH`, `ENCFF981FXV`, `ENCFF440GRZ`
- Peaks: `ENCFF748UZH` (pseudoreplicated peaks)
- Command:
  `nohup bash -c "set -euo pipefail; cd /data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/GM12878; wget -c https://www.encodeproject.org/files/ENCFF962FMH/@@download/ENCFF962FMH.bam -O rep1.bam; wget -c https://www.encodeproject.org/files/ENCFF981FXV/@@download/ENCFF981FXV.bam -O rep2.bam; wget -c https://www.encodeproject.org/files/ENCFF440GRZ/@@download/ENCFF440GRZ.bam -O rep3.bam; wget -c https://www.encodeproject.org/files/ENCFF748UZH/@@download/ENCFF748UZH.bed.gz -O overlap.bed.gz; ln -sf /data1/zhoujiazhen/bylw_atac/hg38.fa hg38.fa; ln -sf /data1/zhoujiazhen/bylw_atac/hg38.fa.fai hg38.fa.fai; ln -sf /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes hg38.chrom.sizes; ln -sf /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/blacklist.bed.gz blacklist.bed.gz" > /data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/logs/GM12878_download.log 2>&1 &`
- PID: `126343`
- Started: `2026-01-22 11:24:52 CST`
- Log: `/data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/logs/GM12878_download.log`

### K562 (ENCSR859USB)
- BAMs: `ENCFF557RSA`, `ENCFF852RRY`
- Peaks: `ENCFF444DLQ` (pseudoreplicated peaks)
- Command:
  `nohup bash -c "set -euo pipefail; cd /data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/K562; wget -c https://www.encodeproject.org/files/ENCFF557RSA/@@download/ENCFF557RSA.bam -O rep1.bam; wget -c https://www.encodeproject.org/files/ENCFF852RRY/@@download/ENCFF852RRY.bam -O rep2.bam; wget -c https://www.encodeproject.org/files/ENCFF444DLQ/@@download/ENCFF444DLQ.bed.gz -O overlap.bed.gz; ln -sf /data1/zhoujiazhen/bylw_atac/hg38.fa hg38.fa; ln -sf /data1/zhoujiazhen/bylw_atac/hg38.fa.fai hg38.fa.fai; ln -sf /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes hg38.chrom.sizes; ln -sf /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/blacklist.bed.gz blacklist.bed.gz" > /data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/logs/K562_download.log 2>&1 &`
- PID: `126443`
- Started: `2026-01-22 11:25:03 CST`
- Log: `/data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/logs/K562_download.log`

### HeLa-S3
- ENCODE ATAC-seq has released bigWig signal tracks but no BAM alignments/peaks found.
- Switched to GEO/SRA source (HeLa-S3 ATAC-seq; BioSample SAMN10391739 confirms HeLa-S3).

### HeLa-S3 (GEO/SRA: SRP166944)
- Runs: `SRR8171302`, `SRR8171303`, `SRR8171304`, `SRR8171305`
- Command:
  `nohup bash -c "set -euo pipefail; prefetch --max-size 200G -O /data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/HeLa-S3 SRR8171302 SRR8171303 SRR8171304 SRR8171305" > /data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/logs/HeLa-S3_download.log 2>&1 &`
- PID: `131179`
- Started: `2026-01-22 11:37:26 CST`
- Log: `/data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/logs/HeLa-S3_download.log`

## Next Steps (after step 1 completes)
- Run long steps with `nohup ... > /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/logs/<step>.log 2>&1 &` and track via `tail -f`.
- Generate chromosome splits:
  `/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin/chrombpnet prep splits -op /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/folds -c /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes -tcr chr1 chr2 -vcr chr3`
- Generate GC-matched nonpeaks:
  `/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin/chrombpnet prep nonpeaks -g /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.fa -o /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/nonpeaks -p /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/overlap.bed.gz -c /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes -fl /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/folds.json -br /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/blacklist.bed.gz`
- Train bias model:
  `/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin/chrombpnet bias pipeline -g /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.fa -c /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes -ibam /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged.bam -o /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/bias -d ATAC -p /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/overlap.bed.gz -n /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/nonpeaks_negatives.bed -fl /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/folds.json -b 0.5`
- Train ChromBPNet with bias model:
  `/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet/bin/chrombpnet pipeline -g /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.fa -c /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes -ibam /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged.bam -o /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/chrombpnet -d ATAC -p /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/overlap.bed.gz -n /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/nonpeaks_negatives.bed -fl /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/folds.json -b /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/bias/models/bias.h5`
