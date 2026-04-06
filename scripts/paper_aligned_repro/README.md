# ChromBPNet 论文对齐复现脚本

本目录提供“论文对齐”版本的复现工具，目标是把教程的单 fold 流程扩展为：
- 5-fold 染色体划分（non-overlapping）
- 多随机种子重复（建议至少 3 个）
- 自动汇总 counts/profile 指标，输出 mean/std

另外，2026-03-27 起本目录也承接 strict-compare 的官方侧工具：
- `run_tutorial_strict_compare_official.sh`
  - 一键起 tutorial `official fidelity` 或 `official controlled`
  - `official controlled` 会内部开启 `--multi-gpu-train`，避免单 fold 时退化成单卡训练
  - 默认从 `CHROMBPNET_OFFICIAL_ROOT` 解析外置 official root，而不是依赖本仓存在官方 `chrombpnet/`
  - 2026-03-29 起支持 `--run-suffix` / `--shared-peaks` / `--shared-nonpeaks` / `--shared-bigwig`
- `select_best_epoch.py`
  - 对 `chrombpnet.epoch_*.h5` 逐个做轻量评估，并按外部统一 metric 选 best epoch
  - 默认在 `valid` split 上选 best，避免 test 泄漏
  - 默认从 `CHROMBPNET_OFFICIAL_ROOT` 调用外置 official `predict.py`
  - 2026-03-30 起支持 `--gpus 0,1` 按 checkpoint 分片；不改指标口径，只减少 wall-clock

本目录的 operator contract 已改为 external official root：
- 6000 canonical official root：`/data1/zhoujiazhen/bylw_atac/chrombpnet_official`
- 本主仓不再 vendoring 官方 `ChromBPNet` 源码；official compare / official metrics / GC helper 都走外置官方仓
- 仍在使用的官方文件族与补丁来历见 [reports/chrombpnet_official_patch_ledger_20260406.md](/home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-chrombpnet-externalization/reports/chrombpnet_official_patch_ledger_20260406.md)

## 1. 环境准备（6000）

```bash
export CHROMBPNET_ENV=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet
export CHROMBPNET_OFFICIAL_ROOT=/data1/zhoujiazhen/bylw_atac/chrombpnet_official
export PATH="$CHROMBPNET_ENV/bin:$PATH"
export LD_LIBRARY_PATH="$CHROMBPNET_ENV/lib:$LD_LIBRARY_PATH"
export PYTHONPATH=/data1/zhoujiazhen/bylw_atac/chromBPNet/vendor/transchrombp:$PYTHONPATH
```

说明：
- `CHROMBPNET_OFFICIAL_ROOT` 是 strict-compare launcher 和 `select_best_epoch.py` 的默认 official root；不显式传 `--official-root` 时，就按这个环境变量解析。
- `PYTHONPATH` 只指向当前主仓里的 `vendor/transchrombp` snapshot，不再把 `/data1/zhoujiazhen/bylw_atac/chromBPNet` 整仓当成官方 `ChromBPNet` 源码根。
- 如果你要核对 official `predict.py` / `metrics.py` / GC helper 现在为什么仍然需要，直接看 [reports/chrombpnet_official_patch_ledger_20260406.md](/home/zhengwei/.config/superpowers/worktrees/chromBPNet/autonomy-20260406-chrombpnet-externalization/reports/chrombpnet_official_patch_ledger_20260406.md)。

## 2. 生成 5-fold 划分

默认输出 5 个 fold，测试/验证对分别为：
- fold_0: test chr1, valid chr2
- fold_1: test chr3, valid chr4
- fold_2: test chr5, valid chr6
- fold_3: test chr7, valid chr8
- fold_4: test chr9, valid chr10

```bash
python scripts/paper_aligned_repro/generate_fold_set.py \
  --chrom-sizes /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes \
  --output-dir /data1/zhoujiazhen/bylw_atac/chrombpnet_paper_repro/folds
```

如果你后续拿到论文精确 fold 定义，可用 `--pairs` 覆盖默认划分。

## 3. 启动批量复现（后台）

下面示例是 K562 ATAC 的 5-fold × 3-seed：

```bash
nohup bash scripts/paper_aligned_repro/run_paper_aligned_matrix.sh \
  --name K562_ATAC \
  --genome /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.fa \
  --chrom-sizes /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes \
  --bam /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/merged.bam \
  --peaks /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/overlap.bed.gz \
  --blacklist /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/blacklist.bed.gz \
  --fold-dir /data1/zhoujiazhen/bylw_atac/chrombpnet_paper_repro/folds \
  --work-root /data1/zhoujiazhen/bylw_atac/chrombpnet_paper_repro \
  --seeds "1234 2345 3456" \
  --gpu 0,1 \
  > /data1/zhoujiazhen/bylw_atac/chrombpnet_paper_repro/logs/k562_matrix.log 2>&1 &
```

## 4. 查看进度

```bash
tail -f /data1/zhoujiazhen/bylw_atac/chrombpnet_paper_repro/logs/k562_matrix.log
grep -RIn "Completed execution\\|failed with exit code" /data1/zhoujiazhen/bylw_atac/chrombpnet_paper_repro/runs/K562_ATAC | head
```

## 5. 汇总指标

```bash
python scripts/paper_aligned_repro/summarize_metrics.py \
  --run-root /data1/zhoujiazhen/bylw_atac/chrombpnet_paper_repro/runs/K562_ATAC \
  --output-dir /data1/zhoujiazhen/bylw_atac/chrombpnet_paper_repro/runs/K562_ATAC/summary
```

输出文件：
- `run_metrics.csv`：每个 fold/seed 一行
- `summary_by_seed.csv`：按 seed 汇总
- `summary_by_fold.csv`：按 fold 汇总
- `summary_overall.json`：总体均值/标准差

## 6. Strict Compare（tutorial 官方侧）

### 启 tutorial 官方 fidelity

```bash
bash scripts/paper_aligned_repro/run_tutorial_strict_compare_official.sh \
  --mode fidelity \
  --work-root /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare \
  --seed 42 \
  --folds "fold_0" \
  --gpus 0,1
```

### 启 tutorial 官方 controlled

```bash
bash scripts/paper_aligned_repro/run_tutorial_strict_compare_official.sh \
  --mode controlled \
  --work-root /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare \
  --seed 42 \
  --folds "fold_0" \
  --gpus 0,1
```

### 启 tutorial 官方 controlled L3（shared-region）

```bash
SHARED_DIR=/data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/shared_regions/tutorial

bash scripts/paper_aligned_repro/run_tutorial_strict_compare_official.sh \
  --mode controlled \
  --work-root /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare \
  --run-suffix "_L3" \
  --seed 42 \
  --folds "fold_0" \
  --gpus 0,1 \
  --shared-peaks "${SHARED_DIR}/shared_filtered_peaks.bed" \
  --shared-nonpeaks "${SHARED_DIR}/shared_filtered_nonpeaks.bed" \
  --shared-bigwig "${SHARED_DIR}/shared_unstranded.bw"
```

说明：
- 传入 `--shared-nonpeaks` 后，底层脚本会跳过 `chrombpnet prep nonpeaks`
- 传入 `--shared-bigwig` 后，official 侧评估会复用这份共享 bigWig
- `run_tutorial_strict_compare_official.sh` 与下游 `run_paper_aligned_fast_1seed.sh` 默认都依赖 `CHROMBPNET_OFFICIAL_ROOT`

### 用外部规则选 best epoch

```bash
python scripts/paper_aligned_repro/select_best_epoch.py \
  --model-glob '/data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/runs/tutorial_official_controlled_s42/fold_0/seed_42/chrombpnet/models/chrombpnet.epoch_*.h5' \
  --genome /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.fa \
  --bigwig /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/runs/tutorial_official_controlled_s42/fold_0/seed_42/chrombpnet/auxiliary/data_unstranded.bw \
  --peaks /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/runs/tutorial_official_controlled_s42/fold_0/seed_42/chrombpnet/auxiliary/filtered.peaks.bed \
  --nonpeaks None \
  --fold-json /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/folds/fold_0.json \
  --output-dir /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/runs/tutorial_official_controlled_s42/fold_0/seed_42/chrombpnet/external_best \
  --split valid \
  --gpus 0,1 \
  --metric profile_metrics.peaks.median_jsd \
  --mode min
```

说明：
- 多 checkpoint selector 可以传 `--gpus 0,1` 做分片并行；这不会改变每个 checkpoint 的评估逻辑，只是把不同 checkpoint 分给不同 GPU。
- 单 checkpoint 的 held-out test 没有天然的双卡收益；这类任务更适合和另一条独立 test/selector 分别占两张卡并行。
- 不传 `--official-root` 时，`select_best_epoch.py` 默认使用 `CHROMBPNET_OFFICIAL_ROOT=/data1/zhoujiazhen/bylw_atac/chrombpnet_official`。
