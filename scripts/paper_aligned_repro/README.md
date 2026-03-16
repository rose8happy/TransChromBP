# ChromBPNet 论文对齐复现脚本

本目录提供“论文对齐”版本的复现工具，目标是把教程的单 fold 流程扩展为：
- 5-fold 染色体划分（non-overlapping）
- 多随机种子重复（建议至少 3 个）
- 自动汇总 counts/profile 指标，输出 mean/std

## 1. 环境准备（6000）

```bash
export CHROMBPNET_ENV=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet
export PATH="$CHROMBPNET_ENV/bin:$PATH"
export LD_LIBRARY_PATH="$CHROMBPNET_ENV/lib:$LD_LIBRARY_PATH"
export PYTHONPATH=/data1/zhoujiazhen/bylw_atac/chromBPNet:$PYTHONPATH
```

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

