# AlphaGenome Pilot

这个目录放的是一个最小 AlphaGenome 对比试验脚手架，目标是：

- 不改我们当前 tutorial 主 benchmark；
- 先在少量位点上做 `AlphaGenome vs ChromBPNet vs TransChromBP` 的 black-box 对比；
- 不把 AlphaGenome API 输出用于训练任何学生模型。

## 默认假设

- tutorial 数据在我们项目里按 `K562` 处理。
- 对应的 AlphaGenome ontology 默认使用 `EFO:0002067`。
- 目前默认只跑 `ATAC` 输出。

依据：

- tutorial 数据在现有脚本里被当作 `K562_ATAC`：
  [scripts/paper_aligned_repro/launch_k562_tutorial_baseline.sh](/home/zhengwei/project/python/chromBPNet/scripts/paper_aligned_repro/launch_k562_tutorial_baseline.sh)
- AlphaGenome 官方 notebook 中示例给出了 `K562 -> EFO:0002067`。

## 文件说明

- `regions_k562_tutorial_selected_loci.csv`
  4 个代表性 tutorial test 位点，来自我们已有的 locus 分析。
- `run_alphagenome_pilot.py`
  调 AlphaGenome API，输出每个位点的 profile 和总量统计。
- `merge_locus_totals.py`
  把 AlphaGenome 结果和我们已有的 `selected_loci_prediction_totals.csv` 合并。

## 运行前提

1. 你的 key 已放到仓库外：

```bash
~/.config/alphagenome/api.env
```

文件里建议至少包含下面任意一种：

```bash
export ALPHAGENOME_API_KEY='...'
```

或

```bash
export ALPHA_GENOME_API_KEY='...'
```

2. 运行脚本的 Python 需要满足：

- Python `>=3.10`
- 已安装 `alphagenome`

注意：

- 项目当前在 6000 机器上默认使用的 `chrombpnet` 环境是
  `/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet`
- 该环境当前实际是 `Python 3.9.x`
- 因此不建议直接把 `alphagenome` 装进这个现有环境里
- 更稳妥的做法是：保留 `chrombpnet` 环境给现有训练/评估流程，另外单独建一个
  `Python >= 3.10` 的 AlphaGenome 专用环境

官方 SDK 仓库说明：

- https://github.com/google-deepmind/alphagenome
- https://www.alphagenomedocs.com/installation.html

## 最小运行示例

```bash
/path/to/python310_or_newer/bin/python scripts/alphagenome_pilot/run_alphagenome_pilot.py \
  --regions-csv scripts/alphagenome_pilot/regions_k562_tutorial_selected_loci.csv \
  --output-dir outputs/alphagenome_pilot/tutorial_selected_loci
```

跑完后会得到：

- `summary.csv`
- `region_metadata.jsonl`
- `profiles/*.npz`
- `run_meta.json`

## 和现有 locus 总量表合并

```bash
/path/to/python310_or_newer/bin/python scripts/alphagenome_pilot/merge_locus_totals.py \
  --alpha-summary outputs/alphagenome_pilot/tutorial_selected_loci/summary.csv \
  --local-totals reports/assets/transchrombp_tutorial_test_20260315/selected_loci_prediction_totals.csv \
  --output-csv outputs/alphagenome_pilot/tutorial_selected_loci/merged_locus_totals.csv
```

## 当前边界

- 这是小规模 pilot，不是正式 benchmark。
- 目前默认只处理少量位点。
- 还没有接入大批量 held-out 评测。
- 不允许把 AlphaGenome API 输出拿来训练学生模型。
