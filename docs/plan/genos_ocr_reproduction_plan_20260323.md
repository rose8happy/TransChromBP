# Genos OCR 官方基准旁线执行计划

更新日期：2026-03-23

基础材料：`docs/research/genos复现ocr.md`

## 0. 定位

这不是新的 Genos 主线计划，而是一条有边界的旁线：

- 目标：复现官方 `human_ocr_ensembl` OCR benchmark，验证本地 `Genos-1.2B` 的 embedding 提取与 pooling 链路
- 作用：作为 `Genos` 环境、权重、layer 选择的 sanity check
- 不做：不把结果写成当前 `cached-fusion` 主线的 gate，不把 RNA-seq 案例直接升级成新的 ATAC 主线

## 1. 当前策略

### 1.1 只跑最小闭环

第一轮只做：

- 模型：`Genos-1.2B`
- 数据集：`human_ocr_ensembl`
- 分类器：`MLP`
- 层：先 `layer 6` 和 `layer 12`
- 流程：`mini smoke -> 全量 6/12 层 -> 需要时再补全层扫描`

### 1.2 默认复用现有环境

优先复用当前已经验证过的本地模型目录：

```bash
MODEL_DIR=/data1/zhoujiazhen/bylw_atac/foundation_models/genos/Genos-1.2B
```

默认不重新下载模型，不新开一条和主线平行的重环境。

只有在现有环境缺依赖时，才补装 benchmark 所需的轻量依赖。

## 2. 执行前提

### 2.1 工作目录

建议单独放到一个旁线目录，不污染主仓库：

```bash
WORK_DIR=/data1/zhoujiazhen/bylw_atac/genos_ocr_sidecar
mkdir -p "$WORK_DIR"/{logs,datasets,results,embeddings}
```

### 2.2 官方 benchmark code

```bash
cd "$WORK_DIR"
if [ ! -d Genos ]; then
  git clone https://github.com/BGI-HangzhouAI/Genos.git
fi

cd Genos
git rev-parse HEAD > "$WORK_DIR"/genos_commit.txt
cd Technical_Notes/benchmarks-code
```

### 2.3 依赖

如果当前 `Genos` 环境里没有这些包，再补装：

```bash
pip install genomic-benchmarks==1.0.0 pandas scikit-learn matplotlib PyYAML xgboost
```

`benchmarks-code` 自己的 `requirements.txt` 很轻，可以直接对照：

```bash
pip install -r requirements.txt
```

如果沿用官方 benchmark code 里的 `torch.use_deterministic_algorithms(True)`，在 CUDA 10.2+ 上还要补：

```bash
export CUBLAS_WORKSPACE_CONFIG=:4096:8
```

## 3. 数据准备

### 3.1 下载 `human_ocr_ensembl`

`genomic-benchmarks` 的真实接口可以直接把 interval-list 数据转成 full-seq 目录结构：

```bash
python - <<'PY'
from genomic_benchmarks.loc2seq import download_dataset
download_dataset(
    "human_ocr_ensembl",
    version=0,
    dest_path="/data1/zhoujiazhen/bylw_atac/genos_ocr_sidecar/datasets/fullseq",
    use_cloud_cache=True,
)
PY
```

预期产物目录类似：

```text
.../datasets/fullseq/human_ocr_ensembl/
  train/
    class_a/
    class_b/
  test/
    class_a/
    class_b/
```

### 3.2 转成官方 benchmark code 需要的 JSONL

仓库内已补一个最小转换脚本：

- `scripts/prepare_genos_ocr_jsonl.py`

用法：

```bash
cd /home/zhengwei/project/python/TransChromBP

python scripts/prepare_genos_ocr_jsonl.py \
  --source-dir /data1/zhoujiazhen/bylw_atac/genos_ocr_sidecar/datasets/fullseq/human_ocr_ensembl \
  --out-dir /data1/zhoujiazhen/bylw_atac/genos_ocr_sidecar/datasets/jsonl/human_ocr_ensembl
```

脚本会：

- 扫描 `train/` 与 `test/` 下的类别子目录
- 以目录名排序建立稳定的 `label -> id` 映射
- 生成：
  - `train.jsonl`
  - `test.jsonl`
  - `label_map.json`
  - `summary.json`

如果后续确实需要强制指定哪个目录是正类，可以再给脚本加 `--positive-class`，但第一轮 OCR benchmark 不需要靠人工标签语义来跑通。

### 3.3 数据完整性检查

```bash
DATASET_JSONL_DIR=/data1/zhoujiazhen/bylw_atac/genos_ocr_sidecar/datasets/jsonl

wc -l \
  "$DATASET_JSONL_DIR"/human_ocr_ensembl/train.jsonl \
  "$DATASET_JSONL_DIR"/human_ocr_ensembl/test.jsonl

python - <<'PY'
import json
base = "/data1/zhoujiazhen/bylw_atac/genos_ocr_sidecar/datasets/jsonl/human_ocr_ensembl"
for split in ["train", "test"]:
    with open(f"{base}/{split}.jsonl") as f:
        rows = [json.loads(next(f)) for _ in range(3)]
    print(split, rows)
with open(f"{base}/label_map.json") as f:
    print("label_map:", json.load(f))
PY
```

## 4. 配置文件

### 4.1 mini smoke

先复制官方 config，生成一个最小版：

```bash
cd "$WORK_DIR"/Genos/Technical_Notes/benchmarks-code
cp config.yaml config_ocr_mini.yaml
```

`config_ocr_mini.yaml` 建议至少改成：

```yaml
model_path: "/data1/zhoujiazhen/bylw_atac/foundation_models/genos/Genos-1.2B"

datasets_feature_path: "benchmarks/datasets_info.yaml"
dataset_path: "/data1/zhoujiazhen/bylw_atac/genos_ocr_sidecar/datasets/jsonl_mini"
eval_datasets:
  - human_ocr_ensembl

embedding_output_dir: "/data1/zhoujiazhen/bylw_atac/genos_ocr_sidecar/embeddings_mini"
eval_result_path: "/data1/zhoujiazhen/bylw_atac/genos_ocr_sidecar/results_mini"

gpu_list: [1]
layer_to_eval: [6]
batch_size: 8
embedding_extract_split: 200000

classifer_type: "MLP"
pooled_embeddings_cat_dim: 1
mlp_dropout: 0.2
mlp_lr: 0.0001
mlp_epochs: 100
save_last_epoch_model: False

process_power_per_gpu: 5
embedding_process_power: 4
eval_process_power: 1

target_best_table:
  - "roc_auc"
  - ["roc_auc", "accuracy"]

target_layer_line_figure:
  - ["roc_auc", "accuracy"]
```

先从全量 JSONL 裁一个 mini 子集，但目录名仍保留 `human_ocr_ensembl`，避免额外改 `datasets_info.yaml`。不要直接 `head`，否则容易裁成单类样本；应按标签均衡抽样，例如：

```bash
MINI_ROOT=/data1/zhoujiazhen/bylw_atac/genos_ocr_sidecar/datasets/jsonl_mini
MINI_DIR="$MINI_ROOT"/human_ocr_ensembl
mkdir -p "$MINI_DIR"
python - <<'PY'
import json
import shutil
from pathlib import Path

src = Path("/data1/zhoujiazhen/bylw_atac/genos_ocr_sidecar/datasets/jsonl/human_ocr_ensembl")
dst = Path("/data1/zhoujiazhen/bylw_atac/genos_ocr_sidecar/datasets/jsonl_mini/human_ocr_ensembl")
limits = {"train": {0: 500, 1: 500}, "test": {0: 100, 1: 100}}

for split, split_limits in limits.items():
    counts = {label: 0 for label in split_limits}
    selected = []
    with (src / f"{split}.jsonl").open() as handle:
        for line in handle:
            item = json.loads(line)
            label = int(item["label"])
            if label in split_limits and counts[label] < split_limits[label]:
                selected.append(line.rstrip("\n"))
                counts[label] += 1
            if all(counts[label] >= split_limits[label] for label in split_limits):
                break
    (dst / f"{split}.jsonl").write_text("\n".join(selected) + "\n")

for name in ("label_map.json", "summary.json"):
    shutil.copy2(src / name, dst / name)
PY
```

### 4.2 全量 6/12 层

mini smoke 通过后，再建全量配置：

```yaml
model_path: "/data1/zhoujiazhen/bylw_atac/foundation_models/genos/Genos-1.2B"

datasets_feature_path: "benchmarks/datasets_info.yaml"
dataset_path: "/data1/zhoujiazhen/bylw_atac/genos_ocr_sidecar/datasets/jsonl"
eval_datasets:
  - human_ocr_ensembl

embedding_output_dir: "/data1/zhoujiazhen/bylw_atac/genos_ocr_sidecar/embeddings"
eval_result_path: "/data1/zhoujiazhen/bylw_atac/genos_ocr_sidecar/results"

gpu_list: [1]
layer_to_eval: [6, 12]
batch_size: 8
embedding_extract_split: 200000

classifer_type: "MLP"
pooled_embeddings_cat_dim: 1
mlp_dropout: 0.2
mlp_lr: 0.0001
mlp_epochs: 100
save_last_epoch_model: False

process_power_per_gpu: 5
embedding_process_power: 4
eval_process_power: 1

target_best_table:
  - "roc_auc"
  - ["roc_auc", "accuracy"]

target_layer_line_figure:
  - ["roc_auc", "accuracy"]
```

## 5. 运行顺序

### 5.1 mini smoke

```bash
cd "$WORK_DIR"/Genos/Technical_Notes/benchmarks-code
python benchmarks.py --config config_ocr_mini.yaml
```

通过标准：

- 能正常加载模型
- 能写出 mini `.pt`
- 能写出 mini `.tsv`
- `roc_auc` 明显高于随机

### 5.2 全量 6/12 层

按仓库的长任务约定，用后台运行，不要交互等：

```bash
cd "$WORK_DIR"/Genos/Technical_Notes/benchmarks-code
nohup python benchmarks.py --config config_ocr_full.yaml \
  > "$WORK_DIR"/logs/genos_ocr_full.log 2>&1 &
echo $! > "$WORK_DIR"/logs/genos_ocr_full.pid
```

查看进度：

```bash
tail -f "$WORK_DIR"/logs/genos_ocr_full.log
```

### 5.3 需要时再补全层扫描

只有在 `layer 6/12` 都跑干净后，才考虑把 `layer_to_eval` 改成 `0..12`。

## 6. 结果判读

### 6.1 这条旁线的成功标准

满足以下条件就够了：

- mini smoke 稳定
- 全量 `layer 6` 或 `layer 12` 的 `roc_auc` 与官方 `0.7569` 同量级
  - 实用上看，大致落在 `0.73-0.78`
- 结果没有明显崩坏
  - 例如 `< 0.65`
  - 或 train/test 维度、标签、embedding 文件不一致

### 6.2 这条旁线不能被怎么用

即便 OCR benchmark 跑得很好，也不能直接推出：

- `cached summary` 一定能提升 `TransChromBP`
- `online G1/G2` 的失败只是工程实现问题
- 现在就值得启动官方 README 里的 RNA-seq full fine-tuning 路线

### 6.3 这条旁线可以怎么用

如果结果正常，可以把它当成：

- `Genos-1.2B` embedding 链路的外部 sanity check
- `layer 6/12` 选择的额外旁证
- 后续 `cached-fusion` 复盘时的背景材料

## 7. 当前不做的事

本计划明确排除：

- `Genos-10B`
- OCR 结果作为 `cached-fusion` gate
- Phase B ATAC 轨迹回归
- 为 OCR 旁线抢占正在运行的主线资源

## 8. 全量 OCR 结束后的固定顺序

为了避免 OCR 旁线和 `cached-fusion` 主线在执行上互相挤占，后续动作固定为：

1. 等当前全量 `layer 6/12` OCR run 完成，并确认 `gpu1` 真正空闲。
2. 先收 OCR 结果，只回答两个问题：
   - `layer 6` 和 `layer 12` 哪个更好
   - 结果是否仍在可接受区间，不是明显崩坏
3. OCR 结果正常后，先跑 `run_genos_summary_probe.py`，不要直接进 `P0/P1/P2`。
4. `probe` 的用途只是判断 `genos_global_mean` 是否对 `count residual` 和 `peak/nonpeak complement` 有独立信号。
5. 只有当 `probe` 结果支持继续时，才启动 cached-fusion pilot。
6. `P0/P1/P2` 不要一把串行全开；在当前 `probe` 结果下，建议按 `P0 -> 看验证 -> P2 -> 看验证 -> P1` 的顺序分步发。
7. 如果 OCR 全量结果异常，或 `probe` 没看到明显增益，则停在复盘，不继续补全层扫描，也不启动 `P2/P1`。
8. 全过程都不改变当前 `cached-fusion` 主线优先级，OCR 旁线不作为主线 gate。

### 8.1 OCR 结束后的第一步

全量 OCR run 结束后，先检查：

- `results/Genos-1.2B/human_ocr_ensembl.tsv`
- `results/Genos-1.2B/reports/`
- `logs/genos_ocr_full_*.log`

若 `layer 6/12` 中至少一层仍落在“没有明显崩坏”的区间，再进入 `probe`。

### 8.2 Probe 阶段

`probe` 固定依赖：

- `G0 best.pt`
- `configs/model/v2fix_baseline.yaml`
- `configs/train/train_genos_cached_short10.yaml`
- `configs/data/data_tutorial_canonical_v1.yaml`
- `outputs/genos_cache/tutorial_canonical_v1/`

`probe` 输出应至少回答：

- `Probe A`: `genos_global_mean -> count residual` 是否有非零解释力
- `Probe B`: `AUC(genos_only)` 是否高于随机，`AUC(concat)` 是否高于 `AUC(encoded_only)`

### 8.3 Pilot 阶段

只有 `probe` 通过后，才进入：

1. `P0 baseline_short10`
2. `P2 global_count_only`
3. `P1 global_late_film`

其中：

- `P0` 是校准基线，不是可省略步骤
- 当前 `probe` 更支持先做 `P2`，因为 `genos_global_mean` 对 `count residual` 有小但稳定的解释力，而朴素分类拼接未见增益
- `P0/P2/P1` 后续统一改为 `2-rank DDP` matched recipe：`nproc_per_node=2`、`batch_size_per_gpu=10`、`global_batch=20`
- `P1/P2` 都要求复用同一份 `genos_cache`
- 如任何一步出现明显退化，应先停下来复盘，再决定是否继续下一步
