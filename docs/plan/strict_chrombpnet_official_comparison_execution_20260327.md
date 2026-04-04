# ChromBPNet 分层对比执行稿（2026-03-27，含三层叙事框架）

## 0. 三层对比框架

本轮对比不再用单一的”严格 vs 不严格”二分法，而是拆成三层，每层回答不同强度的问题：

| Layer | 名称 | 回答的问题 | 锁死的变量 | 未锁定的变量 | 论文中允许的表述 |
|---|---|---|---|---|---|
| **L1** | native pipeline compare | 完整系统谁更强？ | 原始 BAM / genome / fold | preprocessing、batch、optimizer、best rule 各自原生 | “system-level comparison” |
| **L2** | matched-budget system compare | 训练预算下谁更强？ | L1 + 硬件、global batch、seed、epoch 上限、外部 best rule、evaluator | peaks/nonpeaks 候选集合、optimizer/schedule/precision | “controlled comparison under matched hardware and training budget” |
| **L3** | shared-region compare | 能否更接近归因到架构/recipe？ | L2 + 统一的 peaks/nonpeaks/bigWig | optimizer/schedule/precision（仍属 recipe 差异） | “controlled comparison with shared training regions” |

### 核心原则

- **L1 是锚点**：证明”官方流程在我们机器上复现到什么程度”。不进公平主表。
- **L2 是第一轮主结果**：领域里大量工作都在这一层做比较，足以支撑 “system-level” claim。
- **L3 是 robustness check**：如果目标是”优势主要来自架构/recipe 而非训练区域差异”，就必须补这一层。
- optimizer / schedule / bf16 差异始终存在。即便到 L3，也只能写成 “shared-region system compare”，不能写成 “architecture-only ablation”。如果要做到纯架构归因，还需要额外一条 `noTF-controlled` 臂（使用完全相同的 optimizer/schedule/precision），这不在本轮范围内。

### 已实测的 region 差异（followup 2026-03-27 确认）

| | 官方臂 | TransChromBP 臂 | 差异 |
|---|---:|---:|---:|
| peaks | 269,800 | 267,175 | -2,625 |
| nonpeaks | 508,429 | 267,175 | **-241,254** |

nonpeaks 差了将近一半。这不是 fold 划分差异，是训练候选区域本身不同。因此当前工程件只能支撑 L2，要做 L3 必须先统一 region。

---

## 1. 这轮已经落地的资产

本轮已把 L2 = matched-budget system compare 的关键工程件补齐：

1. 官方 `ChromBPNet` 训练新增 `--save-all-checkpoints`
   - 文件：
     - `chrombpnet/parsers.py`
     - `chrombpnet/training/utils/argmanager.py`
     - `chrombpnet/training/train.py`
   - 作用：
     - 在保留 `best.h5` 的同时，额外输出 `*.epoch_001.h5` 这类逐 epoch checkpoint
     - 让官方侧可以走“外部统一选 best”

2. 官方 `ChromBPNet` 评估新增 `--metrics-only`，并补出 `mean_jsd`
   - 文件：
     - `chrombpnet/training/predict.py`
     - `scripts/paper_aligned_repro/summarize_metrics.py`
   - 作用：
     - 外部选 best 时不再为每个 epoch 额外写一大堆预测 HDF5 和图
     - 同时输出 `peak profile mean_jsd / median_jsd / count pearsonr`

3. 官方 tutorial L2 wrapper 已就位
   - 文件：
     - `scripts/paper_aligned_repro/run_tutorial_strict_compare_official.sh`
     - `scripts/paper_aligned_repro/select_best_epoch.py`
   - 作用：
     - 一条命令起 `tutorial official fidelity`
     - 一条命令起 `tutorial official controlled`
     - `controlled` 模式内部会显式启用 `--multi-gpu-train`，让官方训练真正占用 2×A6000，而不是单 fold 误退化成单卡
     - 一条命令按外部 metric 选 best epoch

4. `TransChromBP corrected-B` 的 L2 配置与 selector 已就位
   - 文件：
     - `vendor/transchrombp/transchrombp/configs/train/train_tutorial_corrected_b_strict_compare_6000.yaml`
     - `vendor/transchrombp/transchrombp/configs/train/train_gm12878_corrected_b_strict_compare_6000.yaml`
     - `vendor/transchrombp/scripts/select_best_epoch.py`
     - `scripts/deploy_strict_compare_staging_to_6000.sh`
   - 作用：
     - 双卡 A6000、`global_batch=32`、`max_epochs=50`、`checkpoint_every_epochs=1`
     - 关闭内部 early stop，改用外部统一选 best

5. 机器可读资产已生成
   - 运行矩阵：
     - `reports/assets/strict_chrombpnet_compare_run_matrix_20260327.csv`
   - `Dataset-X` shortlist：
     - `reports/assets/strict_compare_dataset_shortlist_20260327.csv`

---

## 2. 当前冻结的 protocol

### 2.1 L2 protocol（matched-budget system compare）

本轮先冻结以下口径，不再边跑边改：

| 项目 | 官方 controlled | TransChromBP controlled |
|---|---|---|
| 机器 | 6000 | 6000 |
| GPU | 2×A6000 | 2×A6000 |
| 第一轮 seed | `42` | `42` |
| 第一轮 fold | `fold_0` | tutorial 官方 split / 对齐 data config |
| global batch | `32` | `32` |
| epoch 上限 | `50` | `50` |
| checkpoint | 每 epoch 落盘 | 每 epoch 落盘 |
| best 选择规则 | 外部 `peak profile JSD median` | 外部 `peak profile JSD median` |
| 强制汇报 | `median_jsd` + `mean_jsd` + `count_r` | `median_jsd` + `mean_jsd` + `count_r` |
| peaks/nonpeaks | 官方 pipeline 自行生成 | TransChromBP pipeline 自行生成 |

- `official fidelity` 锚点线不进公平主表。
- L2 锁死的：原始输入、硬件、global batch、epoch 上限、外部 best rule、evaluator。
- L2 没有保证：最终进入 trainer 的 peaks/nonpeaks/bigWig 一致。因此只能表述为 `matched-budget system compare`。

### 2.2 L3 protocol（shared-region compare）— 待 L2 结果确认后启动

在 L2 的基础上额外锁死：

| 项目 | 官方 controlled | TransChromBP controlled |
|---|---|---|
| peaks | **同一份** `shared_filtered_peaks.bed` | **同一份** |
| nonpeaks | **同一份** `shared_filtered_nonpeaks.bed` | **同一份** |
| bigWig | 同一份 `merged_unstranded.bw`（如 TransChromBP 侧可直接使用） | 同一份 |
| 其余 | 与 L2 完全一致 | 与 L2 完全一致 |

Region 统一方案（二选一）：

1. **方案 A：TransChromBP 直接吃官方臂的产物**
   - 等 L2 的官方 controlled 跑完后，从其 `auxiliary/` 里取 `filtered.peaks.bed` 和 `candidate_nonpeaks.bed`
   - 让 TransChromBP 的 data config 指向这两份文件
   - 优点：零额外预处理；缺点：TransChromBP 训练时的 region 语义可能和之前的 pilot 不同

2. **方案 B：冻结一份共享 region 集合**
   - 取两边 peaks 的交集作为 `shared_filtered_peaks.bed`
   - 用官方 pipeline 基于交集 peaks 重新生成 `shared_nonpeaks.bed`
   - 两边都用这份共享文件重跑
   - 优点：对称、干净；缺点：需要额外生成一轮

推荐：先尝试方案 A（成本最低），如果方案 A 可行则直接用；如果 TransChromBP 的 data loader 不兼容官方格式，再退到方案 B。

---

## 3. Phase 1：tutorial 官方锚点 + Layer 2 双臂

### 3.1 部署官方 repo 到 6000

先把本地对 `chrombpnet/` 的补丁同步到 6000 上的 `chromBPNet` 仓库：

```bash
bash scripts/sync_project.sh deploy
```

### 3.2 起 tutorial 官方 fidelity

```bash
bash scripts/paper_aligned_repro/run_tutorial_strict_compare_official.sh \
  --mode fidelity \
  --work-root /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare \
  --seed 42 \
  --folds "fold_0" \
  --gpus 0,1
```

### 3.3 起 tutorial 官方 controlled

```bash
bash scripts/paper_aligned_repro/run_tutorial_strict_compare_official.sh \
  --mode controlled \
  --work-root /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare \
  --seed 42 \
  --folds "fold_0" \
  --gpus 0,1
```

### 3.4 用统一规则给官方 controlled 选 best

等 `tutorial_official_controlled_s42` 跑完后执行：

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
  --metric profile_metrics.peaks.median_jsd \
  --mode min
```

补充执行规则（2026-03-28 更新）：

- 这类 `selector` 本质上是“逐 checkpoint 的前向评估”，不是训练。
- 2026-03-30 起，本地 `scripts/paper_aligned_repro/select_best_epoch.py` 已支持 `--gpus 0,1` 做 checkpoint 级分片：不同 checkpoint 会被分配到不同 GPU，各自独立调用 `predict.py`，指标口径不变，只缩短 wall-clock。
- 这不等于“单个 checkpoint 用 MirroredStrategy 双卡推理”；单个 checkpoint 的 held-out test 仍然更适合单卡，或与另一条独立 test/selector 分卡并行。
- 当前更合理的调度规则变成：
  - 单条 multi-checkpoint selector：优先传 `--gpus 0,1`
  - 两条独立 selector / held-out test：优先一条绑 `GPU0`、另一条绑 `GPU1`
  - 不要为了评一个单 checkpoint 去强行做多卡同步推理

### 3.5 部署 matched-budget 专用 `TransChromBP` staged 文件到 6000

这一步只同步 `tmp_remote_edit/` 里本轮新增的 matched-budget 配置和 selector：

```bash
bash scripts/deploy_strict_compare_staging_to_6000.sh
```

### 3.6 起 tutorial `TransChromBP corrected-B controlled`

在 6000 上执行：

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 '
export CHROMBPNET_ENV=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet
export PATH="$CHROMBPNET_ENV/bin:$PATH"
export LD_LIBRARY_PATH="$CHROMBPNET_ENV/lib:$LD_LIBRARY_PATH"
export PYTHONPATH=/data1/zhoujiazhen/bylw_atac/TransChromBP/src:$PYTHONPATH
cd /data1/zhoujiazhen/bylw_atac/TransChromBP
NPROC_PER_NODE=2 MASTER_PORT=29611 \
RUN_NAME=tutorial_corrected_b_strict_compare_s42 \
LOG_FILE=/data1/zhoujiazhen/bylw_atac/logs/tutorial_corrected_b_strict_compare_s42.log \
bash scripts/launch_v2fix_6000_single_dataset.sh \
  /data1/zhoujiazhen/bylw_atac/TransChromBP/configs/model/transchrombp_teacher_v2_center_pool.yaml \
  /data1/zhoujiazhen/bylw_atac/TransChromBP/configs/train/train_tutorial_corrected_b_strict_compare_6000.yaml
'
```

### 3.7 用统一规则给 `TransChromBP controlled` 选 best

等上面的 run 跑完后，在 6000 上执行：

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 '
export CHROMBPNET_ENV=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet
export PATH="$CHROMBPNET_ENV/bin:$PATH"
export LD_LIBRARY_PATH="$CHROMBPNET_ENV/lib:$LD_LIBRARY_PATH"
export PYTHONPATH=/data1/zhoujiazhen/bylw_atac/TransChromBP/src:$PYTHONPATH
cd /data1/zhoujiazhen/bylw_atac/TransChromBP
python scripts/select_best_epoch.py \
  --checkpoint-glob "/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/checkpoints/tutorial_corrected_b_strict_compare_s42/epoch_*.pt" \
  --split valid \
  --metric results.peak.profile_target_jsd_full_median \
  --mode min \
  --output-dir /data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/tutorial_corrected_b_strict_compare_s42_external_best
'
```

补充执行规则（2026-03-28 更新）：

- `TransChromBP` selector 当前接口是单 `--device`，本质也是单 checkpoint 的评估前向。
- 如果 future run 里 official / ours 都需要外部统一选 best，优先考虑把两条 selector 分配到不同 GPU 并行，而不是继续串行。
- 只有在后续确认 selector wall-clock 成为稳定瓶颈时，才值得再做代码层的评估并行改造。

### 3.8 当前结果快照（2026-03-28 18:35 CST）

6000 上 tutorial `Layer 2` 的统一 `valid` selector 队列 `strict_compare_tutorial_selector_queue_20260328_1452` 已于 `2026-03-28 16:49:20 CST` 收口，official / `TransChromBP corrected-B` 两侧都已写出 `epoch_metrics.csv` 与 `best_epoch.json`。

| Arm | best epoch | peak median_jsd(valid) | peak mean_jsd(valid) | peak count_r(valid) | 备注 |
|---|---:|---:|---:|---:|---|
| official controlled | 35 | 0.35783 | 0.35234 | 0.65805 | selector metric=`profile_metrics.peaks.median_jsd` |
| TransChromBP corrected-B | 44 | 0.33538 | 0.32944 | 0.08428 | selector metric=`results.peak.profile_target_jsd_full_median` |

`TransChromBP corrected-B` 的 `peak count_r(valid)` 轨迹进一步显示，问题更像 late-epoch instability，而不是单纯“整条 run 只有坏结果”：

| epoch 区间 | `peak count_r(valid)` 范围 | 判读 |
|---|---:|---|
| 1-28 | `0.73069 ~ 0.81279` | 稳定、与 official 同档甚至更高 |
| 29-31 | `0.40757 ~ 0.70319` | 开始失稳 |
| 32-50 | `0.04403 ~ 0.71969` | 明显震荡，profile-only selector 容易选到坏点 |

同时，同一条 run 中仍存在 profile 几乎不差、但 count 正常的 checkpoint：

| checkpoint | peak median_jsd(valid) | peak mean_jsd(valid) | peak count_r(valid) | 备注 |
|---|---:|---:|---:|---|
| epoch 35 | 0.33641 | 0.33061 | 0.69681 | 与 official best epoch 同编号，profile 仅比 epoch 44 差 `0.00103` |
| epoch 44 | 0.33538 | 0.32944 | 0.08428 | 当前 profile-only selector 选中的 best |

当前 interpretation：

- 在统一 `valid` selector 口径下，`TransChromBP corrected-B` 的 peak profile JSD 明显优于 official。
- 但当前 profile-only selector 最终落在 count 已塌的 late checkpoint，因此当前不能把 tutorial `Layer 2` 结果直接写成“整体优于 official”。
- 更稳妥的阶段结论应写成：profile 侧已有正向信号，但训练后期存在明显的 count instability；当前主风险是 checkpoint selection 与 late-epoch stability，而不是“模型从头到尾都没有可用 checkpoint”。
- 后续优先级应是：先把 `epoch 35` / `epoch 44` 这组对照写入统一 follow-up，再决定是补独立 test / 调整 selector，还是直接进入 `Layer 3 shared-region compare` / `GM12878 Layer 2` 做进一步归因。

---

## 4. Phase 2：GM12878 Layer 2 双臂

### 4.1 GM12878 当前本地资产状态

6000 上已确认可直接进入当前 `Layer 2`：

- 官方训练输入：
  - `/data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/GM12878/merged.bam`
  - `/data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/GM12878/overlap.bed.gz`
- `TransChromBP` data config 所需输入：
  - `/data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/GM12878/prep_v1/inputs/overlap.bed`
  - `/data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/GM12878/prep_v1/background/tutorial_folds/nonpeaks_tutorial_folds.bed`
  - `/data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/GM12878/prep_v1/bigwig/merged_unstranded.bw`

### 4.2 起 GM12878 官方 controlled

```bash
bash scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh \
  --name gm12878_official_controlled_s42 \
  --genome /data1/zhoujiazhen/bylw_atac/hg38.fa \
  --chrom-sizes /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes \
  --bam /data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/GM12878/merged.bam \
  --peaks /data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/GM12878/overlap.bed.gz \
  --blacklist /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/blacklist.bed.gz \
  --fold-dir /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/folds \
  --work-root /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare \
  --seed 42 \
  --gpus 0,1 \
  --max-parallel 2 \
  --batch-size 32 \
  --predict-batch-size 512 \
  --epochs 50 \
  --early-stop 50 \
  --folds "fold_0" \
  --save-all-checkpoints
```

### 4.3 给 GM12878 官方 controlled 选 best

```bash
python scripts/paper_aligned_repro/select_best_epoch.py \
  --model-glob '/data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/runs/gm12878_official_controlled_s42/fold_0/seed_42/chrombpnet/models/chrombpnet.epoch_*.h5' \
  --genome /data1/zhoujiazhen/bylw_atac/hg38.fa \
  --bigwig /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/runs/gm12878_official_controlled_s42/fold_0/seed_42/chrombpnet/auxiliary/data_unstranded.bw \
  --peaks /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/runs/gm12878_official_controlled_s42/fold_0/seed_42/chrombpnet/auxiliary/filtered.peaks.bed \
  --nonpeaks None \
  --fold-json /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/folds/fold_0.json \
  --output-dir /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/runs/gm12878_official_controlled_s42/fold_0/seed_42/chrombpnet/external_best \
  --split valid \
  --metric profile_metrics.peaks.median_jsd \
  --mode min
```

### 4.4 起 GM12878 `TransChromBP corrected-B controlled`

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1 '
export CHROMBPNET_ENV=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet
export PATH="$CHROMBPNET_ENV/bin:$PATH"
export LD_LIBRARY_PATH="$CHROMBPNET_ENV/lib:$LD_LIBRARY_PATH"
export PYTHONPATH=/data1/zhoujiazhen/bylw_atac/TransChromBP/src:$PYTHONPATH
cd /data1/zhoujiazhen/bylw_atac/TransChromBP
NPROC_PER_NODE=2 MASTER_PORT=29621 \
RUN_NAME=gm12878_corrected_b_strict_compare_s42 \
LOG_FILE=/data1/zhoujiazhen/bylw_atac/logs/gm12878_corrected_b_strict_compare_s42.log \
bash scripts/launch_v2fix_6000_single_dataset.sh \
  /data1/zhoujiazhen/bylw_atac/TransChromBP/configs/model/transchrombp_teacher_v2_center_pool.yaml \
  /data1/zhoujiazhen/bylw_atac/TransChromBP/configs/train/train_gm12878_corrected_b_strict_compare_6000.yaml
'
```

选 best 的命令与 tutorial 相同，只需把 `run_name` 改成 `gm12878_corrected_b_strict_compare_s42`。

---

## 5. Phase 3：tutorial L3 shared-region 双臂（待 L2 结果确认后启动）

### 5.1 启动条件

- tutorial L2 双臂已完整闭环（fidelity + controlled + select best + 指标汇总）
- L2 结果显示 TransChromBP 有正向增益，值得进一步归因
- 如果 L2 结果已经显示没有实质差异，则 L3 不必启动

### 5.2 Region 统一（方案 A 路径）

从 L2 官方 controlled 的产物中提取共享 region：

```bash
# 在 6000 上执行
OFFICIAL_AUX=/data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/runs/tutorial_official_controlled_s42/fold_0/seed_42/chrombpnet/auxiliary
SHARED_DIR=/data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/shared_regions/tutorial

mkdir -p $SHARED_DIR
cp $OFFICIAL_AUX/filtered.peaks.bed $SHARED_DIR/shared_filtered_peaks.bed
cp $OFFICIAL_AUX/filtered.nonpeaks.bed $SHARED_DIR/shared_filtered_nonpeaks.bed
cp $OFFICIAL_AUX/data_unstranded.bw $SHARED_DIR/shared_unstranded.bw

# 验证
wc -l $SHARED_DIR/shared_filtered_peaks.bed
wc -l $SHARED_DIR/shared_filtered_nonpeaks.bed
```

### 5.3 起 tutorial 官方 controlled L3

```bash
bash scripts/paper_aligned_repro/run_tutorial_strict_compare_official.sh \
  --mode controlled \
  --work-root /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare \
  --run-suffix "_L3" \
  --seed 42 \
  --folds "fold_0" \
  --gpus 0,1 \
  --shared-peaks $SHARED_DIR/shared_filtered_peaks.bed \
  --shared-nonpeaks $SHARED_DIR/shared_filtered_nonpeaks.bed \
  --shared-bigwig $SHARED_DIR/shared_unstranded.bw
```

更新（2026-03-29，Codex）：
- `run_tutorial_strict_compare_official.sh` 已补 `--run-suffix` / `--shared-peaks` / `--shared-nonpeaks` / `--shared-bigwig`
- `run_paper_aligned_fast_1seed.sh` 已支持在 shared-region 模式下跳过内部 `prep nonpeaks`，并在评估阶段复用输入 peaks 与共享 bigWig
- 因此官方 L3 的剩余工作不再是“补接口”，而是“实际起 run 并收口 selector / held-out test”

### 5.4 准备 TransChromBP L3 data config

需要新建一份 data config，将 peaks / nonpeaks 指向 `$SHARED_DIR/` 下的共享文件：

```yaml
# train_tutorial_corrected_b_strict_compare_L3_6000.yaml
# 与 L2 配置完全一致，只改 peaks / nonpeaks 路径
peaks_bed: /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/shared_regions/tutorial/shared_filtered_peaks.bed
nonpeaks_bed: /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/shared_regions/tutorial/shared_filtered_nonpeaks.bed
```

### 5.5 起 tutorial TransChromBP controlled L3

与 Phase 1 的 3.6 相同，只是把 config 换成 L3 版本，`RUN_NAME` 改成 `tutorial_corrected_b_strict_compare_L3_s42`。

### 5.6 L3 选 best + 汇总

与 L2 完全相同的 selector 和 evaluator。

---

## 6. Dataset-X 当前 gate

本轮不再空谈 `Dataset-X`，当前 gate 已经写死：

1. 必须 `非 tutorial/K562`
2. 必须 `非 GM12878`
3. 必须有完整 `BAM + peaks`
4. 必须能产出 `merged_unstranded.bw + nonpeaks`
5. 最终才允许进入 controlled 双臂

当前本地 shortlist 见：

- `reports/assets/strict_compare_dataset_shortlist_20260327.csv`

基于当前本地文件状态，最值得推进的 `Dataset-X` 候选不是独立 K562，而是：

- `HeLa-S3 (SRP166944)`：
  - 已有 4 个 SRA payload
  - 但还没有 BAM、peaks、bigWig
  - 因此现在是 `P1 local candidate`，不是可立即开跑的数据线

也就是说，`Dataset-X` 这条线目前仍然 **blocked on dataset freeze**，不允许抢到 `GM12878` 前面。

---

## 7. 完整执行顺序

### 第一阶段：L1 + L2（当前可直接启动）

| Step | 动作 | 产出 | Layer |
|---:|---|---|---|
| 1 | `bash scripts/sync_project.sh deploy` | 6000 代码同步 | — |
| 2 | 跑 `tutorial official fidelity` | 官方锚点指标 | L1 |
| 3 | 跑 `tutorial official controlled` (L2) | 官方 L2 逐 epoch checkpoint | L2 |
| 4 | `bash scripts/deploy_strict_compare_staging_to_6000.sh` | TransChromBP 配置同步 | — |
| 5 | 跑 `tutorial corrected-B controlled` (L2) | TransChromBP L2 逐 epoch checkpoint | L2 |
| 6 | 两侧 selector 按 `median_jsd` 选 best | L2 双臂主指标 | L2 |
| 7 | 汇总 L2 结果，判断是否有正向增益 | **L2 go/no-go 决策** | — |

只有 tutorial L2 流程完整闭环后，才平移到 GM12878。

### 第二阶段：L3 shared-region（L2 有正向增益后启动）

| Step | 动作 | 产出 | Layer |
|---:|---|---|---|
| 8 | 从 L2 官方 controlled 的 `auxiliary/` 提取 shared region | `shared_filtered_peaks.bed` + `shared_filtered_nonpeaks.bed` | — |
| 9 | 跑 `tutorial official controlled L3`（吃 shared region） | 官方 L3 逐 epoch checkpoint | L3 |
| 10 | 准备 TransChromBP L3 data config，跑 `tutorial corrected-B controlled L3` | TransChromBP L3 逐 epoch checkpoint | L3 |
| 11 | 两侧 selector 选 best，汇总 L3 结果 | L3 双臂主指标 | L3 |
| 12 | 比较 L2 vs L3 结果，判断 region 差异是否影响结论方向 | **L3 归因判断** | — |

### 第三阶段：扩展（L2/L3 tutorial 闭环后）

| Step | 动作 | 条件 |
|---:|---|---|
| 13 | GM12878 L2 双臂 | tutorial L2 已闭环 |
| 14 | GM12878 L3 双臂（如需要） | tutorial L3 已闭环且 region 差异确实影响结论 |
| 15 | 3 seed 扩展 | 单 seed 方向一致 |
| 16 | Dataset-X | GM12878 已闭环 + Dataset-X 预处理就绪 |

### 论文主表落点建议

- **如果 L2 和 L3 结论方向一致**：主表放 L2（更自然、工程成本低），L3 作为 supplementary robustness check
- **如果 L2 有增益但 L3 缩小/消失**：主表必须放 L3，并在正文讨论 region 差异对结果的影响
- **如果 L2 就没有实质增益**：停止扩展，重新审视模型价值

无论哪种情况，论文中不允许出现以下表述：
- "strictly controlled variable comparison" — 除非 optimizer/schedule/precision 也统一了
- "architecture-only attribution" — 除非补了 `noTF-controlled`（同 optimizer/schedule/precision）
- "same frozen dataset" 指代 region 不同的 L2 — 必须加 "shared raw inputs, native preprocessing" 限定词
