# Genos Cached Fusion Final Execution Plan

更新日期：2026-03-23

本文档用于合并并收口以下三份材料：

- `G0/G1/G2` 的真实运行代码与当前训练口径
- Claude 的实现型方案：`genos_cached_fusion_experiment_design_20260323.md`
- 前一版策略型约束文档：`genos_cached_fusion_refined_experiment_design_20260323.md`

目标不是继续为当前 `online full-seq + pre-transformer gated add` recipe 补证据，而是用一轮低成本、短反馈周期、语义清楚的实验回答：

> `Genos-1.2B` 在当前 `TransChromBP` 主线上，是否还值得继续投入；如果值得，最先应该保留的是哪一种信息和哪一种融合方式？

---

## 1. 已确认的事实

### 1.1 当前 `G0/G1/G2` 的真实代码路径

入口：

- launcher：`vendor/transchrombp/transchrombp/scripts/run_genos_pilot.sh`
- train config：`vendor/transchrombp/transchrombp/configs/train/train_v2fix_genos_profile_select.yaml`
- model config：
  - `G0`：`v2fix_baseline.yaml`
  - `G1`：`v2fix_genos_gate.yaml`
  - `G2`：`v2fix_genos_mean.yaml`

当前 pilot 口径：

- `batch_size_per_gpu = 20`
- `max_epochs = 20`
- `validate_every_epochs = 5`
- `best_metric = peak.profile_target_jsd_full_mean`
- 单卡 run，不用 2-rank DDP 扩大全局 batch

### 1.2 当前 `Genos` 的接法

训练和验证时，每个 step 都会：

1. 从 dataset 取 `seq`
2. 在 `train_ddp.py` 调用 `genos_runtime.extract(seq)`
3. `GenosFeatureExtractor` 做：
   - one-hot -> DNA string
   - tokenizer
   - Genos 正向
   - reverse-complement 再跑一次
   - hidden state 平均
4. 得到 `genos_feat [B, L, 1024]`
5. 在 `local_tower` 之后、Transformer 之前做一次 gated residual add

### 1.3 `G1` 与 `G2` 的差异

- `G1`
  - `pool_mode = none`
  - 用逐位置 `Genos` 特征
- `G2`
  - `pool_mode = mean`
  - 先做全局均值，再广播到每个位置

当前 `G2 >= G1` 的含义应解释为：

- 目前更像是 region-level / global summary signal 在起作用
- 不是“逐位置 `Genos` 已被成功利用”

### 1.4 数据语义约束

真实数据集实现：`vendor/transchrombp/transchrombp/data/real_data.py`

当前训练语义：

- train split：
  - `peak_max_jitter = 500`
  - `random_revcomp = true`
  - `revcomp_prob = 0.5`
  - nonpeak 通过 `PeakNonpeakResamplingSampler` 每个 epoch 重采样子集
- valid split：
  - `jitter = 0`
  - `random_revcomp = false`

因此，缓存系统必须区分两类语义：

- `valid exact cache`
- `train canonical cache`

不能把训练缓存表述成“完全等价于当前增强语义”。

### 1.5 一个必须补写清楚的规模约束

`train_ddp.py` 会打印两种计数：

- `train_regions = len(train_ds)`
- `train_epoch_regions = train_sampler.epoch_size`

两者并不总是相等。

在当前实现里，train split 为了支持 nonpeak 重采样，dataset 可能保留“完整 nonpeak 池”，每个 epoch 再抽一个子集。因此：

- cache 的真实基数由 `train_regions` 决定
- 训练时每个 epoch 实际消费的样本量由 `train_epoch_regions` 决定

这直接影响：

- cache 存储预算
- 是否适合“顺手把 `bins4` 也一起存了”

所以最终方案不能默认套用“236k records”这一数字，必须在 cache 构建前先打印并记录：

- `train_regions`
- `train_epoch_regions`
- `valid_regions`

---

## 2. 合并后的最终设计决策

### 决策 1：保留 `valid exact cache` 与 `train canonical cache` 的命名区分

这是必要的，不只是表述清晰。

- valid：缓存与验证输入完全一致
- train：缓存只是 canonical center 的 record-level side information

这能避免后面把 negative result 错误解释成“cached summary 本身无效”，而忽略“训练语义与 cache 语义不完全对齐”。

### 决策 2：第一轮实验只消费 `global_mean`

第一轮 pilot 只使用：

- `global_mean [1024]`

不直接消费：

- `bins4_mean [4, 1024]`

原因：

- 当前唯一可信信号是 `G2(mean) >= G1(gate)`
- 第一轮的目标是先判断 `global summary` 本身是否值得继续
- 不应同时扩大融合复杂度和 cache 复杂度

### 决策 3：预计算阶段支持同时写出 `global_mean` 与 `bins4_mean`，但受预算门控

这是对 Claude 建议的吸收，但加上一个必须的预算条件。

原因：

- 一旦 full hidden state 已在 GPU 上，额外求 `bins4_mean` 的计算边际开销几乎为零
- 但它的存储成本不是固定“约 1.5 GB”，而取决于真实 `train_regions`

因此最终规则是：

- cache builder 支持一次前向同时导出：
  - `global_mean`
  - `bins4_mean`
- 但是否落盘 `bins4_mean`，由估算大小决定

建议默认门槛：

- 若 `estimated_bins4_bytes <= 4 GiB`，第一轮顺手存出 `bins4_mean`
- 否则只存 `global_mean`，把 `bins4` 推迟到第二轮条件触发

换句话说：

- 实验面第一轮仍然是 `global-only`
- 预计算面则允许“预算内顺手存 bins4”

### 决策 4：主 probe 改为 `count residual ridge`

主 probe：

- `count residual ridge probe`

次 probe：

- `complement logistic probe`

这样 probe 才更贴近当前下游任务，而不只是复刻 Phase 0 的可分性判断。

### 决策 5：缓存后端用 `memmap + manifest.json`

采用：

- `np.memmap` / `np.load(..., mmap_mode="r")`
- `manifest.json`
- `record_sha1`

原因：

- 更适合 DataLoader 多 worker 只读访问
- 更容易做数据配置与缓存对齐校验
- 可以显式打印和追踪 storage 预算

### 决策 6：把 ROI 写成硬性判据

任何 cached summary 方案都必须满足：

- 相对 `P0` 的 step-time overhead `<= 10%`

如果超过这个量级，说明：

- 实现有问题
- 或 cached path 仍然引入了不可接受的工程复杂度

---

## 3. 最终 cache 方案

### 3.1 第一轮必须有的 cache

必须构建：

- `train_global_mean.f16.npy`
- `valid_global_mean.f16.npy`
- `manifest_train.json`
- `manifest_valid.json`

### 3.2 第一轮可选顺手落盘的 cache

在预算允许时额外构建：

- `train_bins4_mean.f16.npy`
- `valid_bins4_mean.f16.npy`

注意：

- 第一轮 pilot 不读取这组文件
- 这只是为第二轮可能的 `P3` 避免重跑 `Genos` 提取

### 3.3 manifest 最少字段

每个 split 的 `manifest.json` 至少记录：

- `data_config_path`
- `split`
- `input_len`
- `n_records`
- `train_epoch_regions`（train split 需要）
- `genos_model_path`
- `layer`
- `bidirectional`
- `features`
- `dtype`
- `record_sha1`

其中：

- `record_sha1` 对 `[(chrom, center, source), ...]` 的顺序列表做哈希
- train split 还要额外记录 `supports_epoch_resampling`

### 3.4 cache builder 的预检查输出

在真正提特征前，脚本先打印：

- `train_regions`
- `train_epoch_regions`
- `valid_regions`
- `global_mean` 预计大小
- `bins4_mean` 预计大小
- 当前输出目录剩余空间

只有预算通过后才开始正式提特征。

---

## 4. probe 设计

### 4.1 Probe A：count residual ridge

目的：

- 判断 `Genos global summary` 是否解释了 baseline 在 count 维度上的残差

输入：

- `G0 best.pt`
- validation split
- `genos_global`

步骤：

1. 用 `G0 best.pt` 在 validation 上前向
2. 记录：
   - `pred_logcount_full`
   - `true_logcount`
   - `genos_global`
3. 构造：
   - `residual = true_logcount - pred_logcount_full`
4. 用 `genos_global [1024]` 做 ridge regression 预测 residual
5. 在 validation 上报告：
   - `R^2`
   - `Pearson r`

判读：

- 若 `R^2 ≈ 0`，说明 `global summary` 很难补 baseline 的 count 残差
- 若 `R^2 > 0` 且稳定，说明 `P2` 至少有继续试的价值

### 4.2 Probe B：complement logistic

目的：

- 判断 `Genos` 是否在已有 `TransChromBP` 表征上提供额外 region-level 可分信息

特征：

- `A`：`G0` pooled encoded
- `B`：`genos_global`
- `C`：`[G0 pooled encoded; genos_global]`

标签：

- `peak = 1`
- `nonpeak = 0`

报告：

- `AUC(A)`
- `AUC(B)`
- `AUC(C)`

判读：

- 若 `AUC(C) ≈ AUC(A)`，说明 `Genos` 的补充信息很弱
- 若 `AUC(C) > AUC(A)`，说明至少存在 summary-level 互补性

Probe B 是辅证，不替代 Probe A。

---

## 5. 第一轮 pilot 设计

统一设置：

- `max_epochs = 10`
- `validate_every_epochs = 1`
- `checkpoint_every_epochs = 1`
- `early_stop_patience = 3`
- pilot 默认改为 `2-rank DDP`
- `batch_size_per_gpu = 10`
- `global_batch = 20`
- `num_workers = 4`（每个 rank）
- 在当前 `profile_bias_stop_gradient` 修补后，双卡默认 `find_unused_parameters = false`
- `P0/P2/P1` 必须复用同一套 `nproc_per_node / batch_size_per_gpu / lr` 口径，避免把训练拓扑差异混进实验对比

### `P0 baseline_short10`

- 不用 `Genos`
- 作用：matched control

### `P1 global_late_film`

输入：

- `genos_global [B, 1024]`

融合位置：

- Transformer 输出之后
- profile/count head 之前

融合方式：

- 零初始化的 `FiLM`

目标：

- 检查 global summary 是否能同时改善 profile 和 count

### `P2 global_count_only`

输入：

- `genos_global [B, 1024]`

融合位置：

- 只注入 `count_signal_head` 中间层

目标：

- 诊断 `Genos` 是否只对 count / occupancy 有帮助

### 第一轮不运行

- `bins4_late_film`
- `offset-aware cache`
- `cross-attention`
- `LoRA`

---

## 6. 第二轮条件触发项

### `P3 bins4_late_film`

前提：

- `P1` 有明确正信号
- 且预算内已顺手存了 `bins4_mean`，或愿意补跑一轮 cache 构建

输入：

- `genos_bins4 [B, 4, 1024]`

目标：

- 判断 coarse positional summary 是否优于纯 global summary

### `P4 offset_aware_global`

前提：

- `P1` 有轻微正信号但波动较大
- 怀疑 `train canonical cache` 与增强语义偏差过大

做法：

- 为若干离散 offset 构建 `global_mean`
- dataset 按本 step 的 jitter 选择最近 offset 的 summary

目标：

- 判断 train-time augmentation mismatch 是否压制了 summary 收益

---

## 7. 模型与训练代码改动

### 7.1 新增文件

- `scripts/build_genos_summary_cache.py`
- `scripts/run_genos_summary_probe.py`
- `scripts/run_genos_cached_pilot.sh`
- `configs/model/v2fix_genos_global_late_film.yaml`
- `configs/model/v2fix_genos_global_count.yaml`
- `configs/train/train_genos_cached_short10.yaml`

### 7.2 修改文件

- `vendor/transchrombp/transchrombp/models/genos_adapter.py`
- `vendor/transchrombp/transchrombp/models/transchrombp.py`
- `vendor/transchrombp/transchrombp/data/real_data.py`
- `vendor/transchrombp/transchrombp/training/train_ddp.py`

### 7.3 `genos_adapter.py`

新增：

- `GenosSummaryFiLM`

接口：

- 输入：`encoded [B, L, 256]`, `genos_summary [B, 1024]`
- 输出：调制后的 `encoded [B, L, 256]`

初始化原则：

- `gamma` 初始接近 `1`
- `beta` 初始接近 `0`

这样训练起点与 baseline 对齐。

### 7.4 `transchrombp.py`

新增 forward 参数：

- `genos_summary: Optional[Tensor] = None`
- `genos_bins4: Optional[Tensor] = None`

第一轮只使用 `genos_summary`。

`P1` 路径：

```python
encoded = self.transformer(tokens)
if genos_summary is not None and self.genos_film is not None:
    encoded = self.genos_film(encoded, genos_summary)
```

`P2` 路径：

- 将 `count_signal_head` 从 `Sequential` 拆成显式层
- 在 `fc1` 后加入 `genos_count_proj(genos_summary)`

### 7.5 `real_data.py`

新增 dataset 参数：

- `genos_cache_dir: str = ""`
- `genos_cache_features: Sequence[str] = ()`

行为：

- init 时读取 `manifest.json`
- 校验：
  - `record_sha1`
  - `input_len`
  - `split`
  - `n_records`
- 用 `np.load(..., mmap_mode="r")` 打开 cache

`__getitem__` 返回：

- `genos_global`
- 若存在且启用，再返回 `genos_bins4`

注意：

- train split 仍然先按当前语义生成 `seq` 和 `profile_counts`
- cached summary 只是附加 side input，不替代原增强

### 7.6 `train_ddp.py`

增加两个能力：

1. 从 batch 读取 cached summary
2. 记录 cached path 的性能与使用统计

第一轮必须新增的统计：

- `step_time_total`
- `step_time_data_to_device`
- `step_time_forward_backward`
- `step_time_genos_cache_fetch`
- `film_gamma_mean/std` 或 `count_proj_norm`

需要显式防止：

- cached 模式意外 fallback 到在线 `GenosFeatureExtractor`

若 `cached_summary` 模式下仍构建了 `genos_runtime`，应直接报错，而不是静默回退。

---

## 8. 下一批数据管线优化候选

以下优化方向有工程价值，但不纳入当前正在跑的 `P0`，也不直接混入本轮 `P0/P2/P1` 对比：

- 预计算 `FASTA + bigWig` 数据管线缓存
- 提高 DataLoader worker / prefetch 配置
- 把 sequence one-hot 延后到 GPU 端完成

这样做的原因不是这些优化“没用”，而是它们会改变当前 real-data backend 的实现路径。若在 `P0` 已经开跑后再引入，容易把“训练拓扑差异”和“数据管线差异”混进同一轮实验。

### 8.1 优先级判断

优先级建议：

1. `FASTA + bigWig` 预计算缓存
2. `num_workers / prefetch_factor / persistent_workers`
3. GPU 端 one-hot

当前真实瓶颈更像是 `ChromBPNetBigWigDataset.__getitem__` 里的在线随机读取，而不是 one-hot 本身：

- `_fetch_onehot()` 里按窗口随机读 `FASTA`
- `_fetch_profile()` 里按窗口随机读 `bigWig`

因此最值得做的是减少在线随机 I/O，而不是优先重写 one-hot。

### 8.2 语义边界

下一批若启用上述优化，必须把它们视为“新数据管线”：

- 当前这条 `P0` 只作为现有 backend 下的 matched control
- 若未来切换到预计算 cache backend 或 GPU one-hot backend，应重跑一条新的 `P0`
- 后续 `P2/P1` 必须与该新 `P0` 共享同一套数据管线口径

换句话说：

- `DDP / batch / workers` 这类训练基础设施修正，不构成实验语义改变
- `FASTA/bigWig` cache、GPU one-hot 这类数据路径重写，应视为新一批实验

### 8.3 建议的落地顺序

若下一批正式尝试，建议顺序为：

1. 先补更准确的 loader profiling，把 `next(dataloader)`、H2D copy、forward/backward 分开记
2. 先实现 `bigWig profile cache`，再考虑 `sequence index cache`
3. 在 cache backend 稳定后，再评估 GPU 端 one-hot 是否仍值得加

缓存介质优先建议：

- 首选 `numpy memmap`
- 次选 `zarr`
- 暂不把单文件 `HDF5` 作为第一实现

原因是当前场景是多 worker、DDP、只读随机访问，`memmap` 更直接，也更接近现有 `genos_cache` 的使用方式。

---

## 9. cache builder 伪代码

```python
def build_cache(split: str, features: list[str]) -> None:
    ds = build_dataset(split=split, disable_train_augmentation=True_for_cache_view)
    records = ds.records

    n_records = len(records)
    epoch_regions = getattr(ds, "target_nonpeak_count", None)
    print_dataset_summary(n_records, epoch_regions)

    est_global = n_records * 1024 * 2
    est_bins4 = n_records * 4 * 1024 * 2
    check_disk_budget(est_global, est_bins4, features)

    global_mem = open_memmap("...global_mean.f16.npy", shape=(n_records, 1024), dtype="float16")
    bins4_mem = maybe_open_memmap("...bins4_mean.f16.npy", shape=(n_records, 4, 1024), dtype="float16")

    extractor = GenosFeatureExtractor(...)

    for batch_indices in batched(range(n_records), batch_size):
        seq = fetch_canonical_onehot(records[batch_indices])
        feat = extractor.extract(seq)           # [B, L, 1024]

        global_mem[batch_indices] = feat.mean(dim=1).cpu().numpy().astype(np.float16)
        if bins4_mem is not None:
            bins4_mem[batch_indices] = pool_bins4(feat).cpu().numpy().astype(np.float16)

    write_manifest(...)
```

关键点：

- 这里的 cache view 必须禁用 jitter/revcomp
- train split 的 `records` 顺序要与训练期 dataset 一致
- cache 构建脚本必须输出 `record_sha1`

---

## 10. Go / No-Go 判据

### 10.1 ROI 硬约束

cached pilot 相对 `P0` 的 step-time overhead：

- 必须 `<= 10%`

超出则先排查实现，而不是继续解释指标。

### 10.2 指标判据

在 `epoch 5-8` 的同 epoch 对齐下，满足以下任一项视为 `Go`：

- `P1` 相对 `P0` 的 `peak.profile_target_jsd_full_mean` 改善 `>= 0.002`
- `P1` 的 `JSD` 改善 `>= 0.0015` 且 `count_pearson_full` 提升 `>= 0.005`
- `P2` 的 `count_pearson_full` 提升 `>= 0.01` 且 profile 不退化

`Conditional Go`：

- `P1` 有连续正趋势但未达硬阈值
- 且 `FiLM` 使用统计明显偏离初始化

`No-Go`：

- `P1/P2` 到 `epoch 4` 仍未超过 `P0`
- usage 统计始终接近初始化
- 或 cached path 的开销超出 ROI 上限

### 10.3 第二轮后的收口条件

若做到 `P3/P4` 后，最佳方案仍只有 `< 0.002` 的稳定 JSD 改善，则 `Genos` 线正式降级。

---

## 11. 最终执行顺序

1. 不再继续扩大当前 `G1/G2` online recipe
2. 先实现 cache builder、probe、`P0/P2/P1`
3. 跑 cache 预检查，记录：
   - `train_regions`
   - `train_epoch_regions`
   - `valid_regions`
   - `global/bins4` 预算
4. 构建：
   - 必选：`global_mean`
   - 可选：预算内顺手构建 `bins4_mean`
5. 先跑 `Probe A`
6. 再跑 `Probe B`
7. 串行跑 `P0 -> P2 -> P1`
8. 根据 `Go / No-Go` 决定是否继续 `P3/P4`

---

## 12. 一句话结论

最终合并后的方案是：

- 保留 Claude 文档里的实现细节和后段融合主线
- 吸收策略文档里的语义边界、task-aligned probe、memmap/manifest 和 ROI 判据
- 再补上一个实现上很关键的约束：
  cache 预算必须按 `train_regions` 而不是 `train_epoch_regions` 估算

这让下一轮实验既不至于继续重复当前 recipe 的浪费，也不至于因为过度简化 cache 规模假设而埋下新的工程风险。
