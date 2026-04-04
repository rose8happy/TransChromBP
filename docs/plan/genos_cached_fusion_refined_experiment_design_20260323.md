# Genos Cached Fusion Refined Experiment Design

更新日期：2026-03-23

本文档是在重新阅读 `G0/G1/G2` 实机运行代码、`Phase 1` 原计划以及 Claude 方案
[`docs/plan/genos_cached_fusion_experiment_design_20260323.md`](./genos_cached_fusion_experiment_design_20260323.md)
之后，对下一轮 `Genos` 实验做的收敛版设计。目标不是继续为当前 `online full-seq + gated add`
recipe 补证据，而是用更低成本回答一个更窄也更关键的问题：

> 如果 `Genos` 真有用，它更可能以什么形式起作用；我们应该用什么最低成本的实验把这件事判清？

---

## 1. 真实代码路径先对齐

### 1.1 `G0/G1/G2` 当前到底在跑什么

- 启动入口是 `vendor/transchrombp/transchrombp/scripts/run_genos_pilot.sh`
- 训练配置是 `vendor/transchrombp/transchrombp/configs/train/train_v2fix_genos_profile_select.yaml`
- `G0` 用 `v2fix_baseline.yaml`
- `G1` 用 `v2fix_genos_gate.yaml`
- `G2` 用 `v2fix_genos_mean.yaml`

当前 pilot 的关键口径：

- `batch_size_per_gpu = 20`
- `max_epochs = 20`
- `validate_every_epochs = 5`
- `best_metric = peak.profile_target_jsd_full_mean`
- 单卡运行，不用 2-rank DDP 扩大单 run 的 global batch

### 1.2 `G1/G2` 的 Genos 用法

真实实现不是“离线特征 + 轻量融合”，而是：

1. `train_ddp.py` 在每个 training / validation step 都调用 `genos_runtime.extract(seq)`
2. `GenosFeatureExtractor.extract()` 做：
   - one-hot -> DNA string
   - tokenizer
   - Genos 正向
   - reverse-complement 再跑一次
   - 两次 hidden state 平均
3. 得到的 `genos_feat [B, L, 1024]` 传给模型
4. 在 `local_tower` 之后、Transformer 之前做一次 gated residual add

### 1.3 `G1` 与 `G2` 的唯一区别

- `G1`：`pool_mode = none`
  - 用逐位置 `Genos` 特征
- `G2`：`pool_mode = mean`
  - 先对整条序列做 mean pooling，再广播回每个位置

这意味着 `G2 >= G1` 的判读应当是：

- 当前真正起作用的更像是 region-level / global summary signal
- 而不是“逐位置 `Genos` 表征已经被这套接法成功利用”

### 1.4 训练数据语义约束

真实数据集实现是 `vendor/transchrombp/transchrombp/data/real_data.py`，其中：

- train split：
  - `peak_max_jitter = 500`
  - `random_revcomp = true`
  - `revcomp_prob = 0.5`
- valid split：
  - `jitter = 0`
  - `random_revcomp = false`

此外，train 还会通过 `PeakNonpeakResamplingSampler` 每个 epoch 重新抽 nonpeak 子集。

这直接决定了缓存策略：

- `valid` 可以做 exact cache
- `train` 可以按 `idx` 对齐缓存，但不能把“canonical cached feature”误写成“与当前增强语义完全等价”

---

## 2. 对 Claude 方案的采纳与修正

## 2.1 我同意的部分

- 在线 `Genos` 推理是当前最大的 wall-clock 浪费
- 下一轮应从 `summary-first` 而不是 `full-token` 开始
- 融合点应后移，不应继续只在 Transformer 前做一次被动 gated add
- 新 pilot 必须改成短周期验证，不能再 `5 epoch` 才看一次验证

## 2.2 需要修正的部分

### 修正 A：`train exact cache` 不能按表面理解直接成立

`idx` 对齐在工程上是可行的，因为 sampler 最终仍然返回 dataset 原始索引；但 `__getitem__`
在取样后还会做 `jitter + revcomp`。因此：

- `idx` 对齐只保证“同一条 canonical record 可以取到同一条 cached summary”
- 不保证“cached feature 与当前这一步实际输入序列完全等价”

所以第一轮应把缓存明确分成：

- `validation exact summary cache`
- `training canonical summary cache`

而不是笼统地写成“训练集离线缓存”。

### 修正 B：第一轮不需要一开始就缓存 `bins4`

Claude 方案建议同时缓存 `global_mean` 和 `bins4_mean`。方向没错，但第一轮更稳妥的做法是：

- 先只落 `global_mean`
- 只有当 `global summary` 先跑出正信号，才继续做 `bins4`

原因：

- 当前唯一有力信号是 `G2(mean) >= G1(gate)`
- 第一轮的首要任务是确认“global summary 是否值得继续”，不是一上来就扩大 cache 面
- 这样可以减少实现、存储和排错面

### 修正 C：probe 要更贴近下游任务，而不只是 peak/nonpeak AUC

`Phase 0` 的 peak/nonpeak probe 有参考意义，但它不是当前任务的主目标。更贴近当前主线的问题应该是：

- `Genos summary` 能否解释 `G0` 在 validation 上的 count residual？
- 在已有 `TransChromBP` 表征上，`Genos` 是否提供额外可分信息？

因此第一轮 probe 应改成：

- 主 probe：`count residual ridge probe`
- 次 probe：`peak/nonpeak complement logistic probe`

### 修正 D：cache 后端不要默认 `torch.save(dict)`

训练期是高频随机读，且 DataLoader 有多 worker。第一轮更稳妥的缓存形态应优先考虑：

- `numpy.memmap`
- 或可 mmap 的单 tensor 存储格式
- 外加 `manifest.json`

这样比一次性 `torch.load` 整个大 dict 更适合 worker 侧只读访问，也更容易做元数据校验。

---

## 3. 新方案的工作假设

### H1：当前可用信号更像 region-level summary，而不是逐位置 token

证据：

- `G2(mean)` 目前不比 `G1(gate)` 差
- 当前 gated add 的表达力很弱，不能说明“逐位置信号一定没用”，但足以说明“不值得继续把 full-token online recipe 当主线”

### H2：如果 `Genos` 有用，最先改善的更可能是 count / occupancy，而不是 profile shape

原因：

- `G2` 本质是 global summary 广播
- 当前模型的 `count_signal_head` 本来就是对 `encoded.mean(dim=1)` 做预测
- 这类全局 signal 更可能先帮 count，再间接影响 profile

### H3：下一轮的评价标准必须包含 ROI，而不只是指标增量

如果新方案相对 baseline 仍然显著放慢 step time，那么即便有小幅收益，也不值得保留。

因此新方案的判断必须同时看：

- 指标增益
- step-time overhead
- 使用统计是否真的偏离初始化

---

## 4. 推荐的总体路线

### Phase A：止损与诊断

目标：

- 不再扩大当前 `online full-seq` recipe
- 把下一轮最需要的诊断信号补齐

动作：

1. 当前 `G1/G2` 不再视为必须跑满 `20 epoch`
2. 若 GPU 资源紧张：
   - `G2` 可直接停
   - `G1` 只保留到下一个 validation 节点，或直接停
3. 新方案统一加入：
   - `step_time_total`
   - `step_time_genos_fetch`
   - `genos_usage_mean/std`
   - `film_gamma_mean/std` 或 `count_proj_norm`

### Phase B：先做 cached global summary

目标：

- 用最低工程复杂度验证“global summary 是否足够有用”

第一轮只做：

- `valid exact global_mean cache`
- `train canonical global_mean cache`

不做：

- full-token cache
- bins4 cache
- offset-aware cache

### Phase C：只有在 global summary 为正时，才增加复杂度

第二轮按信息量从低到高推进：

1. `bins4 summary late fusion`
2. `offset-aware global summary cache`
3. 更强 late fusion（如更深的 summary MLP）

不建议在这之前直接跳：

- LoRA
- cross-attention
- 全量微调 `Genos`

---

## 5. 缓存系统设计

## 5.1 第一轮缓存内容

第一轮只缓存一个特征：

- `global_mean [1024]`

定义：

- 对 `Genos layer_6` hidden state 做 RC-average
- 再对序列长度维做 mean pooling

这与 `G2` 的有效信号最接近，但代价最低。

## 5.2 cache 的两种语义

### `valid exact global cache`

- 以 validation records 的真实窗口为准
- 无 jitter
- 无 revcomp
- 与验证期输入完全一致

### `train canonical global cache`

- 以 `dataset.records[idx]` 的 canonical center 为准
- 无 jitter
- 用 RC-average 使 summary 对方向不敏感
- 训练期只是把这条 canonical summary 当作 record-level side information

这一步是“换问题”，不是“完全等价替换”。

## 5.3 cache 的索引和元数据

建议每个 cache 目录包含：

- `train_global_mean.f16.npy`
- `valid_global_mean.f16.npy`
- `manifest.json`

`manifest.json` 至少记录：

- `data_config_path`
- `split`
- `input_len`
- `genos_model_path`
- `layer`
- `bidirectional`
- `n_records`
- `record_sha1`
- `feature_name`

其中 `record_sha1` 对 `[(chrom, center, source), ...]` 序列做哈希，用来防止 cache 与 data config 悄悄错位。

## 5.4 cache 后端

第一轮推荐：

- 存为 `float16` 的 `numpy` 数组
- 训练时用 `np.load(..., mmap_mode="r")`

原因：

- `global_mean` 只有一条向量，随机读取简单
- worker 只读 mmap 更稳
- 不需要先引入复杂的 chunk/shard 体系

## 5.5 训练集为何仍然可以按 `idx` 读 cache

虽然 train 每个 epoch 会 resample nonpeak 子集，但 sampler 返回的仍然是原始 dataset 索引。
因此 dataset 内部直接按 `idx` 取 `global_mean[idx]` 是可行的。

真正不等价的地方不在“索引错位”，而在：

- 这一步的 `seq` 可能经过 jitter
- 这一步的 `seq` 可能经过 revcomp

也因此，若第二轮需要更接近原增强语义，必须进一步记录并利用：

- 当前 step 的 `jitter`
- 当前 step 的 `is_revcomp`

---

## 6. 先做两个低成本 probe

## 6.1 Probe A：count residual ridge probe（主 probe）

目的：

- 判断 `Genos global summary` 是否解释了 baseline 在 count 维度上的残差

做法：

1. 用 `G0 best.pt` 跑 validation
2. 收集：
   - `pred_logcount_full`
   - `true_logcount`
   - `genos_global`
3. 定义残差：
   - `residual = true_logcount - pred_logcount_full`
4. 用 `genos_global` 训练一个简单的 ridge regression 预测 residual

判读：

- 若验证集 `R^2` 接近 `0`，说明 `Genos` 很难补 baseline 的 count 残差
- 若有稳定正 `R^2`，说明 `count-only` 或 `late-fusion` 至少有继续试的价值

## 6.2 Probe B：complement logistic probe（次 probe）

目的：

- 判断 `Genos` 是否在已有 `TransChromBP` 表征上提供额外可分性

做法：

在 validation 上训练三个 logistic regression：

- `A`：`G0 pooled encoded`
- `B`：`genos_global`
- `C`：拼接 `G0 pooled encoded + genos_global`

标签用 `peak / nonpeak`。

判读：

- 若 `AUC(C) ≈ AUC(A)`，说明 `Genos` 对 region-level discrimination 几乎没有额外信息
- 若 `AUC(C) > AUC(A)`，说明 summary signal 至少提供了额外可分性

Probe B 只是辅证，优先级低于 Probe A。

---

## 7. 第一轮 pilot 矩阵

所有第一轮 pilot 统一：

- `max_epochs = 10`
- `validate_every_epochs = 1`
- `early_stop_patience = 3`
- `batch_size_per_gpu = 20`

`P0/P1/P2` 用同一份短配置，只改 model 分支。

### `P0 baseline_short10`

- 不用 `Genos`
- 作用：matched control

### `P1 global_late_film`

输入：

- `genos_global [B, 1024]`

融合：

- Transformer 输出之后
- 用零初始化的 `FiLM` 调制 `encoded [B, L, 256]`

为什么先试它：

- 它同时影响 profile 和 count
- 比“Transformer 前 gated add”更接近当前假设中的 region-level conditioning
- 零初始化可以保证起点与 baseline 对齐

### `P2 global_count_only`

输入：

- `genos_global [B, 1024]`

融合：

- 仅注入到 `count_signal_head` 的中间层

目的：

- 诊断 `Genos` 是否只对 count / occupancy 有帮助

### 第一轮不跑 `bins4`

原因：

- 当前还没证明 `global summary` 值得继续
- 先把最便宜、最强假设的版本做清楚

---

## 8. 第二轮条件触发项

只有当 `P1` 或 `P2` 给出明确正信号时，才进入第二轮。

### `P3 bins4_late_film`

- 额外缓存 `bins4_mean [4, 1024]`
- 对 `encoded` 做 coarse positional FiLM

启动条件：

- `P1` 明显优于 `P0`
- 且怀疑 global summary 仍然过粗

### `P4 offset_aware_global`

- 仍然用 global summary
- 但为若干离散 offset 预计算 cache，例如 `{-384, -192, 0, 192, 384}`
- dataset 按当前 jitter 选择最近 offset 的 summary

启动条件：

- `P1` 有正信号，但波动较大
- 且怀疑 canonical cache 与 train augmentation 语义偏差过大

`P4` 的意义不是“更强模型”，而是验证：

- 训练期语义不匹配是不是压低了 cached summary 的真实收益

---

## 9. Go / No-Go 判据

## 9.1 ROI 判据

第一轮 cached pilot 相对 `P0` 的 step-time overhead 应满足：

- `<= 10%`

如果连 cached global summary 都明显拖慢训练，这条线的工程价值就已经很低。

## 9.2 指标判据

以同 epoch 对齐判断，在 `epoch 5-8` 内满足以下任一条件可视为 `Go`：

- `P1` 相对 `P0` 的 `peak.profile_target_jsd_full_mean` 改善 `>= 0.002`
- `P1` 的 `JSD` 改善 `>= 0.0015` 且 `count_pearson_full` 提升 `>= 0.005`
- `P2` 的 `count_pearson_full` 提升 `>= 0.01` 且 profile 不退化

`Conditional Go`：

- `P1` 有连续正趋势但未过 `0.002`
- 同时 `FiLM` 使用统计明显偏离初始化

`No-Go`：

- `P1/P2` 到 `epoch 4` 仍未超过 `P0`
- `FiLM` 或 count 注入统计始终接近初始化
- 或 cached path 的 step time 仍显著高于 baseline

## 9.3 第二轮后的收口条件

若做到 `P3/P4` 后，最佳方案仍只有 `< 0.002` 的稳定 JSD 改善，则 `Genos` 线正式降级。

---

## 10. 代码改动建议

### 新增文件

- `scripts/build_genos_summary_cache.py`
- `scripts/run_genos_summary_probe.py`
- `scripts/run_genos_cached_pilot.sh`
- `configs/model/v2fix_genos_global_late_film.yaml`
- `configs/model/v2fix_genos_global_count.yaml`
- `configs/train/train_genos_cached_short10.yaml`

### 修改文件

- `vendor/transchrombp/transchrombp/models/genos_adapter.py`
  - 新增 `GenosSummaryFiLM`
- `vendor/transchrombp/transchrombp/models/transchrombp.py`
  - `forward(..., genos_summary=None)`
  - 新增 late FiLM 与 count-only 分支
- `vendor/transchrombp/transchrombp/data/real_data.py`
  - 支持 mmap 读取 cached summary
  - `__getitem__` 返回 `genos_global`
- `vendor/transchrombp/transchrombp/training/train_ddp.py`
  - 支持 `cached_summary` 模式
  - 增加 step-time 与 usage 统计

### 第一轮刻意不改的部分

- `bias_branch`
- loss 口径
- evaluation 主流程
- `Genos` 主体参数

---

## 11. 推荐执行顺序

1. 停止把当前 `G1/G2` 当成主线扩展
2. 实现 `global summary cache + P0/P1/P2`
3. 先构建 `valid/train global_mean cache`
4. 跑 `Probe A/B`
5. 串行跑 `P0 -> P1 -> P2`
6. 按 `Go / No-Go` 决定是否继续 `P3/P4`

若只看“到第一次明确决策”为止的预算，目标应控制在：

- 一次性 cache 构建：约数小时级
- 第一轮 pilot：约半天级 GPU 时间

这要明显优于当前 `G1/G2` 单 run 就需要数十小时的状态。

---

## 12. 最终建议

新的主线不应再是：

- `online full-seq Genos`
- `Transformer 前 gated add`
- `5 epoch` 才看一次验证

而应改成：

- `cached global summary`
- `late fusion`
- `task-aligned probe`
- `10 epoch 内收口的短 pilot`

如果这套更合理、更便宜的 summary-first 方案仍然只能给出噪声级收益，那么就可以更干净地结束
`Genos-1.2B` 这条线，把 A6000 让回给多 seed、数据扩容或别的 foundation model 方向。
