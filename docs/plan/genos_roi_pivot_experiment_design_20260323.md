# Genos ROI Pivot Experiment Design

更新日期：2026-03-23

本文档用于替代“继续把当前 `G1/G2` 最小 frozen pilot 硬跑完”的默认思路，转而回答一个更实际的问题：

> 在当前 `TransChromBP` 主线里，怎样使用 `Genos-1.2B` 才有机会带来**足够收益**，并且不把 wall-clock 成本放大到不可接受？

核心判断：

- 当前 `Genos` 线的问题，不只是“收益小”，而是“**当前用法的 ROI 很差**”。
- 下一步不应继续扩大“在线 full-seq frozen extractor + 单点 gated add”的投入，而应切换到**缓存/摘要化/后段融合**的方案，用更短反馈周期回答：
  1. `Genos` 是否只提供粗粒度 global signal？
  2. 若是，这个 signal 能否以足够低的代价稳定带来净增益？

---

## 1. 当前背景（截至 2026-03-23 00:13 CST）

当前 `6000 / A6000` 上的 3 条 run：

- `G0 baseline`
  - 已在 `epoch 20/20` 收口
  - `best peak.profile_target_jsd_full_mean = 0.3336565`
  - `best epoch = 20`
- `G1 genos_gate`
  - 仍在运行
  - 标准输出推进到 `epoch 9 step 5160/11809`
  - `epoch_metrics.jsonl` 仍只写到 `epoch 8`
  - 目前唯一 validation 节点仍是 `epoch 5`
  - `epoch 5 peak.profile_target_jsd_full_mean = 0.3419630`
- `G2 genos_mean`
  - 仍在运行
  - 标准输出推进到 `epoch 6 step 5660/11809`
  - `epoch_metrics.jsonl` 仍只写到 `epoch 5`
  - 目前唯一 validation 节点仍是 `epoch 5`
  - `epoch 5 peak.profile_target_jsd_full_mean = 0.3414136`

直接观察：

- `G1/G2` 相对 `G0` 的早期收益只有很小量级。
- `G2(mean)` 目前不比 `G1(gate)` 差，反而略好。
- 当前 wall-clock 代价极大：
  - `G0` 完整 `20 epoch` 已收口；
  - `G1/G2` 到午夜仍远未到下一次 validation。

这已经足够说明：

- 继续把当前 recipe 当成主方案扩下去，不合理。
- 下一步应转为**用法诊断**，不是继续堆算力。

---

## 2. 对当前用法的诊断

### 2.1 不是“Genos 没用”，而是“当前接法没把它用对”

当前实现是：

- `GenosFeatureExtractor` 作为冻结外部 runtime，在 trainer 侧运行；
- 每个 step 都执行：
  - one-hot -> DNA string
  - tokenizer
  - Genos 正向
  - 反向互补再跑一次
  - 两者平均
- 然后把得到的 `genos_feat` 在 `local_tower` 之后、Transformer 之前做一次 gated residual add。

这条链路对“最小可行性验证”是合理的，但对“争取足够收益”并不合理，因为：

- 计算成本主要花在**重复的 frozen 特征提取**上；
- 融合方式只有一个小 adapter，表达能力很弱；
- `gate_bias_init=-2.0` 让模型从一开始就偏向“先别用 Genos”；
- 当前结果显示 `G2(mean) >= G1(gate)`，说明逐位置特征的价值并没有被体现出来。

### 2.2 在线提特征确实是主要浪费，但“训练集 exact cache”没有表面那么简单

当前 `ChromBPNetBigWigDataset` 在 `train` split 上启用了：

- peak `jitter`：`peak_max_jitter = 500`
- `random_revcomp = true`
- `revcomp_prob = 0.5`

而 `val` split 上：

- `jitter = 0`
- `random_revcomp = false`

这意味着：

- **validation exact cache 很容易做**
  - 同一条 valid region 每次取到的序列窗口固定不变；
  - 可以直接缓存 full-token 或 summary 特征。
- **training exact cache 不能按当前语义“直接一键替换”**
  - 同一个 region 在不同 step / epoch 会因为 jitter 和 revcomp 产生不同输入；
  - 如果要完全等价缓存，需要覆盖大量 offset / orientation 组合；
  - 当前 `11809` 是 `steps/epoch`，不是 `samples`。在 `bs=20` 下，对应约 `236k` windows / epoch；再叠加 augmentation，训练窗口远多于一个静态 region 列表。

因此，“把当前训练集逐样本 full-token 特征离线缓存掉”这个说法，方向是对的，但实现路径必须收缩：

- 先做 `validation exact cache`
- 再做 `training canonical summary cache`
- 若 summary 版本有明显信号，再考虑更重的 offset-aware cache

---

## 3. 当前结果真正告诉了我们什么

### 3.1 `G2 >= G1` 是最重要的信号

当前对照组设计本来就是为了回答：

- `G1`：逐位置 Genos 注入
- `G2`：全局均值广播

若 `G1 > G2 > G0`，才说明逐位置 Genos 特征真的有用。

但现在看到的是：

- `G2(mean)` 不比 `G1(gate)` 差；
- 说明当前收益更像来自**粗粒度 global summary**；
- 而不是来自高成本的逐位置 Genos 表征。

这直接决定了 pivot 方向：

- 不要继续优先做 “full-length per-position online fusion”
- 应先做 “summary-first, late-fusion, cached”

### 3.2 当前 feedback cycle 太慢

当前训练配置：

- `validate_every_epochs = 5`
- `max_epochs = 20`

这会导致：

- 一个坏方案要烧很久，才能看到第二个验证点；
- 不适合做“用法重设计”阶段的诊断。

新的方案必须把反馈周期改短到：

- `validate_every_epochs = 1` 或 `2`
- 第一轮 pilot 先做 `6-10 epoch`

---

## 4. Pivot 目标

新的 Genos 实验不再优先回答：

> “逐位置 frozen Genos feature 能否直接改善主线？”

而改为优先回答：

> “如果只保留 Genos 最可能有效的那部分信息（粗粒度 summary），能否用明显更低的代价带来足够收益？”

为此，新的实验设计需要同时满足 3 个目标：

1. 把 wall-clock 成本压回接近 baseline 的量级
2. 单独检验 coarse/global signal 是否是真正有效来源
3. 让每个 pilot 在 `1-2` 次 validation 内就能给出 go/no-go 信号

---

## 5. 新方案总览

建议把下一轮 Genos 实验改成三阶段：

### Phase A：止损 + 基础诊断

目标：

- 不再扩大当前 recipe 的 sunk cost
- 补齐“到底有没有在用 Genos”的诊断信息

动作：

1. 停止把 `G1/G2` 当作需要跑满 `20 epoch` 的主线
2. 为后续所有 Genos run 补两类日志：
   - `genos_extract_time`
   - `gate_mean / gate_std / gate_p90 / proj_norm`
3. 新 pilot 把 `validate_every_epochs` 改成 `1` 或 `2`

当前 run 的建议处理：

- `G2`：可直接停
  - 它离下一次 validation 还远；
  - 当前已经说明 `mean` 口径是有效诊断方向，不需要继续烧到更晚。
- `G1`：只保留到下一次 `epoch 10` validation，或在 GPU 紧张时直接停
  - 不建议再以“跑满 20 epoch”为目标。

### Phase B：缓存化摘要特征（主推荐）

目标：

- 保留 Genos 最可能有用的 coarse signal
- 避免训练期在线 full-token 提特征

核心思路：

- **validation**
  - 做 exact cache
  - 每个 valid region 的 Genos 特征固定缓存
- **training**
  - 不做 exact full-token cache
  - 只做 **canonical summary cache**
    - 基于 record 的 canonical center（无 jitter）
    - 使用 RC-average 后的 Genos summary
    - 训练期不随 jitter 改变

这不是对当前 full-token 方案的语义等价替换，而是刻意回答另一个更窄的问题：

> 若 Genos 真正有用的只是 region-level 粗粒度信息，那么把它做成稳定 summary 后，是否还能带来净增益？

### Phase C：只有在 Phase B 明显为正时，才做 richer 版本

只有当 summary 方案先跑出明确正信号，才继续做：

- `offset-aware summary cache`
- `4/8 bin summary token`
- `late-block FiLM / cross-token conditioning`

不建议在 Phase B 之前直接跳：

- LoRA
- cross-attention
- full-token cache over jitter grid
- Genos 主体微调

---

## 6. 具体实验矩阵

第一轮建议只跑 3 个最便宜、信息量最高的对照：

| 组别 | 名称 | Genos 表示 | 融合位置 | 目的 |
|---|---|---|---|---|
| `P0` | `baseline_short10` | 无 | 无 | 匹配短 pilot 基线 |
| `P1` | `genos_summary_global_count` | 全局 summary（1 向量） | count head / late shared MLP | 验证粗粒度 signal 是否足以改善计数与整体 JSD |
| `P2` | `genos_summary_bins4_late` | 4-bin summary（中心窗口分桶） | 最后 1-2 个 block 后的 late fusion / FiLM | 验证 coarse positional summary 是否优于纯 global summary |

若 `P2 > P1 > P0`，再进入第二轮：

| 组别 | 名称 | Genos 表示 | 融合位置 | 目的 |
|---|---|---|---|---|
| `P3` | `genos_summary_bins8_late` | 8-bin summary | late fusion | 检查更细 summary 是否继续有效 |
| `P4` | `genos_offset_summary` | 离散 offset summary cache | late fusion | 检查 jitter-aware 摘要是否比 canonical summary 更稳 |

### 6.1 `P0 baseline_short10`

要求：

- 完全沿用当前主线 baseline 语义
- 不开 Genos
- `max_epochs = 10`
- `validate_every_epochs = 1`
- 作为所有 summary pilot 的 matched control

### 6.2 `P1 genos_summary_global_count`

输入：

- 对每个 region 预计算一个全局 summary 向量
- 推荐先对 canonical `2114bp` 或中心 `1000bp` 做 mean summary
- 继续使用 RC-average，确保方向不敏感

融合方式：

- 先不要进入 backbone
- 先接到：
  - count head
  - 或 count/profile 共享的 late MLP

原因：

- 当前 `G2` 的正信号更像 global summary；
- 先让最有希望、最便宜的路径回答“值不值得继续”。

### 6.3 `P2 genos_summary_bins4_late`

输入：

- 将 canonical center window 切成 `4` 个 coarse bins
- 每个 bin 各取一个 summary
- 拼成 `[4, H]` 或 flatten 后送入 late adapter

融合方式：

- 不再在 `local_tower` 后做逐位置 add
- 而是在：
  - Transformer 输出后
  - 或最后 `1-2` 个 block 之后
  做 late conditioning / FiLM / gated MLP

理由：

- 当前证据不支持 full-length per-position fusion
- 但也不能过早下结论说“位置结构完全没用”
- `4-bin summary` 是 global 与 full-token 之间更合理的中间点

---

## 7. 训练与评估口径

### 7.1 训练配置

新的 pivot pilot 建议统一：

- `max_epochs = 10`
- `validate_every_epochs = 1`
- `early_stop_patience = 3`
- `best_metric = peak.profile_target_jsd_full_mean`
- `best_metric_mode = min`

理由：

- 目标是缩短反馈周期，而不是追求单次长跑的最终最好值
- 若一个用法在 `epoch 4-6` 仍没有显著优于 baseline，就不该继续吃卡

### 7.2 指标与判断优先级

主指标：

1. `peak.profile_target_jsd_full_mean`
2. `peak.count_pearson_full`

辅助诊断：

- `overall.profile_target_jsd_full_mean`
- `overall.count_pearson_full`
- `gate_mean / gate_std / gate_p90`
- `genos_extract_time`
- `step_time`

### 7.3 Go / No-Go 阈值

第一轮 pilot 建议采用以下规则：

#### 继续推进（Go）

满足以下任一条件：

- 到 `epoch 4-6`，相对 matched baseline 有
  - `peak.profile_target_jsd_full_mean` 绝对改善 `>= 0.003`
- 或同时满足：
  - `peak.profile_target_jsd_full_mean` 改善 `>= 0.002`
  - `peak.count_pearson_full` 提升 `>= 0.005`

#### 直接收口（No-Go）

满足以下任一条件：

- 所有 Genos 变体在前 `2` 次 validation 中都与 baseline 接近：
  - `peak JSD` 差异 `< 0.0015`
  - 且 `count_r` 差异 `< 0.003`
- `P2` 仍不优于 `P1`
  - 说明 coarse positional summary 也没有额外价值
- `gate` 长期接近关闭
  - 说明模型实际忽略 Genos

---

## 8. 工程改动建议

### 8.1 先做的改动

1. 新增缓存脚本
   - `scripts/build_genos_summary_cache.py`
2. 新增摘要 cache 数据结构
   - 先覆盖 validation exact cache
   - 再覆盖 training canonical summary cache
3. dataset / collate 增加稳定 `record_id`
   - 便于根据 region 读取 cache
4. trainer 支持从 cache 读取 `genos_summary`
   - 仅在启用 summary 模式时生效
5. logging 增加：
   - `genos_extract_time`
   - `gate` 统计

### 8.2 暂缓的改动

- 不先做 train split 的 full-token exact cache
- 不先做 LoRA
- 不先做 cross-attention
- 不先让 Genos 进入 optimizer

---

## 9. 推荐执行顺序

### 第 1 步：止损

- 不再新启动任何“在线 full-token Genos + 20 epoch” run
- 当前 `G2` 停止
- `G1` 只保留到下一次 validation，或在 GPU 紧张时直接停止

### 第 2 步：基础设施

- 完成 `validation exact cache`
- 完成 `training canonical summary cache`
- 完成 `record_id` / `summary loader` / timing & gate logging

### 第 3 步：第一轮短 pilot

- `P0 baseline_short10`
- `P1 genos_summary_global_count`
- `P2 genos_summary_bins4_late`

### 第 4 步：按结果决定

- 若 `P1/P2` 没有跨过 go 阈值
  - Genos 线立即降级
  - A6000 释放给多 seed / 数据扩容 / 非 Genos 主线
- 若 `P2` 给出明确正信号
  - 再进入 richer summary / offset-aware 版本

---

## 10. 与 Claude 报告的一致与修正

本方案与 Claude 报告一致的地方：

- 同意“在线 frozen extractor 是当前最大浪费”
- 同意“不应继续把当前 recipe 当主线长期烧卡”
- 同意“若差距继续停留在很小量级，Genos 线应收口”

需要修正的地方：

- 当前日志中的 `11809` 是 `steps/epoch`，不是 `samples`
- 训练集因为 `jitter + revcomp`，**不能**按当前语义直接做一把梭的 exact full-token cache
- 因此更现实的 pivot 不是“先做训练集 full exact cache”，而是：
  - 先做 `validation exact cache`
  - 再做 `training canonical summary cache`

---

## 11. 最终建议

建议把下一轮 Genos 方案定为：

> **“cached summary-first Genos”**

而不是继续沿用：

> **“online full-seq frozen extractor + per-position gated fusion”**

若 summary-first 方案在短 pilot 中仍然没有给出足够收益，则应把 Genos 线正式降级，不再继续消耗 A6000 主训练窗口。
