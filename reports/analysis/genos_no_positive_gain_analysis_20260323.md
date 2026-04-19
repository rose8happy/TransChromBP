# 2026-03-23 Genos 未带来正收益的原因分析

## 1. 要回答的问题

当前仓库里已经有两组看起来“方向相反”的证据：

1. 官方最接近染色质可及性的旁线 benchmark `human_ocr_ensembl` 已被本地成功复现。
2. `TransChromBP` 主线里，不论是 online Genos (`G1/G2`) 还是 cached-fusion (`P2`) 都没有带来正收益，其中 `P2` 甚至在 held-out test 上出现 count 分支明显塌陷。

这份报告的目标不是重复结果表，而是回答：

> 为什么 `Genos-1.2B` 本体看起来没问题，但引入到当前 `TransChromBP` 主线里仍然没有正收益？

## 2. 先排除一种错误解释：不是 Genos 本体或环境坏了

### 2.1 OCR benchmark 已证明本地 Genos 推理链路基本正常

OCR 旁线结果：

- mini smoke (`layer 6`)：`roc_auc=0.5958`
- full `layer 6`：`roc_auc=0.7242`
- full `layer 12`：`roc_auc=0.7535`

其中 `layer 12 = 0.7535` 已非常接近文档中整理的官方 `Genos-1.2B` 参考值 `0.7569`。

证据：

- OCR 摘要与定位：[genos复现ocr.md](/home/zhengwei/project/python/TransChromBP/docs/research/genos复现ocr.md)
- OCR 执行计划：[genos_ocr_reproduction_plan_20260323.md](/home/zhengwei/project/python/TransChromBP/docs/plan/archive/genos/genos_ocr_reproduction_plan_20260323.md)
- full OCR 结果 TSV：`/data1/zhoujiazhen/bylw_atac/genos_ocr_sidecar/results/Genos-1.2B/human_ocr_ensembl.tsv`
- mini OCR 结果 TSV：`/data1/zhoujiazhen/bylw_atac/genos_ocr_sidecar/results_mini/Genos-1.2B/human_ocr_ensembl.tsv`

这说明：

- 权重、tokenizer、`trust_remote_code`、hidden-state 提取、mask mean pooling、MLP benchmark 链路基本都正常
- 因此当前主线失败，不能优先归咎为“Genos 没跑通”或“本地环境坏了”

### 2.2 主线失败不是随机波动，而是跨 recipe 一致出现

online Genos：

- `G1 gate`：test peak `0.3388 / 0.4689`
- `G2 mean`：test peak `0.3387 / 0.6360`

cached-fusion：

- `P0`：test peak `0.3179 / 0.8384`
- `P2`：test peak `0.3174 / 0.6037`

相比之下，baseline `G0` 为 `0.3163 / 0.8410`。

证据：

- online Genos 汇总：[genos_and_6002_run_analysis_20260323.md](/home/zhengwei/project/python/TransChromBP/reports/analysis/genos_and_6002_run_analysis_20260323.md#L3)
- cached-fusion 汇总：[genos_cached_fusion_status_20260323.md](/home/zhengwei/project/python/TransChromBP/reports/analysis/genos_cached_fusion_status_20260323.md#L20)

这说明：

- 失败不是单一配置的偶发问题
- Genos 引入主线的难点是“如何让它对当前 profile/count 目标产生稳定增益”，而不是“如何把模型跑起来”

## 3. 现有证据真正说明了什么

### 3.1 Genos 里确实有信号，但信号弱，而且不是天然可加

cached probe 结果：

- `count residual ridge`：`R^2=0.0167`, `Pearson r=0.1987`
- `genos_only`：`AUC=0.6997`
- `encoded_only`：`AUC=0.9219`
- `concat(encoded, genos)`：`AUC=0.9018`

证据：

- probe JSON：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/genos_cached_probes_20260323/probe_results.json`
- probe 日志：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/logs/genos_summary_probe_20260323_1111.log`

这组数值给出的信息非常关键：

1. `genos_global_mean` 不是纯噪声。
2. 但它对主任务只有“小而稳定”的相关性，不是强解释变量。
3. 更重要的是，`encoded + genos` 反而低于 `encoded_only`，说明 Genos signal 并不是“无脑拼上去就增益”的类型。

一句话说：

- Genos 有信息
- 但这个信息弱、粗，而且和当前 backbone 已提取的信息存在较强重叠或冲突

### 3.2 当前第一轮 cached 方案故意只用了 `global_mean`

执行计划里明确规定，第一轮 pilot 只消费 `global_mean [1024]`，不直接用 `bins4_mean [4, 1024]` 或更细粒度特征。

证据：

- 设计决策：[genos_cached_fusion_final_execution_plan_20260323.md](/home/zhengwei/project/python/TransChromBP/docs/plan/archive/genos/genos_cached_fusion_final_execution_plan_20260323.md#L127)

这意味着当前 cached 线测试的其实不是“Genos 全部能力”，而是：

> 一个非常粗的 region-level summary，是否已经足够给当前 profile/count 主线带来稳定收益

`P2` 的失败说明至少在这个层级上，答案是否定的。

## 4. 最可能的失败原因排序

下面按“解释力”从高到低排序。

### 原因 1：特征粒度不对，`global_mean` 对 profile/count 任务太粗

OCR benchmark 之所以能跑通，是因为它本来就是：

- sequence-level classification
- hidden-state pooling
- MLP classifier

也就是说，OCR 任务天然适合 `global summary`。

但当前 `TransChromBP` 主线不是 sequence-level 分类，而是：

- profile 轨迹预测
- count 标量预测
- 并且两者都和局部 motif 组合、位置结构、bias decomposition 强相关

对这种任务来说，`global_mean` 很可能过于粗糙，只能保留“这个 region 大概像不像开放区域”的弱 summary，却保不住真正决定 profile/count 的局部结构细节。

因此：

- OCR 成功只能说明 `Genos` summary 对 OCR 分类有效
- 不能推出 `global_mean` 对 ATAC-style profile/count 预测也有效

### 原因 2：训练语义与 cached summary 语义存在系统性错位

执行计划里已经明确写过：

- `valid` 是 exact cache
- `train` 只是 canonical center side information

证据：

- 语义约束：[genos_cached_fusion_final_execution_plan_20260323.md](/home/zhengwei/project/python/TransChromBP/docs/plan/archive/genos/genos_cached_fusion_final_execution_plan_20260323.md#L118)

而当前训练语义同时还包含：

- `peak_max_jitter = 500`
- `random_revcomp = true`
- nonpeak resampling

因此，训练时模型实际看到的增强输入，和 cached summary 对应的 canonical record 不是完全同一个对象。

如果 Genos signal 本来就弱，这种错位足以进一步降低它的可学习性：

- model 学到的是“summary 对某个 canonical region 的统计”
- 但 supervised target 来自 jitter/revcomp 后的训练样本语义

这不是致命 bug，但会显著压缩外部 summary 能带来的净收益。

### 原因 3：`P2 count_only` 的融合位置太激进，直接破坏了 count head

`P2` 的实现不是对整个 backbone 轻量调制，而是把 `genos_summary` 直接加到 count head 的隐藏层中间：

- 先 `LayerNorm -> Linear -> GELU`
- 然后 `h = genos_count_proj(h, genos_summary)`
- 再 `Dropout -> Linear -> count_signal`

证据：

- 代码位置：[transchrombp.py](/home/zhengwei/project/python/TransChromBP/vendor/transchrombp/transchrombp/models/transchrombp.py#L370)

这带来的问题是：

1. 注入点离最终 count 输出太近，几乎没有缓冲层来重新吸收错误信号。
2. 注入的是 `global_mean` 这种粗 summary，不是高保真局部特征。
3. 训练目标里 count 分支本来就更敏感，容易被外来 summary 直接改坏。

`P2` 的 held-out 结果正好符合这个机制：

- peak `JSD=0.3174`，说明 profile 路径没有明显被毁掉
- peak `count_r=0.6037`，说明被直接注入的 count 路径严重失真

所以 `P2` 的失败更像：

> “把弱而粗的外部 signal 直接打进了最脆弱的标量读出头”

而不是“Genos 对 count 完全没有任何信息”。

### 原因 4：当前 backbone 已经很强，Genos 的增量空间本来就很小

当前最强 baseline：

- `G0 test = 0.3163 / 0.8410`

连 matched `P0` 也只是：

- `P0 test = 0.3179 / 0.8384`

两者已经非常接近。

证据：

- 汇总表：[genos_cached_fusion_status_20260323.md](/home/zhengwei/project/python/TransChromBP/reports/analysis/genos_cached_fusion_status_20260323.md#L20)

这意味着当前 `TransChromBP` backbone + bias decomposition 已经吃掉了大部分容易拿到的信号。对于这样一个强 baseline，外部 foundation summary 如果不是强增益，就很容易只增加复杂度、不增加效果。

换句话说：

- 当 baseline 很弱时，Genos summary 可能还能刷出表面提升
- 但在当前 baseline 已经稳定、强、偏任务专用的前提下，小信号很难越过噪声门槛

### 原因 5：online 与 cached 两条失败形态不同，但都指向“融合不自然”

online `G1/G2` 的失败表现是：

- validation 看起来不算太差
- held-out test 上尤其是 count 分支塌得很厉害

cached `P2` 的失败表现是：

- validation 没那么难看
- held-out profile 还能追回
- 但 count 路径依然崩

证据：

- online 结果：[genos_and_6002_run_analysis_20260323.md](/home/zhengwei/project/python/TransChromBP/reports/analysis/genos_and_6002_run_analysis_20260323.md#L11)
- cached 结果：[genos_cached_fusion_status_20260323.md](/home/zhengwei/project/python/TransChromBP/reports/analysis/genos_cached_fusion_status_20260323.md#L29)

这说明问题更像：

- “当前监督目标与 Genos signal 的耦合方式不自然”

而不是：

- “只有 online 坏，cached 就一定好”
- “只有 cached 坏，online 就一定值得救”

## 5. 现阶段最合理的结论

基于现有证据，更稳的结论应该是：

1. `Genos-1.2B` 本体与本地推理链路没有明显问题。
2. `Genos` 确实携带和开放染色质有关的序列表征，但当前可验证的主要是 sequence-level OCR 分类能力。
3. 在当前 `TransChromBP` 主线里，`Genos` 没能形成正收益，最主要不是因为模型坏了，而是因为：
   - 当前使用的 summary 太粗
   - train cache 语义和真实训练语义不完全对齐
   - `P2` 的 count-only 融合位置过于激进
   - backbone baseline 已经很强，增量空间很小

因此，当前负结果更应该解释为：

> “当前这套任务定义 + 特征粒度 + 融合位置组合不对”

而不是：

> “Genos 完全无用”

## 6. 对后续动作的直接影响

### 6.1 当前不该做的

- 不继续顺推 `P1 global_late_film`
- 不再追加 `G1/G2` 或 `P2` 同类 recipe 的新 seed
- 不把 OCR benchmark 继续当作主线 gate

### 6.2 如果以后要重开 Genos，优先顺序应该变

如果未来还要重开 Genos 线，优先级应当是：

1. **先换特征粒度**
   - 从 `global_mean` 转到 `bins4_mean` 或其他局部 summary
   - 核心目标是保留一定位置结构，而不是继续押单个全局向量
2. **再换融合位置**
   - 避免像 `P2` 那样直接把 summary 打进 count head 中间层
   - 更偏向 backbone 中后段的轻调制，而不是末端读出头硬注入
3. **最后才考虑是否继续扩大训练规模**
   - 在没有新的更合理 recipe 前，继续烧更多 epoch / 更多 seed 的信息增益很低

## 7. 一句话总结

官方 OCR 复现已经证明 `Genos` 本身是通的；当前没有正收益，最像是“弱而粗的 global summary 被以不合适的方式接进了一个已经很强的任务专用 backbone”，而不是 `Genos` 完全没有可及性相关信息。
