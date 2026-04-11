# 2026-04-11 AlphaGenome-like Factor Ladder 设计

## 1. 设计目标

本设计不把“重新看 U-Net”定义成一次单点架构翻盘尝试，而定义成一条**拆因实验线**。

要回答的核心问题是：

1. 在本项目当前任务上，`长上下文` 是否本身就有价值。
2. 如果只做 `readout-only` 的 U-Net-like patch 不够，`全模型 encoder-decoder` 是否才是关键。
3. 如果结构本身不够，`teacher/distillation` 是否是更接近 AlphaGenome 成功机制的必要因素。

本设计的成功标准不是“必须立刻超过 `corrected B`”，而是：

- 即使没有出现新的最优模型，只要能把 `长上下文 / 全模型 encoder-decoder / distillation` 三者里谁最值得继续押注明确拆出来，这轮就算成功。

## 2. 当前上下文与为什么要这样设计

仓内现有证据已经足够支持以下结论：

- `6002` 的 `U-Net-lite v1` 有效 run `r3/r4` 不具备晋级资格。
- `6000` 的 `multiscale/local-skip decoder v2` formal gate 也明确失败。
- 因此，当前不应再把“在 `corrected B` 上继续追加 readout-only 的 U-Net-like decoder”当作默认高价值动作。

但这些负结果只否掉了一个较窄的假设：

> 在当前 `2114 bp -> 中心 1000 bp`、ChromBPNet-style bias-safe 语义下，readout-only 的 U-Net-like 局部升级会自然带来 clean gain。

它们并不足以否掉 AlphaGenome 那种更大范式：

- 更长上下文
- 整模型 encoder-decoder
- 多任务/多模态监督
- teacher/distillation

因此，下一轮设计必须避免两个错误：

1. 继续给 `unet_lite_v1` / `msdls_v2` 做同配方 cheap rerun。
2. 直接做一个混入太多变量的 “AlphaGenome-lite 大杂烩”，最后无法解释谁起作用。

## 3. 方案选择

本轮曾考虑 3 类路线：

### 3.1 方案 A：因子阶梯

- `E1`：只放大上下文
- `E2`：在同样长上下文上改成全模型 encoder-decoder
- `E3`：在 `E2` 上叠加 teacher/distillation

优点：

- 因果最干净
- 结果正负都能解释
- 最符合“先拆因，再决定是否转主线”的目标

缺点：

- 需要先定义一个过渡任务
- 前期会有少量基础设施搭建成本

### 3.2 方案 B：AlphaGenome-lite 一步到位

- 直接做一个长上下文 + 多尺度主干 + distill 的新模型

优点：

- 最像 AlphaGenome

缺点：

- 好了也不知道是谁起作用
- 坏了也不知道死在结构、上下文还是训练范式

### 3.3 方案 C：只看 distillation

- 保持主体结构接近当前模型，只测试 teacher/distill

优点：

- 起步最快

缺点：

- 几乎不能回答 “encoder-decoder 是否关键”
- 也不能回答 “长上下文本身是否有价值”

### 3.4 选型结论

本设计采用 **方案 A：因子阶梯**。

原因：

- 用户已明确这轮优先目标是“高价值、能最快改变路线判断”，同时成功标准偏向“拆因”
- 因子阶梯最符合这个目标
- 它允许后续根据结果决定是继续向 AlphaGenome-like 范式推进，还是及时停表

## 4. 固定边界

### 4.1 这轮固定不动的东西

1. 主要评价问题仍回到当前 ATAC `profile + count`
2. `output_len=1000`
3. `supervised_bp=1000`
4. `count head` 保持独立
5. `bias-safe` 语义不取消
6. `full/debiased` 口径继续保留

### 4.2 这轮允许变化的东西

1. `input_len`
2. 主干是否从当前结构改成全模型 encoder-decoder
3. teacher/distillation 训练路径

### 4.3 这轮明确不做的事

1. 不继续重跑当前 `unet_lite_v1` / `msdls_v2` family
2. 不在第一轮就引入多模态大任务
3. 不把 `count` 分支并入 decoder
4. 不一次性把长上下文、全模型 encoder-decoder、distill、任务变更全部混在一个 run 里

## 5. 过渡任务定义

本设计采用一个更接近 AlphaGenome、但仍与当前主任务可比较的过渡任务：

- `long-context centered ATAC`

定义如下：

- 输入：更长 DNA 窗口
- 输出：中心 `1000 bp` 的 ATAC `profile + count`
- 评价：继续使用当前主任务评价指标

这样做的目的不是复刻 AlphaGenome，而是把“更大感受野”和“base-resolution 输出形成”这两个问题抽出来单独回答。

## 6. 最小实验矩阵

### 6.1 共同上下文长度

第一轮统一使用：

- `input_len=max_seq_len=4096`

理由：

1. `4096` 已足以把任务从当前 `2114` 变成新的上下文实验线
2. 仓内现有数据与训练管线本身支持改 `input_len / supervised_bp`
3. 当前 `RoPE` 实现支持动态扩 cache，不会被位置编码卡住
4. 第一轮先用 `4096` 拆因更干净；若有清晰正信号，再把胜者升到 `8192`

### 6.2 `E1`: Long-context only

定义：

- 基座保持 `corrected B`
- 只把 `input_len=max_seq_len` 从 `2114` 改到 `4096`
- `output_len=1000`、`supervised_bp=1000` 不变
- profile 仍走当前直读式主路径

要回答的问题：

> 长上下文本身是否就值得押注

### 6.3 `E2`: Full-model encoder-decoder

定义：

- 在 `4096` 输入上，引入真正的层级式 `encoder-decoder`
- 不是只换 profile readout
- 需要在主干中显式包含：
  - 2-3 级 downsample/压缩
  - transformer bottleneck
  - 多级 upsample + skip 恢复
- `count head` 仍独立
- `bias fusion` 仍发生在 debiased profile 之后

要回答的问题：

> 之前失败的关键原因，是不是因为我们只改了 readout，而没有把 encoder-decoder 放进主干

### 6.4 `E3`: `E2 + teacher/distillation`

定义：

- student 结构保持 `E2`
- 额外引入 model-teacher distillation
- teacher 不是旧 NT foundation proxy，也不是 AlphaGenome API 输出
- teacher 默认来自更长 budget 的 `E2` 模型本身，或其小规模 ensemble

要回答的问题：

> teacher/distillation 是否是接近 AlphaGenome 成功机制的重要边际因素

## 7. 评估与判读规则

### 7.1 统一指标

主指标：

1. `peak profile JSD`
2. `peak count Pearson`

补指标：

1. `overall profile JSD`
2. `profile Pearson`
3. `count Spearman`
4. `count RMSE / MAE(logcount)`
5. `full/debiased gap`

判读顺序固定为：

1. 先看主指标
2. 再看补指标是否解释主指标变化
3. 最后看 `gap` 是否恶化

### 7.2 单步结果分类

`positive`：

- 满足以下任一条：
  - `peak JSD` 改善 `>= 0.003`，且 `peak count Pearson` 不下降超过 `0.005`
  - `peak JSD` 变化 `<= 0.001`，且 `peak count Pearson` 提升 `>= 0.01`
- 同时要求：
  - `full/debiased gap` 不明显恶化
  - `count RMSE/MAE` 不出现明显反向塌陷

`flat / weak`：

- `peak JSD` 变化 `< 0.003`
- `peak count Pearson` 变化 `< 0.01`
- 补指标也不给出一致方向

`negative`：

- 满足任一条：
  - `peak JSD` 恶化 `>= 0.003`
  - `peak count Pearson` 下降 `>= 0.01`
  - `gap / calibration` 明显变差

### 7.3 相邻比较关系

`E1 vs corrected B`

- 回答长上下文本身值不值

`E2 vs E1`

- 回答全模型 encoder-decoder 是否才是关键

`E3 vs E2`

- 回答 distillation 是否提供额外边际价值

### 7.4 路线判读

- `E1 positive`：长上下文本身值得继续押
- `E1 flat/negative, E2 positive`：关键不在上下文单独，而在全模型 encoder-decoder
- `E2 positive, E3 positive`：结构和训练范式都有贡献
- `E2 positive, E3 flat`：结构更关键，distill 不是主要矛盾
- `E2 negative`：AlphaGenome-like 路线在当前任务上的短期价值明显下降

## 8. 执行顺序

### 8.1 Phase 0：基础 smoke

`S0-1`：`E1-4096` 数据/前向 smoke

- 验证 `input_len=4096` 的 loader、jitter、中心 `1000bp` 监督是否正常

`S0-2`：`E2-4096` 结构 smoke

- 只跑几十到几百 step
- 验证 shape、显存、吞吐、loss 是否稳定

`S0-3`：`E3` teacher cache spike

- 用极小子集验证：
  - model prediction
  - record alignment
  - teacher cache write/read

`Phase 0` 只排基础设施，不做科学判读。

### 8.2 Phase 1：拆因 short10

- `E1` 与 `E2` 都跑 `short10`
- 默认固定 `seed=42`
- 两者可以并跑

理由：

- 这轮目标是拆因，不是层层晋级
- 即使 `E1` 负，`E2` 也仍然可能正

### 8.3 Phase 2：只在 `E2` 有解释价值时进入 distill

只有以下情况才启动 `E3`：

1. `E2 positive`
2. `E2 flat` 但训练稳定、指标接近噪声带上缘，仍值得看 distill 边际

如果 `E2` 明显 `negative`，`E3` 不做。

## 9. Stop-Rule

### 9.1 基础设施 stop

若 `Phase 0` 任何一步出现以下问题，则先修基础设施，不进入科学判读：

1. OOM
2. shape 不闭合
3. teacher cache 与 record 无法对齐
4. loss 明显病态

### 9.2 科学 stop

`E1`

- 若明显 `negative`，可以停在 `E1`

`E2`

- 即使 `E1 negative` 也照跑
- 若出现强负，则整条 AlphaGenome-like 拆因线直接停在 `E2`

强负定义：

1. `peak JSD` 恶化 `>= 0.01`
2. 且 `peak count Pearson` 下降 `>= 0.02`
3. 或 `gap/calibration` 明显坏掉

`E3`

- 只在 `E2 positive` 或 `E2 flat-but-informative` 时启动

### 9.3 Promotion 规则

1. `E1` 只有在其明显正、而 `E2` 不正时，才值得补更长 budget
2. `E2` 是唯一默认允许升到 teacher budget 的结构线
3. `E3` 只回答 distillation 的边际价值，不承担“证明整条路线最强”的任务

## 10. 实现切口

### 10.1 `E1`

优先走 config-level 新线：

- 新建 `4096` 数据/训练 config
- 保持当前 `corrected B` 模型配置语义

### 10.2 `E2`

需要新增一个真正的分层主干，而不是继续扩当前 readout family：

- 新的 hierarchical encoder-decoder model config
- 新的主干实现
- 尽量复用当前 bias branch、count head、训练/评估框架

### 10.3 `E3`

复用现有 distill 接口：

- `data.teacher_cache_dir`
- `data.teacher_target_names`
- `loss.distill_profile_weight`
- `loss.distill_count_weight`

但需要新增新的 model-teacher export 脚本，避免继续复用当前 `teacher_cache_export.py` 的 feature->ridge 代理路径。

## 11. 明确不做的实现

1. 第一轮不直接上 `8192/16384 + E2 + distill`
2. 第一轮不做 multimodal 大任务
3. 第一轮不改写当前 `bias-safe` contract
4. 第一轮不把 `count head` 和 decoder 强耦合

## 12. 设计结论

本设计把下一轮高价值实验明确限定为一条 **AlphaGenome-like factor ladder**：

- `E1`: 只测长上下文
- `E2`: 在同一长上下文上测全模型 encoder-decoder
- `E3`: 在 `E2` 上测 model-teacher distillation

若这条阶梯线给出清晰信号，后续再决定是否把胜者升到更长上下文或更强 budget；若在 `E2` 前后就连续给出强负结果，则应直接停止把 AlphaGenome-like 迁移当作当前主线高优先级方向。
