# ChromBPNet 官方基线 vs TransChromBP 最优模型 分层对比方案（2026-03-27）

## 1. 为什么现在值得做这件事

当前项目已经基本完成两类内部判断：

1. `TransChromBP` 在我们自己的 tutorial / held-out 体系里优于 `noTF` 对照；
2. clean-matrix 已经表明“Transformer 特有 shortcut”这条强叙事不成立。

因此，下一步最值钱的问题不再是继续扩内部支线，而是：

> 在 **ChromBPNet 官方数据 + 官方代码** 的口径下，若把原始输入、硬件、训练预算与 checkpoint 选择规则尽量统一，并把 `native preprocessing` 与 `shared-region compare` 分层处理，我们当前最优模型相对官方基线到底还有多少真实增益？

这件事的价值也应该分三层看：

- `Layer 1`：把“内部自洽”升级成“对外可对照”的官方锚点。
- `Layer 2`：把“我们模型比自己的历史版本强”升级成“在共享原始输入与匹配训练预算下，完整系统到底强多少”。
- `Layer 3`：在需要更强归因时，再进一步回答“优势能否在共享 candidate regions 下仍然成立”。

---

## 2. 这项对比要回答的核心问题

### Q1. 官方复现 fidelity

在不引入我们自己的额外改造时，使用 ChromBPNet 官方数据和官方代码，能否在我们机器上复现出与官方 tutorial 接近的结果？

### Q2. 公平硬件对比

在 **同一批数据文件、同一台机器、同一 GPU 拓扑、同一 global batch、同一 seed 集合** 下，官方 ChromBPNet 和我们当前最优的 `corrected B` 到底谁更强？

### Q3. 差异来源定位

如果我们的模型更强，这个优势究竟来自：

- Transformer + local tower 的 backbone；
- `center pool` readout；
- `full/debiased` 训练/诊断设计；
- 还是仅仅来自不公平的训练预算/数据预处理差异？

---

## 3. 总体原则：先分清比较层级，再安排执行顺序

这件事不能一步做成一个实验，因为 `official fidelity`、`matched-budget system compare`、`shared-region compare` 回答的是三种不同问题。

### 3.1 目标 A：忠实复现（fidelity）

这里的目标是：

- 尽量少改官方流程；
- 确认官方 tutorial / paper-aligned 流程在我们机器上是通的；
- 给出一个“外部参考锚点”。

这一步允许和后续 `Layer 2/3` 对比在部分变量上不同，例如：

- 官方默认 batch size；
- 官方默认 early stop / epoch 上限；
- 官方默认评估脚本。

因为它的目的不是公平，而是确认“官方口径在我们这里能复现到什么程度”。

### 3.2 目标 B：matched-budget system compare

这里的目标是：

- 固定原始输入与硬件；
- 固定训练预算和 checkpoint 选择规则；
- 让“完整系统在相同预算下谁更强”成为主要问题。

这一步不是“region-identical 的严格归因”，而是“在共享原始输入上做 matched-budget 的 apples-to-apples system compare”。

在这一层里，可以允许：

- 两边保留各自的 native preprocessing；
- 两边保留各自的 optimizer / scheduler / precision；
- 但必须把这些差异写清楚，不能再表述成“只差 architecture”。

### 3.3 目标 C：shared-region compare

这里的目标是：

- 在 `Layer 2` 之上，再统一最终进入 trainer 的 candidate peaks / nonpeaks；
- 尽量统一 bigWig 或其生成链路；
- 把可疑混杂进一步压缩到更小范围。

这一层才更接近“强归因”实验。  
它不是当前 tutorial 主线必须先完成的前置条件，但如果后续想把结论写成“更接近架构/recipe 归因”，就必须进入这一层。

---

## 4. 推荐的实验分层

我建议拆成三层，而不是一口气上最终矩阵。

### Phase 0：资产冻结（不训练）

先把后续所有 run 共用的输入冻结下来，避免边跑边改。

固定资产分两层：

- `Layer 1/2` 共享原始资产：
  - `hg38.fa`
  - `hg38.chrom.sizes`
  - `blacklist.bed.gz`
  - `merged.bam`
  - `overlap.bed.gz`
- 固定 split：
  - 第一轮优先用官方 tutorial 单 fold
  - 第二轮如有必要扩到 `paper_aligned_repro` 的 5-fold
- `Layer 2` 统一 protocol：
  - 统一的 `folds.json` 或 fold set
  - 统一的 evaluator / selector
  - 统一的 run matrix、命名和日志落盘规范
- `Layer 3` 额外冻结的共享 candidate regions：
  - identical `peaks.bed`
  - identical `nonpeaks.bed`
  - 若可能，identical `merged_unstranded.bw`

这一步的关键不是“快”，而是避免把 `Layer 2` 的 system-level compare 误写成 `Layer 3` 的 shared-region compare。

### Phase 1：官方 ChromBPNet 忠实复现

用官方代码 + 官方数据先做一轮 fidelity baseline。

建议先只做：

- 单 fold
- 单 seed
- 官方默认/近官方默认训练配置

目的：

- 看我们当前环境能否稳定复现官方 tutorial 的数量级；
- 确认官方评估口径输出链路无误；
- 不要一开始就上 3 seed × 5 fold，把问题复杂化。

### Phase 2：matched-budget system compare

在共享原始输入、相同硬件和相同训练预算下，做三臂对照：

1. `ChromBPNet-controlled`
2. `TransChromBP-corrected-B-controlled`
3. 可选：`TransChromBP-noTF-controlled`

其中前两条是硬主线，第三条只在你想进一步定位“收益来自 Transformer 还是其他 recipe”时再加。  
这一 phase 默认仍属于 `Layer 2`，不自动声称 `peaks/nonpeaks` 已完全一致。

### Phase 3：必要时扩到 3 seed / 5 fold

只有在 Phase 2 的单 fold 结果有意义时，再扩：

- 3 seed
- 5 fold

否则容易在大矩阵上浪费大量算力，却没有先证明比较口径本身是干净的。

---

## 5. 不同层级下，哪些变量必须锁死

下面这些变量需要按层级区分。

| 维度 | `Layer 2` 要求 | `Layer 3` 要求 | 说明 |
|---|---|---|---|
| 原始数据 | 完全相同 | 完全相同 | 同一份 BAM / peaks / blacklist / genome |
| split | 完全相同 | 完全相同 | 同一份 folds/fold json |
| candidate peaks | 允许 native preprocessing，但必须显式记录 | 完全相同 | 当前 tutorial 主线已确认两侧 peaks 不同 |
| candidate nonpeaks | 允许 native preprocessing，但必须显式记录 | 完全相同 | 当前 tutorial 主线已确认两侧 nonpeaks 差异更大 |
| bigWig / signal 资产 | 最好同源，至少要记录生成链路 | 尽量完全相同 | 这是 preprocessing 混杂的重要来源 |
| 硬件 | 完全相同 | 完全相同 | 都在 6000，2×A6000 |
| GPU 拓扑 | 完全相同 | 完全相同 | 都走双卡，避免单卡/双卡差异混入 |
| global batch | 完全相同 | 完全相同 | 推荐统一成双方都稳定可跑的值 |
| seed 集合 | 完全相同 | 完全相同 | 例如 `42 / 1234 / 2024` |
| epoch 上限 | 完全相同 | 完全相同 | 例如都给 `50` |
| optimizer / scheduler / precision | 可以不同，但必须显式标注为 system-level compare | 若要做更强归因，需统一或加 matched noTF 对照 | 否则不能把结果写成“只差架构” |
| checkpoint 选择规则 | 完全相同 | 完全相同 | 这是最容易被忽略、但最重要的变量之一 |
| 最终评估脚本 | 完全相同 | 完全相同 | 不允许 baseline 用一套 evaluator，我们模型用另一套 |

### 5.1 我最强调的一点：checkpoint 选择规则必须统一

如果一边按：

- 官方默认 early-stop 逻辑选 best

另一边按：

- 我们自己的 `peak.profile_target_jsd_full_mean` 选 best

那么这个比较天然不干净。

更严格的做法是：

1. 两边都保存每个 epoch checkpoint；
2. 用同一个外部 evaluator 在相同 validation 集上逐 epoch 打分；
3. 用同一条规则选 best epoch。

推荐统一的 best rule：

- 主规则：`peak profile JSD median`
- 同时报告：`peak profile JSD mean`、`peak count Pearson r`

这样更接近官方 tutorial 的传统口径，同时仍保留我们主线需要的 `mean` 与 `count r`。

---

## 6. 推荐的三臂实验矩阵

### Arm A：ChromBPNet-official-fidelity

目的：

- 只回答“官方流程在我们这里是否复现”

特点：

- 官方代码
- 官方数据
- 尽量少改默认训练流程
- 不把它直接拿来和我们模型做最终公平比较

### Arm B：ChromBPNet matched-budget system baseline

目的：

- 作为公平对照基线

特点：

- 仍用官方代码
- 但数据 split / nonpeak / batch / seed / 硬件 / best 规则都对齐到 controlled protocol

### Arm C：TransChromBP matched-budget system candidate

目的：

- 我们当前最优模型的公平主结果

当前候选：

- `corrected B = teacher_v2_center_pool + sg=true + deb2`

特点：

- 共享原始输入
- 同一 global batch
- 同一 seed
- 同一 checkpoint 选择规则
- 同一 evaluator
- 当前默认仍允许 native preprocessing

### 可选 Arm D：TransChromBP-noTF-controlled

目的：

- 如果 Arm C 比 Arm B 好，想进一步确认优势主要来自 Transformer 还是其他 recipe
- 这条臂应尽量沿用与 Arm C 相同的 optimizer / scheduler / precision，才有更强归因价值

这条不是第一优先级，但如果你后面准备写 paper 的核心 comparative table，它会很有解释价值。

---

## 7. 指标设计：主指标和诊断指标分开

这件事不能只报一个 JSD。

### 7.1 主比较指标

建议统一为：

1. `peak profile JSD mean`
2. `peak profile JSD median`
3. `peak count Pearson r`

原因：

- `median JSD` 方便和官方 tutorial 口径对齐；
- `mean JSD` 方便和我们现有论文主线对齐；
- `count r` 是当前我们模型最明显的优势来源之一。

### 7.2 诊断指标

对 TransChromBP 还应额外保留：

1. `profile_full_debiased_jsd`
2. `effective_profile_scale`
3. `profile_bias_rms_over_signal_rms`

这三项不一定能直接作为与 ChromBPNet 的一对一公平主指标，但它们能回答：

- 我们的提升是不是建立在更高 bias reliance 上

如果后续能从 ChromBPNet 的 `chrombpnet.h5` / `chrombpnet_nobias.h5` 提出对应 full/debiased 对比，则这部分还可以升级成跨模型诊断比较；但在 Phase 1/2 不必把这件事当硬 gate。

---

## 8. 对“相同 batch size”的具体建议

你的方向是对的，但这里要注意一句：

> 真正应该统一的是 `global batch semantics`，而不是只写“batch size 相同”。

因为两边可能存在：

- 单卡 batch
- 双卡 batch
- grad accumulation

如果只看 `batch_size_per_gpu`，但一边有 `grad_accum=2` 另一边没有，最终 global batch 仍然不同。

### 我建议

统一成：

- 双卡
- 同一个 `global batch`
- 尽量不用 grad accumulation

如果双方都能稳跑，优先考虑：

- `global batch = 32` 或 `64`

选择标准不是“谁更快”，而是：

1. 两边都不 OOM；
2. 两边都能稳定训练；
3. 这个 batch 不明显偏袒某一方。

---

## 9. 我最推荐的执行顺序

### Step 1. 先跑一个最小 fidelity baseline

只做：

- ChromBPNet 官方代码
- 官方 tutorial 数据
- 单 fold
- 单 seed

目标：

- 拿到一个“官方在我机上”的参考点

### Step 2. 再做 controlled single-fold 对比

只做：

- `ChromBPNet-controlled`
- `TransChromBP-controlled`
- 同一 fold
- 同一 seed（例如 `42`）
- 同一 global batch
- 同一 best rule

目标：

- 证明这套 controlled protocol 本身是可跑的、可解释的

### Step 3. 决定是否扩 3 seed

如果 Step 2 已经显示：

- 我们模型显著更好，或者
- 差距比预想小很多

再扩到 `42 / 1234 / 2024`。

### Step 4. 最后才考虑 5 fold

5 fold 是加固，不是第一轮 gate。

如果第一轮单 fold 都还没说明白，先上 5 fold 会把调试和归因成本拉得过高。

---

## 10. 风险点与预防

### 风险 1：表面同数据，实际 preprocessing 不同

例如：

- ChromBPNet 直接从 BAM 取数
- 我们从另一版 bigWig 取数

这时即使原始 BAM 相同，训练信号也可能已有差异。

预防：

- 记录 shared preprocessing 资产的生成方式和路径；
- 尽量让 bigWig 来自官方 tutorial BAM 的单一 canonical 生成链路。

### 风险 2：best checkpoint 选法不同

这是最常见也最致命的“不公平”来源。

预防：

- 用统一 evaluator 逐 epoch 打 validation；
- 再统一选 best。

### 风险 3：官方 baseline 只报 median，我们主线报 mean

预防：

- 一次性同时输出 `mean + median + count_r`；
- 不再允许正文里不同口径混写。

### 风险 4：硬件相同，但吞吐拓扑不同

例如：

- 一边双卡 DDP
- 一边单卡或 DataParallel

预防：

- 所有 controlled run 强制双卡 DDP；
- 明确记录 world size、bs/gpu、global batch、grad accumulation。

---

## 11. 为什么还需要两条额外数据集线

只做官方 tutorial 仍然有一个明显风险：

> 即使我们的模型在官方数据上赢了，也不能完全排除这只是某个特定数据集、特定细胞系、特定预处理口径下的偏好。

因此，建议把最终比较结构从“单一官方 benchmark”升级成：

1. **官方 tutorial 锚点**
   用来证明我们没有偏离社区熟悉的参考口径。
2. **两条额外大数据集的 controlled compare**
   用来证明结论不依赖某一个 benchmark 的偶然性。

这两条额外数据集的任务不是“再讲一个泛化故事”，而是：

- 检查 `ChromBPNet-controlled` vs `TransChromBP-controlled` 的优劣方向是否稳定；
- 看我们在 tutorial 上观察到的优势，是否在独立数据上仍然成立；
- 避免 paper 中被质疑“你们只是挑了一个对自己有利的官方数据”。

### 11.1 这两条额外数据集怎么选

我的建议不是“随便再加两个”，而是按下面规则选：

1. **优先大数据集**
   尽量选择 merged 后测序深度高、峰数充足、重复数足够的数据集，避免比较被小数据高方差主导。
2. **不要两条都和 tutorial 太相似**
   如果当前 tutorial 主线本质上接近 `K562`，那两条额外数据集里不应再把第二条也放成 K562-like。
3. **至少覆盖两个不同生物背景**
   例如一条血液/免疫系（`GM12878`），一条非血液系（如肝/肺/上皮系的高深度数据）。
4. **优先选我们已有缓存或容易补齐的**
   这样能把比较重点放在模型，而不是把大量时间花在新数据清洗上。

### 11.2 当前最推荐的候选

按仓库当前资产，最自然的组合是：

1. **GM12878**
   这是当前最稳的第一条外部大数据集候选，因为本地已缓存，且与 tutorial/K562 不同。
2. **再新增一条高深度、非 K562 的 ENCODE ATAC 数据**
   候选方向可从 `HepG2 / A549 / H1-hESC` 这类与 GM12878 差异较大的体系中挑一条，最终按：
   - read depth
   - replicate 数
   - peaks / BAM 完整性
   - 我们当前预处理链路是否容易复用
   来定。

### 11.3 我不建议的做法

- **不要把“独立 K562”当成两条额外数据集之一的主候选**
  因为 tutorial 本身已经很接近 K562 语义，这样会削弱“避免官方数据偏好”的价值。
- **不要一开始就 3 个数据集同时 3-seed × 5-fold**
  这会在比较口径尚未冻结前，把算力和调试成本拉爆。

---

## 12. 扩展后的实验结构

### Track A：官方锚点

- `ChromBPNet-official-fidelity` on tutorial
- `ChromBPNet-controlled` vs `TransChromBP-controlled` on tutorial

### Track B：外部数据集 1

- `ChromBPNet-controlled` vs `TransChromBP-controlled` on `GM12878`

### Track C：外部数据集 2

- `ChromBPNet-controlled` vs `TransChromBP-controlled` on `Dataset-X`
  - `Dataset-X` = 一条高深度、非 K562、非 GM12878 的 ATAC 数据

这三条 track 里：

- Track A 给“社区可对照锚点”
- Track B/C 给“benchmark 偏好之外的稳定性”

---

## 13. 数据集扩展后的执行顺序

### Step 1. 只在 tutorial 上冻结 protocol

先把下面几件事定死：

- shared raw assets
- unified evaluator
- unified best rule
- unified global batch

### Step 2. tutorial 单 seed 双臂对比

如果 tutorial 的 controlled compare 都还没跑通，就不要急着扩外部数据。

### Step 3. 先加 GM12878

原因：

- 本地已有缓存；
- 与 tutorial 差异足够大；
- 预处理成本最低。

### Step 4. 最后再加 Dataset-X

等前两步已经说明 protocol 稳定后，再补第二条额外大数据集。

这样可以避免出现：

- 数据集选了 3 条；
- 每条的 preprocessing / evaluator / best rule 都还有细节没冻住；
- 最后谁赢谁输反而讲不清。

---

## 14. 扩展后的结果汇总原则

有了额外两条数据集后，主表不能直接横着比较绝对 JSD，因为不同数据集难度不同。

更稳的汇总方式是：

1. **先做数据集内比较**
   每个数据集内部比较 `ChromBPNet-controlled` vs `TransChromBP-controlled`
2. **再做跨数据集方向一致性总结**
   看优势方向是否稳定，而不是看绝对值谁最小

推荐把结果汇总成两层：

### 层 1：每个数据集单独表

| Dataset | Model | Peak JSD mean | Peak JSD median | Peak count r | Notes |
|---|---|---:|---:|---:|---|

### 层 2：跨数据集 delta 表

| Dataset | Δ JSD mean (ours - baseline) | Δ JSD median | Δ count r | 方向是否一致 |
|---|---:|---:|---:|---|

这样更能回答：

> 我们的模型是否在多个独立数据集上，稳定优于受控 ChromBPNet 基线？

---

## 15. 我认为最稳的第一轮最小方案

如果你现在只想做一个“先别跑大矩阵、但足够干净”的版本，我建议：

### Round 1A：官方 fidelity

- `ChromBPNet official`
- 官方 tutorial 数据
- 单 fold
- seed=42

### Round 1B：matched-budget system compare

- `ChromBPNet-controlled`
- `TransChromBP corrected-B-controlled`
- 同一 fold
- 同一 seed=42
- 同一 `global batch`
- 同一双卡 A6000
- 同一 validation 选 best 规则
- 同一 test evaluator
- 允许 native preprocessing，但必须在表注中写清楚

如果这轮结果清晰，再扩到 3 seed。

### Round 1C：外部数据集加固

在 Round 1B 跑通后，优先补：

1. `GM12878` 单 seed controlled compare
2. `Dataset-X` 单 seed controlled compare

只有当 tutorial + 这两条外部数据集的方向基本一致时，才值得进一步扩到 3 seed 或 5 fold。

### Round 2：shared-region compare（按需启动）

如果 paper 需要更强的归因强度，再单独追加：

- 统一 tutorial 的最终 `peaks/nonpeaks`
- 尽量统一 bigWig 或其生成链路
- 复跑 `ChromBPNet` 与 `TransChromBP corrected-B`
- 需要时再补 `TransChromBP-noTF-controlled`

---

## 16. 这件事完成后，最有价值的产出应该是什么

不是只多两个 checkpoint，而是多出一张可以直接进 paper / rebuttal 的表。

理想产出表应包含：

| Model | Code base | Data source | Hardware | Global batch | Selection rule | Peak JSD mean | Peak JSD median | Peak count r |
|---|---|---|---|---|---|---:|---:|---:|
| ChromBPNet official fidelity | official | official tutorial | 2xA6000 | official-like | official-like | ... | ... | ... |
| ChromBPNet matched-budget | official | shared raw inputs + native preprocessing | 2xA6000 | fixed | unified external rule | ... | ... | ... |
| TransChromBP corrected B matched-budget | ours | shared raw inputs + native preprocessing | 2xA6000 | fixed | unified external rule | ... | ... | ... |
| Optional shared-region rows | official / ours | identical peaks + nonpeaks + split | 2xA6000 | fixed | unified external rule | ... | ... | ... |

如果后续再加：

- `TransChromBP-noTF-controlled`

就能把“是不是 Transformer 带来的优势”也顺手回答掉。

---

## 17. 当前建议

1. 先不要一上来就跑 5-fold × 3-seed。
2. 先把“官方 fidelity”与“公平 controlled comparison”拆开。
3. 把最终结构理解成三层：
   - `Layer 1`：官方 tutorial 锚点
   - `Layer 2`：shared raw inputs + matched-budget system compare
   - `Layer 3`：shared-region compare
4. `Layer 2` 里最关键的不是“都是 A6000”，而是：
   - 同一原始输入
   - 同一 global batch
   - 同一 checkpoint 选择规则
   - 同一 evaluator
   - 同时诚实记录 native preprocessing 与优化差异
5. 第一轮最稳妥的起点是：
   - 官方 baseline 单 seed
   - controlled 双臂单 seed
6. 在 tutorial protocol 冻住后，优先把同一 protocol 平移到 `GM12878 + Dataset-X`。
7. 如果 tutorial + 两条外部数据集的单 seed 方向都一致，再扩种子和 fold。

这套顺序的好处是：先证明比较口径干净，再花大算力做统计加固。
