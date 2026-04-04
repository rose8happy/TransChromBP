# 咨询包 03：现有 foundation model 资源状态与已经暴露出来的硬约束

## 这份文件的用途

外部 AI 在讨论“下一轮应该换什么模型、怎么接、要不要直接让 foundation model 做主干”之前，需要先知道：

1. 我们手上到底有哪些环境、模型和脚手架已经准备到什么程度；
2. 哪些约束是之前实验已经逼出来的，不能再当作可忽略细节。

---

## 1. 当前可用硬件与环境

### 1.1 6000 机器

- GPU：`NVIDIA RTX A6000 49140 MiB x2`
- 这是当前最主要的训练与本地 foundation model 推理资源

### 1.2 已明确分开的环境

| 环境 | 路径 | 主要用途 | 当前状态 |
|---|---|---|---|
| `chrombpnet` | `/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet` | 官方 ChromBPNet / TensorFlow / `bedGraphToBigWig` / MoDISco | 可用 |
| `transchrombp` | `/data1/zhoujiazhen/bylw_atac/.mamba/envs/transchrombp` | 当前 PyTorch 主训练、评估、benchmark | 可用 |
| `genos-1.2b` | `/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b` | Genos-1.2B 加载、GPU 推理、embedding 提取 | 已验证可用 |
| `alphagenome` | `/data1/zhoujiazhen/bylw_atac/.mamba/envs/alphagenome` | AlphaGenome SDK / API pilot | 已创建 |

一个重要现实是：

- 我们没有把这些依赖强行揉进一个公共环境
- 之前的经验也不支持这样做

---

## 2. foundation model 相关资源的真实状态

### 2.1 Genos-1.2B

当前状态是 **最完整** 的：

- 本地模型目录存在：`/data1/zhoujiazhen/bylw_atac/foundation_models/genos/Genos-1.2B`
- 已验证：
  - `AutoTokenizer / AutoModelForCausalLM.from_pretrained(..., local_files_only=True)`
  - `flash_attention_2` GPU 最小前向
  - hidden-state 提取
  - next-token 预测
  - `generate`
  - OCR benchmark reproduction

这意味着：

- Genos 线的负结果不是“下载都没下好”
- 外部分析应把它视为一条已经真正做过的方向

### 2.2 Nucleotide Transformer v2-500M multi-species

当前状态是 **部分就绪**：

- 目录：`/data1/zhoujiazhen/bylw_atac/foundation_models/nucleotide_transformer/nt-v2-500m-multi-species/`
- 已完成：
  - 小文件到位
  - `model.safetensors` 已落盘
  - 本地 Hugging Face smoke test 已通过
- 当前未补齐：
  - `pytorch_model.bin`
  - `jax_model/pytree_ckpt.joblib`

这说明：

- 它已经不是“纯想法”，而是接近可用
- 但严格说还没达到像 Genos 那样“从环境到任务 probe 都走通”的成熟度

### 2.3 AlphaGenome

当前状态是 **只适合 black-box / zero-shot pilot**：

- 代码脚手架已存在：`scripts/alphagenome_pilot/`
- 专用环境已建
- 已知边界：
  - 适合外部对比
  - 不适合把 API 输出直接拿来训练学生模型
  - 当前只有极小规模 pilot，不是大规模 held-out benchmark

已知 pilot 级事实：

- 在 4 个代表位点上，AlphaGenome 对绝对总量明显低估，典型案例约 `1168 vs 5346`，约 `4.6x` 偏低
- 但相对强弱排序是对的

这意味着：

- AlphaGenome 更像“上界参考 / 外部理念对照”
- 不是现成可并入当前训练线的监督来源

### 2.4 Gengram-10B

当前状态是 **仅有目录、未验证**：

- 路径：`/data1/zhoujiazhen/bylw_atac/foundation_models/gengram/`
- 目前没有本地 smoke、任务 probe 或训练接入记录

### 2.5 其它候选模型

在这份仓库里，下面这些方向目前更多停留在“讨论对象”，而不是“已经准备好的本地资产”：

- `DNABERT-2`
- `HyenaDNA`
- `Caduceus`
- `Evo`

也就是说，如果外部分析建议优先转向它们，需要把“本地可达性 / 许可证 / 推理成本 / 输入粒度”一起纳入判断，而不是只看论文名气。

---

## 3. 之前实验已经暴露出来的硬约束

### 3.1 当前 baseline 很强

这不是在一个弱 baseline 上“随便加点 foundation model 就可能涨点分”的局面。

当前至少有三层强约束：

1. `V2-full` 已经稳定优于 `V2-noTF`
2. `corrected B` 两个 held-out seed 已经稳定
3. shared-region `L3` 上已经优于 official ChromBPNet

因此，新路线必须面对的是一个成熟任务专用系统，而不是空白地带。

### 3.2 这个任务对位置分辨率要求高

之前最清晰的教训就是：

- 序列级 summary 可以对 OCR 分类有用
- 但对 base-resolution ATAC profile/count 未必有用

所以任何新 foundation model 方向，如果只能稳定给出 sequence-level embedding，而给不出足够自然的 token/base-level 表征，就有很大风险重演 Genos 的问题。

### 3.3 count head 特别脆弱

`P2` 给出的事实已经很明确：

- profile 没坏很多，不代表路线可用
- count 往往是最先暴露结构不匹配的分支

这对外部分析很重要，因为很多“看起来轻量优雅”的接法，可能首先在 count 分支上崩。

### 3.4 validation 不能当最终 gate

Genos `G1/G2` 已经说明：

- validation 可能看起来接近 baseline
- held-out test 才暴露严重的 count collapse

所以任何 foundation model 方向，如果只基于 validation 讲话，都不够可靠。

### 3.5 训练拓扑不能成为暗变量

我们已经明确收口的一条规则是：

- 如果 `P0/P2/P1` 设计成 matched `2-rank DDP`
- 那就不能中途切成单卡 fallback 再拿来横比

这条规则不只适用于 Genos，也适用于后续任何 foundation model 方向。

### 3.6 cache 语义与训练增强存在错位

只要训练里继续存在：

- `jitter`
- `revcomp`

那么离线 cache 和真实训练样本语义就不是完全一致。  
这个问题在外部信号本来就弱时会被放大。

### 3.7 causal vs bidirectional 不是小差别

Genos 作为 causal LM，已经让我们吃过一次亏：

- 即使做 `forward + RC` 平均，仍然只是一种近似补偿

对于当前位置相关的 profile 任务，这个差异可能是结构性的，不是实现细节。

---

## 4. 对下一轮候选模型最实际的筛选标准

这不是实验计划，而是基于已知事实形成的现实筛选框架。

### 4.1 表征层面

优先问：

1. 模型能否提供稳定的 **位置分辨率特征**？
2. 特征是否天然更偏 **双向上下文**，而不是单向语言建模？
3. token/base 粒度是否适合 2114bp 窗口内的 motif 组合与局部结构？

### 4.2 工程层面

优先问：

1. 能否在单张或两张 A6000 上稳定推理或训练？
2. 本地下载、权重许可、依赖链是否可控？
3. 是否需要 API，还是可完全本地运行？

### 4.3 任务层面

优先问：

1. 它是否已经在接近 ATAC / chromatin accessibility 的任务上显示过非分类价值？
2. 它是否只能做 zero-shot black-box 对比，还是可以作为真正可训练的主干 / 特征源？
3. 它的“强”到底是体现在分类、长程上下文、还是真实 profile/count 建模上？

---

## 5. 这份文件希望外部 AI 真正吸收的结论

1. 我们不是没有 foundation model 资源，而是已经至少把 Genos 做到了足以出负结果的成熟度。
2. `NT v2` 和 `AlphaGenome` 都不是空白词条，但准备度和可用场景完全不同。
3. 这个任务真正难的不是“把大模型塞进来”，而是让它在 **强 baseline + 位置分辨率 + count 敏感 + bias-safe** 的约束下形成净收益。
4. 如果外部 AI 建议换模型或换路线，它需要同时回答“为什么这个候选能避开上述约束”，否则建议就不够落地。
