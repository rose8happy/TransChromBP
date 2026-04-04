# 咨询包 02：Genos 接入到底做了什么、跑成什么样、为何判定为高质量负结果

## 这份文件的用途

这份文件只回答一件事：

> `Genos-1.2B` 在我们这里不是“没跑通”，而是“已经做过一轮相对完整的接入、训练和评估，但没有形成正收益”。

外部分析应该基于这个事实继续往前推，而不是默认“也许只是环境坏了”。

---

## 1. 我们实际做过的 Genos 工作

### 1.1 独立 sanity check

先复现官方最接近开放染色质的 side benchmark：`human_ocr_ensembl`

结果：

| 配置 | 指标 | 数值 |
|---|---|---:|
| mini smoke, layer 6 | ROC AUC | `0.5958` |
| full, layer 6 | ROC AUC | `0.7242` |
| full, layer 12 | ROC AUC | `0.7535` |

官方参考值约为 `0.7569`，因此本地 `Genos-1.2B` 的权重、tokenizer、hidden-state 提取链路基本是正常的。

### 1.2 online 融合

核心实现文件：`vendor/transchrombp/transchrombp/models/genos_adapter.py`

已经实现的 online 方案：

| 方案 | 作用方式 | 关键配置 / 脚本 |
|---|---|---|
| `G0` | 无 Genos 的 matched baseline | `configs/model/v2fix_baseline.yaml` |
| `G1` | position-wise gated fusion | `configs/model/v2fix_genos_gate.yaml` |
| `G2` | global-mean gated fusion | `configs/model/v2fix_genos_mean.yaml` |

训练入口：

- `vendor/transchrombp/transchrombp/scripts/run_genos_pilot.sh`
- `vendor/transchrombp/transchrombp/configs/train/train_v2fix_genos_profile_select.yaml`

重要语义：

- online Genos 跑在 `genos-1.2b` 独立环境下
- `GenosFeatureExtractor` 冻结，不参与反传
- 双向特征来自 `forward + reverse-complement` 平均

### 1.3 cached fusion

已经实现的 cached 方案：

| 方案 | 作用方式 | 关键配置 / 脚本 |
|---|---|---|
| `P0` | cached recipe 的 matched baseline | `scripts/run_genos_cached_p0.sh` |
| `P2` | 把 `genos_global_mean` 注入 count hidden layer | `configs/model/v2fix_genos_global_count.yaml` |

缓存与 probe 入口：

- `vendor/transchrombp/transchrombp/scripts/build_genos_summary_cache.py`
- `vendor/transchrombp/transchrombp/scripts/run_genos_summary_probe.py`
- `vendor/transchrombp/transchrombp/scripts/run_genos_cached_pilot.sh`
- `vendor/transchrombp/transchrombp/configs/train/train_genos_cached_short10.yaml`

缓存中真正消耗的主特征是：

- `global_mean [1024]`

虽然也预留了 `bins4_mean`，但第一轮 matched cached 主线没有把它作为主实验推进。

### 1.4 probe

做过两类 probe：

1. `count residual ridge`
   - 问题：`genos_global_mean` 能否解释 baseline 剩余误差
2. `peak vs nonpeak logistic`
   - 问题：Genos 特征本身有无区分度，以及和现有 `encoded` 是否互补

---

## 2. 统一结果表

### 2.1 OCR sanity：说明 Genos 本体不是坏的

| 任务 | 结果 | 结论 |
|---|---:|---|
| `human_ocr_ensembl`, layer 12 | `AUC = 0.7535` | 接近官方 `0.7569`，本地链路正常 |

### 2.2 probe：说明 Genos 不是纯噪声，但也不是天然互补

| probe | 数值 | 解释 |
|---|---:|---|
| residual ridge `R^2` | `0.0167` | 只能解释很小一部分 count 残差 |
| residual ridge `Pearson r` | `0.1987` | 有弱相关，但不强 |
| `genos_only` AUC | `0.6997` | 有独立弱区分能力 |
| `encoded_only` AUC | `0.9219` | 当前主干已经很强 |
| `concat(encoded, genos)` AUC | `0.9018` | 拼接后反而下降，提示冲突或冗余 |

### 2.3 integration：说明目前没有正收益

| run | 方案 | best epoch | val peak JSD | val peak count_r | test peak JSD | test peak count_r | 结论 |
|---|---|---:|---:|---:|---:|---:|---|
| `genos_20260321_baseline_s42` | `G0` baseline | `20` | `0.3337` | `0.7998` | `0.3163` | `0.8410` | 稳定 baseline |
| `genos_20260321_gate_s42` | `G1` online gate | `10` | `0.3360` | `0.8010` | `0.3388` | `0.4689` | validation 误导，held-out count 崩 |
| `genos_20260322_mean_s42` | `G2` online mean | `5` | `0.3414` | `0.7955` | `0.3387` | `0.6360` | 比 `G1` 稍好，但仍明显失败 |
| `genos_cached_P0_baseline_ddp2_20260323_121650` | `P0` cached baseline | — | `0.3353` | `0.8033` | `0.3179` | `0.8384` | 稳定但未超过 `G0` |
| `genos_cached_P2_count_only_ddp2_20260323_173318` | `P2` cached count-only | `9` | `0.3344` | `0.7984` | `0.3174` | `0.6037` | profile 基本追回，但 count 断崖下滑 |

一句话总结：

- `G1/G2`：online 融合没有带来稳定收益
- `P0`：cached baseline 可复现，但没有超过 matched no-Genos
- `P2`：说明“把粗 summary 直接打进 count head”会明显破坏主任务

---

## 3. 目前最稳的解释

### 3.1 不是“Genos 没跑通”

我们已经排除：

1. 权重或 tokenizer 有问题
2. 本地推理链路坏了
3. 单次随机波动
4. 只有 online 坏、cached 一定能救

OCR reproduction 和 probe 都说明：Genos 本体有信号，但这个信号在当前任务里没有自然转化成净增益。

### 3.2 最可能的问题是粒度错位

当前 cached 主线主要消费的是：

- `global_mean [1024]`

这类序列级 summary 更像是给分类任务准备的。  
但我们做的是：

- base-resolution profile
- count regression
- 还要和 bias factorization 共存

因此最稳的解释是：

> `global_mean` 这种 sequence-level summary 与 base-resolution ATAC profile/count 的任务粒度不匹配。

### 3.3 当前 backbone 已经很强，边际空间很小

对 Genos 线最公平的 matched baseline `G0` 已经有：

- `0.3163 / 0.8410`

而全局更强的默认模型还在其之上。  
在这种前提下，一个只带来弱相关、且不天然互补的外部 summary，很难越过噪声门槛。

### 3.4 count head 对外部注入非常敏感

`P2` 的失败很典型：

- profile `JSD = 0.3174`，看起来还行
- count `r = 0.6037`，直接断崖

这说明当前任务里：

- profile 分支对外部弱信号还有一定鲁棒性
- count 分支会先被不自然的融合方式破坏

### 3.5 training semantics 与 cache semantics 有错位

训练数据里存在：

- `peak_max_jitter = 500`
- `train_revcomp = true`

而 cache 基本对应 canonical center。  
当外部 summary 本来就弱时，这种错位会进一步压缩净收益。

### 3.6 Genos 的 causal 架构可能也是问题之一

即使我们用了 `forward + reverse-complement` 平均，它仍只是对单向模型的近似补偿。  
对当前位置相关的局部调控特征，真正双向的模型可能更自然。

---

## 4. 这组负结果真正说明了什么

截至 2026-03-31，我们认为下面这些判断是稳的：

### 可以强说的

1. `Genos-1.2B` 本体和本地链路没有明显问题。
2. `Genos` hidden states 确实包含和开放染色质相关的信号。
3. 在当前 `TransChromBP + current fusion recipes` 下，Genos 没有带来正收益。
4. count 分支对粗粒度外部 summary 注入极度敏感。

### 只能弱说的

1. `task granularity mismatch` 是当前最合理解释。
2. 真正双向的 genomic foundation model 可能比 causal LM 更适合这个任务。
3. 如果改用位置分辨率特征而不是 `global_mean`，结果可能不同。

### 目前仍未测试的

1. foundation model 自身直接作为主干，而不是接入当前结构
2. 位置分辨率更强的 cross-attention / token-level fusion
3. truly bidirectional genomic foundation models
4. 除 Genos 外的其它 foundation models 在本任务上的真实对比

---

## 5. 外部分析时最容易犯的误读

### 误读 1

“既然 OCR 能跑到 `0.7535`，那接到 ATAC 主任务里应该也会有收益。”

问题：OCR 是序列级分类；当前任务是 base-resolution profile/count。

### 误读 2

“`P2` 的 profile JSD 没坏，所以只是 count head 小 bug。”

问题：`P2` 失败不是小 bug，而是非常像“粗外部 summary + 脆弱读出头”的结构性不匹配。

### 误读 3

“只要继续给 Genos 线加 seed，就可能翻案。”

问题：当前失败已经跨 `G1/G2/P0/P2` 多种 recipe，问题更像方向错位，不像单 seed 波动。

---

## 6. 关键源文件

- `reports/genos_integration_full_report_for_review_20260324.md`
- `reports/genos_no_positive_gain_analysis_20260323.md`
- `reports/genos_cached_fusion_status_20260323.md`
- `reports/genos_and_6002_run_analysis_20260323.md`
- `vendor/transchrombp/transchrombp/models/genos_adapter.py`
- `vendor/transchrombp/transchrombp/scripts/run_genos_pilot.sh`
- `vendor/transchrombp/transchrombp/scripts/run_genos_cached_pilot.sh`
- `vendor/transchrombp/transchrombp/scripts/run_genos_summary_probe.py`

---

## 7. 读完这份文件后你应该带着的问题

1. 如果继续保留当前 `TransChromBP` 结构，什么接法才可能比 `G1/G2/P2` 更自然？
2. 如果当前结构本身就不适合吃 foundation model，是否应该反过来让 foundation model 做主干？
3. 如果换模型，哪类模型最值得优先考虑：bidirectional、long-context、ATAC-specific、还是其它？
