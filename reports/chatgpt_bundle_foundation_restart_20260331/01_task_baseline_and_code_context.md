# 咨询包 01：任务、当前强基线与代码入口

## 这份文件的用途

这不是论文定位材料，而是给外部 AI 做下一轮 `foundation model for ATAC` 研究判断的事实底稿。

你需要先理解三件事：

1. 我们要做的不是普通序列分类，而是 **base-resolution ATAC profile/count prediction**。
2. 当前内部 baseline 并不弱，任何新的 foundation model 方向都必须面对一个已经很强的任务专用模型。
3. 这份仓库更像“项目总档案”，实际 `TransChromBP` 代码镜像主要在 `vendor/transchrombp/transchrombp/`。

---

## 1. 任务定义

- 输入：`2114 bp` DNA 序列，one-hot 编码，形状约为 `[B, 2114, 4]`
- 输出 1：中心 `1000 bp` 的 base-resolution ATAC profile
- 输出 2：该区域的总 read count（logcount / count）
- 主要损失：
  - profile：multinomial NLL
  - count：MSE
- 主要指标：
  - profile：`JSD`，越低越好
  - count：`Pearson r`，越高越好

这意味着它和 `OCR / peak-vs-nonpeak` 这类序列级分类任务不是一回事。  
对这个任务有效的特征必须足够接近 **位置分辨率**，并且不能只在分类上有效。

---

## 2. 当前强基线有多强

### 2.1 Backbone 本身已经被验证有效

`V2-full` vs `V2-noTF` 是 matched ablation（3 seeds）：

| 配置 | peak JSD | peak count_r | 结论 |
|---|---:|---:|---|
| `V2-full` | `0.31419 ± 0.00012` | `0.83676 ± 0.00531` | 当前 backbone 主线 |
| `V2-noTF` | `0.32393 ± 0.00021` | `0.82503 ± 0.00489` | 稳定弱于 `V2-full` |

结论：当前任务里，Transformer 主干带来的收益是真实存在的，不是偶然波动。

### 2.2 当前默认模型是 `corrected B = center pool + sg=true + deb2`

这是当前 paper-facing 默认模型，也是本轮 foundation-model 咨询时真正要面对的内部强基线。

代表性 held-out 结果：

| run | best epoch | peak JSD | peak count_r | 结论 |
|---|---:|---:|---:|---|
| `B_s42 = v2fix_20260320_cpool_s42` | `34` | `0.3147` | `0.8503` | 当前单 seed 最优 |
| `B_s1234 = v2fix_20260324_005834_cpool_s1234_gpu0` | `39` | `0.3145` | `0.8488` | 与 `B_s42` 同档 |

clean-matrix 里对应的安全性检查：

| 配置 | peak JSD full | peak JSD debiased | gap | peak count_r full/deb |
|---|---:|---:|---:|---:|
| `corrected B = TF + center pool + sg=true + deb2` | `0.42544` | `0.42556` | `0.00268` | `0.8439 / 0.8420` |

结论：当前默认模型不仅强，而且没有显示出明显的 bias reliance。

### 2.3 外部 shared-region 对比也已经成立

在 tutorial `L3 shared-region` 的 controlled compare 上：

| Arm | best epoch | peak mean JSD | peak count_r | AUROC | 结论 |
|---|---:|---:|---:|---:|---|
| Official ChromBPNet controlled L3 | `29` | `0.33853` | `0.69958` | `0.82295` | 官方对照 |
| TransChromBP corrected-B controlled L3 | `47` | `0.31319` | `0.84016` | `0.87861` | 我们当前系统更强 |

这说明 foundation model 方向想要成立，不能只做到“比 Genos line baseline 好一点”，而要面对一个已经在 shared-region external compare 中表现很强的系统。

### 2.4 Genos 线自己的 `G0` 只是 matched baseline，不是全局最强 baseline

Genos integration 报告里常出现：

- `G0 = 0.3163 / 0.8410`

它的作用是作为 **Genos 线内部 matched recipe baseline**。  
它是有效基线，但不是当前全项目意义上的最强默认模型。  
外部分析时不要把 `G0` 误当作我们唯一 baseline。

---

## 3. 当前默认模型最小配置

核心模型配置：`vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool.yaml`

```yaml
sequence_encoder:
  type: transformer
  enabled: true
  d_model: 256
  n_heads: 8
  n_layers: 6

local_tower:
  n_dil_layers: 8

bias_branch:
  enabled: true
  profile_pool_factor: 32

fusion:
  profile_fusion: add
  count_fusion: logsumexp
  profile_bias_stop_gradient: true

heads:
  profile_output_len: 1000
  count_pool_mode: center
```

代表性训练配置：`vendor/transchrombp/transchrombp/configs/train/train_v2fix_cpool_s1234_6000_single.yaml`

```yaml
seed: 1234
max_epochs: 40

loss:
  profile_weight: 1.0
  count_weight: 0.1
  debiased_profile_weight: 2.0
  debiased_count_weight: 0.0

trainer:
  best_metric: peak.profile_target_jsd_full_mean
  best_metric_mode: min

data:
  input_len: 2114
  output_len: 1000
  batch_size_per_gpu: 16
  peak_max_jitter: 500
  train_revcomp: true
```

对 foundation model 方向最重要的语义只有四个：

1. 当前模型不是纯 Transformer，而是 `conv stem + local tower + transformer + bias branch`
2. `debiased_profile_weight = 2.0` 是当前稳定性的关键变量之一
3. `count_pool_mode = center` 是当前最稳的 count readout
4. `best checkpoint` 默认按 `peak.profile_target_jsd_full_mean` 选

---

## 4. 关键代码入口

如果你想判断“新 foundation model 应该插在哪里”或“是否值得推翻现有结构”，下面这些入口最关键：

- `vendor/transchrombp/transchrombp/models/transchrombp.py`
  - 主模型定义
  - `full / debiased` 输出语义
  - bias 融合
  - count head 的 `full / center / attention` pooling
- `vendor/transchrombp/transchrombp/models/transformer_encoder.py`
  - 当前序列编码器本体
- `vendor/transchrombp/transchrombp/models/bias_branch.py`
  - ChromBPNet-style bias branch
- `vendor/transchrombp/transchrombp/training/train_ddp.py`
  - 主训练循环
  - DDP 语义
  - validation / checkpoint / logging
- `vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py`
  - held-out test 评估入口
- `vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool.yaml`
  - 当前默认模型配置
- `vendor/transchrombp/transchrombp/configs/train/train_v2fix_cpool_s1234_6000_single.yaml`
  - 代表性训练配置

---

## 5. 当前训练与评估语义

### 5.1 `full` 与 `debiased`

当前模型同时保留：

- `full`：`signal + scale * bias`
- `debiased`：`signal` only

这不是装饰性指标。  
它直接决定我们是否认为某个新分支是在学调控信号，还是在引入新的 shortcut / bias reliance。

### 5.2 held-out test 比 validation 更重要

Genos 的 `G1/G2` 已经证明：

- validation 看起来接近 baseline
- held-out test 上 count 可以直接崩掉

所以对于 foundation model 方向，不能只看 validation。

### 5.3 matched recipe 很重要

我们已经踩过一个坑：训练拓扑一变，比较就不再干净。

因此外部分析时请默认：

- 训练拓扑不能偷偷从 `2-rank DDP` 改成单卡再拿来横比
- global batch 不能随意漂
- cached / online / no-FM / with-FM 最好维持 matched recipe

---

## 6. 对外部分析最关键的约束

如果你后面要判断“下一轮 foundation model 线应该怎么走”，请把下面这些约束当真：

1. 我们已经有一个很强的任务专用 baseline，不是从零开始。
2. 这个任务不是 sequence-level classification，而是 base-resolution profile/count。
3. count 分支很敏感，外部特征接法不自然时会先在 count 上崩。
4. 任何新路线都要看 held-out test，而不是只看 validation。
5. 如果一种 foundation model 只在分类 benchmark 上强，但不能提供位置相关特征，它对这个任务未必有用。

---

## 7. 推荐的阅读顺序

读完本文件后，建议继续看：

1. `02_genos_attempts_training_and_eval.md`
2. `03_foundation_model_resources_and_constraints.md`
3. `05_foundation_model_run_matrix.csv`
4. `04_external_analysis_brief.md`
