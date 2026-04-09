# 任务定义、当前最佳模型与代码地图

## 1. 任务是什么

这是一个 base-resolution ATAC-seq 预测任务，不是普通分类任务。

- 输入：`2114 bp` DNA 序列
- 输出：
  - 中心 `1000 bp` 的 profile logits
  - 区域总 count/logcount
- 损失：
  - profile：multinomial NLL
  - count：MSE
- 评价：
  - profile `JSD`
  - count `Pearson r`

因此，一个 foundation model 或外部特征即使对 OCR classification 很有用，也不必然能帮助当前主任务。

---

## 2. 当前最佳模型是什么

当前默认最好模型是 `corrected B`，可以概括成：

- `conv stem + local dilated tower + transformer encoder`
- `ChromBPNet-style bias branch`
- 同时保留 `full` 与 `debiased` 两组输出
- count head 使用 `center pool`
- 训练时显式使用 `debiased_profile_weight=2.0`

代表性 held-out 结果：

- `B_s42`: peak `JSD=0.3147`, `count_r=0.8503`
- `B_s1234`: peak `JSD=0.3145`, `count_r=0.8488`

---

## 3. 这个模型最重要的 5 个语义

1. `full` 和 `debiased` 是两套都要看的输出，不是只保留一个报告口径。
2. `profile_bias_stop_gradient=true` 时，profile 融合路径会对 bias 分支做 detach。
3. count 融合默认是 `logsumexp`，不是简单相加。
4. 当前默认最好 readout 是 `center pool`，不是 attention pool。
5. 最终默认模型是否安全，不能只看主指标，必须同时看 `full/debiased gap`。

---

## 4. 最小配置摘要

当前默认模型的关键配置来自：

- `vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool.yaml`

可压缩成下面这些关键项：

```yaml
sequence_encoder:
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

训练侧最关键的 loss 语义是：

```yaml
loss:
  profile_weight: 1.0
  count_weight: 0.1
  debiased_profile_weight: 2.0
```

---

## 5. 代码地图

### 最关键的源码入口

- `06_transchrombp_current_model.py`
  - 这是当前 `TransChromBP` 顶层模型的逐字快照
  - 原始路径：`vendor/transchrombp/transchrombp/models/transchrombp.py`
- `07_foundation_adapter_current.py`
  - 这是当前 foundation adapter 的逐字快照
  - 原始路径：`vendor/transchrombp/transchrombp/models/foundation_adapter.py`

### 即使这次没外发，也值得知道的关键文件

- `vendor/transchrombp/transchrombp/training/train_ddp.py`
  - 训练主循环、validation、selector 指标、cached feature/teacher target plumbing
- `vendor/transchrombp/transchrombp/models/genos_adapter.py`
  - 早期 Genos online/cached 接入的主要实现
- `vendor/transchrombp/transchrombp/models/caduceus_adapter.py`
  - Caduceus-PS online token fusion 的实现

---

## 6. 读 06 时你应重点看什么

`06_transchrombp_current_model.py` 里最值得关注的是：

1. 顶层结构并不是“纯 Transformer”，而是 conv/local tower/transformer 的组合。
2. bias branch 是模型内显式分支，不是后处理修补。
3. `full` / `debiased` 输出是架构级语义，不只是日志里额外算一个指标。
4. Genos、Caduceus、foundation adapter 都是后来挂载进去的分支，不是主干最初就天然需要的部分。

---

## 7. 读 07 时你应重点看什么

`07_foundation_adapter_current.py` 里最值得看的是：

1. `FoundationCrossAttentionAdapter`
  - 是 residual-gated、零初始化输出投影的接法
  - 目的是让训练起点尽量接近 baseline，而不是一开始就强行改写主干
2. `FoundationResidualHead`
  - 不是直接替换主干，而是只预测 coarse residual correction
  - 这是当前 restart v3 最关键的新假设之一

---

## 8. 为什么这次只外发 2 个源码文件

因为当前文件数被限制在 8 以内，而外部模型第一次读包时最需要的是真相而不是噪声：

- `06` 负责告诉它“当前主模型到底长什么样”
- `07` 负责告诉它“foundation 新主线现在到底准备怎么接”

如果它进一步要求训练语义或 selector 细节，再单独补 `train_ddp.py` 会比一开始把整个仓库堆上去更有效。
