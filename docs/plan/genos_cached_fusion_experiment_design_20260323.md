# Genos 缓存融合实验设计

更新日期：2026-03-23

---

## 1. 来自 G0/G1/G2 实机数据的关键发现

### 1.1 训练速度

| 实验 | 单 epoch 耗时 | 20 epoch 总时 | 每步 Genos 开销 |
|------|-----------|-----------|-----------|
| G0 baseline | ~1457s (24min) | ~8h (已收口) | 0s |
| G1 genos_gate | ~14006s (233min) | ~84-90h (进行中) | ~1.3s |
| G2 genos_mean | 同量级 | ~80-90h (进行中) | ~1.3s |

根因：Genos 使用单碱基 tokenization（vocab=18），2114bp 输入产生 2114 tokens。双向平均需要 2 次 1.2B 模型前向传播。每步额外开销 ~1.3s，占总步时的 ~92%。

### 1.2 验证指标（epoch 5，唯一可比点）

| 实验 | peak JSD_full | count_r | 相对 G0 ΔJSD |
|------|-------------|---------|-------------|
| G0 baseline | 0.3429 | 0.7890 | — |
| G1 genos_gate | 0.3420 | 0.7898 | -0.0009 |
| G2 genos_mean | 0.3414 | 0.7955 | -0.0015 |

G0 最终（epoch 20）：JSD=0.3337, count_r=0.7998。

### 1.3 关键判读

1. **G2(mean) ≥ G1(gate)**：逐位置特征相比全局均值没有额外价值。当前收益（如果有的话）来自粗粒度 global signal。
2. **ΔJSD 在 0.001 量级**：可能是噪声，也可能是真实但微弱的信号被当前融合方式压制。
3. **ROI 极差**：花 10 倍算力换来不确定的 0.001 改善，不可接受。

### 1.4 对照：历史 V2 基线

V2full_s42（不带 Genos，每 epoch 验证）的 JSD 轨迹：

| epoch | JSD_full | 相对 epoch 5 改善 |
|-------|---------|----------------|
| 5 | 0.3420 | — |
| 10 | 0.3371 | 0.0049 |
| 20 | 0.3328 | 0.0092 |
| 30 | 0.3313 | 0.0107 |

这意味着模型自然训练 5→10 epoch 的改善量（0.005）已经远大于 Genos 在 epoch 5 带来的差异（0.001）。

---

## 2. 问题诊断

### 2.1 在线推理是核心浪费

当前每步执行链路（`train_ddp.py:1544-1546`）：

```
seq [B, 2114, 4]
→ onehot_to_dna()          # Python 循环，逐样本字符串拼接
→ tokenizer(fwd_strings)    # 2114 tokens × B
→ model(**fwd_inputs)        # 1.2B 前向，2114 tokens
→ tokenizer(rc_strings)     # 反向互补
→ model(**rc_inputs)         # 1.2B 前向，2114 tokens
→ (fwd + rc) * 0.5          # 平均
→ genos_feat [B, 2114, 1024]
```

Genos 是 frozen 的。同一个 (chrom, center) 的 canonical 窗口，Genos 输出完全确定。但当前实现在 20 epoch × 11809 steps = ~236k 步中重复计算了同样的特征。

**训练集有 ~236k 条 records。每条 record 在整个训练过程中被 Genos 提取了 ~20 次（每 epoch 1 次），总共 ~4.7M 次 1.2B 前向传播，其中绝大部分是纯重复。**

### 2.2 全位置特征缓存不可行

一个直觉是"把 Genos 全位置特征预计算缓存"。但实际计算：

- 训练集 236k records × 2114 positions × 1024 dims × 2 bytes (fp16) = **968 GB**
- 不可能存储

因此缓存只能走**摘要路径**。

### 2.3 jitter 对摘要的影响有限

训练时 peak 有 ±500bp jitter，但这不阻碍 canonical summary cache：

- 序列长度 2114bp，最大偏移 500bp
- canonical 和 jittered 窗口最小重叠 = (2114 - 500) / 2114 = **76.3%**
- 对于 mean-pooled 1024-dim 向量，76%+ 重叠下的变化量微小
- 反向互补：summary 使用 RC-average，方向不敏感
- 实际上 G2 已经证明了 mean pooling 的信号不比 per-position 差

因此：用 canonical center 的 mean summary 近似所有 jittered 版本是合理的首选。

### 2.4 融合位置不合理

当前融合点在 `local_tower` 之后、`Transformer` 之前（`transchrombp.py:351-352`）：

```python
tokens = local_feat.transpose(1, 2)          # [B, L, 256]
if genos_feat is not None:
    tokens = self.genos_adapter(tokens, genos_feat)  # ← 这里
encoded = self.transformer(tokens)            # [B, L, 256]
```

问题：
- Genos 特征先被融合，然后经过 6 层 Transformer。Transformer 有足够的表达能力"抹掉"微弱的辅助信号
- `gate_bias_init=-2.0`（sigmoid ≈ 0.12）让模型初始就倾向忽略 Genos
- 如果 Genos 的价值在于提供"这是什么类型的区域"的粗粒度上下文，那它应该在 Transformer 已经完成序列建模之后再注入——即**后段融合 (late fusion)**

### 2.5 Phase 0 的信号本身很弱

Phase 0 结论：layer_6 的 hold-out AUC = 0.654（peak/nonpeak 二分类）。这说明 Genos 表征对"区域是否开放"有微弱但非零的判别力。

但 TransChromBP 的 ConvStem + LocalTower + Transformer 本身已经是一个强 peak/nonpeak 分类器（模型能准确预测 profile shape 和 count）。Genos 特征能提供多少**互补信息**是未知的。

---

## 3. 核心设计变更

相比当前 G0/G1/G2 方案和 ChatGPT 的 pivot 方案，本设计有三个关键变更：

### 3.1 离线预计算取代在线推理

- 一次性在 A6000 上对所有 records 提取 Genos summary（~3-4h）
- 缓存为 `.pt` 文件，训练时 dataset 直接读取
- 训练速度回到接近 G0 baseline 的水平

### 3.2 后段融合取代前段注入

- 不再在 Transformer 之前注入（当前方案）
- 改为在 Transformer 之后、prediction heads 之前注入
- 使用 FiLM (Feature-wise Linear Modulation) 作为主要融合方式

### 3.3 增加互补性预诊断

- 在启动任何训练 pilot 之前，先用线性探测量化 Genos 特征的互补性
- 这是一个 ~30min 的 CPU 计算，可以避免在注定无效的方向上投入训练算力

---

## 4. 离线预计算缓存系统

### 4.1 缓存内容

对每条 record（由 `idx` 或 `(chrom, center)` 唯一标识），预计算：

```
global_mean:  [1024]   # 全序列 mean pooling 后的 RC-averaged 向量
bins4_mean:   [4, 1024] # 将 2114 tokens 分成 4 段（每段 ~528 tokens），每段取 mean
```

**存储开销：**

| 缓存类型 | 训练集 (236k) | 验证集 (23k) | 总计 |
|---------|------------|-----------|------|
| global_mean (fp16) | 484 MB | 48 MB | 532 MB |
| bins4_mean (fp16) | 1.9 GB | 192 MB | 2.1 GB |

两者都完全可以放在 CPU 内存中。

### 4.2 预计算脚本设计

新增 `scripts/build_genos_cache.py`，核心逻辑：

```python
"""
用法：
  python scripts/build_genos_cache.py \
    --data-config configs/data/data_tutorial_canonical_v1.yaml \
    --folds-json .../folds.json \
    --split train \
    --genos-model-path /data1/.../Genos-1.2B \
    --layer 6 \
    --batch-size 20 \
    --output-dir /data1/.../genos_cache/tutorial_canonical_v1/
"""

def build_cache(args):
    # 1. 加载 records（与训练时相同的 peak + nonpeak 列表）
    records = load_records(args.data_config, args.split)

    # 2. 加载 Genos extractor
    extractor = GenosFeatureExtractor(
        model_path=args.genos_model_path,
        layer=args.layer,
        device="cuda",
        bidirectional=True,
    )

    # 3. 逐批次提取 canonical 特征
    #    对每条 record，使用 canonical center（无 jitter），无 revcomp
    global_means = []  # 累积 [1024] 向量
    bins4_means = []   # 累积 [4, 1024] 向量

    for batch_records in batched(records, args.batch_size):
        # 取 canonical 窗口序列（center ± 1057）
        seqs = [fetch_onehot(r.chrom, r.center - 1057, 2114) for r in batch_records]
        seq_tensor = torch.stack(seqs).to("cuda")  # [B, 2114, 4]

        # Genos 提取（含 RC-average）
        feat = extractor.extract(seq_tensor)  # [B, 2114, 1024]

        # 全局 mean
        global_means.append(feat.mean(dim=1).half().cpu())  # [B, 1024]

        # 4-bin mean
        L = feat.shape[1]  # 2114
        bin_size = L // 4   # 528
        bins = []
        for i in range(4):
            start = i * bin_size
            end = start + bin_size if i < 3 else L
            bins.append(feat[:, start:end, :].mean(dim=1))  # [B, 1024]
        bins4_means.append(torch.stack(bins, dim=1).half().cpu())  # [B, 4, 1024]

    # 4. 保存
    torch.save({
        "global_mean": torch.cat(global_means),   # [N, 1024] fp16
        "bins4_mean": torch.cat(bins4_means),      # [N, 4, 1024] fp16
        "record_keys": [(r.chrom, r.center, r.source) for r in records],
    }, output_path)
```

**预计算耗时估算：**

| split | records | batches (bs=20) | Genos 提取 (~1.3s/batch) | 总计 |
|-------|---------|-----------------|----------------------|------|
| train | 236,188 | 11,810 | ~15,353s | ~4.3h |
| valid | ~23,417 | ~1,171 | ~1,522s | ~25min |

总计约 **5 小时**，一次性投入。作为对比，G1 一个 20 epoch 的 run 需要 ~85 小时。

### 4.3 数据集集成

修改 `real_data.py` 的 `ChromBPNetBigWigDataset`，增加可选的 cache 加载：

```python
class ChromBPNetBigWigDataset(Dataset):
    def __init__(self, ..., genos_cache_path: str = ""):
        # ... 现有初始化 ...
        self._genos_cache = None
        if genos_cache_path:
            cache = torch.load(genos_cache_path, map_location="cpu")
            self._genos_cache = {
                "global_mean": cache["global_mean"],    # [N, 1024] fp16
                "bins4_mean": cache.get("bins4_mean"),   # [N, 4, 1024] fp16 (可选)
            }

    def __getitem__(self, idx: int) -> dict[str, Any]:
        # ... 现有逻辑 (jitter, revcomp) ...
        result = {
            "seq": torch.from_numpy(seq),
            "profile_counts": torch.from_numpy(profile_counts),
            "region_source": record.source,
        }
        if self._genos_cache is not None:
            result["genos_global"] = self._genos_cache["global_mean"][idx].float()
            if self._genos_cache.get("bins4_mean") is not None:
                result["genos_bins4"] = self._genos_cache["bins4_mean"][idx].float()
        return result
```

**关键决策：cache 按 `idx` 索引，无需引入 `record_id` 字段。** 因为 records 列表在同一 data config 下的构建顺序是确定的（同样的 peaks_bed + nonpeaks_bed + folds_json + split → 同样的 records 顺序）。预计算脚本和训练脚本使用相同的数据配置即可保证对齐。

### 4.4 训练循环修改

修改 `train_ddp.py` 的训练/验证循环，支持从 batch 读取 cached summary：

```python
# 当前（在线提取）：
genos_feat = None
if genos_runtime is not None:
    genos_feat = genos_runtime.extract(seq)
outputs = model(seq, genos_feat=genos_feat)

# 改为（缓存读取）：
genos_summary = None
if "genos_global" in batch:
    genos_summary = batch["genos_global"].to(device, non_blocking=True)  # [B, 1024]
outputs = model(seq, genos_summary=genos_summary)
```

**在线提取代码路径保留但不再作为 Genos pilot 的默认方式。** 新的 config 通过 `genos_branch.mode: cached_summary` vs `genos_branch.mode: online` 区分。

---

## 5. 融合架构

### 5.1 当前架构的数据流（复习）

```
seq [B, 2114, 4]
  → ConvStem     → [B, 256, 2114]
  → LocalTower   → [B, 256, 2114]
  → transpose    → [B, 2114, 256]  (tokens)
  → Transformer  → [B, 2114, 256]  (encoded)
  ├→ Profile Head: Linear(256→1) + center_crop → [B, 1000]  (profile_signal)
  └→ Count Head:  mean_pool → LayerNorm → Linear(256→128) → GELU → Dropout → Linear(128→1)  (count_signal)
  → Bias Branch  → profile_bias [B, 1000], count_bias [B, 1]
  → Fusion       → profile_full, count_full
```

### 5.2 方案 P1: Global FiLM (主推荐)

**在 Transformer 输出之后、两个 prediction head 之前，用全局 Genos summary 做 FiLM conditioning。**

新增模块 `GenosSummaryFiLM`：

```python
class GenosSummaryFiLM(nn.Module):
    """FiLM conditioning on Transformer output using cached Genos summary."""

    def __init__(self, d_model: int = 256, genos_dim: int = 1024):
        super().__init__()
        self.proj = nn.Linear(genos_dim, d_model * 2)  # 生成 gamma 和 beta
        # 初始化：gamma 接近 1，beta 接近 0 → 初始行为等价于无 Genos
        nn.init.zeros_(self.proj.weight)
        nn.init.zeros_(self.proj.bias)

    def forward(self, encoded: Tensor, genos_summary: Tensor) -> Tensor:
        """
        Args:
            encoded:       [B, L, d_model]   Transformer 输出
            genos_summary: [B, genos_dim]    全局 Genos summary
        Returns:
            [B, L, d_model]  调制后的特征
        """
        params = self.proj(genos_summary)              # [B, d_model * 2]
        gamma, beta = params.chunk(2, dim=-1)          # 各 [B, d_model]
        gamma = 1.0 + gamma.unsqueeze(1)               # [B, 1, d_model]
        beta = beta.unsqueeze(1)                       # [B, 1, d_model]
        return gamma * encoded + beta
```

修改 `TransChromBP.forward()`：

```python
def forward(self, seq_onehot, ..., genos_summary=None):
    # ... ConvStem → LocalTower → transpose ...
    encoded = self.transformer(tokens)                   # [B, L, 256]

    # ── 新增：Genos FiLM conditioning ──
    if genos_summary is not None and self.genos_film is not None:
        encoded = self.genos_film(encoded, genos_summary)

    profile_signal = self.profile_signal_head(encoded).squeeze(-1)
    profile_signal = center_crop_1d(profile_signal, self.output_len)
    pooled = encoded.mean(dim=1)
    count_signal = self.count_signal_head(pooled)
    # ... bias, fusion ...
```

**特性：**
- 参数量：Linear(1024, 512) = 524,800 参数，约占模型总参数量的 <1%
- 初始化为恒等变换（gamma=1, beta=0），训练开始时行为与 baseline 完全一致
- 同时影响 profile 和 count 预测
- FiLM 的 scale/shift 是位置无关的（global summary），但通过调制 Transformer 特征的通道幅度和偏移来影响下游头的预测

**为什么选 FiLM 而不是 gated add：**
- Gated add 在当前实验中表现不如 mean pooling
- FiLM 不直接修改特征向量的方向，而是调制通道的幅度——这更适合"微调"型的辅助信号
- FiLM 的零初始化保证了完美的 baseline 兼容性

### 5.3 方案 P2: Global Summary → Count Head 注入

**单独验证 Genos summary 对 count estimation 的帮助。**

将 count_signal_head 从 Sequential 拆成显式步骤，在中间层注入 Genos：

```python
# 原始：
self.count_signal_head = nn.Sequential(
    nn.LayerNorm(d_model),              # 256
    nn.Linear(d_model, d_model // 2),   # 256 → 128
    nn.GELU(),
    nn.Dropout(dropout),
    nn.Linear(d_model // 2, 1),         # 128 → 1
)

# 修改为：
self.count_ln = nn.LayerNorm(d_model)
self.count_fc1 = nn.Linear(d_model, d_model // 2)  # 256 → 128
self.count_fc2 = nn.Linear(d_model // 2, 1)        # 128 → 1
self.count_dropout = nn.Dropout(dropout)

# 新增：
self.genos_count_proj = nn.Linear(genos_dim, d_model // 2)  # 1024 → 128
nn.init.zeros_(self.genos_count_proj.weight)
nn.init.zeros_(self.genos_count_proj.bias)
```

前向：

```python
pooled = encoded.mean(dim=1)                   # [B, 256]
h = self.count_ln(pooled)
h = self.count_fc1(h)                          # [B, 128]
if genos_summary is not None:
    h = h + self.genos_count_proj(genos_summary)  # [B, 128]
h = F.gelu(h)
h = self.count_dropout(h)
count_signal = self.count_fc2(h)               # [B, 1]
```

**诊断价值：** P2 只影响 count，不影响 profile。如果 P2 改善了 count_pearson_full 但 JSD 没变，说明 Genos 的粗粒度信号对 count estimation 有用，但不足以改善 profile shape。

### 5.4 方案 P3: 4-bin FiLM (条件触发)

**用 4-bin 摘要替代 global summary，引入粗粒度位置信息。**

```python
class GenosBinnedFiLM(nn.Module):
    """FiLM conditioning with binned Genos summary (position-varying)."""

    def __init__(self, d_model: int = 256, genos_dim: int = 1024, n_bins: int = 4):
        super().__init__()
        self.n_bins = n_bins
        self.proj = nn.Linear(genos_dim, d_model * 2)
        nn.init.zeros_(self.proj.weight)
        nn.init.zeros_(self.proj.bias)

    def forward(self, encoded: Tensor, genos_bins: Tensor) -> Tensor:
        """
        Args:
            encoded:    [B, L, d_model]
            genos_bins: [B, n_bins, genos_dim]
        """
        B, L, D = encoded.shape
        params = self.proj(genos_bins)            # [B, n_bins, d_model * 2]
        gamma, beta = params.chunk(2, dim=-1)     # 各 [B, n_bins, d_model]
        gamma = 1.0 + gamma                       # [B, n_bins, d_model]

        # 上采样到序列长度
        bin_size = L // self.n_bins
        gamma_up = gamma.repeat_interleave(bin_size, dim=1)[:, :L, :]  # [B, L, d_model]
        beta_up = beta.repeat_interleave(bin_size, dim=1)[:, :L, :]    # [B, L, d_model]

        return gamma_up * encoded + beta_up
```

**只在 P1 给出正信号且 P1 > P2 时启动 P3。** 如果 P1 ≈ P2，说明 global summary 已经足够，位置信息没有额外价值。

---

## 6. 可选预诊断：特征互补性线性探测

**目的：** 在投入训练算力之前，用 ~30min CPU 计算判断 Genos summary 是否与 TransChromBP 特征互补。

**步骤：**

1. 用 G0 best.pt 对验证集提取 Transformer 最后一层的 pooled 特征 `[N_val, 256]`
2. 用预计算缓存中的 global_mean `[N_val, 1024]`
3. 标签：peak=1, nonpeak=0
4. 训练三个 logistic regression：
   - A: TransChromBP-only `[256]` → AUC_A
   - B: Genos-only `[1024]` → AUC_B
   - C: 拼接 `[1280]` → AUC_C
5. 判读：
   - 若 AUC_C > max(AUC_A, AUC_B) + 0.01 → 有互补信号，继续 pilot
   - 若 AUC_C ≈ AUC_A → Genos 没有额外信息，直接降级

**说明：这一步是可选的，不是必须的。** 如果时间紧迫或更倾向直接用训练实验验证，可以跳过。但作为一个几乎零成本的预筛，建议先做。

---

## 7. 实验矩阵

### 第一轮（~30h 总 GPU 时间）

| 组别 | 名称 | Genos 输入 | 融合方式 | 训练速度 | 目的 |
|------|------|-----------|---------|---------|------|
| P0 | baseline_short10 | 无 | 无 | ~baseline | matched control |
| P1 | global_film | global_mean [1024] | FiLM on Transformer output | ~baseline (+5%) | 主实验：粗粒度 FiLM 能否改善 profile + count |
| P2 | global_count | global_mean [1024] | Count head 中间层注入 | ~baseline (+2%) | 诊断：Genos 是否只能帮 count |

三组用同一张 A6000 串行跑（每组 10 epoch × ~1500s/epoch ≈ 4.2h），总计 ~13h 训练 + 5h 预计算 = ~18h。

### 第二轮（条件触发）

仅当第一轮满足 Go 条件时：

| 组别 | 名称 | Genos 输入 | 融合方式 | 触发条件 |
|------|------|-----------|---------|---------|
| P3 | bins4_film | bins4_mean [4, 1024] | Binned FiLM | P1 达到 Go 阈值 |
| P4 | global_film_deep | global_mean [1024] | FiLM + 2 层 MLP | P1 有信号但未跨 Go |

### P0 与 G0 的关系

P0 **不等于** G0：
- P0: max_epochs=10, validate_every_epochs=1, early_stop_patience=3, batch_size=20
- G0: max_epochs=20, validate_every_epochs=5, early_stop_patience=8, batch_size=20

P0 是一个更短、更频繁验证的 matched control，专为快速对比设计。G0 的最终指标可以作为参考，但 P0 的曲线才是直接可比对象。

---

## 8. 训练与评估配置

### 8.1 统一训练参数

```yaml
# train_genos_pivot_short10.yaml
seed: 42
max_epochs: 10
trainer:
  # ... 沿用现有 trainer 配置 ...
  validate_every_epochs: 1
  checkpoint_every_epochs: 1
  best_metric: peak.profile_target_jsd_full_mean
  best_metric_mode: min
  early_stop_patience: 3
  early_stop_min_delta: 0.0

optimizer:
  name: adamw
  learning_rate: 5.0e-4
  weight_decay: 0.01
  betas: [0.9, 0.95]

loss:
  profile_weight: 1.0
  count_weight: 0.1
  debiased_profile_weight: 2.0
  debiased_count_weight: 0.0

schedule:
  name: cosine
  warmup_ratio: 0.06
  min_lr_ratio: 0.1

data:
  # ... 沿用 data_tutorial_canonical_v1.yaml ...
  batch_size_per_gpu: 20
  genos_cache_path: ""  # P0 留空；P1/P2/P3 指向预计算文件
```

### 8.2 额外日志

所有 Genos pivot run 增加以下 per-epoch 日志字段（附加到 `epoch_metrics.jsonl`）：

| 字段 | 说明 |
|------|------|
| `genos_film_gamma_mean` | FiLM gamma 均值（应接近 1.0） |
| `genos_film_gamma_std` | FiLM gamma 标准差 |
| `genos_film_beta_norm` | FiLM beta L2 范数 |
| `step_time_mean` | 平均步时间（验证缓存加速效果） |

### 8.3 Go / No-Go 判断标准

**Go（继续到第二轮）：**

满足以下任一条件（在 epoch 5-8 时点）：

- P1 的 `peak.profile_target_jsd_full_mean` 相对 P0 同 epoch 改善 >= 0.002
- P1 同时满足：ΔJSD >= 0.001 **且** Δcount_r >= 0.005
- P2 的 `peak.count_pearson_full` 相对 P0 同 epoch 改善 >= 0.01

**Conditional Go（启动 P4 但不启动 P3）：**

- P1 有 0.001-0.002 的 JSD 改善趋势，但未达到 Go 阈值
- 且 FiLM 参数没有塌到零（gamma_std > 0.01）

**No-Go（降级）：**

满足以下任一条件：

- P1、P2 在 epoch 4 之前都没有超过 P0
- FiLM 参数在 epoch 3 后仍接近初始化值（gamma ≈ 1.0 ± 0.005, beta ≈ 0）
- P1 的 nonpeak JSD 显著恶化（> 0.01 绝对退化）

**阈值说明：** 参考 V2full_s42 从 epoch 5→10 的自然改善量（ΔJSD = 0.005），要求 Genos 贡献 0.002 相当于约 40% 的自然训练改善量。这比 ChatGPT 方案的 0.003 阈值更宽松，但仍然要求收益是可感知的，而不是噪声级的 0.001。

---

## 9. 代码实现清单

### 9.1 新增文件

| 文件 | 说明 |
|------|------|
| `scripts/build_genos_cache.py` | 离线预计算脚本 |
| `scripts/run_genos_probe.py` | （可选）互补性线性探测脚本 |
| `configs/model/v2fix_genos_global_film.yaml` | P1 模型配置 |
| `configs/model/v2fix_genos_global_count.yaml` | P2 模型配置 |
| `configs/model/v2fix_genos_bins4_film.yaml` | P3 模型配置（第二轮） |
| `configs/train/train_genos_pivot_short10.yaml` | pivot 统一训练配置 |
| `scripts/run_genos_pivot.sh` | 新 pilot 启动脚本 |

### 9.2 修改文件

| 文件 | 改动说明 |
|------|---------|
| `models/genos_adapter.py` | 新增 `GenosSummaryFiLM` 和 `GenosBinnedFiLM` 类 |
| `models/transchrombp.py` | forward() 增加 `genos_summary` 参数；增加可选的 `self.genos_film` 模块；P2 需拆 count_signal_head |
| `data/real_data.py` | `__init__` 增加 `genos_cache_path` 参数；`__getitem__` 返回缓存的 summary |
| `training/train_ddp.py` | 训练/验证循环支持从 batch 读取 genos_summary；增加 FiLM 统计日志 |

### 9.3 不改动的文件

- `models/bias_branch.py` — Genos 不进入 bias branch
- `losses/` — 损失函数不变
- `evaluation/` — 评估流程不变

### 9.4 代码改动量估算

| 文件 | 新增行 | 修改行 |
|------|-------|-------|
| `genos_adapter.py` | ~60 (FiLM + BinnedFiLM) | 0 |
| `transchrombp.py` | ~30 (FiLM 初始化 + forward 分支) | ~5 (forward 签名) |
| `real_data.py` | ~15 (cache 加载 + getitem) | ~2 (init 签名) |
| `train_ddp.py` | ~20 (cache 模式分支 + FiLM 日志) | ~10 (训练/验证循环) |
| `build_genos_cache.py` | ~120 | — |
| `run_genos_pivot.sh` | ~80 | — |
| 配置文件 | ~60 | — |
| **总计** | **~385** | **~17** |

---

## 10. 执行排程

### 阶段 0：止损（立即）

- G2：停止（已经提供了"mean ≥ gate"的结论）
- G1：等到 epoch 10 validation 出结果后停止
  - 如果 epoch 10 的 JSD 相对 G0 epoch 10（0.3370）没有改善 → 立即停止
  - 如果有 >= 0.003 的改善 → 记录后仍停止（因为 ROI 依然差），但纳入后续分析

### 阶段 1：基础设施（~1 天）

在本地完成所有代码改动和配置编写：

1. 实现 `GenosSummaryFiLM` 和 count head 注入
2. 实现 `build_genos_cache.py`
3. 修改 dataset 支持 cache 读取
4. 修改 trainer 支持 cached summary 模式
5. 编写 P0/P1/P2 的模型和训练配置
6. `py_compile` 全部新增/修改文件
7. 编写 launcher 脚本

### 阶段 2：预计算（~5h，一次性）

在 6000 A6000 上：

```bash
# 停掉 G1/G2 释放 GPU 后
python scripts/build_genos_cache.py \
  --data-config configs/data/data_tutorial_canonical_v1.yaml \
  --split train --output-dir .../genos_cache/ --gpu 0

python scripts/build_genos_cache.py \
  --data-config configs/data/data_tutorial_canonical_v1.yaml \
  --split valid --output-dir .../genos_cache/ --gpu 0
```

### 阶段 2.5：（可选）互补性探测（~30min CPU）

```bash
python scripts/run_genos_probe.py \
  --g0-checkpoint .../genos_20260321_baseline_s42/best.pt \
  --genos-cache .../genos_cache/valid_global_mean.pt \
  --val-data-config configs/data/data_tutorial_canonical_v1.yaml
```

### 阶段 3：第一轮 pilot（~13h）

```bash
# P0 → P1 → P2 串行
bash scripts/run_genos_pivot.sh --groups P0,P1,P2 --gpu 0
```

### 阶段 4：判读与决策（~1h）

- 汇总 P0/P1/P2 的 epoch 1-10 验证曲线
- 按 Go/No-Go 标准判断
- 若 Go → 进入第二轮 P3/P4
- 若 No-Go → Genos 线正式降级，A6000 释放

### 总时间预算

| 阶段 | 耗时 | GPU 占用 |
|------|------|---------|
| 止损 + 基础设施 | ~1 天（主要是代码工作） | 不占 GPU |
| 预计算缓存 | ~5h | 1× A6000 |
| 互补性探测 | ~30min | CPU only |
| 第一轮 pilot (P0+P1+P2) | ~13h | 1× A6000 |
| 判读 | ~1h | — |
| **总计到 Go/No-Go 决策** | **~2 天** | **~18h GPU** |

作为对比：当前 G1 一个 run 就需要 ~85h GPU 时间且还没跑完。

---

## 11. 降级条件

如果出现以下任一情况，Genos 线正式降级：

1. 互补性探测（若执行）显示 AUC_C ≈ AUC_A
2. 第一轮 P1/P2 全部落入 No-Go
3. 两轮 pilot 做完后，最好的 Genos 方案相对 P0 的 JSD 改善 < 0.002

降级后：
- Genos 相关代码保留但不再作为训练主线
- A6000 优先分配给：V2fix 多 seed、GM12878/K562 数据扩容实验、更长 epoch 训练
- 保留 Phase 0 + Phase 1 的结论文档，作为"Genos-1.2B 在 ATAC-seq profile 建模中收益有限"的实验证据

---

## 12. 与 ChatGPT 方案的关系

本方案保留了 ChatGPT 方案的核心框架（缓存化、summary-first、短 pilot、Go/No-Go），在以下方面做了调整：

| 方面 | ChatGPT 方案 | 本方案 | 原因 |
|------|------------|--------|------|
| 融合方式 | P1 接 count head / late MLP | P1 用 FiLM on Transformer output | FiLM 同时影响 profile + count，且零初始化保证 baseline 兼容 |
| Count head 实验 | 合并在 P1 中 | 拆为独立的 P2 | 分开诊断 count 和 profile 维度的收益 |
| P2 (ChatGPT 的 4-bin FiLM) | 第一轮 | 推到第二轮 (P3) | 先验证 global signal 是否有效再增加复杂度 |
| Go 阈值 | ΔJSD >= 0.003 | ΔJSD >= 0.002 | 参照历史训练的自然改善量，0.003 过于严格 |
| 互补性预诊断 | 未包含 | 可选步骤 | 几乎零成本但可以避免无效投入 |
| 缓存 record_id | 需新增 record_id 字段 | 直接用 idx 索引 | records 列表顺序在同一 config 下确定，无需额外 ID |
| 预计算耗时 | 未估算 | ~5h (一次性) | 让用户知道这是一个可接受的一次性投入 |
