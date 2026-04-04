# TransChromBP × Genos-1.2B 基因组基座模型接入完整技术报告

> **目的**：供外部 AI 评审判断——当前负面结果是否是 Genos 特有的问题，还是基因组基座模型与本任务之间的系统性不匹配。如果换用其他基座模型（如 Nucleotide Transformer、DNABERT-2、HyenaDNA、Caduceus、Evo 等），是否可能获得正收益。

---

## 1. 任务定义与主模型架构

### 1.1 任务

**染色质可及性碱基分辨率建模（ATAC-seq）**

- **输入**：2114 bp DNA 序列（one-hot 编码，`[B, 2114, 4]`）
- **输出 1 — Profile**：中心 1000 bp 的碱基分辨率 read count 分布（softmax 概率，`[B, 1000]`）
- **输出 2 — Count**：该区域的总 read count 对数值（标量，`[B, 1]`）
- **损失函数**：Profile 用多项分布负对数似然（multinomial NLL），Count 用 MSE
- **评估指标**：
  - Profile 用 Jensen-Shannon Divergence (JSD)，越低越好
  - Count 用 Pearson correlation (count_r)，越高越好

### 1.2 TransChromBP 模型架构

```
输入 seq_onehot [B, 2114, 4]
    │
    ▼
Conv Stem (1D Conv → BN → ReLU)          # 将 4 通道映射到 d_model=256
    │
    ▼
Local Tower (多层 dilated conv)           # 捕捉局部 motif 组合
    │
    ▼
tokens [B, L, 256]
    │
    ├── [可选] Genos Gated Adapter        # ← 在线 Genos 注入点 (G1/G2)
    │
    ▼
Transformer Encoder (4 层, 4 heads)       # 全局上下文建模
    │
    ├── [可选] Genos FiLM Modulation      # ← 缓存 Genos 注入点 (P1)
    │
    ▼
encoded [B, L, 256]
    │
    ├── Profile Head: Linear → [B, 1000]  # 碱基分辨率 profile
    │
    └── Count Head:
        │  mean_pool → LayerNorm → Linear(256→128) → GELU
        │     ├── [可选] Genos Count Proj # ← 缓存 Genos 注入点 (P2)
        │     ▼
        │  Dropout → Linear(128→1)
        └── → count [B, 1]

+ Bias Branch (独立的偏差模型，处理酶切偏好，预训练后冻结)
最终输出 = signal + scale * bias
```

### 1.3 偏置分解机制

TransChromBP 继承 ChromBPNet 的核心思想：先训练一个 Bias 模型捕捉酶切偏好（Tn5 transposase bias），再用主模型回归去偏后的真实 TF 信号。最终输出 = signal_branch + learnable_scale × bias_branch。Bias branch 预训练后冻结。

---

## 2. Genos-1.2B 模型概述

- **来源**：基因组基座大模型，1.2B 参数
- **架构**：Causal Language Model（类 GPT），基于 Transformers 库加载
- **Tokenizer**：单碱基 tokenizer（1 token = 1 base），输入 2114 bp → 2114 tokens
- **Hidden size**：1024
- **Layers**：12 层 Transformer
- **预训练目标**：DNA 序列自回归生成（next-token prediction）
- **官方 benchmark**：在 OCR（Open Chromatin Region）二分类任务上表现良好

### 2.1 我们验证的 OCR benchmark 复现结果

为排除本地环境/权重问题，我们复现了 Genos 官方的 `human_ocr_ensembl` 基准：

| 配置 | ROC AUC | 说明 |
|------|---------|------|
| Mini smoke (layer 6) | 0.5958 | 小样本冒烟测试 |
| Full (layer 6) | 0.7242 | 中间层 |
| **Full (layer 12)** | **0.7535** | **接近官方参考值 0.7569** |

**结论**：本地 Genos 推理链路正常，权重/tokenizer/hidden-state 提取均无问题。

---

## 3. Genos 特征提取方式

### 3.1 在线特征提取器（用于 G1/G2 实验）

冻结 Genos-1.2B，在训练时实时提取 hidden states，**不参与反向传播**。

```python
class GenosFeatureExtractor:
    """冻结 Genos-1.2B wrapper。NOT an nn.Module — 活在 DDP 外面。"""

    def __init__(self, model_path, layer=6, device="cuda",
                 attn_implementation="flash_attention_2", bidirectional=True):
        from transformers import AutoModelForCausalLM, AutoTokenizer

        self.layer = layer
        self.device = torch.device(device)
        self.bidirectional = bidirectional

        self.tokenizer = AutoTokenizer.from_pretrained(model_path, local_files_only=True)
        self.model = AutoModelForCausalLM.from_pretrained(
            model_path, local_files_only=True,
            dtype=torch.bfloat16,
            attn_implementation=attn_implementation,
        ).to(self.device).eval()

        for p in self.model.parameters():
            p.requires_grad_(False)

        self.hidden_size = self.model.config.hidden_size  # 1024

    @torch.no_grad()
    def extract(self, seq_onehot):
        """
        输入: [B, L, 4] one-hot DNA
        输出: [B, L, 1024] float32 hidden states

        双向模式: forward hidden + reverse-complement hidden 取平均
        """
        # one-hot → DNA 字符串
        dna_strings = onehot_to_dna(seq_onehot)

        fwd_inputs = self._tokenize(dna_strings)
        fwd_hidden = self._forward_hidden(fwd_inputs)  # [B, L, 1024]

        if not self.bidirectional:
            return fwd_hidden

        # 反向互补
        rc_strings = [reverse_complement(s) for s in dna_strings]
        rc_inputs = self._tokenize(rc_strings)
        rc_hidden = self._forward_hidden(rc_inputs)
        rc_hidden = rc_hidden.flip(dims=[1])  # 翻转位置轴

        return (fwd_hidden + rc_hidden) * 0.5  # 双向平均

    def _forward_hidden(self, token_inputs):
        outputs = self.model(**token_inputs, output_hidden_states=True)
        hidden = outputs.hidden_states[self.layer]  # [B, L, 1024]
        return hidden.float()
```

**关键设计**：
- Genos 是 causal LM，天然只看左侧上下文。通过双向（forward + reverse-complement 平均）弥补这一缺陷
- 提取 layer 6（中间层），因为 feasibility 测试表明中间层在 peak/nonpeak 区分上表现最佳
- 冻结权重，`torch.no_grad()`，不参与梯度计算

### 3.2 离线缓存特征（用于 P0/P1/P2 实验）

预计算所有数据集区域的 Genos 特征，存为 numpy memmap，训练时直接查表：

```python
# build_genos_summary_cache.py 核心逻辑

def pool_bins4(feat, n_bins=4):
    """将 [B, L, H] 池化为 [B, 4, H]，分成 4 个空间 bin 取均值"""
    B, L, H = feat.shape
    bin_size = L // n_bins
    feat = feat[:, :bin_size * n_bins, :]
    return feat.reshape(B, n_bins, bin_size, H).mean(dim=2)

# 对每个区域提取两种粒度的特征:
for batch_idx in range(n_batches):
    feat = extractor.extract(seq_batch)  # [B, L, 1024]

    # 特征 1: 全局均值 — 将整个序列压缩为单个 1024 维向量
    global_mean = feat.mean(dim=1)  # [B, 1024]

    # 特征 2: 4-bin 均值 — 将序列分成 4 段，每段取均值
    bins4_mean = pool_bins4(feat)   # [B, 4, 1024]
```

**存储格式**：
- `train_global_mean.f16.npy` — shape: `[n_train, 1024]`，float16
- `valid_global_mean.f16.npy` — shape: `[n_valid, 1024]`，float16
- `manifest_{split}.json` — 元数据 + SHA1 校验

**缓存语义的重要限制**：
- **Valid split**：exact cache，与验证数据完全一致（jitter=0, revcomp=false）
- **Train split**：只是 canonical center 的 side info。实际训练使用 `peak_max_jitter=500, random_revcomp=true`，即训练样本看到的序列和缓存对应的 canonical 序列不同

---

## 4. 融合方式详解（所有已实现方案）

### 4.1 方案 G1 — 逐位置门控融合（online, 注入 Transformer 之前）

```python
class GenosGatedAdapter(nn.Module):
    """
    公式: fused = tokens + gate * proj(genos_feat)
    gate 从 local tokens 计算，初始化接近 0 (sigmoid(-2) ≈ 0.12)
    """
    def __init__(self, d_model=256, genos_hidden_size=1024,
                 gate_bias_init=-2.0, pool_mode="none"):
        super().__init__()
        self.pool_mode = pool_mode  # "none" = 逐位置, "mean" = 全局广播
        self.proj = nn.Linear(genos_hidden_size, d_model)   # 1024 → 256
        self.ln = nn.LayerNorm(d_model)
        self.gate_proj = nn.Linear(d_model, d_model)        # 门控
        nn.init.constant_(self.gate_proj.bias, gate_bias_init)

    def forward(self, tokens, genos_feat):
        """
        tokens:     [B, L, 256]   来自 conv_stem + local_tower
        genos_feat: [B, L, 1024]  来自 GenosFeatureExtractor
        """
        if self.pool_mode == "mean":
            # G2 变体: 先全局平均再广播到每个位置
            genos_feat = genos_feat.mean(dim=1, keepdim=True).expand_as(genos_feat)

        projected = self.ln(self.proj(genos_feat))       # [B, L, 256]
        gate = torch.sigmoid(self.gate_proj(tokens))     # [B, L, 256]
        return tokens + gate * projected
```

**注入位置**：Local Tower 输出之后、Transformer 之前

```python
# transchrombp.py forward() 中的调用:
tokens = local_feat.transpose(1, 2)  # [B, L, 256]
if genos_feat is not None and self.genos_adapter is not None:
    tokens = self.genos_adapter(tokens, genos_feat)  # ← G1/G2 注入
encoded = self.transformer(tokens)
```

- **G1 (pool_mode="none")**：逐位置注入，每个碱基位置有独立的 Genos 特征
- **G2 (pool_mode="mean")**：先对 Genos 特征做全局平均池化，再广播到每个位置

### 4.2 方案 P1 — 后置 FiLM 调制（cached, Transformer 输出之后）

```python
class GenosSummaryFiLM(nn.Module):
    """
    Feature-wise Linear Modulation:
    output = (1 + gamma) * encoded + beta

    gamma, beta 由 genos_summary 通过线性映射得到
    零初始化: 训练开始时 gamma=0, beta=0，即 output = encoded（与 baseline 一致）
    """
    def __init__(self, d_model=256, genos_dim=1024):
        super().__init__()
        self.proj = nn.Linear(genos_dim, d_model * 2)  # 输出 gamma + beta
        nn.init.zeros_(self.proj.weight)   # 零初始化
        nn.init.zeros_(self.proj.bias)

    def forward(self, encoded, genos_summary):
        """
        encoded:       [B, L, 256]  Transformer 输出
        genos_summary: [B, 1024]    全局均值 Genos 特征
        """
        params = self.proj(genos_summary)             # [B, 512]
        gamma_delta, beta = params.chunk(2, dim=-1)   # [B, 256] each
        gamma = 1.0 + gamma_delta.unsqueeze(1)        # [B, 1, 256]
        beta = beta.unsqueeze(1)                      # [B, 1, 256]
        return gamma * encoded + beta
```

**注入位置**：Transformer 输出之后、Profile/Count Head 之前

```python
# transchrombp.py forward():
encoded = self.transformer(tokens)
if genos_summary is not None and self.genos_film is not None:
    encoded = self.genos_film(encoded, genos_summary)  # ← P1 注入
```

> **注意**：P1 方案因 P2 已明确失败而被取消执行，代码已实现但未跑实验。

### 4.3 方案 P2 — Count Head 直接注入（cached, 最激进的方案）

```python
class GenosCountProj(nn.Module):
    """
    在 count head 中间层做加法注入:
    count_hidden = count_hidden + proj(genos_summary)

    零初始化: 训练开始时注入量为 0
    """
    def __init__(self, d_count_hidden=128, genos_dim=1024):
        super().__init__()
        self.proj = nn.Linear(genos_dim, d_count_hidden)  # 1024 → 128
        nn.init.zeros_(self.proj.weight)
        nn.init.zeros_(self.proj.bias)

    def forward(self, count_hidden, genos_summary):
        """
        count_hidden:  [B, 128]   count head 中间层激活
        genos_summary: [B, 1024]  全局均值 Genos 特征
        """
        delta = self.proj(genos_summary)
        return count_hidden + delta
```

**注入位置**：Count Head 的 GELU 之后、最终 Linear 之前

```python
# transchrombp.py forward():
if genos_summary is not None and self.genos_count_proj is not None:
    count_layers = list(self.count_signal_head)
    h = count_layers[0](pooled)   # LayerNorm
    h = count_layers[1](h)        # Linear(256→128)
    h = count_layers[2](h)        # GELU
    h = self.genos_count_proj(h, genos_summary)  # ← P2 加法注入
    h = count_layers[3](h)        # Dropout
    count_signal = count_layers[4](h)  # Linear(128→1)
```

---

## 5. Probe 实验：独立验证 Genos 特征的信息量

在进行端到端训练实验之前/之后，我们通过独立的线性 probe 来评估 Genos 特征的信息含量：

### 5.1 Probe A — Count 残差岭回归

验证 Genos 特征能否解释 baseline 模型的 count 预测残差：

```python
# 用 G0 (baseline) 模型预测 logcount
# residual = true_logcount - pred_logcount
# 用 genos_global_mean 对 residual 做 5-fold OOF 岭回归

residual = true_logcount - pred_logcount
r2, pearson_r, pred_residual = ridge_oof_metrics(
    genos_global, residual, n_splits=5, alpha=1.0
)
```

### 5.2 Probe B — 互补性逻辑回归

验证 Genos 特征与 TransChromBP 已学特征的互补性：

```python
# 分别用三种输入做 peak vs nonpeak 分类:
# 1. 仅 TransChromBP encoded pooled 特征
# 2. 仅 Genos global_mean 特征
# 3. 两者拼接

auc_encoded = logistic_oof_auc(pooled, labels)
auc_genos = logistic_oof_auc(genos_global, labels)
auc_concat = logistic_oof_auc(
    np.concatenate([pooled, genos_global], axis=1), labels
)
```

### 5.3 Probe 结果

| Probe | 指标 | 值 | 含义 |
|-------|------|----|------|
| A: Count 残差岭回归 | R² | **0.0167** | Genos 仅能解释 1.67% 的 count 残差方差 |
| A: Count 残差岭回归 | Pearson r | **0.1987** | 弱但统计显著的相关性 |
| B: 仅 encoded | AUC | **0.9219** | TransChromBP 自身特征已高度饱和 |
| B: 仅 genos | AUC | **0.6997** | Genos 有独立的弱区分能力 |
| B: encoded + genos | AUC | **0.9018** | **拼接后反而下降！说明两者存在冲突而非互补** |

### 5.4 Phase 0 可行性测试（更早期的探索）

在 Genos 接入前，我们还做了更基础的可行性验证——多层双向 embedding 对 peak/nonpeak 的线性可分性：

```python
# genos_feasibility.py: 对 layer 0, 3, 6, 9, 12 分别做:
# 1. 提取双向 (fwd + rc_flipped) 平均 hidden states
# 2. 对 global_mean 和 center_1000bp_mean 做逻辑回归探针
# 3. 测量 center vs flank 位置方差比
# 4. 不同 batch size 下的 VRAM 和吞吐 benchmark
```

这个测试确认了 Genos embeddings 在 peak/nonpeak 二分类上有一定区分度（center CV AUC > 0.7），因此决定推进后续融合实验。但后来发现，**二分类的区分度不等于碱基分辨率 profile/count 的预测增益**。

---

## 6. 端到端训练实验结果

### 6.1 完整结果表

| 实验 | 类型 | 融合方式 | Val Peak JSD↓ | Val Count_r↑ | **Test Peak JSD↓** | **Test Count_r↑** | 判定 |
|------|------|----------|---------------|--------------|---------------------|---------------------|------|
| G0 (baseline) | online | 无 | 0.3337 | 0.7998 | **0.3163** | **0.8410** | 最强基线 |
| G1 (gate) | online Genos | 逐位置门控 | 0.3360 | 0.8010 | 0.3388 | 0.4689 | **count 崩溃** |
| G2 (mean) | online Genos | 全局均值门控 | 0.3414 | 0.7955 | 0.3387 | 0.6360 | **count 严重下降** |
| P0 (cached baseline) | cached | 无 | 0.3353 | 0.8033 | **0.3179** | **0.8384** | 可复现基线 |
| P2 (count-only) | cached Genos | count head 注入 | 0.3344 | 0.7984 | 0.3174 | 0.6037 | **count 断崖下滑** |

### 6.2 关键观察

1. **Validation 具有误导性**：G1/G2 在 val 上看起来接近 baseline，但 held-out test 上 count 严重崩塌。这提示 Genos 特征可能导致过拟合验证集
2. **Profile 相对 robust**：所有方案的 profile JSD 都没有大幅恶化，说明 profile head 对 Genos 注入有一定鲁棒性
3. **Count 极度敏感**：所有引入 Genos 的方案都出现 count 分支下降，P2 最严重（0.8384 → 0.6037，下降 28%）
4. **在线 vs 缓存都失败**：排除了"缓存语义错位"作为唯一原因——在线提取 (G1/G2) 同样失败

### 6.3 训练配置

```yaml
# 共同配置:
seed: 42
max_epochs: 10        # cached pilot 用 10 epoch 短训
validate_every_epochs: 1
early_stop_patience: 3

# 在线 Genos (G1/G2):
genos_branch:
  enabled: true
  model_path: /path/to/Genos-1.2B
  layer: 6
  hidden_size: 1024
  gate_bias_init: -2.0
  bidirectional: true
  freeze: true                    # Genos 权重完全冻结
  pool_mode: none                 # G1: 逐位置 / G2: mean

# 缓存 Genos (P2):
genos_cached:
  enabled: true
  genos_dim: 1024
  fusion_mode: count_only         # 直接注入 count head

# 数据增强 (所有实验共用):
data:
  input_len: 2114
  supervised_bp: 1000
  peak_max_jitter: 500            # ← 训练时有 jitter
  random_revcomp: true            # ← 训练时有反向互补增强
```

---

## 7. 失败原因分析

### 7.1 原因排序（按解释力从高到低）

#### 原因 1: 特征粒度太粗 — `global_mean` 对碱基分辨率任务不合适

当前所有 cached 实验只用了 `global_mean [1024]`，即将整个 2114bp 序列压缩为一个 1024 维向量。

- OCR benchmark 本质上是 **序列级二分类**，天然适合全局池化
- 但 ATAC-seq profile 是 **碱基分辨率轨迹预测**，需要局部 motif 组合和位置信息
- `global_mean` 丢失了所有位置结构，只保留"这个区域大概像不像开放区域"

#### 原因 2: 训练增强与缓存语义不对齐

- 缓存特征对应 canonical center（jitter=0, revcomp=false）
- 实际训练使用 `peak_max_jitter=500, random_revcomp=true`
- 模型在训练时看到的增强后序列，和缓存的 Genos 特征对应的不是同一个位置

#### 原因 3: P2 融合位置过于激进

将弱信号直接注入 count head 的倒数第二层，距离最终输出太近，没有缓冲层来吸收错误：

```
LayerNorm → Linear(256→128) → GELU → [注入点] → Dropout → Linear(128→1)
                                         ↑
                                    这里注入弱而粗的 global_mean
                                    太接近输出，count head 直接被破坏
```

#### 原因 4: 基线已经很强，增量空间极小

- G0 baseline: JSD=0.3163, count_r=0.8410
- TransChromBP 的 conv + dilated conv + Transformer + bias 分解已经吃掉了大部分可学习信号
- 对于如此强的基线，弱外部信号很难越过噪声门槛

#### 原因 5: Genos 的 causal 架构可能不适合此任务

Genos 是 causal LM（类 GPT），每个位置只能看到左侧上下文。虽然我们通过双向平均（forward + reverse-complement）来缓解，但这只是一种近似。**真正的双向模型**（如 BERT-style masked LM）可能更适合提取位置相关的局部特征。

---

## 8. 已排除的假设

| 假设 | 证据 | 结论 |
|------|------|------|
| Genos 模型/权重有问题 | OCR 复现 AUC=0.7535 ≈ 官方 0.7569 | **排除** |
| 只是随机波动 | 4 个独立实验 (G1/G2/P0/P2) 一致失败 | **排除** |
| 只是融合代码有 bug | 零初始化确保训练开始时等价于 baseline；probe 验证了特征本身 | **排除** |
| 只是缓存方案的问题 | 在线 (G1/G2) 同样失败 | **排除** |

---

## 9. 对"换用其他基座模型"的评估框架

基于以上实验，如果要评估换用其他基因组基座模型的潜力，需要回答以下问题：

### 9.1 当前实验已确认的约束

1. **任务特征**：碱基分辨率 profile 预测 + count 回归，不是序列级分类
2. **基线很强**：TransChromBP 已经是 Transformer + bias 分解的强模型，增量空间有限
3. **Count head 对外部注入极度敏感**
4. **全局池化特征（global_mean）信息量不足**（R² = 0.0167）
5. **特征拼接产生冲突而非互补**（concat AUC < encoded_only AUC）

### 9.2 换用其他基座模型需要满足的条件

1. **必须提供位置分辨率特征**，不能只有序列级 embedding
2. **最好是双向架构**（BERT-style / bidirectional），而非 causal LM
3. **hidden states 需与碱基分辨率任务相关**，而非仅在序列级分类上表现好
4. **必须用 probe 先验证互补性**：concat(encoded, new_model) AUC > encoded_only AUC
5. **融合方式需要调整**：
   - 避免直接注入 count head
   - 优先尝试 backbone 中段的轻量调制
   - 考虑 cross-attention 而非简单加法/门控

### 9.3 候选模型对比维度

| 维度 | Genos-1.2B（已测试） | 理想候选模型 |
|------|---------------------|-------------|
| 架构 | Causal LM (单向) | 双向更佳 (MLM / bidirectional) |
| 预训练目标 | Next-token prediction | MLM, span corruption, 或 ATAC-specific |
| 位置分辨率 | 有（但 global_mean 后丢失） | 需保留并使用 |
| 与 ATAC 任务关联 | 弱（OCR 分类 OK，profile 预测无增益） | 需在类似任务上验证 |
| 参数量 | 1.2B（推理成本高） | 轻量级可能更实用 |
| Token 粒度 | 单碱基 | k-mer 可能更好匹配 motif |

---

## 10. 附录：完整文件清单

### 核心模型代码
- `transchrombp/models/genos_adapter.py` — GenosFeatureExtractor, GenosGatedAdapter, GenosSummaryFiLM, GenosCountProj
- `transchrombp/models/transchrombp.py` — 主模型 forward()，包含 3 个 Genos 注入点

### 训练集成
- `transchrombp/training/train_ddp.py` — DDP 训练循环中的 Genos 特征提取与前向传播

### 数据管线
- `transchrombp/data/real_data.py` — ChromBPNetBigWigDataset 的 Genos cache 加载与验证
- `transchrombp/scripts/build_genos_summary_cache.py` — 离线缓存构建（memmap + manifest）

### 验证脚本
- `transchrombp/scripts/run_genos_summary_probe.py` — Count 残差 ridge + 互补性 logistic probe
- `scripts/genos_feasibility.py` — Phase 0 可行性验证（多层 embedding + throughput benchmark）

### 分析报告
- `reports/genos_no_positive_gain_analysis_20260323.md` — 详细失败原因分析
- `reports/genos_cached_fusion_status_20260323.md` — 完整结果表与 Go/No-Go 判定
- `reports/genos_and_6002_run_analysis_20260323.md` — 跨服务器实验指标对比

---

*报告日期: 2026-03-24*
*TransChromBP 项目组*
