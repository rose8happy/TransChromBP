# TransChromBP V2 读出头改进与训练优化消融实验

> 文档版本：v2.1（2026-03-20，含 Codex 审查修复）
> 实验代号：v2fix
> 实验状态：**准备就绪，待部署启动**

---

## 一、为什么要做这轮实验

### 1.1 背景：Transformer 已验证有效，但收益有限

在 `ablation_tf_20260318` 消融实验中，我们首次验证了 Transformer 对 TransChromBP 模型的贡献。
核心结论是 **Transformer 有效但收益偏小**：

| 指标 | V2-full (有 TF) | V2-noTF (无 TF) | 差值 | 改善方向 |
|---|---|---|---|---|
| peak JSD full ↓ | **0.31419±0.00012** | 0.32393±0.00021 | -0.00974 | ✅ 好 |
| peak JSD debiased ↓ | **0.31425±0.00012** | 0.32403±0.00021 | -0.00978 | ✅ 好 |
| count pearson ↑ | **0.83676±0.00531** | 0.82503±0.00489 | +0.01173 | ✅ 好 |
| full-debiased gap ↓ | **0.00006±0.00001** | 0.00010±0.00000 | -0.00004 | ✅ 更少 shortcut |

- Transformer 在 profile JSD 上带来了约 **3%** 的相对改善（(0.32393-0.31419)/0.32393 ≈ 3.0%）
- Count pearson 提升约 **1.4%**
- Bias shortcut 没有被重新引入（gap 极小）

这个收益是**真实的**但**不够显著**，在统计上三个 seed 的效应一致但绝对值偏小。

### 1.2 核心问题：收益为什么偏小？

ChatGPT Pro 对此做了深入的四层分析，结合我们自己的 scale 轨迹数据和训练曲线，
核心判断是：

> **问题不像是 Transformer 没有学到东西，而更像是"学到的东西没有被最后的 head 和训练设置充分兑现成指标"。**

支撑这个判断的三个关键证据：

#### 证据 1：训练未充分收敛

六个 run（3 full + 3 noTF）的 best epoch 分布：

| run | best_epoch | 训练总 epoch |
|---|---|---|
| full_s42 | 29 | 30 |
| full_s1234 | 29 | 30 |
| full_s2024 | 28 | 30 |
| noTF_s42 | 30 | 30 |
| noTF_s1234 | 29 | 30 |
| noTF_s2024 | 27 | 30 |

Best epoch 集中在 27-30，且验证曲线到 epoch 30 仍在下降。这意味着模型（尤其是 Transformer）
**还没完全展开学习**。Transformer 模块通常比纯 conv 更依赖足够的训练长度和 LR schedule。

#### 证据 2：Count head 的全长 pooling 稀释了 Transformer 信号

当前 count 预测流程：
```
encoded = transformer(tokens)           # [B, 2114, 256]
pooled = encoded.mean(dim=1)            # [B, 256]  ← 全长 2114 token 平均
count = count_signal_head(pooled)       # [B, 1]
```

但实际监督的 count 对应**中心 1000bp 窗口**。全长 2114 token 的 mean pooling 会把
**两端各 557 token 的上下文**一起平均进来，稀释中心信号。Profile head 已用
`center_crop_1d(output_len=1000)` 做了裁剪，但 count head 没有对齐。

更重要的是：Transformer 的 long-range context 往往只对**某些特定位置**有帮助，
不是整条序列都同等重要。全局 mean pooling 会把"少数真正有信息的位置"冲淡。

#### 证据 3：Profile head 过于简单

当前 profile 预测：
```
profile_signal = Linear(256, 1)(encoded)    # 每个位置独立映射
profile_signal = center_crop_1d(profile_signal, 1000)
```

这是一个逐位置的线性映射——Transformer 输出的每个 token 独立通过 `Linear(d_model, 1)`。
对于 profile 这种高分辨率峰形（sharp peaks），这有点"太平"。

Transformer 负责建模远程依赖（比如"500bp 之外有一个增强子，所以这里的可及性应该更高"），
但最终的峰形精修（"这个 peak 的具体形状应该更尖"）需要**局部卷积**来完成。
当前 Linear head 不具备这种局部精修能力。

### 1.3 本轮实验的核心假设

基于以上分析，本轮实验的核心假设是：

> **Transformer 已经学到了有用的长程上下文信息，但当前的读出头（count head 和 profile head）
> 以及训练设置无法充分将这些信息兑现为指标改善。**

如果这个假设正确：
- 改进 count head pooling → count_pearson 提升
- 改进 profile head → peak JSD 提升
- 延长训练 → 两者同时提升

如果这个假设不正确（即 Transformer 确实没学到什么）：
- 改进读出头不会带来改善
- 延长训练也不会有显著变化
- 这时需要回到 Transformer 本体进行结构改动（降采样 TF、gated residual 等）

---

## 二、实验设计总览

### 2.1 设计原则

| 原则 | 说明 |
|---|---|
| **单变量** | 每个实验只改一个因素，其余与 baseline 完全一致 |
| **共享 bias pretrain** | 所有实验复用 `ablation_tf_20260318` 的 shared bias checkpoint |
| **3 seeds** | 42, 1234, 2024（与现有 ablation 一致）|
| **40 epochs** | 从 30 提高（修复训练不充分问题）|
| **默认关闭** | 所有新功能通过 config 参数控制，默认 off，完全向后兼容 |

### 2.2 实验清单

```
第零批 ─── 零训练成本 ───────────────────
  E: Checkpoint Soup（对已有 checkpoint 做 weight average）

第一批 ─── 核心实验（40 epoch）──────────
  A: Freeze 修复 + 40ep（正确性修复 + 新 baseline）
  B: Count Center Pooling + 40ep（缩小 pooling 范围到中心 1000 token）
  F: Count Attention Pooling + 40ep（让模型学习 pooling 权重）    [NEW]
  G: Profile Refinement Conv + 40ep（轻量 depthwise conv 精修峰形）[NEW]

第二批 ─── 可选（视第一批结果决定）──────
  C: profile_scale_init=0.1（已被 scale 轨迹否定，效果可忽略）
  D: Bias BN→GN（前提假设已被数据否定）
```

### 2.3 为什么新增 F 和 G

ChatGPT Pro 在第二层建议中提出了"改读出头，不要先堆更大的 Transformer"的策略。
我们对比分析了这些建议与已有计划的重叠度：

| ChatGPT Pro 建议 | 已有计划 | 动作 |
|---|---|---|
| 40 epoch 复试 | ✅ 已有 | 确认 |
| Checkpoint soup / EMA | ✅ 实验 E | 确认，EMA 暂缓 |
| Count head: center pooling | ✅ 实验 B | 确认 |
| Count head: **attention pooling** | ❌ 新增 | → **实验 F** |
| Profile head: **轻量 conv refinement** | ❌ 新增 | → **实验 G** |
| profile_scale_init 0.1 | ✅ 实验 C (降级) | 确认降级 |
| Frozen bias eval 语义 | ✅ 实验 A | 确认 |
| 降采样 TF / gated residual / stratified training | 第 4 层 | **暂缓** |

F 和 G 的价值在于：它们直接作用于"Transformer 输出到最终指标"的瓶颈环节，
改动量小、风险低、与已有实验正交，是"低成本高回报"的改动。

---

## 三、各实验详细设计与预期

### 实验 E：Checkpoint Soup（零训练成本）

**原理**：对训练末尾连续 K 个 epoch 的模型权重做算术平均（Stochastic Weight Averaging 的后处理版本）。
当模型在训练末尾处于 loss landscape 的平坦区域时，相邻 epoch 的权重虽然在 weight space 中位置不同，
但性能接近。对它们取平均通常能找到更好的泛化点。

**操作**：
```bash
# 对每个 baseline run 做两种 soup
python3 scripts/checkpoint_soup.py --ckpt-dir outputs/checkpoints/full_s42/ --mode last_k --k 5 --output soup_last5.pt
python3 scripts/checkpoint_soup.py --ckpt-dir outputs/checkpoints/full_s42/ --mode top_k_jsd --k 3 --output soup_top3.pt
```

说明：
- 新版 orchestrator 优先按 `peak.profile_target_jsd_full_mean` 选 `soup_top3`
- 对旧 run 若该字段不存在，会自动回退到 `peak.profile_target_jsd_debiased_mean`，再回退到 `peak.loss_profile`

**预期**：
- soup 通常比单点 best 提升 0.5–2%
- 如果 soup 反而变差 → 模型末尾在震荡，weight space 不平坦

---

### 实验 A：Freeze 修复 + 40 epoch

**问题**：当前 `freeze_core()` 只设 `requires_grad=False`，但 bias branch 里的
`BatchNorm1d` 和 `Dropout` 在 `model.train()` 下仍会更新 running stats / 随机 drop。

**修复措施**（已在代码中实现）：
1. 新增 `_bias_branch_frozen` 状态标记
2. `freeze_bias_core()` 同时调用 `bias_branch.eval()` 强制关闭 BN/Dropout
3. Override `train()` 方法：当 bias 被冻结时，保持 `bias_branch.eval()`
4. Forward 中冻结时用 `torch.no_grad()` 包裹 bias 调用

**为什么仍然做**（尽管 gap 只有 0.00006）：
- 这是一个**正确性修复**，不管效果大小都应该合入
- 40 epoch 的结果可以作为新的干净 baseline
- 如果效果不大 → 证实当前 4 层 dilated conv + BN 在 eval vs train 差异有限

**对照**：

| 组别 | freeze 语义 | max_epochs |
|---|---|---|
| baseline (已完成) | BN/Dropout 仍 train | 30 |
| A (本实验) | 强制 eval + no_grad | 40 |

**预期**：
- debiased 指标的 seed 间方差可能减小
- profile/count 指标持平或略好
- 主要价值是提供一个干净的 40 epoch baseline

---

### 实验 B：Count Head Center Pooling + 40 epoch

**问题**：当前 count pooling 对全长 2114 token 做 mean，但监督的 count 对应中心 1000bp。

**改动**：
```python
# 原始
pooled = encoded.mean(dim=1)                              # [B, 2114, D] -> [B, D]

# 改为
start = (encoded.size(1) - self.output_len) // 2          # (2114-1000)//2 = 557
pooled = encoded[:, start:start+self.output_len, :].mean(dim=1)  # [B, 1000, D] -> [B, D]
```

**为什么 center 而不是其他方式**：
- 最简单、最符合直觉：count 监督的是中心 1000bp，pooling 也对齐到中心 1000bp
- Profile head 已经做了同样的 center crop，count head 理应对齐
- 不引入额外参数，不增加过拟合风险

**对照**：

| 组别 | pooling 范围 | token 数 |
|---|---|---|
| baseline | 全长 | 2114 |
| B (center pool) | 中心 1000 | 1000 |

**预期**：
- `count_pearson_full` 提升（中心信号更聚焦）
- `count_mae_full` 下降
- profile 指标不应受影响（profile head 不变）
- 如果 count 反而变差 → Transformer 学会了利用全长上下文做 count，这本身是重要发现

---

### 实验 F：Count Head Attention Pooling + 40 epoch [NEW]

**动机与 B 的区别**：
实验 B 使用"硬"截取（固定取中心 1000 token），假设所有有用信息都在中心。
但实际上：
- 不同 region 的信息分布可能不同
- Transformer 可能让**特定位置**学到了跨区域上下文信息
- 全局 mean pooling 把"少数真正有信息的位置"冲淡

Attention pooling 让模型**自己学习**每个位置对 count 预测的贡献权重：

```python
# 学习的 attention pooling
attn_logits = self.count_pool_proj(encoded).squeeze(-1)    # [B, L] — 每个位置一个标量 logit
attn_weights = torch.softmax(attn_logits, dim=1)           # [B, L] — 归一化为概率分布
pooled = torch.sum(encoded * attn_weights.unsqueeze(-1), dim=1)  # [B, D] — 加权求和
```

**新增参数**：`count_pool_proj = nn.Linear(d_model, 1)` — 只有 257 个参数（256+1）。

**与 B 的实验逻辑关系**：

```
baseline (full mean pooling)
    ├── B (center pool): 限定范围但权重均匀
    └── F (attention pool): 不限范围但权重学习

如果 B 好、F 也好 → 信号集中在中心，且 attention 可以学到
如果 B 好、F 差 → 信号确实在中心，attention 的自由度反而有害
如果 B 差、F 好 → 信号不只在中心，attention 能捕获更灵活的分布
如果 B 差、F 差 → 全长上下文对 count 确实有正面贡献
```

**诊断工具**：训练完成后可可视化 attention weights 的空间分布，
观察模型是否学会聚焦中心、或者有其他有趣的 pattern。

**Config 改动**：`heads.count_pool_mode: "attention"`

---

### 实验 G：Profile Head 轻量 Conv Refinement + 40 epoch [NEW]

**动机**：

当前 profile head 是 `Linear(d_model, 1)`——每个位置的 transformer 输出独立映射到一个标量。
这等价于一个 1×1 的"卷积"，没有任何**局部空间交互**。

对于 profile 预测而言，最终的峰形（peak shape）是一个高分辨率的空间信号。
Transformer 建模的是 long-range dependency（"远处的增强子影响这里的可及性"），
但峰形的精细结构（"peak 顶部该多尖、边缘该多陡"）更多取决于**局部 motif 组合**。

最佳设计是 **分工合作**：
- Transformer 负责 long-range context → 确定"这个位置大概应该有多高"
- Conv 负责 local refinement → 精修"这个 peak 的具体形状"

**代码实现**：

```python
# 当 use_profile_refine=True 时
self.profile_refine = nn.Sequential(
    nn.Conv1d(d_model, d_model, kernel_size=5, padding=2, groups=d_model),  # depthwise
    nn.GELU(),
    nn.Conv1d(d_model, d_model, kernel_size=1),  # pointwise
    nn.GELU(),
)
self.profile_signal_head = nn.Conv1d(d_model, 1, kernel_size=1)

# forward 中
feat = encoded.transpose(1, 2)                          # [B, D, L]
feat = feat + self.profile_refine(feat)                  # ← 残差连接
profile_signal = self.profile_signal_head(feat).squeeze(1)  # [B, L]
```

**设计选择详解**：

| 设计点 | 选择 | 理由 |
|---|---|---|
| Conv 类型 | **Depthwise separable** | 参数量极小：`d_model×5 + d_model² ≈ 67K`（d_model=256 时） |
| Kernel size | **5** | 覆盖 ±2bp 局部上下文，足够精修峰形但不会太大 |
| 连接方式 | **残差** (`feat + refine(feat)`) | 保证模型至少不比 baseline 差——如果 refine 无用，可以学为零 |
| 最终 head | **Conv1d(d,1,1)** 替代 Linear(d,1) | 保持 channel-first 一致性，数学上等价 |
| 激活函数 | **GELU** | 与 Transformer 和 Conv stem 一致 |

**参数预算**：

| 组件 | 参数量 | 计算 |
|---|---|---|
| Depthwise Conv1d(256,256,5,groups=256) | 1,536 | 256 × (5+1) |
| Pointwise Conv1d(256,256,1) | 65,792 | 256 × 256 + 256 |
| Conv1d(256,1,1) head | 257 | 256 + 1 |
| **总新增** | **~67,585** | 占模型总参数 ~1% |

对比：模型总参数约 ~6.7M（d_model=256, n_layers=6），新增 67K 约占 1%，不会显著增加过拟合风险。

**Config 改动**：`heads.use_profile_refine: true`

**预期**：
- `peak.profile_target_jsd` 改善（峰形更尖锐，JSD 更低）
- 如果效果不大 → Linear head 已经足够，瓶颈不在局部峰形
- 如果效果明显 → 说明 Transformer 的全局信息确实需要局部精修才能兑现

---

## 四、代码改动总览

所有代码改动都在 `tmp_remote_edit/` 目录下，待部署到远端。

### 4.1 模型代码

**文件**：`src/transchrombp/models/transchrombp.py`

| 改动位置 | 内容 | 影响 |
|---|---|---|
| `__init__` 参数 | 新增 `count_pool_mode="attention"` 选项和 `use_profile_refine` 参数 | 向后兼容（默认 off） |
| `__init__` 模块创建 | 新增 `profile_refine` (depthwise+pointwise conv)、`count_pool_proj` (attention projection) | 仅当启用时创建 |
| `_pool_for_count()` | 新增 attention pooling 分支 | 3 行新代码 |
| `forward()` | Profile head 分流：有 refine 时走 conv+残差，否则走 Linear | if/else 分支 |
| `build_transchrombp_from_config()` | 传递 `use_profile_refine` 和更新的 `count_pool_mode` | config 接口 |

**关键设计决策**：
- 两个新功能**默认关闭**，现有训练不受任何影响
- Profile refinement 用**残差连接**，最坏情况退化为 baseline
- Attention pooling 只新增 257 个参数，不影响训练速度

### 4.2 Config 文件

| 文件 | 用途 | 相对 baseline 的改动 |
|---|---|---|
| `v2fix_baseline.yaml` | 实验 A (freeze fix only) | 无改动（与 teacher_v2 相同） |
| `v2fix_center_pool.yaml` | 实验 B | `heads.count_pool_mode: center` |
| `v2fix_attn_pool.yaml` **[NEW]** | 实验 F | `heads.count_pool_mode: attention` |
| `v2fix_profile_refine.yaml` **[NEW]** | 实验 G | `heads.use_profile_refine: true` |
| `train_ablation_v2_main.yaml` | 旧主训练 config | 保留 `peak.loss_total` 选模，仅用于回溯 |
| `train_ablation_v2_main_profile_select.yaml` **[NEW]** | **推荐主训练 config** | `max_epochs: 40`, `early_stop_patience: 10`, `best_metric: peak.profile_target_jsd_full_mean` |
| `train_ablation_v2_main_throughput_smoke.yaml` **[NEW]** | 吞吐/短 smoke 预设 | 在主训练策略基础上改 `batch_size_per_gpu: 24`，**不用于严格主结论** |

### 4.3 Orchestrator 脚本

**文件**：`TransChromBP/scripts/run_v2fix_ablations.sh`

改动：
- 第一批从 A/B/C → **A/B/F/G**
- 第二批从 D → **C/D**（均可选）
- 新增 `MODEL_F` 和 `MODEL_G` config 引用
- Port 分配扩展到 12 个（4 组 × 3 seeds）
- 新增 `TRAIN_STRATEGY`：`scientific`（默认）/`throughput_smoke`/`legacy_loss_total`
- 新增 `--groups` / `--seeds`，支持按组、按 seed 分阶段启动

**使用方式**：
```bash
# DRY RUN：只打印命令不执行
DRY_RUN=1 bash scripts/run_v2fix_ablations.sh --batch 1

# 正式执行第一批（A/B/F/G × 3 seeds = 12 runs）
bash scripts/run_v2fix_ablations.sh --batch 1

# 吞吐/短 smoke 预设（基于 benchmark 的 bs/gpu=24 建议）
TRAIN_STRATEGY=throughput_smoke bash scripts/run_v2fix_ablations.sh --batch 1

# Soup 评估（零成本）
bash scripts/run_v2fix_ablations.sh --soup-only

# 推荐首次正式启动：只跑 A/B 的单 seed 预热
TRAIN_STRATEGY=scientific bash scripts/run_v2fix_ablations.sh --batch 1 --groups A,B --seeds 42

# A/B 单 seed 稳定后，再扩成 3 seeds
TRAIN_STRATEGY=scientific bash scripts/run_v2fix_ablations.sh --batch 1 --groups A,B --seeds 42,1234,2024

# A/B 跑完后，再补 F/G
TRAIN_STRATEGY=scientific bash scripts/run_v2fix_ablations.sh --batch 1 --groups F,G --seeds 42,1234,2024
```

**推荐启动顺序**：
1. 先在真实训练环境里做一次 `evaluate_checkpoint.py --split valid` smoke
2. 跑 `--soup-only`，只用 `valid` 做零成本筛查
3. 跑 `A/B + seed42`
4. 若稳定，再把 `A/B` 扩成 3 seeds
5. 最后再跑 `F/G`

---

## 五、训练配置详解

### 5.1 关键超参数

| 参数 | baseline (30ep) | 本轮 (40ep) | 理由 |
|---|---|---|---|
| max_epochs | 30 | **40** | Best epoch 全贴 27-30，给更多收敛空间 |
| early_stop_patience | 6 | **10** | 与更长训练匹配，避免过早停止 |
| best_metric | peak.loss_total | **peak.profile_target_jsd_full_mean** | JSD 直接反映 profile 质量（Codex 审查修复） |
| learning_rate | 5e-4 | 5e-4 | 不变（cosine schedule 会自动调整） |
| batch_size_per_gpu | 16 | 16 | 不变 |
| debiased_profile_weight | 2.0 | 2.0 | Signal branch 的独立 profile 监督 |
| debiased_count_weight | 0.0 | 0.0 | 保持不变，避免引入新变量（见 §6.4）|
| compile | false | false | `compile` 单独做 throughput benchmark，不混入主 ablation |

### 5.2 训练策略分层

为了吸收 2×A6000 benchmark 的结论，同时避免把“系统优化”和“科学结论”混为一谈，
本轮后续实验明确拆成两条策略线：

| 策略 | 用途 | 关键参数 | 是否用于主结论 |
|---|---|---|---|
| `scientific` | **默认主 ablation** | 2GPU, bf16, `bs/gpu=16`, `best_metric=peak.profile_target_jsd_full_mean` | **是** |
| `throughput_smoke` | 短 smoke / pilot | 2GPU, bf16, `bs/gpu=24`, 同样的 JSD 选模 | 否 |
| `legacy_loss_total` | 回溯旧口径 | 保留 `peak.loss_total` | 否 |

这样做的原因是：
- benchmark 已证明 `bs/gpu=24/32` 更适合作为**吞吐类短实验**默认点；
- 但 global batch 从 32 提到 48/64 会改变训练动力学，不应直接混进严格的 A/B/F/G 对比。

### 5.3 LR Schedule 变化

使用 cosine schedule + warmup_ratio=0.06：
- 30 epoch：warmup ≈ 1.8 epoch → 从 epoch 2 开始衰减 → epoch 30 到 min_lr=5e-5
- 40 epoch：warmup ≈ 2.4 epoch → 从 epoch 3 开始衰减 → epoch 40 到 min_lr=5e-5

更长的 schedule 意味着：
- **中段 LR 更高**（衰减更慢）→ Transformer 有更多"高 LR 探索阶段"
- **末段 LR 衰减更缓**→ 末尾 fine-tuning 更充分
这对 Transformer 的学习特别有利，因为 attention 模式通常需要更多训练步才能稳定。

---

## 六、预期结果与判读逻辑

### 6.1 各实验预期收益

| 实验 | 主要改善指标 | 预期幅度 | 信心 |
|---|---|---|---|
| E (soup) | peak JSD | 0.5-2% 相对改善 | 高（成熟技术） |
| A (freeze+40ep) | 稳定性 + 整体 | seed 方差减小 + 末尾收敛 | 中（正确性修复） |
| B (center pool) | **count_pearson_debiased** | 明确改善 | 高（逻辑清晰） |
| F (attn pool) | **count_pearson_debiased** | 可能 ≥B | 中（依赖 TF 信号分布） |
| G (profile refine) | peak JSD | 局部峰形精修 | 中高（成熟设计模式） |

### 6.2 判读决策树

```
40ep 的 best_epoch 是否还贴尾部（38-40）？
├── 是 → 再补 50ep 单 seed 探针
└── 否 → 训练已充分，关注指标

B 的 count_pearson 是否提升？
├── 是 → center pool 有效，与 F 对比选更优
│        ├── F > B → attention 更灵活，用 F
│        ├── F ≈ B → centeer pool 足够简单，用 B
│        └── F < B → center 的硬约束更好，用 B
└── 否 → 全长上下文对 count 有正面贡献
          → 尝试 center+global 拼接或保持 full

G 的 peak JSD 是否改善？
├── 是 → profile refine 有效
│        → 与最佳 B/F 叠加测试组合效果
└── 否 → 瓶颈不在局部峰形
          → 考虑第 4 层改动（gated residual 等）
```

### 6.3 合并决策规则

每个实验跑完 3 seeds 后：

1. **必要条件**：`peak.profile_target_jsd_debiased_mean` 均值**不劣于** baseline
   - 判定：paired t-test（同 seed 对比），p < 0.1 为显著劣化
2. **充分条件**（满足任一即合入）：
   - peak JSD debiased 均值显著更好（p < 0.1）
   - count_pearson 显著提升 **且** peak JSD 不劣
   - seed 间标准差显著减小
3. **组合叠加顺序**：A → B/F → G → 确认每步组合不冲突

### 6.4 间接监督说明（B/F 特别注意）

当前 `debiased_count_weight: 0.0`，即 count_signal 没有被**直接**监督。
梯度需要通过 logsumexp 融合反传。这意味着：

- **B/F 的主判读指标应为 `count_pearson_debiased`**（而非 `count_pearson_full`），
  因为前者直接反映 signal branch 的 count 质量，不受 bias 融合稀释
- 如果 B/F 提升 `count_pearson_debiased` 但 `count_pearson_full` 不变 → pooling 改进有效，
  但需要加 `debiased_count_weight > 0` 来放大效果
- 如果 B/F 均无提升 → **不能直接下结论说 pooling 没用**，可能是 signal count 路径
  本身没被充分约束。此时应追加一轮加 `debiased_count_weight: 0.1` 的对照实验

---

## 七、资源与时间规划

### 7.1 硬件

| 服务器 | GPU | 显存 | 估计单 run 时间 |
|---|---|---|---|
| 6000 | 2×NVIDIA A6000 | 48GB each | ~2.3h (40ep, DDP) |
| 6002 | 1×NVIDIA RTX 3080 | 10GB | ~19h (40ep, 单卡) |

### 7.2 时间线

| 阶段 | 操作 | 预计耗时 |
|---|---|---|
| 第零批 | E: soup 评估 (6 runs × ~70s) | ~30min |
| 第一批 | A+B+F+G (4 组 × 3 seeds = 12 runs) | 6000: ~28h 串行; 或分流到 6002 |
| 分析 | 汇总 + 对比 + 统计检验 | ~2h |
| 第二批（可选）| C/D | 仅在异常情况补做 |

### 7.3 当前 6002 占用

6002 正在运行 bias-safe ablation (B3→B4)，预计还需 8-27h。
可以等 B3/B4 完成后再分流部分 v2fix 实验到 6002。

---

## 八、风险与应急措施

| 风险 | 影响 | 应对 |
|---|---|---|
| 40 epoch 仍不够（best 仍贴尾） | 训练上限未探明 | 补 50 epoch 单 seed 探针 |
| B count pooling 反而变差 | TF 利用了全长上下文做 count | 改用 F (attention pooling) |
| F attention pooling 过拟合 | 自由度太大 | 加 dropout 或约束 attention 分布 |
| G profile refine 无效 | 瓶颈不在局部峰形 | 试更大 kernel_size (7/9)，或加第二层 conv |
| A freeze 修复无差异 | BN drift 确实可忽略 | 仍合入（正确性修复）|
| F 和 G 叠加时相互干扰 | 无法判断各自贡献 | 先单独验证，再做叠加 |
| 代码部署出错 | 远端 import 失败 | 部署后先 smoke test: `python -c "from transchrombp.models.transchrombp import TransChromBP"` |

---

## 九、与后续工作的关系

本轮实验属于 ChatGPT Pro 四层优先级策略中的**第 1-2 层**（训练预算 + 读出头改进）。

如果本轮实验成功放大了 Transformer 收益 → 直接用于论文，不需要动 Transformer 本体。

如果本轮实验收益仍然有限 → 进入第 3-4 层：

| 后续方向 | 描述 | 工程量 |
|---|---|---|
| Gated residual | `encoded = tokens + gate * (TF(tokens) - tokens)` | 中 |
| 降采样 TF | conv 保留高分辨率，TF 在粗 token 上建模 | 大 |
| Stratified training | 按"需要长程上下文"的难度分层 region | 前置分析量大 |
| EMA | 训练中维护 EMA 参数 | 改训练循环 |

这些改动的优先级严格低于本轮的读出头改进。
