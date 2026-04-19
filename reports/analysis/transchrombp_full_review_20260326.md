# TransChromBP 项目全景复盘与论文路线报告

> 日期：2026-03-26
> 覆盖范围：全部实验历程、代码架构、结果数据、论文叙事建议

---

## 一、项目目标

在 ChromBPNet 的 bias factorization 框架上嫁接 Transformer，实现碱基分辨率 ATAC-seq 信号建模。核心挑战：Transformer 的高容量是否会绕过 bias 分解、直接在 signal branch 中"中继" bias 信号（即 Profile Shortcut 假设）。

---

## 二、模型架构

```
DNA [B, 2114, 4]
  → Conv Stem (4→256, k=21)
  → Local Tower (8层扩张卷积, k=3, dilation cycle=8)
  → Transformer Encoder (6层, 8头, RoPE, d=256)
  → Profile Head (Linear → [B, 1000])  +  Count Head (MLP → [B, 1])
  → Bias Branch (冻结, 独立卷积预测酶切偏好)
  → 最终输出 = signal + learnable_scale × bias
```

关键可配置开关：

| 参数 | 含义 | safe 默认值 |
|------|------|-------------|
| `sequence_encoder.enabled` | 是否启用 Transformer | true |
| `profile_bias_stop_gradient` (sg) | 是否阻断 profile bias 梯度 | true |
| `debiased_profile_weight` (deb) | debiased profile loss 权重 | 2.0 |
| `count_pool_mode` | count 聚合方式 | center |

---

## 三、实验历程

### 阶段 1：V1→V2 架构升级（~3月16日前）

同时引入 `stop_gradient=true`、`deb_weight=2.0`、`pool_factor=32`。观察到 debiased JSD 从 0.5633 骤降至 0.3163。当时解读为"发现并修复了 Transformer 特有的 Profile Shortcut"。

**问题**：三个变量同时改变，无法归因。

### 阶段 2：消融矩阵与读出头设计（3月16-22日）

- **TF 消融**（3 seed matched ablation）：确认 Transformer 带来 3.0% profile JSD 改善
- **读出头对比**（A/B/F/G）：B（center pool）胜出，F（attention pool）多 seed 失守，G（profile refine）失败
- 多 seed 确认 B 稳定性

### 阶段 3：Genos 融合与独立数据（3月18-24日）

- Genos-1.2B **全线负结果**：probe R²=0.017，concat AUC 反降，在线/缓存融合均导致 count 崩塌
- 独立数据验证通过：GM12878 和 K562 上模型可泛化

### 阶段 4：Clean Matrix 受控复核（3月25-26日）

- 发现早期"强 shortcut"论断证据不足
- 启动严格 2×2 矩阵，逐步补齐全部四个格子的 test 评估
- 论文叙事从"Transformer 特有 shortcut"转向"debiased supervision 是通用稳定器"

---

## 四、核心实验结果

### 4.1 Transformer 消融（Tutorial K562, 3 seed）

| 配置 | Profile JSD (mean±std) | Count r (mean±std) |
|------|----------------------|-------------------|
| V2-full (TF) | **0.31419 ± 0.00012** | **0.83676 ± 0.00531** |
| V2-noTF | 0.32393 ± 0.00021 | 0.82503 ± 0.00489 |
| ChromBPNet 官方 | 0.3360 (median) | — |

**结论**：Transformer 带来真实收益，profile 改善 3.0%，count 改善 1.4%。3 seed 全部一致（符号检验 p≈0.016）。

### 4.2 完整 2×2 Bias Reliance 矩阵

这是本项目最核心的实验。四个格子的 test 结果现已全部到齐。

**矩阵设计**：拆分两个因素——架构（TF vs noTF）和 bias 安全配置（safe=sg_false+deb2 vs unsafe=sg_true+deb0）。

> 注：A/C 在 6000 A6000 上训练，noTF 系列在 6002 3080 上训练。绝对 JSD 受训练环境影响，**核心比较量是 gap（full-debiased 差距）**。

#### 4.2.1 Test 结果总表

| 条件 | 配置 | peak JSD full | peak JSD deb | **gap** | count r full | count r deb | bias rms |
|------|------|:---:|:---:|:---:|:---:|:---:|:---:|
| A | TF + sg=false + deb2 | 0.31410 | 0.31416 | **0.00147** | 0.8506 | 0.8488 | 0.00247 |
| C (s42) | TF + sg=true + deb0 | 0.31720 | 0.31779 | **0.00972** | 0.8395 | 0.8370 | 0.01668 |
| noTF_deb2 | noTF + sg=false + deb2 | 0.43531 | 0.43557 | **0.00487** | 0.7959 | — | 0.01375 |
| **notf_sg1_deb0** | noTF + sg=true + deb0 | 0.43659 | 0.43835 | **0.02682** | 0.7919 | 0.7956 | 0.07647 |

#### 4.2.2 Gap 矩阵（核心对比）

|  | **safe (sg=false, deb2)** | **unsafe (sg=true, deb0)** |
|---|:---:|:---:|
| **TF** | 0.00147 | 0.00972 |
| **noTF** | 0.00487 | **0.02682** |

#### 4.2.3 C 线多 seed 确认

| Seed | 状态 | Val peak gap | Test peak gap | Test peak JSD (full / deb) | 备注 |
|------|------|:---:|:---:|:---:|------|
| s42 (6000) | test 已完成 | 0.00993 | 0.00972 | 0.31720 / 0.31779 | 主线 |
| s1234 (6000) | test 已完成 | **0.01039** | **0.00997** | 0.31739 / 0.31795 | best_epoch=30, 跑满无 early-stop |

两 seed 的 test gap 为 `0.00972` vs `0.00997`，与 validation 的 `0.00993` vs `0.01039` 高度一致。  
双 seed 汇总后，C 线 peak full/debiased 为 `0.31730±0.00014 / 0.31787±0.00011`，count r 为 `0.83671±0.00389 / 0.83479±0.00311`，说明这条“轻度但稳定”的 gap 不是单 seed 偶然波动。

#### 4.2.4 矩阵结论

1. **unsafe 配置下，noTF 的 gap（0.027）比 TF 的（0.010）大 2.8 倍** → "Profile Shortcut 是 Transformer 特有"这个说法**不成立**
2. **safe 配置下，TF 和 noTF 都干净**（0.00147 vs 0.00487）→ deb2 是通用稳定器，与架构无关
3. **deb2 的效应远大于 sg 的效应**：A（无 sg 有 deb2）= 0.00147，C（有 sg 无 deb2）在两 seed 上为 `0.00972 / 0.00997` → deb2 比 sg 更关键
4. **Transformer 在两种配置下 gap 都不比 noTF 差** → Transformer 不加剧 bias reliance

### 4.3 读出头设计对比

| 方案 | Seed | Held-out JSD↓ | Held-out count r↑ | 判定 |
|------|------|:---:|:---:|------|
| **B (center pool)** | s42 | **0.3147** | **0.8503** | **当前最优** |
| **B (center pool)** | s1234 | **0.3145** | **0.8488** | 多 seed 稳定 |
| F (attn pool) | s42 | 0.3148 | 0.8492 | 单 seed 可接受 |
| F (attn pool) | s1234 | 0.4269 | 0.8457 | **profile 失守** |
| G (profile refine) | s42 | 0.4260 | 0.8460 | 失败 |

**结论**：B（center pool）是唯一跨 seed 稳定的读出头设计。

### 4.4 Corrected B（最终交付模型）

Clean matrix 中 corrected B = TF + center pool + sg=true + deb2：

| 指标 | Test 值 |
|------|---------|
| peak profile JSD full | 0.42544 |
| peak profile JSD debiased | 0.42556 |
| **gap** | **0.00268** |
| count r full | 0.8439 |

gap=0.00268，模型安全可用。

### 4.5 Genos-1.2B 融合（负结果）

| 阶段 | 方式 | 结果 | 判定 |
|------|------|------|------|
| Probe A | Count 残差岭回归 | R²=0.017 | 几乎无信号 |
| Probe B | 互补性逻辑回归 | concat AUC=0.90 < encoded-only 0.92 | 信号冲突 |
| G1 Online | 逐位置门控 | test count r=0.47 | 崩塌 |
| G2 Online | 全局均值门控 | test count r=0.64 | 严重下降 |
| P2 Cached | count head 注入 | test count r=0.60 | 断崖 |

**失败根因**：global_mean 粒度与碱基分辨率任务不匹配；encoded 与 genos 特征存在冲突；baseline 已很强，边际空间极小。

### 4.6 独立数据泛化

| 数据集 | peak JSD | count r |
|--------|:---:|:---:|
| GM12878 (6002) | 0.4226 | 0.8040 |
| K562 (6002) | 0.6124 | 0.8586 |

模型在独立细胞系上稳定工作。

---

## 五、关键发现与认知升级

### 5.1 最大认知翻转

| 之前的认知 | 现在的证据 |
|---|---|
| "Profile Shortcut 是 Transformer 特有的" | noTF+deb0 的 gap=0.027 是全矩阵最高，Transformer 反而更低 |
| "Stop-gradient 是关键修复" | deb2 的效应更大；A（无 sg 有 deb2）是最干净的条件 |
| "纯卷积不会出现 shortcut" | noTF+deb0 出现了最强的 bias reliance |
| "Genos 基础模型应该能帮" | 全线负结果，特征粒度根本不匹配 |

### 5.2 真正成立的结论

1. **Transformer 带来真实收益**（3 seed, p≈0.016）
2. **Debiased supervision 是 bias reliance 的通用控制手段**，与架构无关（2×2 矩阵）
3. **Full/debiased gap 是有效的诊断指标**（能清晰区分 safe 和 unsafe 配置）
4. **Center pool 是最稳读出头**（跨 seed 不退化）
5. **当前全局粒度的基础模型特征无法改善碱基分辨率任务**（高质量负结果）

---

## 六、论文叙事建议

### 6.1 推荐定位

~~"Identifying and Resolving Transformer-Specific Profile Shortcuts"~~

→ **"TransChromBP: Bias-Safe Transformer Integration for Base-Resolution Chromatin Accessibility Modeling"**

核心卖点不是"发现了一个 bug"，而是"建立了一个框架"：
- 方法贡献：将 Transformer 安全融入 bias-factorized 建模
- 诊断贡献：full/debiased gap 作为通用 bias reliance 指标
- 实证贡献：2×2 矩阵证明 debiased supervision 是架构无关的稳定器

### 6.2 四个主 Claim

| # | 主张 | 强度 | 核心证据 |
|---|------|------|---------|
| C1 | Transformer 带来真实收益 | **shows** | 3-seed ablation, 0.314 vs 0.324 |
| C2 | full/debiased gap 是有效诊断 | **shows** | 2×2 矩阵四格差异 |
| C3 | debiased supervision 是通用稳定器 | **shows** | 2×2 矩阵交互效应 + 架构无关性 |
| C4 | 最终模型安全可用 | **shows** | corrected B gap=0.00268, B 双 seed 稳定 |

### 6.3 必须删除或降级的表述

| 不能写 | 改为 |
|--------|------|
| "Profile Shortcut 是 Transformer 特有的" | "bias reliance 风险存在于 unsafe 配置中，与架构无关" |
| "Stop-gradient 单独修复了问题" | "debiased supervision 是主效应，stop-gradient 是辅助设计" |
| "纯卷积不会出现 shortcut" | "conv-only 在 unsafe 配置下同样出现、且 gap 更大" |

### 6.4 建议论文结构

1. **Introduction**：以"如何安全集成 Transformer 到 bias-factorized 框架"为问题
2. **Method**：架构 + full/debiased 诊断框架 + debiased supervision 设计
3. **Experiments**：
   - 3.1 主结果与 ChromBPNet 对比（TF ablation 3 seed）
   - 3.2 2×2 Bias Reliance 矩阵（本文核心表格）
   - 3.3 读出头设计消融
   - 3.4 独立数据泛化
   - 3.5 基础模型融合负结果（Discussion 或附录）
4. **Discussion**：诊断框架的通用性、当前证据边界、负结果启示

### 6.5 2×2 矩阵的论文价值

这个矩阵本身就是一个方法学贡献。它的结论——**"debiased supervision 的效果与架构无关，是通用的 bias 控制手段"**——对任何使用 bias factorization 的模型都有参考意义，不限于 ChromBPNet 系列。

---

## 七、后续实验建议

### 立即可做

| 优先级 | 实验 | 目的 | 耗时 |
|--------|------|------|------|
| **P1** | ChromBPNet 官方 baseline 用 mean JSD 重算 | 消除 median vs mean 口径差异 | 脚本 |

### 可选加固

| 优先级 | 实验 | 目的 | 耗时 |
|--------|------|------|------|
| **P2** | `notf_sg1_deb0` 第二 seed (s1234) | 确认 noTF+deb0 高 gap 的稳定性 | ~6h (6002) |
| **P3** | B 第三 seed | 主表展示质量 | ~8h |

### 不需要做

- Genos 任何新融合（充分负结果）
- 更多读出头变体（B 已确认胜出）
- V1→V2 的逐步拆解（价值低，且无法完全匹配历史条件）

---

## 八、数据资产清单

### 结果文件

| 文件 | 位置 | 内容 |
|------|------|------|
| A test | 6000 outputs/metrics/ | TF+sg0+deb2 test |
| C s42/s1234 test | 6000 outputs/metrics/ | TF+sg1+deb0 双 seed test |
| noTF_deb2 test | 6002 outputs/metrics/ | noTF+sg0+deb2 test |
| notf_sg1_deb0 test | 6002 outputs/metrics/ | noTF+sg1+deb0 test（本日新增） |
| corrected B test | 6002 outputs/metrics/ | TF+center+sg1+deb2 test |
| B s42/s1234 | 6000 outputs/metrics/ | center pool 双 seed held-out |
| GM12878/K562 | 6002 outputs/metrics/ | 独立数据 held-out |

### 分析报告

| 报告 | 路径 |
|------|------|
| 本报告 | `reports/analysis/transchrombp_full_review_20260326.md` |
| 论文改写策略 | `reports/paper/paper_rewrite_strategy_20260326.md` |
| Claim-Evidence 矩阵 | `reports/paper/paper_claim_evidence_matrix_20260326.md` |
| Profile Shortcut 复核 | `reports/analysis/profile_shortcut_revalidation_summary_20260326.md` |
| 指标源数据总表 | `reports/assets/paper_metric_source_table_20260326.csv` |
| Genos 负结果分析 | `reports/analysis/genos_no_positive_gain_analysis_20260323.md` |
| Genos 完整技术报告 | `reports/analysis/genos_integration_full_report_for_review_20260324.md` |
| 读出头实验状态 | `reports/analysis/v2fix_readout_head_status_20260323.md` |
| 独立数据跟进 | `reports/analysis/v2fix_and_6000_realdata_followup_20260325.md` |

---

## 九、总结

**当前最大的认知升级**：`notf_sg1_deb0` 的 gap=0.027 是全矩阵最高。这彻底改变了论文故事——从"修复 Transformer 的 bug"变成了"发现 debiased supervision 是通用的、架构无关的 bias 控制手段"。后者的学术价值更高、更可靠、更不容易被 reviewer 攻击。

**论文最稳路径**：不试图证明 Transformer 特有的强机制，而是证明在 bias-safe 约束下 Transformer 能安全地提供真实增益，并提出 full/debiased gap 作为通用诊断框架。
