# ChatGPT 咨询包 02：实验历史、关键证据与当前已知事实

## 一分钟版本
- 我们最初以为发现了一个“Transformer 特有的 Profile Shortcut”，但后来的 clean-matrix 复核基本推翻了这个强叙事。
- 现在最稳的结论是：`debiased profile supervision` 是通用稳定器；`full/debiased gap` 是有价值的诊断口径；Transformer 的收益来自真实表征提升，而不是更强地偷 bias。
- 当前最优默认 readout 是 `B = center pool`，它在 held-out 上跨 seed 稳定；`F = attention pool` 在 seed=1234 明显失守。
- Genos-1.2B 融合目前是高质量负结果，几乎可以视为当前阶段应降级的支线。

## 1. 项目主线是什么
目标：在 ChromBPNet 的 bias factorization 框架上接入 Transformer，做 base-resolution ATAC-seq 的 profile/count 建模，同时避免模型在 signal branch 中“中继” bias 信号。

最初问题意识：
- Transformer 容量更强，是否更容易绕过 bias factorization？
- 如果是，模型表面指标提高可能只是更依赖 bias branch，而不是真学到调控信号。

## 2. 实验历程（按时间）

### 阶段 1：V1 -> V2 架构升级，产生早期误判
早期版本同时改了三件事：
- `stop_gradient=true`
- `debiased_profile_weight=2.0`
- `bias profile pool_factor=32`

当时观察到 debiased JSD 从大约 `0.56` 降到 `0.316` 左右，于是曾解读成：
- “发现并修复了 Transformer 特有的 Profile Shortcut”

后来的问题：这三个变量是同时改的，所以完全不能归因。

### 阶段 2：matched ablation + readout 设计
#### 2.1 Transformer 是否真的带来收益
三 seed matched ablation：

| 配置 | Profile JSD (mean±std) | Count r (mean±std) | 结论 |
|---|---:|---:|---|
| V2-full (TF) | **0.31419 ± 0.00012** | **0.83676 ± 0.00531** | 当前最佳 backbone |
| V2-noTF | 0.32393 ± 0.00021 | 0.82503 ± 0.00489 | 稳定弱于 TF |

最稳结论：Transformer 有真实收益，这条是成立的。

#### 2.2 Readout 设计：B vs F vs G
| 方案 | Seed | Held-out JSD | Held-out count r | 结论 |
|---|---|---:|---:|---|
| B (center pool) | s42 | **0.3147** | **0.8503** | 胜出 |
| B (center pool) | s1234 | **0.3145** | **0.8488** | 多 seed 稳定 |
| F (attention pool) | s42 | 0.3148 | 0.8492 | 单 seed 看着可接受 |
| F (attention pool) | s1234 | 0.4269 | 0.8457 | profile 明显失守 |
| G (profile refine) | s42 | 0.4260 | 0.8460 | 失败 |

最稳结论：`B=center pool` 是当前唯一跨 seed 稳定的默认 readout。

### 阶段 3：Genos-1.2B 融合与独立数据
#### 3.1 Genos 融合：基本是负结果
| 阶段 | 方法 | 结果 | 判断 |
|---|---|---|---|
| Probe A | count 残差岭回归 | `R²=0.017` | 几乎没新信息 |
| Probe B | encoded + genos 逻辑回归 | concat AUC 反而下降 | 特征冲突/冗余 |
| G1/G2/P2 | 在线或缓存融合 | count 明显崩塌 | 任务粒度不匹配 |

当前判断：Genos summary 这种全局粒度特征不适合当前 base-resolution profile/count 任务。

#### 3.2 独立数据泛化
| 数据集 | Held-out peak JSD | Held-out count r |
|---|---:|---:|
| GM12878 | 0.42265 | 0.80396 |
| independent K562 | 0.61235 | 0.85857 |

最稳结论：模型不是只会 tutorial benchmark，跨数据集仍能稳定工作。

### 阶段 4：Clean Matrix 受控复核（最关键）
这是当前最重要的一组证据。我们把两个因素拆开：
- backbone：`TF` vs `noTF`
- bias 安全配置：`safe = sg=false + deb2` vs `unsafe = sg=true + deb0`

#### 4.1 判读门槛
- `gap <= 0.005`：基本无实质 shortcut / bias reliance
- `0.005 < gap <= 0.02`：轻度 bias reliance
- `gap > 0.02`：明显 shortcut

这里的 `gap` 指：`profile_full_debiased_jsd`

#### 4.2 Clean matrix test 结果
| Run | 条件 | peak JSD full | peak JSD deb | gap | count r full/deb |
|---|---|---:|---:|---:|---:|
| A | `TF + sg=false + deb2` | 0.31410 | 0.31416 | **0.00147** | 0.8506 / 0.8488 |
| C (s42) | `TF + sg=true + deb0` | 0.31720 | 0.31779 | **0.00972** | 0.8395 / 0.8370 |
| C (s1234) | `TF + sg=true + deb0` | 0.31739 | 0.31795 | **0.00997** | 0.8340 / 0.8326 |
| noTF_deb2 | `noTF + sg=false + deb2` | 0.43531 | 0.43557 | **0.00487** | 0.7959 / 0.7966 |
| notf_sg1_deb0 (s42) | `noTF + sg=true + deb0` | 0.43659 | 0.43835 | **0.02682** | 0.7919 / 0.7956 |
| notf_sg1_deb0 (s1234) | `noTF + sg=true + deb0` | 0.43731 | 0.43841 | **0.01679** | 0.8287 / 0.8267 |
| corrected B | `TF + center pool + sg=true + deb2` | 0.42544 | 0.42556 | **0.00268** | 0.8439 / 0.8420 |

#### 4.3 这张表真正说明了什么
1. 当前最高风险格不是 TF，而是 conv-only 的 `notf_sg1_deb0 = 0.02682 / 0.01679`；两 seed 都高于 TF 的 `C = 0.00972 / 0.00997`，但第二个 seed 没有越过 `0.02`。
2. `A = 无 sg 但有 deb2` 很干净，说明 stop-gradient 不是唯一关键开关。
3. `C = 有 sg 但无 deb2` 在两个 seed 上稳定出现轻度 gap，说明 `deb2` 更像主效应。
4. `corrected B` gap 只有 `0.00268`，说明最终默认模型仍然是“安全”的。

#### 4.4 因此必须放弃的说法
- “Transformer 特有 shortcut”
- “纯卷积不会出现 shortcut”
- “stop-gradient 单独修复了问题”
- “我们已经证明了强而稳定的 shortcut 机制”

## 3. 现在哪些结论最稳
### 可以强写（shows）
1. Transformer 在当前 factorized ATAC 建模里带来真实收益。
2. `full/debiased gap` 是有价值的 bias reliance 诊断口径。
3. `center pool` 是当前最稳的默认 readout。
4. 最终默认模型 `corrected B` 没有明显 bias reliance。

### 只能弱写（suggests）
1. 在当前受控矩阵下，`debiased profile supervision` 比 stop-gradient 更关键。
2. Genos 失败可能是任务粒度不匹配，而不只是工程实现差。

### 不能再写
1. Transformer 特有 shortcut
2. conv-only 不会 shortcut
3. stop-gradient 单独修复
4. 强机制已被证明

## 4. clean-matrix 最后一个关键加固已完成
- 任务：`notf_sg1_deb0_s1234_6002`
- 最新结果：peak `profile_full_debiased_jsd=0.01679`，`count_r=0.82869 / 0.82669`
- 新判断：conv-only unsafe 仍是 clean matrix 的最高风险格，但幅度存在 seed 级波动；因此不能写成“稳定复现强 shortcut”，但更足以反证“Transformer 特有”

## 5. 当前科学结论的真正边界
### 成立的
- Transformer 有用
- bias-safe 训练/诊断框架有用
- `deb2` 很关键
- `B=center pool` 稳定
- final model 可交付

### 不成立或证据不足的
- 强 shortcut 已被明确证明
- 这是 Transformer 特有问题
- 机制已经被完全解释
- Genos 线还有高优先级投入价值

## 6. 这份文件背后的关键源报告
如果你想核查更细节的证据链，核心源头是：
- `reports/analysis/transchrombp_full_review_20260326.md`
- `reports/paper/paper_rewrite_strategy_20260326.md`
- `reports/analysis/profile_shortcut_revalidation_summary_20260326.md`
- `reports/paper/paper_claim_evidence_matrix_20260326.md`
- `reports/analysis/v2fix_and_6000_realdata_followup_20260325.md`
