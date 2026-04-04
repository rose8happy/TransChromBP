# Tutorial L2 Locked Follow-up 对比表（2026-03-29）

> 前置：[strict_compare_v2_updated_methodology_20260329.md](../docs/plan/strict_compare_v2_updated_methodology_20260329.md)
> Codex 复核：[strict_compare_v2_methodology_codex_review_20260329.md](strict_compare_v2_methodology_codex_review_20260329.md)

## 1. Follow-up 表

对 3 个固定 checkpoint 分别在 `valid` 和 `test` split 上的完整指标。test 口径仅用于描述固定 checkpoint 表现，不参与重新选 best。

| Arm | Checkpoint | split | peak median_jsd | peak mean_jsd | peak count_r |
|---|---|---|---:|---:|---:|
| official controlled | epoch 35 | valid | 0.35783 | 0.35234 | 0.65805 |
| official controlled | epoch 35 | **test** | **0.33983** | **0.33769** | **0.68164** |
| TransChromBP corrected-B | epoch 35 | valid | 0.33641 | 0.33061 | 0.69681 |
| TransChromBP corrected-B | epoch 35 | **test** | **0.31547** | **0.31349** | **0.84234** |
| TransChromBP corrected-B | epoch 44 | valid | 0.33538 | 0.32944 | 0.08428 |
| TransChromBP corrected-B | epoch 44 | **test** | **0.31437** | **0.31219** | **0.83778** |

## 2. 关键发现

### 2.1 epoch 44 在 test 上 count_r 正常，不再支持"训练真实退化"的解释

- epoch 44 valid count_r = **0.084**（严重崩塌）
- epoch 44 test count_r = **0.838**（完全健康，甚至远高于 official 的 0.682）

这意味着 epoch 44 **不应被称为"明显坏 checkpoint"**。之前观察到的 valid count_r 崩塌不是训练真正退化。

但需要注意：count_r 在 valid/test 上的巨大差异（0.084 vs 0.838）的根因尚不清楚。fold_0 的 valid 是 chr2、test 是 chr1，两者都是单条染色体，因此"valid 单染色体导致估计不稳定"这个解释不充分——test 同样是单染色体但 count 完全正常。更准确的问题表述是：**valid 和 test 对 count 指标的一致性有限**，而不是已确定根因。

### 2.2 TransChromBP 在 test 上显著优于 official（两个 checkpoint 均如此）

| 指标 | official epoch 35 (test) | ours epoch 35 (test) | ours epoch 44 (test) | ours 优势（epoch 35 vs official） |
|---|---:|---:|---:|---|
| peak median_jsd | 0.33983 | 0.31547 | 0.31437 | **-0.02436**（更低 = 更好） |
| peak mean_jsd | 0.33769 | 0.31349 | 0.31219 | **-0.02420** |
| peak count_r | 0.68164 | 0.84234 | 0.83778 | **+0.16070**（更高 = 更好） |

TransChromBP 在 profile（JSD 降低约 7%）和 count（pearson r 提升约 0.16）上均明显优于 official。

### 2.3 epoch 35 与 epoch 44 在 test 上几乎无差异

| 指标 | ours epoch 35 (test) | ours epoch 44 (test) | 差异 |
|---|---:|---:|---|
| peak median_jsd | 0.31547 | 0.31437 | -0.00110 |
| peak mean_jsd | 0.31349 | 0.31219 | -0.00130 |
| peak count_r | 0.84234 | 0.83778 | -0.00456 |

两个 checkpoint 在 test 上表现几乎相同。profile-only selector 选中的 epoch 44 在 test 上没有踩空。

## 3. 对决策的影响

### 3.1 可以确认的

1. **epoch 44 不是坏 checkpoint**：test 上 count_r=0.838，profile 和 count 均优于 official
2. **tutorial L2 单 seed matched-budget compare 可以更正面地写**：test 口径下 TransChromBP 在 profile 和 count 上均显著优于 official，且两个候选 checkpoint 表现一致
3. **selector 改造从 blocker 降级**：不再需要先修 selector 才能继续推进
4. **优先推进 L3 shared-region compare 或 GM12878 L2**：这是当前提升主结果说服力的最高价值项

### 3.2 不能确认的

1. **"profile-only selector 已经足够"**：这次 valid-selector 选中的 epoch 44 在 test 上碰巧没问题，但 valid/test 对 count 指标的一致性仍然有限。如果后续在其他 fold / dataset 上出现类似分歧，selector 健壮性问题仍需重新审视
2. **valid count_r 崩塌的根因**：chr2（valid）和 chr1（test）都是单条染色体，但 count 表现截然不同。根因可能是染色体级别的 count 分布差异、模型对特定区域的 count 预测敏感性等，目前没有足够证据确定
3. **selector 改造是否彻底不需要**：当前只是不再是 blocker。如后续实现，仍建议保持"默认关闭 + 显式启用"的 safeguard 形态，并在多 fold / 多 dataset 上观察后再决定是否升格

## 4. 建议的 L2 主表口径

| Arm | representative checkpoint | peak median_jsd (test) | peak mean_jsd (test) | peak count_r (test) |
|---|---|---:|---:|---:|
| ChromBPNet official controlled | epoch 35 (valid-selector best) | 0.33983 | 0.33769 | 0.68164 |
| TransChromBP corrected-B controlled | epoch 44 (valid-selector best) | 0.31437 | 0.31219 | 0.83778 |

口径限定：matched-budget system comparison, single seed (42), tutorial dataset, native preprocessing (region 未统一)。

## 5. 数据来源

- 官方 valid 指标：`/data1/.../external_best_valid/best_epoch.json`
- 官方 test 指标：`/data1/.../followup/tutorial_L2_locked_20260329/official_epoch35_test/best_epoch.json`
- TransChromBP valid 指标：`/data1/.../external_best_valid/epoch_metrics.csv`
- TransChromBP epoch 35 test：`/data1/.../followup/tutorial_L2_locked_20260329/ours_epoch35_test.json`
- TransChromBP epoch 44 test：`/data1/.../followup/tutorial_L2_locked_20260329/ours_epoch44_test.json`
