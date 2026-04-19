# TransChromBP vs ChromBPNet 严格对比 V2 方法论（2026-03-29，修订版）

> 本文档基于 tutorial L2 第一轮对比的经验教训，更新对比方法论，供下一轮执行使用。
> 前置文档：[strict_chrombpnet_official_comparison_execution_20260327.md](strict_chrombpnet_official_comparison_execution_20260327.md)、[strict_compare_code_audit_20260327.md](../../reports/analysis/strict_compare_code_audit_20260327.md)、[strict_compare_followup_assessment_20260327.md](../../reports/analysis/strict_compare_followup_assessment_20260327.md)
> Codex 复核：[strict_compare_v2_methodology_codex_review_20260329.md](../../reports/analysis/strict_compare_v2_methodology_codex_review_20260329.md)

### 修订说明（2026-03-29）

根据 Codex 复核意见（F1/F2/F3），本文档做如下修订：
- **F1 采纳**：`count_r floor=0.5` 不再作为默认 protocol，降级为"默认关闭的可选 robustness 开关"
- **F2 采纳**：执行顺序调整为"先补 locked follow-up 表 → 再决定是否改 selector → 再决定是否进 L3"
- **F3 采纳**：如实施 selector 改动，需补全 NaN 处理、fallback 行为、元数据字段等协议定义

---

## 一、第一轮经验教训总结

### 1.1 已踩坑与已修复

| # | 问题 | 严重度 | 根因 | 修复状态 | 教训 |
|---|---|---|---|---|---|
| E1 | 官方 selector 在 test 集上选 best | P0 | `predict.py` 硬编码 `mode="test"` | 已修复：新增 `--split` 参数，默认 `valid` | 任何 selector 必须先验证其评估 split，不能假定上游默认正确 |
| E2 | Region 集合不一致（peaks 差 2,625，nonpeaks 差 241,254） | P1 | 两侧各自独立 preprocessing | 已识别，L3 方案已设计 | L2 结论只能写 system-level，不能写 architecture-level；要做更强归因必须走 L3 |
| E3 | Late-epoch count_r 崩塌（epoch 29+ 开始失稳） | P1 | 待分析（疑似学习率/loss 权重/训练动态） | 未修复，需独立分析 | 不能只看 profile 选 best；需要对固定 checkpoint 做完整描述后再决策 |
| E4 | Profile-only selector 选中 count 已塌 checkpoint | P1 | `median_jsd` 单指标选 best | 已识别 | 先描述问题（locked follow-up），再讨论 selector 改进 |
| E5 | `--max-parallel` 重复传参 | P2 | wrapper 脚本两处同时写入 | 已修复 | 脚本参数只在一处赋值 |
| E6 | `metrics_only` 仍写 PNG | P3 | `counts_metrics()` 无条件 `plt.savefig` | 已修复 | 逐 epoch 评估时额外 I/O 开销 |
| E7 | MirroredStrategy global batch 语义待验证 | P1 | TF/Keras + Sequence 的分发行为 | 已闭环：synthetic smoke 证实 global=32 | 多卡训练的 batch 语义必须有 runtime 证据，不能仅靠代码推断 |

### 1.2 第一轮 tutorial L2 关键数据

| Arm | best epoch | peak median_jsd (valid) | peak mean_jsd (valid) | peak count_r (valid) | 说明 |
|---|---:|---:|---:|---:|---|
| official controlled | 35 | 0.35783 | 0.35234 | 0.65805 | profile-only selector |
| TransChromBP corrected-B (epoch 44) | 44 | 0.33538 | 0.32944 | 0.08428 | profile-only selector 选中，count 已塌 |
| TransChromBP corrected-B (epoch 35) | 35 | 0.33641 | 0.33061 | 0.69681 | 同 run 中 count 正常的优质 checkpoint |

**关键发现**：TransChromBP 在 profile 侧已有明确正向信号（优于 official），但 profile-only selector 选到了 count 已崩的晚期 checkpoint。同 run 的 epoch 35（count 正常）profile 仅比 epoch 44 差 0.00103，同时 count_r=0.69681 高于 official 的 0.65805。

### 1.3 尚未解决的关键问题

1. **count_r 崩塌根因不明**：TransChromBP 的 count_r 在 epoch 29 后开始失稳，具体原因待独立分析
2. **epoch 35 vs 44 的 test 表现未知**：当前只有 valid 口径，不确定 count 崩塌是 valid 偶然现象还是 held-out 上也成立
3. **Region 差异对结论的影响未定量**：nonpeaks 差了近一半
4. **单 seed 结论不够稳健**：当前只有 seed=42

---

## 二、V2 执行策略（修订后）

### 2.1 核心原则：先描述，再决策

Codex 复核正确指出：当前最值钱的下一步不是改 selector，而是先把 `epoch 35/44` 这组固定 checkpoint 的 `valid/test` 指标描述清楚。只有在这一步完成后，才有足够信息决定是否需要改 selector、是否进 L3。

### 2.2 修订后的执行顺序

```
Step 1（即时，不重训）: Locked Follow-up 表
  → 对 3 个固定 checkpoint 补 test 口径评估
  → 产出 valid + test 双口径描述表
  → test 只用于"固定 checkpoint 描述"，不参与重新选 best

Step 2（Step 1 完成后决策）: 判断 count 崩塌是否影响 test
  → 如果 epoch 44 的 test count_r 也崩 → count 崩塌是真实问题，不是 valid 偶然
  → 如果 epoch 44 的 test count_r 正常 → count 崩塌可能是 valid set 特有现象

Step 3（视 Step 2 结论）: 决定下一步
  → 路径 A：count 崩塌在 test 上确认 → 讨论 selector 改进（可选 robustness 开关）+ floor 灵敏度分析
  → 路径 B：count 崩塌仅 valid 偶然 → 当前 profile-only selector 可能仍足够，优先进 L3
  → 两条路径都可以进入 GM12878 L2，但叙事不同

Step 4（Step 3 之后）: L3 shared-region compare
  → 仍然是提升主结果说服力的硬价值项
  → 复合 selector 不能替代 L3 的作用

Step 5: 多 seed 扩展 / Dataset-X
  → 在 tutorial + GM12878 方向一致后
```

### 2.3 指标汇报规范

每条 arm 的 locked follow-up 必须汇报：

| 指标 | 必报 | 说明 |
|---|---|---|
| peak profile median_jsd | 是 | valid 和 test 双口径 |
| peak profile mean_jsd | 是 | valid 和 test 双口径 |
| peak count_r | 是 | valid 和 test 双口径 |
| 评估 split | 是 | 明确标注 valid / test |

### 2.4 论文表述规范（不变）

| 如果... | 允许写 | 不允许写 |
|---|---|---|
| L2 TransChromBP 优于 official | "matched-budget system comparison" | "strictly controlled variable comparison" |
| L3 进一步确认 | "shared-region system comparison" | "architecture-only attribution" |
| count 存在 late-epoch instability | "profile advantage with noted count instability" | "overall superior" |

---

## 三、Phase 1：Tutorial L2 Locked Follow-up（即时执行）

### 3.1 目标

对以下 3 个**固定** checkpoint 分别产出 `valid` + `test` 双口径指标：

| # | Arm | Checkpoint | valid 指标 | test 指标 |
|---|---|---|---|---|
| 1 | official controlled | epoch 35 | 已有（selector 产物） | **需补** |
| 2 | TransChromBP corrected-B | epoch 35 | 已有（epoch_metrics.csv） | **需补** |
| 3 | TransChromBP corrected-B | epoch 44 | 已有（epoch_metrics.csv） | **需补** |

**关键约束**：test 口径只用于描述这 3 个固定 checkpoint 的表现，不参与重新选 best。

### 3.2 官方侧：对 epoch 35 补 test 口径

官方 selector 已有 valid 口径产物。只需对 epoch 35 单独跑一次 `--split test`：

```bash
# 在 6000 上执行
export CHROMBPNET_ENV=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet
export PATH="$CHROMBPNET_ENV/bin:$PATH"
export LD_LIBRARY_PATH="$CHROMBPNET_ENV/lib:$LD_LIBRARY_PATH"
export PYTHONPATH=/data1/zhoujiazhen/bylw_atac/chromBPNet:$PYTHONPATH

OFFICIAL_RUN=/data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/runs/tutorial_official_controlled_s42/fold_0/seed_42/chrombpnet
FOLLOWUP_DIR=/data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/followup/tutorial_L2_locked_20260329

mkdir -p $FOLLOWUP_DIR

python scripts/paper_aligned_repro/select_best_epoch.py \
  --model-glob "$OFFICIAL_RUN/models/chrombpnet.epoch_035.h5" \
  --genome /data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.fa \
  --bigwig "$OFFICIAL_RUN/auxiliary/data_unstranded.bw" \
  --peaks "$OFFICIAL_RUN/auxiliary/filtered.peaks.bed" \
  --nonpeaks None \
  --fold-json /data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/folds/fold_0.json \
  --output-dir "$FOLLOWUP_DIR/official_epoch35_test" \
  --split test \
  --metric profile_metrics.peaks.median_jsd \
  --mode min
```

注意：这里用 `select_best_epoch.py` 对单个 checkpoint glob 求值，本质上只是借用它的评估链路，不是在"选 best"。

### 3.3 TransChromBP 侧：对 epoch 35 和 epoch 44 补 test 口径

```bash
# 在 6000 上执行
export CHROMBPNET_ENV=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet
export PATH="$CHROMBPNET_ENV/bin:$PATH"
export LD_LIBRARY_PATH="$CHROMBPNET_ENV/lib:$LD_LIBRARY_PATH"
export PYTHONPATH=/data1/zhoujiazhen/bylw_atac/TransChromBP/src:$PYTHONPATH

OURS_RUN=tutorial_corrected_b_strict_compare_s42
CKPT_DIR=/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/checkpoints/$OURS_RUN
FOLLOWUP_DIR=/data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/followup/tutorial_L2_locked_20260329

# Epoch 35 test
python -m transchrombp.evaluation.evaluate_checkpoint \
  --checkpoint "$CKPT_DIR/epoch_035.pt" \
  --split test \
  --device cuda:0 \
  --output "$FOLLOWUP_DIR/ours_epoch35_test.json"

# Epoch 44 test
python -m transchrombp.evaluation.evaluate_checkpoint \
  --checkpoint "$CKPT_DIR/epoch_044.pt" \
  --split test \
  --device cuda:1 \
  --output "$FOLLOWUP_DIR/ours_epoch44_test.json"
```

两个 epoch 可以分配到不同 GPU 并行。

### 3.4 产出：Locked Follow-up 对比表

评估完成后，汇总为以下表格（写入 `reports/` 下独立报告）：

```
| Arm | Checkpoint | split | median_jsd | mean_jsd | count_r |
|---|---|---|---|---|---|
| official | epoch 35 | valid | 0.35783 | 0.35234 | 0.65805 |
| official | epoch 35 | test | ? | ? | ? |
| ours | epoch 35 | valid | 0.33641 | 0.33061 | 0.69681 |
| ours | epoch 35 | test | ? | ? | ? |
| ours | epoch 44 | valid | 0.33538 | 0.32944 | 0.08428 |
| ours | epoch 44 | test | ? | ? | ? |
```

### 3.5 决策规则

1. **如果 ours epoch 44 test count_r 也很低（< 0.3）**：
   - count 崩塌是真实训练问题，不是 valid set 偶然
   - epoch 35 成为更可信的 checkpoint
   - 值得讨论 selector 改进（作为可选 robustness 分析）
   - 论文中 tutorial L2 应以 epoch 35 为主，epoch 44 作为"late-epoch instability"的证据

2. **如果 ours epoch 44 test count_r 正常（> 0.5）**：
   - count 崩塌可能是 valid set 特有的估计噪声
   - 当前 profile-only selector 可能仍足够
   - 优先推进 L3 shared-region 而不是改 selector

3. **两种情况下 ours epoch 35 的 test profile 都有参考价值**：
   - 如果 test profile 仍优于 official → 正向信号在 held-out 上确认
   - 如果 test profile 不如 official → L2 结论需要降级

---

## 四、复合 Selector（可选 Robustness 开关，非默认 protocol）

> 本节基于 Codex 复核意见 F1/F3 修订：selector 改动定位为"默认关闭的可选 robustness 分析工具"，不直接升格为新的默认 protocol。

### 4.1 定位

- 不是默认选择规则，不进主表 protocol
- 是灵敏度分析 / robustness check 工具
- 只有在 locked follow-up 确认 count 崩塌是真实问题后，才讨论是否在特定 run 上启用

### 4.2 如实施，需满足的协议定义（F3 要求）

| 行为 | 定义 |
|---|---|
| 默认值 | `--count-r-floor` 不传 = 关闭，行为与当前完全一致 |
| `count_r` 为 NaN | 视为不通过 floor |
| `metric` 为 NaN 但 count 正常 | 不参与候选（与现有逻辑一致） |
| 多 epoch 满足 floor 且 metric 并列 | 取 epoch 编号更小的（deterministic tie-break） |
| 全部 epoch 低于 floor | 回退到 profile-only selector，`best_epoch.json` 标注 `"fallback_used": true`，论文中不允许写入主表，只作为 warning case 单列 |
| GM12878 / 多 seed 是否沿用同一阈值 | 不默认沿用；如需启用，先在 tutorial 做 floor 灵敏度分析（0.3 / 0.4 / 0.5 / 0.6），确认结论对阈值不敏感后再推广 |

### 4.3 `best_epoch.json` 输出格式扩展

启用 floor 时，`best_epoch.json` 需要包含以下额外字段：

```json
{
  "metric": "...",
  "mode": "...",
  "split": "...",
  "n_checkpoints": 50,
  "count_r_floor": 0.5,
  "count_r_metric": "counts_metrics.peaks.pearsonr",
  "n_excluded_by_floor": 18,
  "excluded_epochs": [29, 30, 32, ...],
  "fallback_used": false,
  "best": { ... }
}
```

未启用 floor 时，这些字段全部缺省（不出现在 JSON 中），产物格式与当前完全兼容。

### 4.4 代码改动清单

| 文件 | 改动 | 优先级 | 说明 |
|---|---|---|---|
| `scripts/paper_aligned_repro/select_best_epoch.py` | 增加 `--count-r-floor`（默认不传 = 关闭）；启用时过滤 + 输出扩展字段 | P2（在 locked follow-up 之后） | 不改变现有默认行为 |
| `vendor/transchrombp/scripts/select_best_epoch.py` | 同步增加，逻辑对齐 | P2 | 同上 |

---

## 五、后续阶段（locked follow-up 之后）

### 5.1 L3 shared-region compare

仍然是提升主结果说服力的**硬价值项**。复合 selector 不能替代 L3 的作用。

执行方案与原执行稿一致（见 [execution_20260327.md 第 5 节](strict_chrombpnet_official_comparison_execution_20260327.md#5-phase-3tutorial-l3-shared-region-双臂待-l2-结果确认后启动)），需要的代码改动：

| 文件 | 改动 | 优先级 |
|---|---|---|
| `scripts/paper_aligned_repro/run_tutorial_strict_compare_official.sh` | 增加 `--shared-peaks` / `--shared-nonpeaks` 参数 | P1 |
| 新建 TransChromBP L3 data config | 指向共享 region | P1 |

### 5.2 GM12878 L2

tutorial locked follow-up 闭环后平移。命令模板与原执行稿一致。

### 5.3 多 seed 扩展

单 seed 方向一致后再启动。

---

## 六、与现有文档的关系

| 文档 | 关系 |
|---|---|
| [strict_chrombpnet_official_comparison_execution_20260327.md](strict_chrombpnet_official_comparison_execution_20260327.md) | 原执行稿，L1/L2/L3 框架、protocol、命令模板仍有效 |
| [strict_compare_v2_methodology_codex_review_20260329.md](../../reports/analysis/strict_compare_v2_methodology_codex_review_20260329.md) | Codex 复核，本文档的修订依据 |
| [paper_claim_evidence_matrix_20260326.md](../../reports/paper/paper_claim_evidence_matrix_20260326.md) | Claim matrix 第 5.1 条与本文档 Phase 1 对齐 |
