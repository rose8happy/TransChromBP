# Codex 审查回应：4 个阻塞点分析与修复方案

> 创建日期：2026-03-20
> 背景：Codex 对 v2fix 实验代码进行了全面审查，识别出 4 个执行链路阻塞点
> 状态：**全部修复完成**

---

## 问题 1：evaluate_checkpoint.py 使用了错误的数据集 API

### 验证

**确认存在。** 对比两处代码：

| 位置 | API 调用方式 |
|---|---|
| `evaluate_checkpoint.py:69-76` | `ChromBPNetBigWigDataset(config=ds_config, split="val", max_seq_len=..., output_len=...)` |
| `train_ddp.py:954-994` | `ChromBPNetBigWigDataset(genome_fasta=..., bigwig_path=..., peaks_bed=..., folds_json=..., split="valid", input_len=..., supervised_bp=...)` |

两个关键差异：
1. **构造函数参数**：evaluate 用的是旧 API (`config=`, `max_seq_len=`)，真实 API 需要 `genome_fasta=`, `bigwig_path=`, `peaks_bed=` 等
2. **split 名称**：evaluate 用 `"val"`，训练器和 folds.json 用 `"valid"`

### 影响范围

- **实验 E (soup)** 的评估阶段会**直接失败**
- 实验 A/B/F/G 训练本身不受影响（训练器 `train_ddp.py` 的 API 是正确的），但**训练后的独立评估**也会失败

### 修复

重写为 `build_eval_loader()`，直接使用 `ChromBPNetBigWigDataset` 的正确 API，支持 `valid`/`test` split 选择。
同时 `main()` 中补齐了 `data_source_cfg` 加载链路，复用 `train_ddp` 的 config 解析函数。

---

## 问题 2：soup_top3 依赖的 JSD 指标可能不存在

### 验证

**部分确认。**

- `train_ddp.py:432-433` 中 **已经实现了** `profile_target_jsd_full_mean` 和 `profile_target_jsd_debiased_mean` 的计算（在 `finalize_selection_metrics()` 中）
- 但 TRACKING.md:17 说"validation 指标对齐改造**尚未写回远端**"

关键问题：**已有 `ablation_tf_20260318` 的 `epoch_metrics.jsonl` 是用哪个版本的 `train_ddp.py` 跑的？**
如果旧版没有 `finalize_selection_metrics` → `epoch_metrics.jsonl` 里没有 JSD 字段 → `soup_top3` 会失败。

### 修复

Orchestrator 中 soup_top3 增加三级 fallback：
1. 优先 `peak.profile_target_jsd_full_mean`
2. 回退到 `peak.profile_target_jsd_debiased_mean`
3. 再回退到 `peak.loss_profile`

---

## 问题 3：best_metric 应该用 JSD 而不是 loss_total

### 验证

**确认存在，且影响重大。**

`peak.loss_total` 是混合指标（profile NLL + count MSE + debiased profile NLL），问题在于：
- NLL 和 JSD **不完全单调**——NLL 最低的 checkpoint 不一定 JSD 最低
- 对 B/F（count 改进实验），`peak.loss_total` 里 count 权重只有 0.1，改善容易被掩盖
- 对 G（profile 改进实验），profile NLL 和 JSD 的发散在训练末期最明显

### 修复

新建 `train_ablation_v2_main_profile_select.yaml` 作为推荐主训练 config：
```yaml
best_metric: peak.profile_target_jsd_full_mean
best_metric_mode: min
early_stop_patience: 10
```

同时保留原 `train_ablation_v2_main.yaml`（`peak.loss_total`）用于回溯对比。

> **前提**：部署到远端的 `train_ddp.py` 必须包含 `finalize_selection_metrics()`，
> 否则 JSD 字段不会出现在 validation metrics 中，best checkpoint 永远不更新。

---

## 问题 4：Count 信号分支缺乏直接监督

### 验证

**确认存在，需要认真对待。**

当前 `debiased_count_weight: 0.0`，即 count_signal 没有被直接监督。
梯度需要通过 logsumexp 融合反传。B/F 改的恰好是 count_signal 路径上的 pooling——
如果 count_signal 本身没被显式约束，改 pooling 的效果会被稀释。

### 修复

**采用方案 A + C**（保持单变量原则）：
1. **这轮先不改 loss**：`debiased_count_weight: 0.0`，与 baseline 对齐
2. **判读时用 `count_pearson_debiased` 做主指标**：直接反映 signal branch 的 count 质量
3. 如果 B/F 不提升，在第二轮加 `debiased_count_weight: 0.1` 做对照

实验文档 v2.1 新增 §6.4 详细说明了这一判读逻辑。

---

## 修复优先级与依赖关系

```
修复任务                              阻塞什么？            状态
──────────────────────────────────────────────────────────────
[1] 修 evaluate_checkpoint.py API     实验 E (soup 评估)     ✅ 已修复
[3] 改 best_metric 为 JSD             实验 A/B/F/G 选模      ✅ 已修复
[2] soup_top3 增加 fallback           实验 E (top3 子项)     ✅ 已修复
[4] 判读指标改为 count_pearson_debiased 结果分析              ✅ 已更新文档
```
