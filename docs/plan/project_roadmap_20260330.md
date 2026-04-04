# TransChromBP 项目长远规划（2026-03-30 起）

## Context

TransChromBP 项目已进入论文收口阶段。核心证据基本齐备（clean matrix、TF ablation、readout design、cross-dataset），论文叙事已从"Transformer-specific shortcut"成功转向"bias-safe Transformer framework + dual-metric diagnostics"。当前瓶颈不是缺实验，而是：
1. L3 shared-region 技术链路已收口，但数字尚未完全并入论文主表、claim matrix 与 supplementary
2. 论文撰写仅 ~25%（结构在、内容缺）
3. 若需更高层级投稿，可能需要 GM12878 第二数据集加固

### 2026-03-30 晚间更新

- tutorial `L3 shared-region` 已完成 selector + 双侧 held-out test 收口，且结论方向与 `L2` 一致：`TransChromBP` 在 shared-region 的 held-out test 上同时优于 official `ChromBPNet`
- Official `classification_metrics` 缺失已定位为 6000 远端 `predict.py` / `metrics.py` 版本落后；同步后对 official `epoch 29` 单 checkpoint 重评估，新 JSON/CSV 已带出 `AUROC/AUPRC/F1`
- 因此 Phase 0 与 Phase 1 的 strict-compare 技术闭环已基本完成；下一步重心应从“补链路”切换到“把 L3 数字并入论文主表与 supplementary”

本计划定义从当前到论文提交的完整路径，含"快速通道"和"完整通道"两条分支。

---

## Phase 0: 收尾当前运行（3/30 - 3/31，已完成）

### 任务

| # | 任务 | 依赖 | 预计耗时 | 执行位置 |
|---|------|------|---------|---------|
| 0.1 | **TransChromBP L3 selector**: 对已完成的 50 epoch 跑 `select_best_epoch.py --split valid`，产出 `best_epoch.json` + `epoch_metrics.csv` | 无（L3 已完成） | 已完成 | 6000 |
| 0.2 | **Official L3 训练 + selector 监控**：确认自动队列正常推进到 `external_best_valid` | 无（被动） | 已完成 | 6000 |
| 0.3 | **TransChromBP L3 held-out test**：用 selector 选中的 checkpoint 跑 test split 评估 | 0.1 | 已完成 | 6000 |
| 0.4 | 开始草拟**架构图**（`fig:arch`）—— 论文最突出的缺失图 | 无 | 本地并行 | 本地 |

### 决策门 DG-0（Day 1 结束）
- Official L3 是否正常跑起来？若队列脚本卡住 → 手动启动
- TransChromBP L3 selector 是否有异常（count 崩塌等）？若有 → 立即补跑 backup epoch test（参照 L2 locked follow-up 模式）

---

## Phase 1: L3 锁定（4/1 - 4/4，技术链路已完成）

### 任务

| # | 任务 | 依赖 | 说明 |
|---|------|------|------|
| 1.1 | **Official L3 完成 + selector** | 0.2 训练完成 | 已完成；best checkpoint 为 `epoch 29` |
| 1.2 | **Official L3 test 评估** | 1.1 | 已完成；held-out test peak `mean_jsd=0.33853`、`median_jsd=0.33989`、`count_r=0.69958` |
| 1.3 | **组装 L3 对比表** | 0.3 + 1.2 | 已完成；shared-region 对比报告见 `reports/tutorial_L3_shared_region_closure_20260330.md` |
| 1.4 | **ATAC 分类指标验证** | 1.2 | 已完成；official 通过同步远端 `predict.py` / `metrics.py` 后补齐 `AUROC/AUPRC/F1` |
| 1.5 | **更新 TRACKING.md + 报告** | 1.3 | 已完成；`TRACKING.md`、source table 与路线图已同步 |

### 决策门 DG-1（Day 5）—— **快速/完整通道分叉点**

| 情景 | L3 结果 | 路径 |
|------|---------|------|
| **L3 确认 L2** | TransChromBP 在 shared-region 下 profile + count 双赢 | → **快速通道** Phase 2 |
| **L3 部分确认** | Profile 赢但 count 接近，或反之 | → 快速通道 + 谨慎措辞 |
| **L3 矛盾** | Official 追平或反超 | → 根因诊断 sprint（额外 2-3 天），然后决策 |

---

## Phase 2: 论文冲刺（4/5 - 4/18，~14 天）

### 2A: 内容填充（Days 6-10）

论文现有 9 个 `\todo{}` 标记，按优先级：

| 优先级 | TODO 位置 | 内容 | 信息来源 |
|--------|-----------|------|---------|
| **P0** | Table 1 | L3 shared-region 数字 | Phase 1 产出 |
| **P0** | Line 641 | 架构图 fig:arch | Phase 0.4 开始的草图 |
| **P1** | Line 656 | 架构参数 (d_model, dilated layers) | `configs/model/transchrombp_teacher_v2.yaml` |
| **P1** | Line 666 | Transformer 规格 (N, H, d_ff, dropout) | 同上 |
| **P1** | Line 743 | 训练超参 (LR, wd, warmup, epochs, batch) | `configs/train/` 系列 yaml |
| **P1** | Line 509 | 6002 单卡 runs 说明 | 重写或移至 supplementary |
| **P2** | Line 33 | 作者单位 | 需用户提供 |
| **P2** | Line 797 | 致谢 | 需用户提供 |
| **P2** | Line 802 | 仓库 URL + 数据编号 | 需用户确认 |

其他内容任务：
- 更新 clean matrix 表格（确保 conv-only 双 seed 数字已反映）
- 补充 ATAC 分类指标表（若 Phase 1.4 验证通过）
- 训练曲线 / gap 轨迹图（可视化 clean matrix 时间维度）

### 2B: 撰写与打磨（Days 11-16）

- 修订 Discussion：围绕 L3 证据收紧系统对比论述
- 修订 Limitations：若 L3 确认，削弱"single benchmark"限制
- 中文稿同步：`transchrombp_paper_cn_v1.tex` 与英文稿保持一致
- **参考文献扩充**：当前仅 6 条，目标 30-50 条（RoPE, BPNet, Enformer, ENCODE, GC-matching, 多项式 NLL, AdamW 等）

### 2C: 内部审核（Days 17-19）

- 逐条核对 `paper_claim_evidence_matrix_20260326.md`，确保文中措辞不超过允许强度
- 数字一致性检查：每个数字追溯到源 JSON/CSV
- 上传 ChatGPT 咨询包获取外部审核意见（用户操作）

### 决策门 DG-2（Day 19）
- 所有 TODO 是否解决？数字是否一致？
- 外部审核是否揭示关键缺口？若是 → 考虑 Phase 3

---

## Phase 3: 可选扩展（4/19 - 5/3，仅在需要时执行）

**进入条件**：用户决定需要更多广度（目标更高层级投稿），或 DG-2 揭示"单数据集"是致命弱点。

### 3A: GM12878 Shared-Region 对比（推荐优先项，若执行 Phase 3）

| 步骤 | 任务 | 耗时 | 硬件 |
|------|------|------|------|
| 3A.1 | 准备 GM12878 shared peaks/nonpeaks/bigWig | ~2h | 6000 |
| 3A.2 | TransChromBP GM12878 L3 训练 | ~6-8h | 6000 2×A6000 |
| 3A.3 | Official ChromBPNet GM12878 L3 训练 | ~6-8h | 6000 (序贯于 3A.2) |
| 3A.4 | 双方评估 + 对比表 | ~1h | 6000 |

### 3B: Genos Phase 1 → **建议归档为 future work**

理由：
- Phase 0 已是 yellow light，250+ GPU 小时产出大概率仍是负结果
- 论文已有充分的负结果描述（Section 5.6, 418-493 行 + probe diagnostics）
- ROI 极低，不值得投入到论文 deadline 紧张的阶段

### 3C: L3 第三 seed（低优先级）
若 6000 空闲，可跑第三 seed 收窄置信区间。机械性工作，不阻塞写作。

### 3D: 6002 独立数据集验证（可并行于 6000 工作）
GM12878/K562 在 6002 上的 smoke test，与 6000 上的 Phase 3A 并行。

---

## Phase 4: 提交准备（快速通道 4/19-4/25，完整通道 5/4-5/10）

| 任务 | 说明 |
|------|------|
| 4.1 出版级图表 | 所有 figure 导出 300dpi PDF/EPS |
| 4.2 补充材料 | S1: L3 全 epoch 指标; S2: ATAC 分类指标; S3: Genos probe; S4: L2 locked follow-up |
| 4.3 仓库清理 | 按 `docs/plan/git_cleanup_and_archive_plan_20260321.md` 执行：移除 .h5/.pt、清理 tmp_remote_edit/、更新 README |
| 4.4 LaTeX 编译 | 解决所有 warning（当前主要是 undefined `fig:arch`） |
| 4.5 合作者审核 + 单位确认 | 用户操作 |
| 4.6 选定目标期刊/会议 + 投稿 | 用户操作 |

---

## 时间线总览

### 快速通道（最小可行论文）

```
3/30 ──── Phase 0 ──── 3/31
              │
4/1  ──── Phase 1 (L3 锁定) ──── 4/4
              │
              ├── DG-1: L3 确认? ── 是 ──→
              │
4/5  ──── Phase 2 (论文冲刺) ──── 4/18
              │
4/19 ──── Phase 4 (提交准备) ──── 4/25
              │
        目标：4 月底提交
```

### 完整通道（含第二数据集）

```
3/30 - 4/4: Phase 0 + 1 (同上)
4/5 - 4/11: Phase 2A 开始（与 Phase 3 并行写作）
4/12 - 5/3: Phase 3 (GM12878 扩展)
5/4 - 5/10: Phase 2 完成（整合 Phase 3 结果）
5/11 - 5/17: Phase 4 (提交准备)

目标：5 月中旬提交
```

---

## 并行化策略

| GPU 密集 (6000) | GPU 密集 (6002) | CPU/本地 |
|----------------|----------------|---------|
| Phase 0: Official L3 训练中 | 空闲，可跑 smoke test | Phase 0.4: 架构图 |
| Phase 1: L3 评估 | Phase 3D: 独立数据集 | Phase 2: 全部写作 |
| Phase 3A: GM12878 训练 | | Phase 4: 图表制作 |

**核心原则：写作永远不等 GPU。**

---

## 风险与缓解

| 风险 | 影响 | 缓解措施 |
|------|------|---------|
| Official L3 训练失败/NaN | 无法完成 L3 对比 | 退回 L2 作为主对比 + 脚注 L3 准备中 |
| L3 count valid/test 分歧（类似 L2） | 主表可信度受损 | 立即执行 locked follow-up（已有成熟 pattern） |
| L3 结果矛盾（official 追平） | 需重新定位论文 | 诊断 sprint + 考虑写成"matched performance with richer diagnostics" |
| 论文写作速度不足 | 延误提交 | 优先填 P0/P1 TODO，P2 留给合作者 |
| Genos Phase 1 诱惑 | 浪费 250+ GPU 小时 | 坚持归档为 future work |

---

## 关键文件索引

| 文件 | 用途 |
|------|------|
| `reports/transchrombp_paper_draft_v1.tex` | 英文论文主稿（839 行，9 TODO） |
| `reports/transchrombp_paper_cn_v1.tex` | 中文论文主稿（833 行） |
| `reports/paper_claim_evidence_matrix_20260326.md` | 论证强度权威指南 |
| `reports/paper_rewrite_strategy_20260326.md` | 论文叙事策略 |
| `reports/tutorial_L2_locked_followup_20260329.md` | L2 follow-up 模板（L3 可复用） |
| `scripts/paper_aligned_repro/select_best_epoch.py` | L3 selector 脚本 |
| `scripts/paper_aligned_repro/summarize_metrics.py` | 指标聚合脚本 |
| `reports/assets/paper_metric_source_table_20260326.csv` | 数字溯源表（需 L3 更新） |
| `docs/plan/git_cleanup_and_archive_plan_20260321.md` | 仓库清理方案 |
| `TRACKING.md` | 实时状态（每个 phase 完成后更新） |

---

## 验证计划

每个 Phase 完成时的检查点：

- **Phase 0 完成**：TransChromBP L3 `best_epoch.json` 存在 + test JSON 存在 + Official L3 epoch 进度正常
- **Phase 1 完成**：L3 对比表完整（双方 test 指标）+ TRACKING.md 已更新 + 新报告已写入 `reports/`
- **Phase 2 完成**：论文 0 个 TODO + 数字全部有源文件 + claim-evidence 逐条核对通过 + LaTeX 编译无错误
- **Phase 4 完成**：仓库无敏感文件 + README 可复现 + 最终 PDF 生成
