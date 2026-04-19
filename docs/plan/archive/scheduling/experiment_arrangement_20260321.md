# 实验训练任务安排方案（Genos 优先级提升版）

> 日期：2026-03-21
> 目标：基于项目当前进度、两台服务器的硬件能力，在不打断已启动任务的前提下，提升 Genos 相关实验优先级
> 原则：Genos 相关实验必须放 A6000；当前已启动任务不主动中断；在下一轮 GPU 窗口里，Genos Phase 1 优先于 V2fix 继续扩 seeds，也优先于数据扩容线

---

## 一、硬件约束总览

| 服务器 | GPU | 显存 | 单 run 时间 (40ep) | 适合的任务类型 |
|---|---|---|---|---|
| 6000 (A6000×2) | 2×NVIDIA A6000 | 48GB each | ~2.3h (DDP) | 重型实验、多 seed 串行、Genos 推理 |
| 6002 (3080×1) | 1×NVIDIA RTX 3080 | 12GB | ~19h (单卡) | 轻量实验、单 seed smoke test、小规模消融 |

> **重要**：3080 单个 40-epoch 训练需要 **~19h**，是 A6000 DDP 的 **~8 倍**。分配到 3080 的实验必须严格控制数量和规模。

---

## 二、当前各任务线状态

### 2.1 TransChromBP V2fix 消融系列（当前在跑）

| 实验 | 优先级 | 当前状态 | 服务器 |
|---|---|---|---|
| **E: Checkpoint Soup** | P0 | ✅ 已完成：soup 零成本评估已出结果（`soup_top3` JSD 基本持平 `best`，`count_pearson_debiased` 更好） | 6000 |
| **A: Freeze 修复 + 40ep** | P1 | ✅ 已完成，held-out peak `JSD_full=0.31479` | 6000 |
| **B: Count Center Pool + 40ep** | P1 | ✅ 已完成，held-out peak `JSD_full=0.31467` | 6000 |
| **F: Attention Pool + 40ep** | P1 | ⏳ `v2fix_20260321_attnpool_s42` 已启动 | 6000 |
| **G: Profile Refine + 40ep** | P1 | ⏳ `v2fix_20260321_profref_s42` 已启动 | 6002 |
| **C: pscale=0.1** | P3 可选 | ⏸ 降级，除非出现意外 | — |
| **D: BN→GN** | P3 可选 | ⏸ 降级 | — |

**V2fix 当前策略**：
1. 不打断已经启动的 `F_s42` / `G_s42`
2. `F/G` 单 seed 收口后，暂不立刻扩多 seed
3. `A/B` 多 seed、`F/G` 多 seed 均暂缓到 Genos Phase 1 得出首轮结论之后

### 2.2 Genos 融合系列（优先级已提升）

| 阶段 | 状态 | 结论 |
|---|---|---|
| Phase 0 可行性（正式 500+500） | ✅ 已完成 | **黄灯**：`layer_6` AUC=0.6717（CV）/ 0.6538（holdout），仅边缘可用 |
| Phase 1 最小 frozen pilot | 🔜 下一优先级 | 在当前 `F_s42/G_s42` 结束后立即接棒 |

**Genos 的资源特点**：
- Phase 0 独立推理粗估：Genos-1.2B 在 `batch=8`、双向设置下约 **21.8 GB 显存**；流式版 Phase 0 的 host RSS 约 **7.64 GB**
- 但当前 Phase 1 实际训练链路的 `G1 batch=20` 50-step quick sweep 峰值仅 **8.1 GiB**
- 旧版脚本曾出现 **~45 GB RAM**，但那是“缓存五层完整 embedding”的旧实现，现已不再是当前约束
- 3080 当前不作为 Genos 首选，不是因为“绝对装不下”，而是因为 wall-clock 更差，且 6002 已有独立任务队列
- A6000 仍是首选，主要因为吞吐和调度优先级

**优先级提升后的解释**：
- 虽然 Phase 0 仍然是黄灯，不是绿灯，但其结果已经足以支持“做一个受控、最小、frozen 的 Phase 1 pilot”。
- 预训练基因组大模型这条线具有更高的潜在上界与论文价值，因此在下一轮 A6000 窗口里，应优先回答它是否对当前主线有净增益。
- 这不改变技术边界：优先级提高，不等于风险降低。Phase 1 仍必须是 matched baseline + 最小接入，不得直接跳 LoRA / cross-attention / 长窗口。

### 2.3 6002 上已完成的实验

| 实验 | 状态 |
|---|---|
| Bias-Safe 2×2 析因（B1-B4） | ✅ 全部完成，报告已出 |
| tutorial 单卡 smoke test | ✅ 已验证 3080 可以跑 TransChromBP 训练 |

---

## 三、推荐任务分配方案（更新后）

### 3.1 A6000 (6000) — 先收当前 run，再切 Genos 主线

```
阶段 1（当前已启动，不主动中断）—————————————————
  ├── `F: attnpool_s42` 在 6000 双卡继续跑完
  └── `G: profref_s42` 在 6002 单卡继续跑完

阶段 2（下一轮最高优先级：Genos 实现 + smoke）————————
  ├── 当前一块 A6000 窗口释放后，立刻切到 Genos Phase 1 代码接入
  ├── 完成 `genos_adapter.py` / `forward(genos_feat=...)` / trainer runtime
  ├── 做 `G0/G1/G2` 三组最小 smoke，确认：
  │     → 当前主线代码接入正常
  │     → matched baseline / genos_gate / genos_mean 均可解析
  │     → validation metrics / gate / full-debiased gap 能正常写出
  ├── 2026-03-21 当前状态：上述工程改造、corrected smoke、`G1` quick batch sweep 已全部完成
  └── 当前已知正式 pilot 口径：单卡、非 DDP、统一 `batch=20`

阶段 3（Genos Phase 1 pilot）——————————————————
  ├── G0 baseline_v2：matched baseline，20ep，单卡 A6000，先跑
  ├── G1 genos_gate：layer_6，20ep，单卡 A6000，`bs=20` 的 50-step quick sweep 为 `71s / 8090 MiB`
  ├── G2 genos_mean：layer_6，20ep，单卡 A6000，预计与 `G1` 同量级
  └── 推荐排法：
        → 先收掉 G0
        → 再把两张 A6000 并行给 G1 / G2
        → 不再沿用旧的 `20–24h/run` 估计；按当前 `G1` step time 粗算，单个 20-epoch run 可能是 `~82–93h`

阶段 4（Genos 之后再决定其它线）———————————————
  ├── 若 `G1 > G0` 且 `G1 > G2`：
  │     → 继续保留 Genos 高优先级，可考虑扩第二个 seed 或做更稳的 pilot
  ├── 若 `G1 ≈ G2 > G0`：
  │     → 说明主要收益更像粗粒度 global signal，再决定是否做更便宜的替代
  └── 若 `G1 ≈ G2 ≈ G0`：
        → Genos 线降级，回到 V2fix 多 seed / 数据扩容
```

> **注意**：文档里同时保留了两类数字，但语义不同。`~21.8 GB` 指的是 Phase 0 独立推理粗估；`8.1 GiB` 指的是当前 Phase 1 真实训练链路在 `G1 batch=20` 下的 50-step quick sweep。对正式 pilot 的调度判断，应以后者为准。当前真正的约束已从“显存装不下”转向“wall-clock 太长”。

### 3.2 RTX 3080 (6002) — 继续承接 Genos 无关的独立实验

**适合 3080 的任务**（单 run ≤ 19h，且不依赖 A6000 权重同步）：

| 任务 | 耗时估计 | 备注 |
|---|---|---|
| ① Genos 无关的 TransChromBP 单 seed 实验 | ~19h/run | 继续消化原计划里适合单卡的实验 |
| ② 独立数据集实验（GM12878/K562） | smoke ~2-5h；40ep ~19-22h | 不依赖 A6000 checkpoint |
| ③ Smoke test / 代码部署验证 | ~1-2h | 快速确认链路可通 |

**推荐给 3080 的具体任务**：

```
推荐方案（Genos 提优后）：

  场景 1：当前
    → 继续让 `G: profref_s42` 在 3080 上自然收口

  场景 2：A6000 全力给 Genos pilot 时
    → 3080 不再承担 Genos
    → 也不默认承担 A6000 模型的从属评估
    → 改为继续做 Genos 无关的独立 TransChromBP 实验

  场景 3：优先顺序
    → 先收掉当前 `G_s42`
    → 然后优先接 `GM12878-only` smoke
    → smoke 正常后，再跑 `GM12878-only` 单 seed baseline
    → `K562` 先做 smoke，不急着直接上完整长跑

  场景 4：Genos 首轮结论出来之后
    → 若 Genos 无收益，再决定 3080 是继续独立数据集实验，还是回补剩余 V2fix 单 seed
```

> **避免**：
> 1. 把 3080 当成 A6000 的默认评估从属机
> 2. 为了评估而频繁同步大 checkpoint
> 3. 在 3080 上串行排 3+ 个长训练 seed → 总耗时 57h+，效率太低

### 3.3 两台机器的协同调度

| 时段 | 6000 (A6000×2) | 6002 (3080×1) |
|---|---|---|
| **当前 ~ F_s42 结束** | 继续跑 `F_s42`；不打断 | 继续跑 `G_s42`；不打断 |
| **F_s42 结束后当天** | 当前已完成：修正 Genos launcher / matched baseline 配置、`G0/G1/G2` corrected smoke、`G1` quick sweep；下一步是 `G0(batch=20)` | 若 `G_s42` 仍在跑，则继续等待其收口 |
| **随后 3-4 天** | `G0` 收口后，两张 A6000 分别单卡跑 `G1 genos_gate` / `G2 genos_mean`（统一 `batch=20`，不启用两卡 DDP） | `G_s42` 结束后切 `GM12878-only` smoke / baseline |
| **Genos 第一批结论后** | 决定继续 Genos 还是回补 V2fix / 数据线；若 `20ep` wall-clock 不可接受，则先缩短 epoch 再继续 | 继续独立数据集实验，或再接剩余的 Genos 无关单 seed |

---

## 四、关键决策点

### 决策 1：Genos Phase 1 是否立刻接棒

当前答案是：`是，但只限最小 frozen pilot`。

原因：
- Phase 0 仍然只是黄灯，但已经足以支持受控 pilot；
- 预训练基因组大模型的潜在收益和研究价值更高；
- 当前 `A/B` 已经通过 held-out test 证明“差距很小”，继续扩多 seed 的边际信息增益不如先回答 Genos 是否有净增益。

因此新的顺序是：
1. 不打断当前 `F_s42/G_s42`
2. 它们一结束，A6000 立即切到 Genos G0/G1/G2
3. 只有 Genos 第一批给出明确负结论后，才回头扩 V2fix 多 seed

### 决策 2：Genos Phase 1 的边界

Phase 0 结果是黄灯（AUC ~0.67，当前候选层位为 `layer_6`），所以即便优先级提高，进入 Phase 1 仍需要满足 checklist（`docs/plan/archive/genos/genos_phase1_rebase_checklist.md`）：
- 代码基于当前主线（`vendor/transchrombp/transchrombp`）
- matched baseline 对照（G0/G1/G2 三组）
- Genos frozen + no_grad + eval()

同时应遵守最小接入方案（`docs/plan/archive/genos/genos_phase1_minimal_integration_plan.md`）：
- trainer 外挂 frozen extractor
- model 内只做 gated fusion
- `G0/G1/G2` corrected smoke 与 `G1` quick sweep 已经完成
- 正式 pilot 当前锁定为单卡 `batch=20`，并通过“两张 A6000 各跑一条单卡 run”利用资源，不用两卡 DDP 改变单个 run 的训练口径

**建议时机**：在当前已经启动的 `F_s42/G_s42` 收口后立即启动，不再等待 V2fix 多 seed。

### 决策 3：3080 参与方式

| 方案 | 优点 | 缺点 | 推荐程度 |
|---|---|---|---|
| 主要做从属评估 | 不需要新设计 | 依赖权重同步，性价比低 | ⭐ |
| 跑独立单 seed / 独立数据集实验 | 真正并行，不依赖 A6000 checkpoint | 单 run 仍较慢 | ⭐⭐⭐（推荐） |
| 跑 3+ 个长训练 seed | 最大化利用 | 占满 3080 3+ 天，且优先级不高 | ⭐ |

---

## 五、总结

| 服务器 | 分配原则 | 主要任务 |
|---|---|---|
| **A6000 (6000)** | 当前 run 收口后，优先切 Genos | `F_s42` 收口后，立即执行 Genos G0/G1/G2 最小 smoke/pilot |
| **3080 (6002)** | 继续承接 Genos 无关的独立实验 | 跑完 `G_s42` 后，优先做 `GM12878/K562` 等独立数据集实验或其它不依赖 A6000 权重的单 seed 实验 |

**核心策略**：当前已启动的 `F_s42/G_s42` 不主动中断；它们结束后，下一轮 A6000 训练窗口优先给 Genos Phase 1。与此同时，3080 不再默认做评估从属机，而是继续推进 Genos 无关、无需频繁同步大 checkpoint 的独立 TransChromBP 实验。Genos 相关一切仍只放 A6000。

---

## 六、相关文档

- 实验详细设计：`docs/plan/archive/v2/v2_code_improvement_ablations.md`
- 读出头实验设计：`docs/plan/archive/v2/v2fix_readout_head_experiment_design.md`
- Genos Phase 1 checklist：`docs/plan/archive/genos/genos_phase1_rebase_checklist.md`
- Genos 最小接入方案：`docs/plan/archive/genos/genos_phase1_minimal_integration_plan.md`
- 6002 独立实验线清单：`docs/plan/archive/scheduling/6002_independent_experiment_queue_20260321.md`
- 项目追踪清单：`TRACKING.md`
