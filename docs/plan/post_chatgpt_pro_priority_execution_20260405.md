# 基于 ChatGPT Pro 外部意见的后续优先级执行计划（2026-04-05）

> 更新（`2026-04-06 03:30 CST`）：本文档定义的 foundation 决策链
> `restart v3 -> matched no-foundation control -> bins16 center-aligned residual`
> 已经全部走完，且第 3 步仍明显落后于 matched control。当前默认动作因此已经从
> “保留一个窄门继续试”收口为“foundation 线不再自动分配 GPU，只保留为受控负结果”。

> 当前权威读法（`2026-04-06`）：
> 1. 本文档仍是 foundation 停表规则的默认权威入口；
> 2. 其中 `3.3-3.5` 保留为“当时允许的一次受控后续链”的历史决策树，不再代表待执行动作；
> 3. 当前真实优先级已经切回“论文主线与文档收口优先，foundation 仅保留为 appendix / future work 边界负结果”。

## 1. 总目标

从这一刻起，后续工作围绕两条主轴展开：

1. 避免 `A6000` 训练窗口空置浪费，用最短反馈周期回答“预训练基因组大模型到底该怎么用，或者当前阶段应不应该再用”。
2. 把 ChatGPT Pro 的回流意见真正转化成行动规则，而不是只停留在“外部意见说得对”。

这份计划默认接受一个前提：

> 论文主线已经足够成立，foundation-model 线不再是新的项目主线，只是有严格停表规则的受控 side quest。

补充说明（`2026-04-05 18:01 CST`）：

- 本文档定义的是“默认自动执行”的停表与扩线规则。
- 若用户后续给出显式授权，允许在 `fail-or-unsafe` 后按受控重入继续占用 A6000，则以
  `docs/plan/a6000_foundation_reentry_candidates_20260405.md` 与 `TRACKING.md` 的最新记录为准。

补充说明（`2026-04-05` 当日晚间结果）：

- `short10 matched no-foundation control` 的 full held-out 已证实 `short10` 预算本身仍可维持健康 baseline；
- `ntv2_residual_short10_gate_smoke_20260405_fix1` 的 64-region bounded smoke 则给出 near-random / negative correlation。

因此，本文档定义的“默认停表”规则现在不只是原则判断，也已经被最新结果实证支持；后续若再进 foundation，只能按显式授权引用
`docs/plan/a6000_foundation_reentry_candidates_20260405.md` 的新假设重入方案。

---

## 2. 当前优先级排序

### P0：把 foundation 停表结论写回主线

当前最高优先级已经不再是继续发车，而是把已经得到的负结果统一写回：

- 论文主线收口
- claim matrix / supporting writeup 改写
- foundation 停表规则固定

### P1：维护论文主线优先，foundation 退回受控 side quest

在最新结果下，默认资源分配应当回到：

- `L3 shared-region` system compare 写回主文与 supplementary
- clean matrix / readout / bias-safe framework 这些 paper-facing 主证据
- 远端结果归档与工作目录整理

### P2：历史上唯一允许的一条 foundation 后续链（已走完，保留为记录）

当前允许的 foundation 后续链只有：

`restart v3` -> `matched no-foundation short10 control`（仅当当前 run 有潜力） -> `中心对齐、高带宽 residual head A/B`（仅当前两步都过线）

除此之外，不允许自动切到：

- `seed1234`
- `cross_attention`
- `Genos`
- `Caduceus`
- foundation 主干重构

---

## 3. 当前 run 结束后的固定决策树

## 3.1 第一步：先做 run 收口，不先脑补

当前 `NT v2 cached residual short10` 结束后，第一批必须读取的产物固定为：

- `best.pt`
- `epoch_metrics.jsonl`
- best checkpoint 对应的 held-out `test-full` 指标

在 held-out 结果出来前，不做“它看起来好像有戏”的扩线决定。

---

## 3.2 第二步：先判定是否直接停表

用当前 paper-facing `corrected B` 参考锚点做第一道硬过滤：

- 参考：
  - peak `JSD ~= 0.3145-0.3147`
  - peak `count_r ~= 0.8488-0.8503`
  - clean gap `~= 0.00268`

若 `restart v3` 的 held-out 结果满足任一条件，直接判为 `fail-or-unsafe`，本阶段停表：

1. peak `profile_target_jsd_full_mean > 0.3197`
2. peak `count_pearson_full < 0.8400`
3. `full/debiased gap` 明显放大到 `> 0.010`
4. 训练过程出现 count 不稳定、gap 异常抬升或其它已知 unsafe 信号

判为 `fail-or-unsafe` 后的动作固定为：

- foundation 线本阶段停止
- 不自动扩 `seed1234`
- 不自动切 `cross_attention`
- 不重开旧 `Genos/Caduceus` 线
- A6000 后续转给 paper-aligned 的更高价值 GPU 工作，而不是继续烧 foundation 试错

---

## 3.3 第三步：只有“看起来有潜力”时，才允许开 matched control

如果当前 run 没有触发上面的 `fail-or-unsafe` 条件，只能说明它“值得继续判定”，还不能说明它已经有 clean gain。

这时下一条唯一允许的 A6000 任务是：

> 跑一个 matched `no-foundation short10 control`

固定语义：

- 训练配置继续用 `configs/train/train_tutorial_foundation_short10.yaml`
- backbone / data semantics / DDP / batch / early stop 全部保持一致
- 只把模型配置改回 `configs/model/transchrombp_teacher_v2_center_pool.yaml`
- seed 固定为 `42`
- 用同样的 held-out `test-full` 口径出结果

这一步的目的不是新开一条研究线，而是给当前 `restart v3` 找到真正可比的 short10 对照。

---

## 3.4 第四步：只有 clean gain 才能进入唯一 follow-up

补充说明（`2026-04-06 03:30 CST`）：

- 这一步的前提在最新结果下并未成立；
- `bins16 center-aligned residual` 没有形成 clean gain；
- 因此下面这段“唯一 follow-up”定义现在只保留为历史决策树，不再代表默认后续动作。

`restart v3` 相对 matched no-foundation short10 control 同时满足下面 3 条，才算 `clean gain`：

1. peak `profile_target_jsd_full_mean` 改善 `>= 0.002`
2. peak `count_pearson_full` 不允许下降超过 `0.005`
3. `full/debiased gap` 绝对值保持 `<= 0.010`，且相对 matched control 不增加超过 `0.002`

只要有任一条不满足，都按 `near-zero / not-worth-expanding` 处理。

`near-zero / not-worth-expanding` 的动作固定为：

- foundation 线本阶段停表
- 不扩第二个 seed
- 不转 `cross_attention`
- 不重开旧 foundation recipe

---

## 3.5 第五步：唯一允许的 follow-up 长什么样

补充说明（`2026-04-06 03:30 CST`）：

- 当前并不存在一个仍然开放的默认 follow-up；
- 若未来还要重开 foundation，必须先写出一个不同于当前 residual short10 family 的新 hypothesis；
- `bins16 late cross_attention` 不再是“自然下一步”，只能在新的独立论证里单独决定。

只有 `clean gain` 成立，才允许一个唯一 follow-up：

> 把当前 `FoundationResidualHead` 改成“中心对齐、较高带宽”的版本，再做 matched A/B

这个 follow-up 必须同时满足：

- 仍然是 residual 路线，不换成 foundation 主干
- 不回头重开 `Genos` / `Caduceus`
- 不直接跳到大规模多 seed
- 目标是回答“粗 summary residual 不够时，中心对齐的更高带宽 residual 是否才是正确使用方法”

在这个 follow-up 跑完前，仍不扩 `seed1234`。

只有当它也给出 clean gain，才允许进入第二个 seed 或更长训练。

---

## 4. 论文与 foundation 线的并行分工

## 4.1 论文主线

从现在开始，论文相关文档默认按下面这句主线改写：

> A bias-safe Transformer framework for base-resolution ATAC profile/count prediction, with explicit full/debiased diagnostics and stable readout design.

主证据固定为：

- matched backbone ablation
- clean matrix
- readout design
- `L3 shared-region` system compare

foundation 只保留为：

- appendix
- future work
- 或带停表规则的 side quest

## 4.2 foundation 线

foundation 线的使命被压缩成一句话：

> 回答“是否存在低干扰、可复现、值得继续的互补信号”，而不是再去争取成为主故事。

---

## 5. 资源使用规则

1. 当前双卡 run 活着时，不拆任务，不抢卡。
2. 当前若有新的 foundation run 被显式授权，也必须在同一工作块内完成“读结果 -> 归类 -> 明确停表或单独论证下一假设”。
3. 若 foundation 线被停表，A6000 不得因为“还没想好下一条 foundation 任务”而空置；此时直接切到 paper-aligned GPU 工作或其它已确认高价值任务。
4. 任何人若想重开 `Genos`、旧 `Caduceus` recipe、默认 `cross_attention`、或 foundation 主干路线，必须先显式说明它为什么不违反本计划的停表规则。

---

## 6. 一句话执行标准

后续所有 foundation 决策都按下面这句检查：

> 这一步是在用最短 GPU 反馈周期验证“低干扰互补信号是否存在”，还是又在没有停表规则地扩一条旧支线？
