# 在“不要吝惜 A6000×2 算力”前提下的 foundation 重开候选（2026-04-05）

> 更新（`2026-04-06 03:30 CST`）：本文档当时推荐的唯一 genuinely new-hypothesis 路线
> `NT v2 bins16 center-aligned residual` 已在 6000 双卡完成 full held-out。结果仍明显落后于
> matched `short10_nofoundation_control`：overall `count_pearson_full=0.7727` vs `0.8457`、
> `profile_target_jsd_full_mean=0.4625` vs `0.4395`；peak `count_pearson_full=0.7516` vs
> `0.8298`、`profile_target_jsd_full_mean=0.3588` vs `0.3193`。因此本文档现在只保留为
> “为何当时允许受控重入”的历史依据，不再构成默认继续 foundation 的许可。

## 1. 这份文档回答什么

当前 `restart v3` 已经在 tutorial canonical 上触发 `fail-or-unsafe`，因此“按原计划默认继续 foundation”这件事已经结束。

但如果用户**明确授权**：

> 可以继续占用 `6000 / A6000 x2`，前提是不是盲目重复旧 recipe，而是真正用算力去回答还没被回答的问题。

那么我们需要明确：

1. 哪些双卡实验仍然值得做；
2. 哪些看起来“还能再试试”，但其实只是重复旧失败；
3. 为什么这些候选不等于推翻当前停表规则。

一句话总判断：

> 不是“已经没有任何 foundation 实验值得做”，而是“只剩极少数 genuinely new-hypothesis 的实验值得做；重跑旧 `Genos/Caduceus/NT v2 short10 coarse-summary` family 不值得”。

---

## 2. 当前证据已经排除了什么

以下路线已经回答得足够清楚，不值得继续烧双卡：

### 2.1 不值得重开的

- 旧 `Genos` recipe
  - 结论已稳定：粗 summary 对 count 很危险，且没有 clean gain。
- 旧 `Caduceus-PS` recipe
  - 结论已稳定：工程上可跑，但只有 `near-null / marginal positive`，不足以支撑继续扩线。
- 当前 `NT v2 cached residual short10` 同配方重跑
  - 当前配方本质是：
    - `layer_07__bins4_mean`
    - `feature_tokens=4`
    - `FoundationResidualHead`
    - `profile_bin_count=16`
    - `encoded.mean + foundation.mean`
  - held-out 已直接掉到 peak `0.3560 / 0.7729`，不是“差一点点”，而是明显出界。
- 在当前配方上直接扩：
  - 第二个 seed
  - 更长训练
  - 同样 `bins4_mean` 的默认 `cross_attention`
  - 旧 `Genos/Caduceus` 再开多 seed

这些动作的共同问题是：

> 它们消耗更多算力，但不测试新的信息假设，只会让我们更昂贵地重复已经很清楚的负结果。

---

## 3. 现在还有哪些双卡实验值得做

只保留一个已完成的校准项，以及两条 genuinely new-hypothesis 候选。

### 3.1 已完成校准项：`short10 matched no-foundation control`

这不是因为它本身更有希望，而是因为它能回答一个仍有价值的问题：

> 当前 `short10` 训练预算本身到底有多差？

当前 `restart v3` 是在 `train_tutorial_foundation_short10.yaml` 这个短程预算下跑的。虽然它已经明显失败，但在本轮结束前，我们还没有一个**完全 matched 的 no-foundation short10 control** 去校准“短训练预算本身会损失多少”。

这条 run 的价值是：

- 给所有后续 foundation short-budget 结果提供真正可比的底座；
- 判断 `0.3560 / 0.7729` 里有多少是 foundation 伤害，有多少是 short10 预算伤害；
- 如果 short10 no-foundation 自身也很差，后续 foundation A/B 就必须重新定义预算和门槛。

它的定位是：

- 不是为了证明 foundation 有用；
- 而是为了防止我们误把“budget 不够”写成“foundation 无用”或反过来。

这条校准项现在已经完成，结果是：

- full held-out `n=62342`
- `count_pearson_full=0.8457`
- `profile_pearson_full_mean=0.7116`
- `profile_target_jsd_full_mean=0.4395`
- 分类指标 `auroc=0.8677`、`auprc=0.8647`、`f1=0.7902`

它回答的问题也已经明确：

> 当前 `short10` 预算本身并不差，至少不足以解释 `NT v2` 这条 coarse-summary residual 线为什么失守。

因此它现在的角色已经从“下一条待跑候选”变成：

> 所有后续 foundation short-budget gate 的固定校准底座。

### 3.2 候选 B（已执行并判负）：NT v2 `center-aligned / higher-bandwidth residual`

这在当时是唯一还值得认真花双卡回答的**新假设**：

> `restart v3` 失败，也许不是“NT v2 完全没有互补信号”，而是“当前 residual head 太粗、太低带宽、太全局化”。

当前实现真正做的是：

- `feature_name=layer_07__bins4_mean`
- `feature_tokens=4`
- `profile_bin_count=16`
- `FoundationResidualHead` 先把 local encoded 和 foundation tokens 都压成 summary，再做 coarse residual

所以它只否证了：

> coarse summary / low-bandwidth residual 没有给出净增益。

它**没有**否证：

- 更高带宽的 cached token summary（例如 `bins16_mean`）
- 更强的中心对齐 profile residual
- 不再把 local/foundation 都先做全局均值再融合的 residual 设计

若要重开，这条线必须满足：

- feature 改成 `bins16_mean`，而不是继续 `bins4_mean`
- residual 预测带宽显著提高，不能继续只做 `profile_bin_count=16`
- 明确强调中心区域对齐，而不是“再做一次全局 summary”
- 仍保持 residual / zero-init / low-disturbance 思路，不直接把 foundation model 升成 backbone 主干

### 3.3 候选 C（不再是默认下一步）：NT v2 `bins16` late cross-attention

如果用户明确说“算力不是问题，想把 token-level 路线也做一次真正像样的 gate”，那只保留一个版本：

- `feature_name=bins16_mean`
- `feature_layout=token`
- late `cross_attention`
- gate 继续保持强保守初始化

这条路线现在不再具有“自然下一步”的资格。若未来还要重开，必须先单独论证：为什么在
`summary / token-fusion / coarse residual / bins16 residual` 都已判负之后，它仍然回答的是
一个足够新的问题，而不是更昂贵地重复当前 family。

但它的优先级低于候选 B，因为：

- 风险更高；
- 更容易把 backbone 扰动和 feature 互补性混在一起；
- 如果连更保守的 higher-bandwidth residual 都没有起色，cross-attention 的解释价值会很差。

---

## 4. 当前不值得做、但容易误判为“可以多烧卡”的路线

以下路线最容易在“算力充足”的氛围下被误开，但当前不建议：

### 4.1 默认 `cross_attention` + 旧 `bins4_mean`

这不是新假设，只是把注入方式从 coarse residual 换成默认 cross-attention，而输入信息量仍然过低。

### 4.2 直接多 seed

在当前 effect size 还停留在“严重负结果或接近零结果”时，先上多 seed 只是更昂贵地确认负结果。

### 4.3 直接上 foundation 主干重构

当前任务是 base-resolution ATAC profile/count，现有 foundation feature 长度和空间语义都不天然对齐。现在直接让 foundation 主干接管模型，只会让“特征价值”和“主干扰动”完全缠在一起。

### 4.4 `distill_only` 立刻开跑

代码接口虽然在，但当前训练配置里：

- `teacher_target_names=[]`
- `distill_profile_weight=0`
- `distill_count_weight=0`
- `distill_rank_weight=0`

这说明 `distill_only` 还不是一个“拿现成配置就能直接跑的成熟路线”，而是要先定义 teacher target 语义和缓存产物。它更像下一阶段方法开发，不是当前最该占双卡的第一顺位。

---

## 5. 重新进入 foundation 线时的推荐顺序

如果用户坚持“可以继续用 A6000×2，把 foundation 线再推进一点”，推荐顺序固定为：

1. `short10 matched no-foundation control` 已完成，作为校准底座保留
2. 跑 NT v2 `bins16` center-aligned / higher-bandwidth residual（已执行，结果判负）
3. 原设想是只有第 2 步出现接近 clean gain 的信号，才跑 NT v2 `bins16` late cross-attention

这意味着：

> 截至 `2026-04-06 03:30 CST`，前两步都已经走完，且第 2 步判负；当前不再存在一个仍可直接顺延执行的默认下一步。

任何顺序里都不建议插入：

- 旧 `Genos`
- 旧 `Caduceus`
- 同配方 `bins4` residual 重跑
- 直接多 seed
- foundation 主干重构

---

## 6. 为什么这不违反当前停表规则

当前停表规则首先否定的是：

> 旧 `restart v3` 所代表的 coarse-summary / low-bandwidth residual family。

在完成 `bins16 center-aligned residual` 重入之后，当前可稳写的更新版结论是：

> 当前 `summary / token-fusion / coarse residual / bins16 center-aligned residual` 这些已实测 recipe family
> 都没有给出 clean gain。

它仍然没有逻辑上否定：

> 所有尚未实测的 token/base-level foundation 路线。

因此，若用户显式授权继续用双卡，唯一合理的做法不是“撤销停表”，而是：

> 在保留停表结论的前提下，只为 genuinely new-hypothesis 的下一条路线开一次受控重入。

---

## 7. 一句话建议

如果目标是“把双卡尽量用在最有信息增益的地方”，当前答案不是“随便再开一个 foundation run”，而是：

> 截至 `2026-04-06 03:30 CST`，这条“最窄、最有解释力的新路线”已经执行完且结果判负；当前默认动作回到论文主线，而不是继续保留一个开放的 foundation 候选。

---

## 8. 当前执行记录与结论（2026-04-05）

用户已明确授权继续使用 `6000 / A6000 x2` 推进 foundation 线，但要求不要吝惜算力的同时，也不要把算力浪费在旧 recipe 上。

本轮已完成推荐顺序中的第 1 步：

- 任务：`short10 matched no-foundation control`
- run name：`short10_nofoundation_control_s42_20260405_dual`
- 机器 / GPU：`6000 / GPU0,1`
- launcher：`scripts/run_short10_no_foundation_control.sh`
- 日志：`/data1/zhoujiazhen/bylw_atac/logs/short10_nofoundation_control_s42_20260405_dual_6000.log`
- 启动时刻：`2026-04-05 18:01 CST`

这条 run 的唯一目的不是追求新的最好指标，而是回答：

> 在和 `restart v3` 同样的 `short10` 预算、同样的双卡 DDP 语义下，不接 foundation 时基线到底能到哪里。

现在这个问题已经被回答：

- matched no-foundation control full held-out：
  - 输出：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/short10_nofoundation_control_s42_20260405_dual_best_test_full_20260405_2100.json`
  - 核心指标：`count_pearson_full=0.8457`、`profile_pearson_full_mean=0.7116`、`profile_target_jsd_full_mean=0.4395`
- `NT v2` bounded smoke held-out：
  - 输出：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/ntv2_residual_short10_gate_smoke_20260405_fix1_best_test_full_20260405_205905.json`
  - 核心指标：`n=64`、`count_pearson_full=-0.0980`、`profile_pearson_full_mean=-0.0142`、`profile_target_jsd_full_mean=0.6668`、`peak_auroc=0.5543`

因此当前结论收口为：

1. `short10` 预算本身不是主要问题；matched no-foundation baseline 仍然健康。
2. 当前 `FoundationResidualHead + layer_07__bins4_mean + feature_tokens=4 + profile_bin_count=16` family 已足够判为 `unsafe / not-worth-expanding`。
3. 当前 6000 两张 A6000 已空闲；若没有新的显式授权，不应继续把 GPU 花在旧 `bins4/coarse-summary residual` family 上。
4. 若用户未来仍明确要求继续 foundation，必须先提出一个不同于当前 `summary / token-fusion / coarse residual / bins16 residual` family 的新 hypothesis；不能把这份文档里的旧候选继续当成默认路线。
