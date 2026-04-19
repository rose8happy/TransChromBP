# ChatGPT Pro 外部分析回流与行动化（2026-04-05）

## 1. 这份文件是什么

这是一份“外部判断回流报告”，不是新的实验结果文件。

它的作用有三个：

1. 固定这轮 ChatGPT Pro 对 `chatgpt_bundle_project_handoff_20260405` 的完整外部判断，避免后续只记得口头摘要。
2. 把原始回流意见转写成当前仓库语境下可执行的高层结论。
3. 给后续 `TRACKING.md`、论文写作和 foundation-model 停表规则提供统一参照。

若本报告与后续新实验事实冲突，以更新的实验事实为准；若只是和旧叙事冲突，则以本报告和最新证据链为准。

---

## 2. 来源与上下文

- 外发材料：`reports/external/chatgpt_bundle_project_handoff_20260405/`
- 时间：`2026-04-05`
- 当前外部回复明显已经读到了：
  - `06_transchrombp_current_model.py`
  - `07_foundation_adapter_current.py`
- 因此这轮意见不是只基于高层摘要，而是已经结合了当前模型本体和 `FoundationResidualHead` / `FoundationCrossAttentionAdapter` 的实现语义。

---

## 3. 原文保留

下面保留用户回传的 ChatGPT Pro 原文，后续如需核对措辞或判断是否误读，应回到这段原文：

```text
我看完之后的判断很明确：

**现在最值得做的是把论文主线收口；最不值得做的是继续追逐“Transformer 特有 shortcut”这种戏剧化叙事。**
你真正已经成立的，不是一个强机制发现，而是一个**把 Transformer 安全接入 ChromBPNet-style bias factorization，并用 full/debiased 双口径做严谨诊断**的框架。这个主线同时有状态总结、实验历史和代码结构三层支撑。  

更具体地说，你现在的“真成果”有五块。第一，Transformer backbone 本身有真实收益：matched backbone ablation 里，TF 版本的 peak JSD 从 0.32393 降到 0.31419，count r 从 0.82503 升到 0.83676，这不是噪声级差异。第二，`full/debiased gap` 确实是必须保留的诊断口径，单看 full 指标会误判。第三，`corrected B` 是稳定默认模型。第四，`center pool` 比 attention pool 更稳。第五，L3 shared-region 的 matched system compare 很强：你们版本是 0.31319 / 0.84016，official controlled L3 是 0.33853 / 0.69958。 

所以如果我要替你重写一句论文主线，我会写成：
**“A bias-safe Transformer framework for base-resolution ATAC profile/count prediction, with explicit full/debiased diagnostics and stable readout design.”**
而不是“我们发现并修复了 Transformer 特有 shortcut”。后者已经被 clean matrix 实质性降级了：conv-only 也会出问题，`stop-gradient` 也不是单独的万灵药。 

这里还有一个我觉得很关键、但容易被你们自己忽略的点：
**你们赢的并不只是 backbone。你们赢的是“架构 + 训练语义 + 诊断口径”的整套系统。**
从 06 的代码看，当前模型本体不是“纯 Transformer”，而是 `conv stem + local dilated tower + transformer + bias branch` 的混合体，而且 `full/debiased` 是结构化输出，不是日志里临时多算两个指标。换句话说，任何把成果简化成“Transformer 机制更强”或“纯 transformer 替代 conv baseline” 的叙事，都会被 reviewer 质疑为过度归因。 

我会把你们当前结论分成三层。

**可以强判断的：**
Transformer 有真实收益；`full/debiased gap` 有必要；`corrected B` 稳；`center pool` 稳；L3 shared-region 对外部系统对比成立。独立 GM12878 / K562 结果也足够支持“不是只会 tutorial”，但只能写成**non-collapse / usable transfer sanity**，不要写成跨数据集 benchmark win。 

**只能弱判断的：**
当前 clean matrix 下，`debiased profile supervision` 可能比 `stop-gradient` 更关键；NT v2 比 Genos 更像一个还没被完全否掉的候选；foundation line 还没死，但只剩“是否存在低干扰互补信号”这一个窄问题值得继续。 

**当前不该再坚持的：**
“Transformer 特有 shortcut”“conv-only 天然安全”“stop-gradient 单独修复问题”“我们已经证明了强 shortcut 机制”。尤其从 06 的实现看，`profile_bias_stop_gradient` 只是在 **profile 融合路径**里对 bias 分支做 detach，count 融合并没有因此获得一个通用的“机制修复器”；这和 clean matrix 对 stop-gradient 叙事的降级是相互吻合的。 

对 foundation model 这条线，我的判断是：

**它不值得再做成新的项目主线，但值得保留一次“有停表规则的最后验证”。**
Genos 这条旧主线基本已经回答完了：模型本体没坏，本地 OCR sanity 也对，但它的信号形态更像粗 summary，不自然补强你们这个 base-resolution profile/count 任务；而且一旦打进 count 路径，count 很容易崩。Caduceus 当前 recipe 说明“在线 token fusion 可以不炸”，但只有 near-null / marginal positive，不够支撑扩线。NT v2 的位置稍微特殊：它比 Genos 更像真有独立信号，但目前只证明了“有信号”，没证明“有互补净增益”。

我会这样概括这三条支线各自“真正说明了什么”：

1. **Genos**：说明“独立信号存在”不等于“对主任务可用”；也说明粗 summary 注入对 count head 很危险。它没有说明“所有 foundation model 都没用”。
2. **Caduceus**：说明 token-level 融合在工程上可以安全，但在当前 recipe 下增益太小，不够过门槛。它没有说明“token-level 路线永远无效”。
3. **NT v2**：说明它比 Genos 更像候选，但简单 concat / simple probe 仍然过不了 complementarity gate。它没有说明“保守 residual/distill 路线也一定无效”。

从代码角度看，restart v3 的方向是对的。
`FoundationCrossAttentionAdapter` 是 residual-gated，而且输出投影零初始化；`FoundationResidualHead` 的 profile/count 输出也零初始化。这个设计的价值很大：它让实验更接近回答“foundation feature 是否能提供**净增益**”，而不是“你是不是一开始就把 backbone 搞坏了”。这比之前的粗暴注入更科学。 

但我也会提醒一个你们现在最可能继续自我误导的地方：
**当前 residual head 实际测试的是“粗全局 summary 的互补性”，不是“真正 token/base 级互补性”。**
因为 07 里的 `FoundationResidualHead` 会先把 `encoded` 做全局均值池化，再把 foundation tokens 也压成 summary，然后只预测很粗的 `profile_bin_count=16` 残差，再线性插值回 1000bp。这个设计非常保守，很适合做“最后一次低风险判定”；但如果它跑不出来，只能否证“粗全局 residual 没用”，不能一口气否证所有更细粒度、中心对齐、token-level 的 foundation 接法。

所以我的路线建议很简单，只留 3 个动作：

1. **立刻把论文主线定稿为 bias-safe framework + diagnostics，不再围绕 dramatic shortcut 写。**
   主证据就用：matched backbone ablation、clean matrix、readout 选择、L3 shared-region compare。foundation 只当 appendix 或 future work。 

2. **把 `NT v2 cached residual short10 restart v3` 当成 foundation 线的最后一个正式 gate。**
   继续条件应该很严：必须在 held-out 上给出干净净增益，而且不能扩大 `full/debiased gap`，也不能再碰 count stability 的红线。否则就停。你们自己对 Caduceus 已经有 `0.002` 级门槛意识了，这里应该沿用甚至更严格。

3. **只有当 v3 真有 clean gain，才做一个唯一的 follow-up：把 07 的 residual head 改成“中心对齐、较高带宽”的版本再做 matched A/B。**
   不是再开 Genos 新坑，也不是扩 Caduceus 旧 recipe，更不是上 foundation 主干。因为当前任务是 base-resolution profile/count，而现有 pretrained 特征长度和分辨率都不天然对齐；在这种前提下，直接让 foundation model 接管主干，风险远大于信息量。 

一句话收束：

**你这个项目已经足够像一篇论文了；现在最危险的不是证据不够，而是叙事写偏。**
主线该收口在“安全接入、严格诊断、稳定设计”，而不是“夸张机制发现”。foundation model 还可以留一个 NT v2 residual 的判定机会，但它现在只配做受控 side quest，不配再当主线。

如果只再看一份真实代码，我会优先看 **07**，重点验证一件事：你现在到底在测试“foundation 有无互补信号”，还是只是在测试“粗 summary 能不能救一点残差”。
```

---

## 4. 归一化后的结论

### 4.1 总判断

这轮外部意见和我们仓库里到 `2026-04-05` 为止最稳的证据是高度一致的：

- 论文主线已经足够成立。
- 真正需要收口的是叙事，不是再去制造新的 dramatic 机制故事。
- foundation-model 这条线不再配当主线，只配保留一个有停表规则的 side quest。

### 4.2 可以强写的

1. Transformer backbone 在当前任务里带来真实收益。
2. `full/debiased gap` 是必要诊断。
3. `corrected B` 是稳定默认模型。
4. `center pool` 比 attention pool 更稳。
5. `L3 shared-region` 的 matched system compare 成立。
6. GM12878 / K562 只能写成 `non-collapse / usable transfer sanity`，不能写成跨数据集 benchmark win。

### 4.3 只能弱写的

1. 当前 clean matrix 下，`debiased profile supervision` 比 `stop-gradient` 更关键。
2. NT v2 比 Genos 更像还未被否掉的候选。
3. foundation 线只剩“是否存在低干扰互补信号”这一个窄问题值得继续。

### 4.4 当前不该再坚持的

1. `Transformer-specific shortcut`
2. `conv-only 天然安全`
3. `stop-gradient 单独修复`
4. `已经证明了强 shortcut 机制`

---

## 5. 这份外部意见与当前仓库证据如何对齐

### 5.1 与 clean matrix 对齐

外部意见把 shortcut 叙事从“强机制发现”降级为“不能再强写”，这和当前 `clean matrix` 的主结论一致：

- `conv-only unsafe` 并不比 TF 更安全
- `stop-gradient` 不是单独的万灵药
- `full/debiased gap` 是必要诊断，而不是可有可无的附加分析

### 5.2 与当前主模型实现对齐

外部意见强调“你们赢的是整套系统，不只是 backbone”，这也和 `06_transchrombp_current_model.py` 一致：

- 当前主模型不是纯 Transformer
- `full/debiased` 是结构化输出
- bias branch 是显式模型语义，不是后处理补丁

### 5.3 与 foundation 支线结果对齐

- Genos：相对完整地失败，支持“粗 summary 不自然补强当前任务”
- Caduceus：安全但增益太小
- NT v2：独立信号更强，但仍未出现互补净增益证据

这和当前 `reports/analysis/nt_v2_probe_20260404.md`、`reports/analysis/genos_cached_p1_restart_20260403.md`、`docs/plan/archive/foundation/caduceus_tutorial_ab_online_fusion_20260331.md` 的阶段结论是一致的。

---

## 6. 对当前工作的直接影响

### 6.1 论文主线

后续论文与 supporting writeup 应显式收口到：

> bias-safe Transformer framework + full/debiased diagnostics + stable readout design

而不是再围绕 shortcut dramatic story 做包装。

### 6.2 foundation line

`NT v2 cached residual short10 restart v3` 当时是当前 foundation 线唯一还值得跑完的正式 gate。

后续落实情况（`2026-04-06 03:30 CST`）：

- 建议中的唯一 follow-up `NT v2 bins16 center-aligned residual` 已完成 full held-out；
- 结果仍明显落后于 matched `short10_nofoundation_control`：
  - overall `count_pearson_full=0.7727` vs `0.8457`
  - overall `profile_target_jsd_full_mean=0.4625` vs `0.4395`
  - peak `count_pearson_full=0.7516` vs `0.8298`
  - peak `profile_target_jsd_full_mean=0.3588` vs `0.3193`
- 因此这份外部意见现在应理解为“主线收口 + foundation 严格停表”的支持证据，而不是仍待执行的后续建议。

若从当前时点回看，这条线的可执行结论已经变成：

- 不默认扩 `seed1234`
- 不默认切 `cross_attention`
- 不重开旧 `Genos` / `Caduceus` recipe

### 6.3 代码层特别提醒

外部意见最值得我们额外重视的一点是：

> 当前 `FoundationResidualHead` 更像在测试“粗 summary residual 是否有互补性”，而不是在测试“更一般的 token/base-level foundation 接法是否有效”。

这意味着未来若 `restart v3` 失败，我们最多只能说：

- 这条“保守 coarse-summary residual”路线没有净收益

而不能直接扩写成：

- 所有更细粒度、中心对齐、token-level foundation 路线都没用

---

## 7. 一句话行动化结论

这轮 ChatGPT Pro 外部意见应被视为：

> 对当前主线“论文收口优先、foundation 线严格停表”的强支持，而不是新的扩线许可。
