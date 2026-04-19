# TransChromBP 内部设计与实验思路长报告（2026-03-31）

> 目的：这不是对外论文，也不是阶段总结 PPT，而是一份给自己回看用的“连续复盘文档”。
>
> 核心用途有两个：
>
> 1. 阶段性回看代码时，知道当前每个关键模块和配置是为什么出现的。
> 2. 阶段性思考论文时，知道哪些结论是后来保住的，哪些结论已经被推翻或降级。

如果这份报告与更晚的结果文件出现冲突，以最新的 source table 和收口报告为准：

- [paper_metric_source_table_20260326.csv](assets/paper_metric_source_table_20260326.csv)
- [paper_claim_evidence_matrix_20260326.md](paper_claim_evidence_matrix_20260326.md)
- [tutorial_L3_shared_region_closure_20260330.md](tutorial_L3_shared_region_closure_20260330.md)

---

## 1. 先给一句当前版本的总判断

如果把整个项目从头到现在压成一句话，当前最稳的表述是：

> `TransChromBP` 最终不是写成“发现了 Transformer 特有 shortcut 的论文”，而是写成“把 Transformer 安全接到 ChromBPNet-style bias factorization 框架里，并用 full/debiased 双口径建立 bias-safe 训练与诊断体系”的论文。

这句话里的每一部分，都是被实验历史一步一步逼出来的，不是一开始就想清楚的：

- `Transformer` 是后来通过 matched ablation 才确认“确实有真实收益”。
- `bias-safe` 是在 clean matrix 之后才真正站稳。
- `full/debiased` 从“看起来有意思的诊断”变成“必须保留的方法学部分”。
- `center pool` 是在读出头多轮对照后才变成当前默认。
- `L3 shared-region` 是 strict compare 技术链路真正收口后，才成为当前最强 external evidence。
- no-bias 最终没有翻案，它只变成 supplementary 的边界结果。

---

## 2. 看这份报告时的使用方法

建议按下面两种方式之一来用。

### 2.1 如果你现在主要想看代码

按阶段读，每一阶段只抓三样东西：

1. 当时要解决什么问题
2. 代码/配置改动落在哪几个入口
3. 这些改动最后留下了什么，哪些又被后续实验推翻

### 2.2 如果你现在主要想想论文

按阶段读，但重点抓：

1. 当时的主 claim 是什么
2. 这个 claim 靠什么证据支撑
3. 后来是被加强、保留、降级还是推翻

---

## 3. 项目起点：我们最初到底想做什么

最初问题其实很直接：

1. `ChromBPNet` 的 bias factorization 思路是对的，但 backbone 相对保守。
2. `Transformer` 在长距离建模上看起来更强，理论上可能更适合 base-resolution ATAC profile/count 预测。
3. 但高容量 backbone 也带来一个担心：它会不会更容易绕开 bias 分解，在 signal branch 里“偷回” bias 信息？

因此项目从一开始就同时有两条线：

- **建模线**：能不能做出比旧 ChromBPNet 更强的 backbone。
- **诊断线**：即使指标变好，怎么证明它不是靠 bias leakage 拿到的。

这也是为什么当前代码结构从一开始就不是“只做一个纯 Transformer 模型”，而是保留了大量 ChromBPNet/BPNet 风格的结构语义。

### 代码锚点

| 入口 | 作用 |
|---|---|
| [transchrombp.py](../vendor/transchrombp/transchrombp/models/transchrombp.py) | 当前主模型定义，能看出项目最终落成的是“局部卷积塔 + 可开关 Transformer + bias branch + full/debiased 输出” |
| [bias_branch.py](../vendor/transchrombp/transchrombp/models/bias_branch.py) | ChromBPNet-style bias branch 在当前项目里的模型内化实现 |
| [BPNet与ChromBPNet模型对比.md](../docs/learning/BPNet与ChromBPNet模型对比.md) | 理解最初要保留 bias factorization 的背景 |

---

## 4. 第一阶段：先把真实训练语义和工程地基打稳

在项目一开始，真正阻碍推进的并不只是“大方向是不是要加 Transformer”，还有很多现实问题：

- `pyBigWig`/FASTA 句柄在 dataloader worker 里的复用问题
- nonpeak 采样是否在每个 epoch 重新抽样
- zero-count profile loss 是否在给背景窗口施加不合理监督
- `max_epochs` 怎么设才不像拍脑袋

这些看起来像工程细节，但它们对实验解释影响很大。因为如果训练语义本身不稳，那么后面所有关于 shortcut、bias reliance、architecture 的讨论都会混进噪声。

这一步的特点是：

- 还不是“论文故事”
- 但它决定了哪些实验后面值得信

### 当时做了什么

- DataLoader worker 按 PID 重新打开 FASTA / bigWig 句柄
- train split 的 nonpeak 改成按 epoch 重采样
- zero-count profile 不再被强行监督成均匀分布
- 在训练器里补了 early stopping 能力

### 这一步为什么重要

它把实验从“可能连训练语义都不稳定”推进到“至少每条 run 的比较是有意义的”。

### 代码锚点

| 入口 | 关键点 |
|---|---|
| [train_ddp.py](../vendor/transchrombp/transchrombp/training/train_ddp.py) | 训练主循环、early stopping、validation、best metric |
| [transchrombp_training_adjustments_20260316.md](transchrombp_training_adjustments_20260316.md) | 这一阶段最完整的工程调整记录 |

### 留下了什么

- “训练语义需要显式留档” 这件事后来一直保留下来。
- 之后所有主线结果更倾向于通过独立报告收口，而不是只看命令行日志。

### 没有在这一步解决什么

- 还没回答 Transformer 是否真的有用
- 还没回答 bias reliance 是不是 Transformer 特有
- 还没收敛 paper story

---

## 5. 第二阶段：V1 -> V2 升级，产生了第一次强叙事

这一步是整个项目里最重要、也最容易误导人的阶段。

当时我们同时引入了几件变化：

- `profile_bias_stop_gradient=true`
- `debiased_profile_weight=2.0`
- bias profile pooling 的 v2 语义

然后观察到一个非常显眼的现象：

- debiased JSD 从很差的水平一下掉到接近 `0.316`

于是项目早期自然形成了一个非常强、也非常诱人的叙事：

> “我们发现并修复了 Transformer 特有的 Profile Shortcut。”

这条叙事当时不是完全空穴来风，因为现象确实大，而且直觉上也说得通：

- Transformer 容量更高
- bias branch 又是显式存在的
- 如果 signal branch 可以借路 bias，确实可能出现“看起来 full 很好，debiased 很差”的情况

### 但问题也埋在这里

这一步是**多变量同时改动**：

- `sg`
- `deb2`
- v2 的若干结构/配置变化

所以这个阶段最大的问题是：**现象很强，但归因不干净**。

### 代码锚点

| 入口 | 为什么看它 |
|---|---|
| [transchrombp_teacher_v2.yaml](../vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2.yaml) | teacher v2 的核心模型语义 |
| [train_tutorial_teacher_v2_main.yaml](../vendor/transchrombp/transchrombp/configs/train/train_tutorial_teacher_v2_main.yaml) | 早期 v2 主训练入口 |
| [transchrombp.py](../vendor/transchrombp/transchrombp/models/transchrombp.py) | `profile_bias_stop_gradient`、bias branch 融合、full/debiased 输出 |

### 这一步留下了什么

- full/debiased 双口径从这时开始变成项目真正的中心诊断口径
- “单看 full 指标不够” 这条意识后来被保留了

### 这一步后来被推翻了什么

- “Transformer 特有 shortcut 已被证明” 这个强断言后来被 clean matrix 明确推翻

---

## 6. 第三阶段：先把 backbone 价值单独问清楚

在早期 shortcut 叙事出现之后，一个更基础的问题其实必须先回答：

> 即使先不谈 shortcut，Transformer 到底有没有带来真实收益？

如果这个问题都答不稳，那么后面关于 bias-safe、diagnostics、external comparison 的故事都会失去基础。

因此才有了 matched ablation 主线：

- `V2-full`
- `V2-noTF`

这里的关键不是只关掉一个开关，而是尽量保持其它因素不动，只让 `sequence_encoder.enabled` 成为 backbone 层面的核心差异。

### 关键结果

三 seed 汇总后：

- `V2-full`: peak profile `0.31419±0.00012`，peak count `0.83676±0.00531`
- `V2-noTF`: peak profile `0.32393±0.00021`，peak count `0.82503±0.00489`

这一步的意义非常大，因为它把一件事稳住了：

> Transformer 的收益是真实存在的，不是只靠某种巧合拿来的。

这也是后来 paper story 为什么还能保留“Transformer 有价值”这一主结论的原因。

### 代码锚点

| 入口 | 关键点 |
|---|---|
| [transchrombp.py](../vendor/transchrombp/transchrombp/models/transchrombp.py) | `sequence_encoder.enabled` 决定是否启用 `SequenceTransformerEncoder` |
| [transformer_encoder.py](../vendor/transchrombp/transchrombp/models/transformer_encoder.py) | 当前 Transformer 编码器本体 |
| [ablation_no_transformer.yaml](../vendor/transchrombp/transchrombp/configs/model/ablations/ablation_no_transformer.yaml) | noTF 对照模型 |
| [train_ablation_v2_main_profile_select.yaml](../vendor/transchrombp/transchrombp/configs/train/train_ablation_v2_main_profile_select.yaml) | 后来 paper-facing 的 ablation 选模口径 |
| [paper_metric_source_table_20260326.csv](assets/paper_metric_source_table_20260326.csv) | 当前 backbone 主数字来源 |

### 这一步留下了什么

- “Transformer 提升真实表征能力”成为当前仍可强写的主 claim
- `noTF` 从“弱基线”变成理解 architecture / bias-safe 关系的关键对照

### 这一步没有回答什么

- 没有回答 Transformer 是否更容易 shortcut
- 没有回答最终默认 readout 是什么
- 没有回答 external comparison 应该怎样讲

---

## 7. 第四阶段：读出头设计，为什么最后是 center pool

当 backbone 价值基本稳住之后，项目进入了一个更细、但对最终交付非常重要的阶段：

> count 头应该怎么聚合？

这里的核心分支后来演化成：

- `B = center pool`
- `F = attention pool`
- `G = profile refine`

早期如果只看某些单 seed 结果，`F` 并没有立即显得很差，所以这条线一度不是一眼就能收口的。

真正把它收住的是 matched seed 的 held-out 结果：

- `B_s42`: `0.3147 / 0.8503`
- `B_s1234`: `0.3145 / 0.8488`
- `F_s1234`: `0.4269 / 0.8457`

因此当前的结论不是“B 看起来不错”，而是：

> `center pool` 是目前唯一跨 seed 稳定、且 paper-facing 最容易解释的默认 readout。

### 代码锚点

| 入口 | 关键点 |
|---|---|
| [transchrombp.py](../vendor/transchrombp/transchrombp/models/transchrombp.py) | `count_pool_mode` 支持 `full/center/attention` |
| [transchrombp_teacher_v2_center_pool.yaml](../vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool.yaml) | 当前默认 paper-facing 模型配置 |
| [v2fix_center_pool.yaml](../vendor/transchrombp/transchrombp/configs/model/v2fix_center_pool.yaml) | center pool 线的模型配置入口 |
| [v2fix_attn_pool.yaml](../vendor/transchrombp/transchrombp/configs/model/v2fix_attn_pool.yaml) | attention pool 对照入口 |
| [train_v2fix_cpool_s1234_6000_single.yaml](../vendor/transchrombp/transchrombp/configs/train/train_v2fix_cpool_s1234_6000_single.yaml) | center pool 第二个 seed 代表训练配置 |
| [v2fix_and_6000_realdata_followup_20260325.md](v2fix_and_6000_realdata_followup_20260325.md) | `B/F` 收口后的代表报告 |

### 这一步留下了什么

- `corrected B = center pool + sg=true + deb2` 成为当前最终默认模型
- 论文里 readout 线已经不该再写成 `B/F` 并列候选

---

## 8. 第五阶段：Genos 线为什么最后成了高质量负结果

项目中期还开过一条很自然、但后来被降级的支线：

> 能不能把 `Genos-1.2B` 这类 foundation model 的信息融进来，进一步提升 profile/count 任务？

这条线一开始有合理性：

- foundation model 在 sequence-level classification 上确实常常有效
- 当前 backbone 已经不弱，如果再加外部 summary，也许能补更多全局表征

后来这条线没有被粗暴丢掉，而是做出了比较完整的一组负结果：

- online 融合 (`G1/G2`) 没有带来稳定收益
- cached-fusion (`P2`) 甚至在 count 上明显塌陷
- probe 表明 Genos 不是纯噪声，但它和当前主干特征强烈重叠，且粒度不对

因此现在对这条线的判断不是“环境坏了”或“实现没跑通”，而是：

> `global_mean` 这类 sequence-level summary 与当前 base-resolution profile/count 任务粒度不匹配；在强 baseline 前提下，它的边际增益小到不足以越过噪声门槛。

### 代码锚点

| 入口 | 关键点 |
|---|---|
| [genos_adapter.py](../vendor/transchrombp/transchrombp/models/genos_adapter.py) | Genos 相关模块入口 |
| [run_genos_pilot.sh](../vendor/transchrombp/transchrombp/scripts/run_genos_pilot.sh) | online Genos 入口 |
| [run_genos_cached_pilot.sh](../vendor/transchrombp/transchrombp/scripts/run_genos_cached_pilot.sh) | cached fusion 入口 |
| [build_genos_summary_cache.py](../vendor/transchrombp/transchrombp/scripts/build_genos_summary_cache.py) | summary cache 构建 |
| [genos_no_positive_gain_analysis_20260323.md](genos_no_positive_gain_analysis_20260323.md) | 当前最完整的失败原因分析 |

### 这一步留下了什么

- 一组可以直接写进 discussion / appendix 的高质量负结果
- “task granularity mismatch” 成为 foundation model 线的当前解释

### 这一步没有留下什么

- 没有留下当前主线中的正向模型组件
- 没有进入 paper 主表

---

## 9. 第六阶段：clean matrix 让整个项目认知发生翻转

如果说前面几步还在“逐渐变强”，那么 clean matrix 是真正改变整个项目写法的一步。

为什么必须做 clean matrix？

因为在经历了 V2 强现象、Transformer ablation、readout 线之后，一个根本矛盾仍然存在：

- 我们手上已经有很多结果
- 但关于“shortcut 到底是不是 Transformer 特有”的解释仍然不干净

因此才有了后来的 clean matrix：

- 架构轴：`TF` vs `noTF`
- 安全配置轴：`safe` vs `unsafe`

当前最关键的几条 test 结果是：

| 条件 | gap（full - debiased） |
|---|---:|
| `A = TF + sg=false + deb2` | `0.00147` |
| `C = TF + sg=true + deb0` | `0.00972 / 0.00997` |
| `noTF + sg=false + deb2` | `0.00487` |
| `noTF + sg=true + deb0` | `0.02682 / 0.01679` |
| `corrected B = TF + center pool + sg=true + deb2` | `0.00268` |

这组结果把很多旧判断直接推翻了。

### clean matrix 真正教会了什么

1. `full/debiased gap` 确实是有用的诊断，不是偶然噪声。
2. `debiased_profile_weight=2.0` 比单独 `stop-gradient` 更像主效应。
3. 当前没有证据支持 “Transformer 特有 shortcut”。
4. conv-only unsafe 条件下的风险甚至更高。
5. 最终默认模型 `corrected B` 仍然是干净的。

这一步之后，项目的主线叙事被迫改变：

- 从“强 shortcut 机制发现”
- 变成“bias-safe 训练/诊断框架”

### 代码锚点

| 入口 | 关键点 |
|---|---|
| [train_ddp.py](../vendor/transchrombp/transchrombp/training/train_ddp.py) | `compute_losses()` 里 `debiased_profile_weight` / `count_weight` / debiased 分支监督语义 |
| [train_ddp.py](../vendor/transchrombp/transchrombp/training/train_ddp.py) | `run_validation()` 如何同时计算 full/debiased 指标 |
| [train_ddp.py](../vendor/transchrombp/transchrombp/training/train_ddp.py) | `resolve_best_metric_name()` 和选模逻辑 |
| [transchrombp.py](../vendor/transchrombp/transchrombp/models/transchrombp.py) | `profile_bias_stop_gradient` 与 full/debiased 融合 |
| [profile_shortcut_revalidation_summary_20260326.md](profile_shortcut_revalidation_summary_20260326.md) | clean matrix 当前最稳的收口报告 |
| [paper_rewrite_strategy_20260326.md](paper_rewrite_strategy_20260326.md) | 论文叙事如何被 clean matrix 改写 |

### 这是整个项目最大的认知翻转

如果只挑一件最重要的事来记，那就是：

> clean matrix 不是“又补了一组实验”，而是它迫使整个项目从机制发现论文，改写成诊断框架论文。

---

## 10. 第七阶段：strict compare 从 L2 走到 L3，外部证据真正收口

clean matrix 之后，内部方法学主线基本稳了。但论文还缺一个问题：

> 对外比较到底有多稳？

这就是 strict compare 线的意义。

它本身也经历了认知升级：

### 10.1 L2 的问题

`Layer 2` 给出的方向已经支持 `TransChromBP` 优于 official，但后来发现它仍有一个难以忽略的混杂：

- candidate region 不完全一致

因此 `L2` 更像 system-level compare，但还不够像 matched-region compare。

### 10.2 为什么必须做 L3

`L3 shared-region` 的目的就是更进一步统一：

- peaks
- nonpeaks
- bigWig

这样 external evidence 才更能站住。

### 10.3 L3 最终收口后的结果

held-out test 上：

- Official controlled L3: `0.33853 / 0.33989 / 0.69958`
- TransChromBP controlled L3: `0.31319 / 0.31547 / 0.84016`

同时分类指标也补齐：

- Official: `0.82295 / 0.83148 / 0.75248`
- TransChromBP: `0.87861 / 0.87530 / 0.80991`

这一步的意义是：

> external comparison 终于可以写成“matched shared-region system compare 仍然成立”。

注意这里仍不能写成 architecture-only attribution，但它已经是当前最强的 external evidence。

### 代码锚点

| 入口 | 关键点 |
|---|---|
| [train_tutorial_corrected_b_strict_compare_L3_6000.yaml](../vendor/transchrombp/transchrombp/configs/train/train_tutorial_corrected_b_strict_compare_L3_6000.yaml) | 当前 L3 训练配置入口 |
| [evaluate_checkpoint.py](../vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py) | held-out test、count/profile/classification 统一评估 |
| [select_best_epoch.py](../vendor/transchrombp/scripts/select_best_epoch.py) | 自研侧外部 selector |
| [select_best_epoch.py](../scripts/paper_aligned_repro/select_best_epoch.py) | official 侧外部 selector |
| [tutorial_L3_shared_region_closure_20260330.md](tutorial_L3_shared_region_closure_20260330.md) | 当前 L3 收口报告 |

### 这一步留下了什么

- 论文主文外部比较最终改成双轨：
  - 历史 tutorial baseline
  - `L3 shared-region`

### 这一步纠正了什么

- 不再把早期 external compare 全都混写成“严格控制变量比较”

---

## 11. 第八阶段：no-bias 补证为什么最后只是 supplementary

在 clean matrix 和 L3 都收口之后，一个自然的追问是：

> 在现在这套 Transformer + center-pool recipe 下，显式 bias factorization 还必要吗？

于是有了 no-bias 单 seed 补证：

- `bias_branch.enabled=false`
- 沿用当前 paper-facing recipe
- 训练、external valid selector、held-out test 全闭环

### 关键结果

- external best valid: `epoch 26`
- held-out test-full peak `mean/median JSD=0.31496/0.31715`
- peak count `r=0.84978`
- count-scale `AUROC/AUPRC/F1=0.88990/0.88418/0.82037`

这组数和当前 corrected-B / center-pool 主配置是同档的。

因此这条线最后没有得出“bias branch 必不可少”的强结论，也没有得出“去掉它反而更强”的强结论。

当前最准确的理解是：

> 在现有 Transformer + center-pool recipe 下，显式 bias factorization 在单 seed 上没有表现出清晰额外收益；但这条边界结果也不能替代 clean matrix 和 `L3 shared-region` 的主证据角色。

### 代码锚点

| 入口 | 关键点 |
|---|---|
| [ablation_no_bias.yaml](../vendor/transchrombp/transchrombp/configs/model/ablations/ablation_no_bias.yaml) | no-bias 模型配置 |
| [train_ablation_v2_main.yaml](../vendor/transchrombp/transchrombp/configs/train/train_ablation_v2_main.yaml) | no-bias 训练配置基座 |
| [select_best_epoch.py](../vendor/transchrombp/scripts/select_best_epoch.py) | no-bias 外部 valid selector |
| [evaluate_checkpoint.py](../vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py) | no-bias held-out test |
| [paper_writeup_flow_20260330.md](paper_writeup_flow_20260330.md) | no-bias 在论文中如何降级为 supplementary |

### 这一步留下了什么

- 一个对 reviewer 友好的 boundary result
- 但没有改写主线

---

## 12. 到今天为止，哪些结论是稳定的

### 12.1 可以保留为当前主结论的

1. Transformer 在当前 factorized ATAC 建模里带来真实收益。
2. full/debiased gap 是必要诊断。
3. `debiased profile supervision` 是当前最关键的 bias-safe 训练语义之一。
4. `center pool` 是当前默认 readout。
5. `corrected B` 是当前 paper-facing 默认模型。
6. `L3 shared-region` 证明 external comparison 在更严格条件下仍成立。

### 12.2 已经被推翻或必须降级的

1. “Transformer 特有 shortcut”
2. “纯卷积不会出现 shortcut”
3. “stop-gradient 单独修复了问题”
4. “显式 bias factorization 是当前 recipe 的绝对必要组件”

### 12.3 当前仍然只是边界结果的

1. no-bias 单 seed
2. Genos 失败的最终原因
3. architecture sensitivity 的更普适结论

---

## 13. 当前论文为什么会写成现在这个样子

如果只从论文角度看，整个项目大致经历了三次改写。

### 第一次版本

- 主角：Transformer 特有 shortcut
- 叙事：发现问题 -> 修复问题

### 第二次版本

- 主角：Transformer 收益 + bias reliance 风险
- 叙事：模型有效，但机制不那么确定

### 当前版本

- 主角：bias-safe Transformer framework
- 叙事：
  1. Transformer 带来真实收益
  2. full/debiased gap 是必要诊断
  3. clean matrix 说明 `deb2` 比单独 sg 更关键
  4. center pool 是最终默认 readout
  5. `L3 shared-region` 保住 external evidence
  6. no-bias 只是 supplementary boundary result

这也是为什么当前论文材料会被拆成：

- [paper_rewrite_strategy_20260326.md](paper_rewrite_strategy_20260326.md)
- [paper_claim_evidence_matrix_20260326.md](paper_claim_evidence_matrix_20260326.md)
- [paper_writeup_flow_20260330.md](paper_writeup_flow_20260330.md)
- [transchrombp_paper_cn_v1.tex](transchrombp_paper_cn_v1.tex)

---

## 14. 如果现在要阶段性读代码，推荐顺序是什么

这是给“想回到代码细节，同时又不想丢掉论文主线”的阅读顺序。

### 第一步：先看模型结构，知道当前默认模型到底是什么

1. [transchrombp.py](../vendor/transchrombp/transchrombp/models/transchrombp.py)
2. [transformer_encoder.py](../vendor/transchrombp/transchrombp/models/transformer_encoder.py)
3. [bias_branch.py](../vendor/transchrombp/transchrombp/models/bias_branch.py)

要重点盯这些问题：

- `ConvStem` 和 `LocalDilatedTower` 为什么还在
- Transformer 是怎么接进去的
- bias branch 怎么和 signal branch 融合
- `full` / `debiased` 是怎么同时输出的
- `count_pool_mode=center` 的实现到底在哪里

### 第二步：看训练语义，知道当前哪些 loss/metric 是真正的主线

1. [train_ddp.py](../vendor/transchrombp/transchrombp/training/train_ddp.py)
2. [train_ablation_v2_main_profile_select.yaml](../vendor/transchrombp/transchrombp/configs/train/train_ablation_v2_main_profile_select.yaml)
3. [transchrombp_teacher_v2_center_pool.yaml](../vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool.yaml)

要重点盯：

- `compute_losses()`
- `run_validation()`
- `resolve_best_metric_name()`
- `apply_training_mode_defaults()`
- `validate_semantics_profile()`

### 第三步：看评估与 selector，知道结果是怎么被选出来的

1. [evaluate_checkpoint.py](../vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py)
2. [select_best_epoch.py](../vendor/transchrombp/scripts/select_best_epoch.py)
3. [select_best_epoch.py](../scripts/paper_aligned_repro/select_best_epoch.py)

要重点盯：

- held-out test 如何同时输出 profile/count/classification
- 自研侧 external selector 的 metric 路径是什么
- official 侧 selector 如何补齐 `classification_metrics`

### 第四步：再回来看报告和论文

1. [paper_claim_evidence_matrix_20260326.md](paper_claim_evidence_matrix_20260326.md)
2. [paper_writeup_flow_20260330.md](paper_writeup_flow_20260330.md)
3. [tutorial_L3_shared_region_closure_20260330.md](tutorial_L3_shared_region_closure_20260330.md)
4. [transchrombp_paper_cn_v1.tex](transchrombp_paper_cn_v1.tex)

---

## 15. 如果现在要阶段性想论文，推荐顺序是什么

### 先问三个问题

1. 当前最强主 claim 是什么？
2. 哪些结论已经不该再写？
3. 哪组结果最适合进入主文，哪组只适合放 supplementary？

### 然后对应去看

| 问题 | 先看哪里 | 再看哪里 |
|---|---|---|
| 当前哪些 claim 还能写 | [paper_claim_evidence_matrix_20260326.md](paper_claim_evidence_matrix_20260326.md) | [paper_metric_source_table_20260326.csv](assets/paper_metric_source_table_20260326.csv) |
| 论文为什么不再写 shortcut 机制发现 | [paper_rewrite_strategy_20260326.md](paper_rewrite_strategy_20260326.md) | [profile_shortcut_revalidation_summary_20260326.md](profile_shortcut_revalidation_summary_20260326.md) |
| 外部比较现在该怎么写 | [tutorial_L3_shared_region_closure_20260330.md](tutorial_L3_shared_region_closure_20260330.md) | [transchrombp_paper_cn_v1.tex](transchrombp_paper_cn_v1.tex) |
| no-bias 该怎么放 | [paper_writeup_flow_20260330.md](paper_writeup_flow_20260330.md) | [TRACKING.md](../TRACKING.md) |

---

## 16. 最后一页记忆卡片

如果以后隔一段时间再回来，只记下面 8 句话就够了：

1. 这个项目的起点是“把 Transformer 接进 ChromBPNet-style bias factorization”。
2. 早期确实看到了很强现象，但当时是多变量一起变，不能干净归因。
3. matched ablation 保住了“Transformer 有真实收益”。
4. `center pool` 是 readout 线最后真正站住的默认。
5. Genos 线是高质量负结果，不是简单没跑通。
6. clean matrix 推翻了“Transformer 特有 shortcut”。
7. `L3 shared-region` 保住了当前最强 external evidence。
8. no-bias 只是一条 supplementary boundary result，不改主线。

这就是为什么今天的论文、代码导读和结论口径，会同时落在：

> bias-safe Transformer framework + dual-metric diagnostics + shared-region external comparison

而不是：

> Transformer shortcut 机制发现论文
