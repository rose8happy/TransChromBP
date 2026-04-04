# 论文收口写作顺序（2026-03-30）

## 1. 当前主线怎么讲

这轮论文不再围绕“Transformer-specific shortcut”组织，而是按下面四层证据写：

1. `Transformer + bias factorization` 能带来真实收益  
   主证据：`V2-full vs V2-noTF` matched ablation。
2. `full/debiased gap` 是必要诊断  
   主证据：clean matrix（A/C/noTF/corrected B）。
3. `corrected B = center pool` 是当前 paper-facing 默认模型  
   主证据：`B/F` readout 对照与双 seed held-out。
4. 对外比较在 `L3` shared-region 下仍成立  
   主证据：official vs TransChromBP 的 matched system compare。

## 2. 主文表格怎么排

主文外部比较采用双轨结构，而不是继续只保留旧 tutorial baseline：

- `Panel A`：历史 tutorial baseline  
  作用：与公开可复用的 ChromBPNet 数字建立延续关系。
- `Panel B`：`L3` shared-region system compare  
  作用：提供当前更强、matched-region 的 external comparison。

分类指标单独做 supplementary 风格附表：

- `AUROC/AUPRC/F1` 只作为辅助比较
- 不进入主表
- 不单独扩写成新的主结果节

## 3. 当前不该再写成什么

- 不写成 “Transformer 特有 shortcut”
- 不写成 “纯卷积不会出现 shortcut”
- 不把 `L3` 写成 architecture-only attribution
- 不把 no-bias 消融写成新的主结论来源

## 4. no-bias 补证怎么定位

`bias_branch.enabled=false` 的单 seed run 只回答一个问题：

> 在已有 Transformer + center-pool recipe 下，显式 bias factorization 是否仍必要？

这条结果的定位是：

- supplementary 的 component-necessity check
- 不替代 clean matrix
- 不替代 `L3` shared-region 主表
- 若结果很弱或训练不稳，也只写成负对照/边界说明，不拖住主稿

实际结果（`2026-03-31`）：

- 外部 valid selector 选中 `epoch 26`
- held-out `test-full` peak `mean/median JSD=0.31496/0.31715`
- held-out peak `count_r=0.84978`
- count-scale 分类 `AUROC/AUPRC/F1=0.88990/0.88418/0.82037`

因此这条补证的最终口径应固定为：

- 单 seed 下 no-bias 与当前 corrected-B / center-pool 主配置基本同档
- 不支持把“显式 bias factorization 是绝对必要组件”写成强结论
- 继续保持 supplementary 的 component-necessity / boundary-result 定位

## 5. 本轮写作顺序

1. 先收中文主稿
   - 双轨主表
   - `L3` 分类附表
   - `fig:arch`
   - methods 参数补齐
2. 再同步英文主稿的对应段落
3. 最后只把 no-bias 结果补进 supplementary 相关位置

一句话口径：

> 这篇论文现在要写成 “bias-safe Transformer framework + dual-metric diagnostics + matched shared-region external comparison”，而不是 “shortcut 机制发现论文”。
