# Profile Shortcut 复核结论（2026-03-26）

## 1. 本轮要回答的问题

围绕 draft 中的 “Profile Shortcut” 叙事，本轮最关键的 3 个问题是：

1. 我们是否真的复现出了明显的 profile shortcut？
2. `debiased_profile_weight=2.0` 和 `profile_bias_stop_gradient=true`，到底谁更关键？
3. “纯卷积不会出现 shortcut / 这是 Transformer 特有问题” 这类表述还能不能保留？

本报告只使用本轮 clean matrix 与 paper-facing final model 的同口径结果：

- `A = TF + sg=false + deb2`
- `C = TF + sg=true + deb0`（s42 / s1234）
- `noTF + sg=false + deb2`
- `noTF + sg=true + deb0`
- `corrected B = TF + center pool + sg=true + deb2`

其中 TF 系列来自 6000，`noTF` 与 `corrected B` 来自 6002；clean matrix 四格与 `corrected B` 现在都已经有 `split=test` 结果，且 C 线与 conv-only unsafe 线都已补齐第二个 seed。

## 2. 判读门槛

沿用 2026-03-25 的复核计划：

- `profile_full_debiased_jsd <= 0.005`
  视为“没有看到实质 shortcut”
- `0.005 < profile_full_debiased_jsd <= 0.02`
  视为“轻度 bias reliance，需要弱化结论”
- `profile_full_debiased_jsd > 0.02`
  视为“明显 shortcut”

辅助参考：

- `effective_profile_scale`
- `profile_bias_rms_over_signal_rms`
- `peak.profile_target_jsd_full_mean / debiased_mean`
- `count_pearson_full / debiased`

## 3. Clean Matrix Test 结果

| Run | 条件 | peak JSD full | peak JSD debiased | `profile_full_debiased_jsd` | `effective_profile_scale` | `profile_bias_rms_over_signal_rms` | peak count `r` full/deb |
|---|---|---:|---:|---:|---:|---:|---:|
| A | `TF + sg=false + deb2` | 0.31410 | 0.31416 | 0.00147 | 0.00977 | 0.00247 | 0.8506 / 0.8488 |
| C (s42) | `TF + sg=true + deb0` | 0.31720 | 0.31779 | 0.00972 | 0.06418 | 0.01668 | 0.8395 / 0.8370 |
| C (s1234) | `TF + sg=true + deb0` | 0.31739 | 0.31795 | 0.00997 | 0.06590 | 0.01737 | 0.8340 / 0.8326 |
| noTF_deb2 | `noTF + sg=false + deb2` | 0.43531 | 0.43557 | 0.00487 | 0.03156 | 0.01375 | 0.7959 / 0.7966 |
| notf_sg1_deb0 (s42) | `noTF + sg=true + deb0` | 0.43659 | 0.43835 | 0.02682 | 0.17193 | 0.07647 | 0.7919 / 0.7956 |
| notf_sg1_deb0 (s1234) | `noTF + sg=true + deb0` | 0.43731 | 0.43841 | 0.01679 | 0.10793 | 0.04826 | 0.8287 / 0.8267 |
| corrected B | `TF + center pool + sg=true + deb2` | 0.42544 | 0.42556 | 0.00268 | 0.01761 | 0.00717 | 0.8439 / 0.8420 |

对应文件：

- [A test](/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/profile_shortcut_20260325_tf_sg0_deb2_s42_6000_best_test.json)
- [C s42 test](/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/profile_shortcut_20260325_tf_sg1_deb0_s42_6000_best_test.json)
- [C s1234 test](/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/profile_shortcut_20260326_tf_sg1_deb0_s1234_6000_best_test.json)
- [noTF_deb2 test](/home/zhengwei/bylw_atac/TransChromBP/outputs/metrics/profile_shortcut_20260325_notf_sg0_deb2_s42_6002_best_test.json)
- [notf_sg1_deb0 s42 test](/home/zhengwei/bylw_atac/TransChromBP/outputs/metrics/profile_shortcut_20260326_notf_sg1_deb0_s42_6002_best_test.json)
- [notf_sg1_deb0 s1234 test](/home/zhengwei/bylw_atac/TransChromBP/outputs/metrics/profile_shortcut_20260326_notf_sg1_deb0_s1234_6002_best_test.json)
- [corrected B test](/home/zhengwei/bylw_atac/TransChromBP/outputs/metrics/profile_shortcut_20260325_tf_center_sg1_deb2_s42_6002_retry1_best_test.json)

## 4. 结论

### 4.1 有没有复现出“明显的 Profile Shortcut”？

没有复现出“Transformer 特有且稳定的强 shortcut”。

按预先定义的门槛，在当前矩阵各条件里：

- A：`0.00147`，没有明显 shortcut
- `noTF + sg=false + deb2`：`0.00487`，仍在“无实质 shortcut”区间内
- `notf_sg1_deb0`：`0.02682 / 0.01679`，两 seed 都落在全矩阵最高档，但只有 s42 超过 `0.02`
- corrected B：`0.00268`，没有明显 shortcut
- C（双 seed）：`0.00972 / 0.00997`，只有轻度且稳定的 bias reliance

因此，本轮复核不支持把 draft 里的现象继续写成：

- “已经稳健复现的明显 shortcut”
- “Transformer 下会系统性出现强 shortcut”

更准确的说法是：

- 当前 clean matrix 中，风险最高的仍是 `noTF + sg=true + deb0`，而不是 TF 条件；
- 但它的 matched 第二个 seed 降到 `0.01679`，说明 conv-only unsafe 的风险大小存在 seed 级波动，而不是“同幅稳定的强 shortcut”；
- `TF + sg=true + deb0` 只表现为轻度且稳定的 bias reliance 增强；
- 因而原来的 “Transformer 特有强 shortcut” 叙事被当前矩阵直接反证。

### 4.2 `deb2` 和 `sg` 谁更关键？

当前证据明显更支持：`deb2` 比 `sg` 更关键。

理由很直接：

- 关掉 `sg` 但保留 `deb2` 的 A，仍然非常干净：`0.00147`
- 保留 `sg` 但关掉 `deb2` 的 C，双 seed gap 上升到 `0.00972 / 0.00997`
- 同时 C 的 `effective_profile_scale` 和 `profile_bias_rms_over_signal_rms` 也同步抬升到
  `0.06418 / 0.01668` 与 `0.06590 / 0.01737`
  明显高于 A 的
  `0.00977 / 0.00247`

因此当前最稳的结论不是：

- “stop-gradient 单独修复了 shortcut”

而是：

- `debiased profile supervision` 看起来是主效应；
- `stop-gradient` 更像辅助稳态/补强项，而不是唯一关键开关。

### 4.3 “纯卷积不会出现 shortcut / Transformer 特有” 还能不能保留？

不能保留原来的强断言。

现在能稳妥说的是：

- 在 `deb2` 存在时，`TF` 和 `noTF` 都可以维持很小的 full/debiased gap；
- 在 `sg=true + deb0` 的 unsafe 配置下，`noTF` 的 gap 为 `0.02682 / 0.01679`，两 seed 都高于 TF 的 `0.00972 / 0.00997`；
- 因而“只要用了 Transformer 就会明显 shortcut”不仅不成立，而且与当前 matched matrix 相矛盾。

所以论文里 architecture-specific 的结论必须进一步收紧为：

- 当前证据**反对** “Transformer 特有” 这类表述；
- 更合理的解释是：bias reliance 的主导因素是 supervision 设计是否安全，而不是 backbone 是否带 Transformer。

### 4.4 corrected B 是否支持“最终最好模型也没有明显 shortcut”？

支持。

`corrected B` 的 test：

- peak `JSD_full/debiased = 0.42544 / 0.42556`
- `profile_full_debiased_jsd = 0.00268`
- `effective_profile_scale = 0.01761`

这说明即使换到 paper-facing 的 center-pool 最终模型口径，当前也没有看到明显 shortcut。

因此：

- “最终默认模型也没有明显 profile shortcut”
  这句可以保留。

## 5. 论文改写建议

### 5.1 应删除/弱化的句子

下面这些写法不再成立：

- “Profile Shortcut 已被明确证明”
- “纯卷积架构中不会出现 Profile Shortcut”
- “stop-gradient 单独修复了该问题”
- “这是 Transformer 特有现象”

### 5.2 推荐替代表述

可以改成：

> 当前受控复核表明，显式的 debiased profile supervision 是抑制 bias reliance 的关键因素；在其存在时，无论 Transformer scaffold 还是 conv-only scaffold，full/debiased gap 都保持在较低水平。相较之下，仅保留 stop-gradient 而去除 debiased profile supervision 会带来可见但仍属轻度的 bias reliance 增强。

如果想再短一点，可以写成：

> 我们没有观察到稳定而强烈的 “Profile Shortcut” 复现；更准确的结论是，去偏 profile 监督缺失时会增加 bias reliance，而这一现象目前不足以支持 “Transformer 特有” 的强断言。

## 6. 还要不要追加实验？

### 6.1 对当前 paper 结论：不必再为 shortcut 主线强行追加

如果接受上面的降级表述，现在已经够写 paper：

- clean matrix 四格 test 已齐，且 C 线与 conv-only unsafe 线都补齐了第二个 seed
- 最终默认模型 `corrected B` 也已经验证
- 现有结论已经能支持“弱化 shortcut 叙事 + 保留最终模型稳健性”

### 6.2 如果仍要继续补，优先级应转向别的问题

现在更值得补的已经不是 clean matrix，而是：

- tutorial `Layer 2` 的 checkpoint / selector follow-up
  当前更大的开放问题是 late-epoch instability 为什么会把 profile-only selector 推到 count 已塌的 checkpoint
- `Layer 3 = shared-region compare`
  如果要把 system-level 主结果进一步往归因层推进，这条比继续加 shortcut seed 更有信息量

### 6.3 当前不推荐的追加项

- `corrected F`
  对 “Profile Shortcut 是否成立” 不是关键路径
- `tf_sg1_deb0` 更多 seed
  两个 seed 的 test 已经足够说明这是一条轻度且稳定的 C 线现象
- 多 seed 全矩阵 repeat
  只有在你想把 `deb2 > sg` 写成更强因果断言时才值得

## 7. 最终判断

当前最稳的 paper 级结论是：

1. 我们没有复现出 “Transformer 特有的强 `Profile Shortcut`”；当前矩阵里风险最高的 gap 仍出现在 conv-only 的 `notf_sg1_deb0`，但其幅度存在 seed 波动（`0.02682 / 0.01679`）。
2. `debiased profile supervision` 比 `stop-gradient` 更像主效应。
3. “Transformer 特有 / 纯卷积不会出现” 这类强断言应删除，并改写为 supervision-driven 的 bias-safe 叙事。
4. 最终默认模型口径 `corrected B` 也没有显示出明显 shortcut，可以继续作为主模型。
