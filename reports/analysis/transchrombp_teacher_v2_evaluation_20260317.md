# TransChromBP Teacher V2 训练与评估摘要

本摘要整理自 `teacher_v2_20260317` 这条线的阶段性评估结论，用于替代根目录散落的临时说明文本。

## 核心结论

- Teacher V2 不只是分数提升，更关键的是修复了早期 teacher 线里的 profile-dimension bias shortcut。
- 在 `epoch_024` 上，`full` 与 `debiased` 的 peak profile JSD 基本重合，说明 signal branch 已能独立学习 profile 形状。
- 相比早期 Bias->Main 和 tutorial baseline，Teacher V2 在 held-out peak 的 `count_r` 与 `JSD` 都明显更优。

## Held-out `test-full` 概览

| 模型 | Overall count r | Overall JSD | Peak count r | Peak JSD | Nonpeak count r | Nonpeak JSD |
|---|---:|---:|---:|---:|---:|---:|
| Teacher V2 `epoch_024` | 0.8623 | 0.4354 | 0.8465 | 0.3141 | 0.2526 | 0.5567 |
| Teacher V2 `best.pt` | 0.8618 | 0.4367 | 0.8466 | 0.3162 | 0.2563 | 0.5573 |
| Bias->Main Best | 0.8416 | 0.4455 | 0.8207 | 0.3310 | 0.2276 | 0.5599 |
| Baseline `epoch_20` | 0.8253 | 0.4536 | 0.8048 | 0.3437 | 0.1875 | 0.5635 |

## 外部锚点对比

### 对比 ChromBPNet tutorial baseline

| 模型 | Peak count r | Peak mean JSD | Peak median JSD |
|---|---:|---:|---:|
| Teacher V2 `epoch_024` | 0.8465 | 0.3141 | 0.3164 |
| ChromBPNet official tutorial | 0.7274 | — | 0.3360 |

### 对比 AlphaGenome 4 位点 pilot

| 模型 | 4 位点平均 profile JSD | 峰位点平均 profile JSD | 平均绝对计数误差 |
|---|---:|---:|---:|
| Teacher V2 `epoch_024` | 0.1643 | 0.0677 | 818.4 |
| AlphaGenome | 0.1687 | 0.0729 | 1355.4 |

## 当前判断

- 当前论文与长期归档更适合把 `epoch_024` 视为 Teacher V2 的代表 checkpoint。
- 后续若还要继续沿 Teacher V2 深挖，优先改 checkpoint 选择准则，而不是回退 architecture。
- 相关代码快照已正式归档到 `vendor/transchrombp/`，不再以根目录临时文本作为代码入口。
