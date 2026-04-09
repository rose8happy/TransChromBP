# NT v2 bins16 center-aligned residual gate（2026-04-05）

## 1. 目标

在用户已显式授权继续使用 `6000 / A6000 x2` 的前提下，执行唯一还值得做的 foundation 受控重入：

- `NT v2`
- `layer_07__bins16_mean`
- center-aligned / token-aware residual
- `short10` 双卡训练
- 自动 held-out `test-full` gate

这条任务的目的不是“继续试试看”，而是回答一个更窄的问题：

> 当前被判负的是 `bins4 + coarse-summary residual` family；如果把 residual 改成真正按中心区域对齐的 `bins16` token residual，是否还能保留低扰动前提下的 clean gain 空间？

---

## 2. 为什么这次不能只改 YAML

当前仓库已经支持：

- cached feature `bins16_mean`
- token layout foundation input
- dual-card cache / training / held-out gate launcher

但还**不**支持我们真正想测的新假设，因为当前 residual path 仍有两个旧限制：

1. `train_ddp.py` 会给 residual head 传 `foundation_summary = foundation_tokens.mean(dim=1)`
2. `FoundationResidualHead` 也会把 `encoded` 和 `foundation` 都压成全局均值，再预测 `profile_bin_count` 个 coarse bins

因此，单纯把 `feature_name` 从 `bins4_mean` 改到 `bins16_mean`，并不能构成真正的新假设，只会变成：

> 更细的输入 token + 同样粗的全局 summary residual

这条路线信息增益不够，不值得继续占两张 A6000。

---

## 3. 推荐实现

### 3.1 最小可行新假设

保留以下不变：

- backbone 仍是 corrected-B / center-pool
- foundation 仍是 cached feature，不接管主干
- residual head 仍零初始化
- `short10` 预算仍用于快速 gate

只改两点：

1. residual head 改成 `center_token_bins` 路径
2. launcher/runtime config 真正切到 `layer_07__bins16_mean`

### 3.2 center_token_bins residual 的具体语义

- 输入：
  - `encoded`: 本地 backbone token `[B, L, D]`
  - `foundation_tokens`: NT v2 cached token `[B, 16, H]`
- 对 `encoded`：
  - 先取中心 `output_len=1000` 对齐区域
  - 再自适应池化到 `16` 个局部 bins
- 对每个 bin：
  - 本地 token bin 与对应 foundation token 各自投影
  - 拼接后过共享 `fuse` 层
- profile residual：
  - 不再从全局 fused summary 一次性预测全部 bins
  - 改成每个 token bin 预测一段局部 residual bins，再按顺序拼接
- count residual：
  - 对 token-wise fused bins 做均值池化，再预测一个 `count_delta`

这条实现仍然是“低扰动 residual”，但已经不再是“全局 summary residual”。

### 3.3 本轮推荐配置

- model config：
  - `feature_name: layer_07__bins16_mean`
  - `feature_tokens: 16`
  - `profile_bin_count: 128`
  - `alignment_mode: center_token_bins`
- training：
  - 继续 `train_tutorial_foundation_short10.yaml`
  - 双卡 DDP
  - 自动 held-out `test-full`
- cache：
  - 仅构建当前 run 真正需要的 `layer_07__bins16_mean`

---

## 4. 执行步骤

### Step 1. 代码

- 修改 `vendor/transchrombp/transchrombp/models/foundation_adapter.py`
  - 为 `FoundationResidualHead` 增加 `alignment_mode`
  - 新增 `center_token_bins` 路径
- 修改 `vendor/transchrombp/transchrombp/models/transchrombp.py`
  - 把 `alignment_mode` 等配置传入 residual head
- 修改 `vendor/transchrombp/transchrombp/scripts/run_ntv2_residual_short10.sh`
  - runtime train config 能同步 `foundation_cache_features`
  - cache feature types 默认可由 `feature_name` 推导或显式覆盖
- 新增 model config
  - `vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_ntv2_bins16_residual.yaml`

### Step 2. 测试

- 新增 pytest 单测
  - residual head 在 `center_token_bins` 路径下 shape 正确
  - `feature_tokens=16, profile_bin_count=128` 可正常前向
- 运行
  - `python -m pytest tests/test_ntv2_bins16_residual_head.py -q`
  - `python -m pytest tests/test_foundation_cache_alignment.py -q`
  - `python -m py_compile ...`
  - `bash -n vendor/transchrombp/transchrombp/scripts/run_ntv2_residual_short10.sh`

### Step 3. 发车

- 机器：`6000`
- GPU：`0,1`
- 建议 run name：
  - `ntv2_bins16_centerres_short10_s42_20260405_dual`
- 日志：
  - `/data1/zhoujiazhen/bylw_atac/logs/ntv2_bins16_centerres_short10_20260405_6000.log`

### Step 4. 文档同步

- 启动后立即回写 `TRACKING.md`
- 若实现细节超过状态粒度，再同步 `reports/`

---

## 5. stop/go 规则

这条 run 的解释前提固定为：

- 若结果仍明显落后于 matched `short10 no-foundation control`
  - 说明“更细 token + center-aligned residual”也没有补出 clean gain
  - foundation 线进一步收口
- 若结果接近或超过 matched control
  - 才允许讨论第二步 `bins16 cross-attention`

在结果出来前，不自动扩：

- 多 seed
- 更长训练
- 旧 `bins4 residual`
- 旧 `Genos/Caduceus`

---

## 6. 当前执行记录（2026-04-05 23:32 CST）

本轮实现与启动状态：

- 已新增 `center_token_bins` residual 路径，并在 `tests/test_ntv2_bins16_residual_head.py` 上用 6000 远端环境验证通过（`2 passed`）
- 已新增配置：
  - `configs/model/transchrombp_teacher_v2_center_pool_ntv2_bins16_residual.yaml`
  - `configs/train/train_tutorial_foundation_short10_ntv2_bins16.yaml`
- 已在 6000 后台启动双卡 smoke：
  - run name：`ntv2_bins16_centerres_gate_smoke_20260405`
  - launcher PID：`3983963`
  - 机器 / GPU：`6000 / GPU0,1`
  - 日志：`/data1/zhoujiazhen/bylw_atac/logs/ntv2_bins16_centerres_gate_smoke_20260405_6000.log`
  - 当前阶段：`Step 0` 双卡构建 `train/valid` 的 `bins16_mean` cache

本次 smoke 的唯一目的：

> 先验证 `center_token_bins + bins16` 新路径能否在真实 6000 双卡环境里完整跑通 cache -> short10 train -> held-out gate。

若 smoke 通过，下一步不是重建 cache，而是直接复用：

- cache：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/foundation_cache/ntv2_tutorial_canonical_bins16_v1`
- model config：`transchrombp_teacher_v2_center_pool_ntv2_bins16_residual.yaml`
- train config：`train_tutorial_foundation_short10_ntv2_bins16.yaml`

直接起 full dual-card short10 run。

---

## 7. 进展更新（2026-04-06 01:42 CST）

### 7.1 smoke 结果

`ntv2_bins16_centerres_gate_smoke_20260405` 已完整结束，并成功写出：

- held-out JSON：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/ntv2_bins16_centerres_gate_smoke_20260405_best_test_full_20260406_012209.json`

关键事实：

- 链路层面：通过
  - `train/valid` cache 正常复用 / 构建
  - `2-step` 双卡训练正常
  - `test` cache 正常构建
  - held-out gate 正常写 JSON
- 指标层面：仍为明显负结果
  - overall `count_pearson_full=-0.1451`
  - overall `profile_target_jsd_full_mean=0.6668`
  - peak `count_pearson_full=-0.2968`
  - peak `profile_target_jsd_full_mean=0.6201`

这里的判断固定为：

> smoke 的目的只是验证 `bins16 center_token_bins` 新路径在真实 6000 上能否完整跑通，不把 `global_step=2` 的 held-out 指标当成 stop/go 依据。

### 7.2 full dual-card run 已启动

基于 smoke 已通过链路验证，当前已直接复用已有 `train/valid` bins16 cache 启动 full short10 训练：

- run name：`ntv2_bins16_centerres_short10_s42_20260406_dual`
- 启动时刻：`2026-04-06 01:42 CST`
- launcher PID：`4020074`
- train `torchrun` PID：`4020091`
- 机器 / GPU：`6000 / GPU0,1`
- 日志：`/data1/zhoujiazhen/bylw_atac/logs/ntv2_bins16_centerres_short10_20260406_6000.log`
- 复用 cache：
  - `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/foundation_cache/ntv2_tutorial_canonical_bins16_v1`

### 7.3 当前运行信号

从日志看：

- `Step 0` 已正确跳过，说明 `train/valid` bins16 cache 已被复用
- `Step 2` 双卡训练已经进入稳定 step 输出
- `foundation_residual_profile_rms` 与 `foundation_residual_count_rms` 都在随 step 增长

这说明当前 full run 不再停留在 smoke 的“近零更新”状态，而是在真实学习。

### 7.4 下一步

当前不再做新的实现改动，只做两件事：

1. 等 full short10 训练收口
2. 读取 final held-out JSON，并与 matched `short10_nofoundation_control` 并排对照

这一步的历史判断条件已经兑现：final held-out 明显落后于 matched baseline，因此
`bins16 center-aligned residual` 这条 genuinely new-hypothesis 也已经判负。

## 8. 最终结果（2026-04-06 03:25 CST）

`ntv2_bins16_centerres_short10_s42_20260406_dual` 已完成 full held-out，并写出：

- `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/metrics/ntv2_bins16_centerres_short10_s42_20260406_dual_best_test_full_20260406_032529.json`

与 matched `short10_nofoundation_control_s42_20260405_dual` 的并排对照如下：

- overall
  - `count_pearson_full=0.7727` vs `0.8457`
  - `profile_target_jsd_full_mean=0.4625` vs `0.4395`
  - `profile_pearson_full_mean=0.6580` vs `0.7116`
  - `peak_auroc=0.8082` vs `0.8677`
- peak
  - `count_pearson_full=0.7516` vs `0.8298`
  - `profile_target_jsd_full_mean=0.3588` vs `0.3193`
  - `profile_pearson_full_mean=0.7837` vs `0.8641`

这条 run 的失败语义需要固定为：

1. 它不是因为单纯的 bias leakage 指标爆炸而失败；
2. 而是因为 full held-out 与 matched control 对照时，profile/count 两条主指标都没有守住；
3. 问题落在“互补性没有转化成净增益”，而不是“链路没跑通”。

最终 stop/go 判断固定为：

> 当前 NT v2 residual short10 family 已经完成两道真实 gate：旧 `bins4/coarse-summary residual`
> 判负，新的 `bins16 center-aligned residual` 也仍判负。当前证据不支持再把这条 family
> 自动扩到 `bins16 late cross_attention`、更多 seed 或更长训练。

因此，本文件现在只保留为这次受控重入的完整执行记录；后续若还要继续 foundation，
必须先提出一个不同于当前 `summary / token-fusion / coarse residual / bins16 residual`
family 的新问题，而不是把这里继续当作默认后续链。
