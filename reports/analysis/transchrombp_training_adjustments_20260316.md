# TransChromBP 训练代码调整记录（2026-03-16）

这份记录只总结本轮高置信修改和暂不硬改的存疑点，避免后续回看时把“已经改了的”和“只是建议实验验证的”混在一起。

## 已直接修改

### 1. DataLoader worker 重新按 PID 打开 FASTA / bigWig 句柄

- 位置：
  `.codex_remote_edit/TransChromBP/src/transchrombp/data/real_data.py`
- 原因：
  这和你们之前真实训练里已经碰到过的 `pyBigWig` 句柄继承问题是同一类根因；当前本地可见代码里还没回写。
- 当前做法：
  `Fasta` 和 `pyBigWig` 缓存句柄现在都会检查当前 PID；一旦 worker 进程变化，就重新打开。

### 2. train split 的 nonpeak 改为按 epoch 重采样

- 位置：
  `.codex_remote_edit/TransChromBP/src/transchrombp/data/real_data.py`
  `.codex_remote_edit/TransChromBP/src/transchrombp/training/train_ddp.py`
- 原因：
  之前 `nonpeak_ratio < 1` 时只在数据集初始化时抽一次 nonpeak 子集，后续整个训练都复用这批背景窗口。
  这和原始 ChromBPNet 的“每个 epoch 重采样 negatives”不一致，也容易让 low-nonpeak 训练更依赖一次抽样的偶然性。
- 当前做法：
  train split 会保留完整 nonpeak 池，由 `PeakNonpeakResamplingSampler` 在每个 epoch 重新抽取目标数量的 nonpeak，再和全部 peaks 混合。

### 3. zero-count profile loss 不再把背景窗口强行监督成均匀 profile

- 位置：
  `.codex_remote_edit/TransChromBP/src/transchrombp/training/train_ddp.py`
- 原因：
  之前 profile loss 对全零窗口会构造近似均匀 target distribution，这会给背景窗口引入额外、且不太合理的 profile 监督压力。
- 当前做法：
  全零窗口在 profile loss 里记为 0 贡献；非零窗口仍按当前归一化版本计算。

### 4. 增加 early stopping 能力，不再靠硬猜“20 epoch 该不该改”

- 位置：
  `.codex_remote_edit/TransChromBP/src/transchrombp/training/train_ddp.py`
- 新配置：
  `trainer.early_stop_patience`
  `trainer.early_stop_min_delta`
- 说明：
  默认仍为关闭，不强行改变历史配置；但后续可以在保留 `max_epochs` 上限的同时，用 `best_metric` 做提前停止。

## 暂不硬改，先留档

### A. signal backbone 仍然偏弱

当前主干还是“浅 conv stem + Transformer + 轻量 profile/count head”。这条线我仍然认为很可能限制相对 ChromBPNet 的上限，但这已经是结构级修改，不属于这轮“高置信、小步快跑”的范畴。

建议后续优先做一个明确对照：

- `BPNet/ChromBPNet-style local tower + current bias branch`
- `BPNet/ChromBPNet-style local tower + Transformer`

如果这组对照成立，再决定是否把 Transformer 放在 local tower 上方，而不是继续沿“浅 stem 直接接 Transformer”细磨。

### B. profile loss 目前只修了 zero-count 语义，尚未完全对齐原始 ChromBPNet 的 multinomial 标度

这轮没有直接把 profile loss 改成完全按原始 counts 权重的版本，因为那会明显改变 loss 标度，可能连带影响现有 `count_weight` 的经验值和对照口径。

当前建议：

- 先用修过 zero-count 的版本重跑一轮主线；
- 如果结果仍显示 bias 分支长期主导，再专门做“normalized profile loss vs count-weighted multinomial-style loss”的小实验。

### C. hybrid 路径里 `track_total_count_target` 与 auto count-weight 仍然耦合

这件事现在还没改。后续如果还要认真分析 `hybrid_data` 为什么更差，建议优先拆成两个实验：

- 只开 `track_total_count_target`
- 只开 `count_weight_strategy=chrombpnet_auto`

不要再继续用当前“两个变量一起动”的 hybrid 口径直接解释方法优劣。

### D. debiased 分支的训练目标语义仍待明确

当前默认仍只优化 `full` 输出；如果未来要把 `debiased_*_weight` 打开，必须先定义清楚目标语义，否则很容易把 `debiased` 分支重新拉回去拟合 full signal。

### E. `count_weight_strategy=chrombpnet_auto` 不能直接迁移到当前 `learnable_scales` 融合

这件事现在已经不只是“要不要继续试”的建议项，而是需要正式留档的训练动力学结论。

问题不在于 `chrombpnet_auto` 本身绝对错误，而在于它来自 ChromBPNet 的刚性 bias 分解语境：ChromBPNet 的 bias 是固定减法，高 `count_weight` 会逼主模型自己把 counts 学好；而当前 TransChromBP 的 `learnable_scales` 属于柔性加性融合，高 `count_weight` 会给优化器一条更容易的捷径，即放大 bias 相关分支而不是认真优化 debiased 主分支。

这意味着：

- `count_weight_strategy=chrombpnet_auto` 不能再被视为“canonical 数据语义”的自然组成部分；
- `hybrid_data` 路径中观测到的退步，也不能再简单归因于数据处理本身，因为其中混入了 `count_weight × learnable_scales` 这个独立混杂因子；
- 后续凡是使用 `learnable_scales` 的主训练，都应把 `count_weight` 作为单独受控变量，而不是直接继承 ChromBPNet 的经验值。

基于这点，已经追加了一条受控对照：保持 canonical 数据、teacher 架构、`local_tower`、低 `nonpeak_ratio` 等全部不变，只把 `count_weight` 固定回 `0.1` 的 `teacher_canonical_v1_fixcw_20260316`。这条对照的定位不是“最终方案”，而是确认根因：

- 如果它恢复到旧主线那类健康 scale 区间，就能基本确认本次回归的主因是 `count_weight`；
- 如果它仍然异常，再回头检查 `local_tower`、canonical 数据预处理或其他结构改动。

## 关于 epoch 数量

当前我的判断不是“直接把 20 改成更大”或“直接砍成更小”，而是：

- 不建议再单独押注一个固定 epoch 数；
- 更稳的做法是保留上限，再配 `early stopping`。

一个保守的下一轮起点是：

- `max_epochs = 30`
- `trainer.early_stop_patience = 6`
- `trainer.early_stop_min_delta = 0.0`

这样做的理由是：

- 你们当前已有 run 的 `best.pt` 明显早于最后一个 epoch；
- 但在修完 nonpeak 重采样和 zero-count loss 之后，最佳 epoch 位置很可能会移动；
- 先给训练留一点上限，再用 patience 收口，比重新拍脑袋定 12、15 或 24 更稳。

## 已补的 rerun 配置

为了避免直接改动原始 tutorial 基线，这轮另外新增了一份专用配置：

- `.codex_remote_edit/TransChromBP/configs/train/train_tutorial_rerun_earlystop.yaml`

当前参数是：

- `max_epochs = 30`
- `trainer.early_stop_patience = 6`
- `trainer.early_stop_min_delta = 0.0`

它的定位是“修完训练语义后的一轮保守 rerun 起点”，不是新的永久默认值。
