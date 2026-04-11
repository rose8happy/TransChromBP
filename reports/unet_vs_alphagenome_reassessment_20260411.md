# U-Net 与 AlphaGenome 关系复核（2026-04-11）

## 一句话结论

当前仓库里的实验结果不支持把“引入 U-Net”继续表述成一个默认值得推进的近端动作，但这条结论的精确含义必须写窄：

- 现有负向证据足以否掉的是：在 `corrected B` 基座上，只替换 `2114 bp -> 中心 1000 bp` 的 debiased profile readout 为 `U-Net-lite / multiscale-local-skip decoder`，就在当前数据、预算和训练语义下带来 clean gain。
- 现有负向证据不足以否掉的是：AlphaGenome 那种 `1 Mb` 长上下文、全模型 `U-Net-style encoder + transformer tower + decoder`、多模态多任务、并带 teacher/distillation 的更大 sequence-to-track 范式。

换句话说，当前证据否掉的是我们这版 drop-in decoder family，不是泛化意义上的 “U-Net 永远不适合基因组 sequence-to-track”。

## 1. 这份复核回答什么

要回答的不是“AlphaGenome 用了 U-Net，所以我们是不是一定错了”，而是更精确的三个问题：

1. 仓库里到底已经实测过哪些 `U-Net / multiscale / local-skip` 家族变体。
2. 这些实验的负结果，究竟能外推到多大范围。
3. AlphaGenome 的 U-Net 设计和本项目当前负结果，是否真在同一个假设空间里。

## 2. 本项目当前主线与当时为什么会试 U-Net

当前 paper-facing 稳定基座仍是 `corrected B`，即：

- `conv stem + local dilated tower + transformer + ChromBPNet-style bias branch`
- 保留 `full/debiased` 双口径
- `count` 头独立，默认 `center-pool`
- `profile_bias_stop_gradient=true`

这条基座在仓内已有较稳定的“bias-safe + matched external compare”证据链。之后外部路线评审把怀疑重点压到了：

- 不是 backbone 本身是否太弱；
- 而是 `bottleneck representation -> 1 bp profile` 的 readout / decoder 机制是否过于直读式。

因此，仓内后续确实有意识地试了两层更接近 `U-Net-like` 的路线：

1. `6002` 上的 `U-Net-lite v1` cheap-screen。
2. `6000` 上更正式的 `multiscale/local-skip decoder v2` A6000 formal gate。

## 3. 仓内已经实测过什么，以及它们各自说明什么

### 3.1 `6002`: `U-Net-lite v1` 是真实有效的 runtime 实验，不是概念空转

本次复核额外确认了一个重要事实：

- `6002` runtime 仓 `/home/zhengwei/bylw_atac/TransChromBP/src/transchrombp/models/profile_decoder.py` 确实包含 `UNetLiteProfileDecoder`。
- 对应 `transchrombp.py` 也确实注册了 `profile_readout_mode='unet_lite_v1'`。
- `r4` 日志头明确写到：
  - `Model config: .../transchrombp_teacher_v2_center_pool_unet_lite_v1.yaml`
  - `Data config: .../data_tutorial_canonical_v1_6002.yaml`
  - `BATCH_SIZE_PER_GPU: 16`
  - `RUN_NAME: teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r4`

因此，这不是“配置名叫 U-Net-lite，但实际没跑到 U-Net-lite”的伪实验。

### 3.2 `6002`: `U-Net-lite v1` 的有效结果仍然不占优

根据 [reports/unet_lite_v1_rigor_review_20260409.md](./unet_lite_v1_rigor_review_20260409.md)：

- `r1/r2` 都是无效启动，不算负样本。
- `r3/r4` 是有效 run。
- `r4` 的 best peak 约为：
  - `profile_target_jsd_full_mean=0.4484`
  - `count_pearson_full=0.7990`

对照仓内现有锚点：

- `skipprobe_v1_wide` best peak 约 `0.4475 / 0.8024`
- matched `short10_nofoundation_control` peak 约 `0.3193 / 0.8298`

因此，`U-Net-lite v1` 的信息量是：

- 它不是完全跑不起来；
- 但它没有把自己推上 shortlist；
- 也不像“只差一次幸运 seed 就能翻正”的边缘候选。

### 3.3 `6000`: 更正式的 `multiscale/local-skip decoder v2` 也没有给出 clean gain

根据 [reports/multiscale_local_skip_v2_30ep_gate_closeout_20260410.md](./multiscale_local_skip_v2_30ep_gate_closeout_20260410.md)：

- `teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1`
- `6000 / A6000 x2`
- `30 epoch` formal gate
- 唯一问题被定义成：
  - 在 `corrected B` 其余语义不动时，只升级 `profile readout` 为 `multiscale/local-skip decoder`，能否带来 clean gain

结果：

- `best peak profile_target_jsd_full_mean=0.3338`
- `best peak count_pearson_full=0.7948`

对比稳定引用的 `corrected B` two-seed mean：

- `peak JSD ≈ 0.3146`
- `peak count_r ≈ 0.8496`

差值约为：

- peak JSD `+0.019`
- peak count `-0.055`

这不是边缘摇摆，而是 formal gate 明确失败。

### 3.4 因此，仓内对“readout-only 引入 U-Net-like 结构”的当前判读已经不是单点负例

截至本次复核，至少有两层独立但方向一致的内部证据：

1. `6002` 轻量 `U-Net-lite v1` cheap-screen 不晋级。
2. `6000` 更正式的 `multiscale/local-skip decoder v2` formal gate 也失败。

据此做出的 family-level 解释是一个推断：

> 在本项目当前这套 `corrected B + 短上下文 base-resolution ATAC profile/count + bias-safe additive fusion` 设定下，把 U-Net-like 结构作为 profile readout 的局部升级项，并没有显示出 clean gain。

这个推断我认为已经够稳。

## 4. AlphaGenome 的 “U-Net” 到底是什么

这里不能只抓 “用了 U-Net” 这一个词。

根据 AlphaGenome 官方论文《Advancing regulatory variant effect prediction with AlphaGenome》（Nature, 2025）：

- AlphaGenome 输入是 `1 Mb` DNA sequence，而不是 `~2 kb` 局部窗口。
- 其核心架构是一个整模型级别的 `U-Net-style encoder + transformer tower + decoder`。
- Encoder 会把序列逐步下采样到 `128 bp` 分辨率，Transformer 在该尺度上做长程建模，再由 Decoder 逐步上采样回更高分辨率。
- 模型同时预测跨 `11` 个 output types、`5930` 条 human 或 `1128` 条 mouse genome tracks。
- 论文 ablation 明确指出：
  - `1 bp` target resolution 对细粒度任务很重要；
  - `1 Mb` 训练和推理上下文给出最佳整体结果；
  - distillation 对 variant-effect 任务有实际帮助；
  - fully multimodal training 整体优于单模态训练。

论文入口：

- Nature article: <https://www.nature.com/articles/s41586-025-10014-0>
- 代码/权重说明页中提到的研究仓：<https://github.com/google-deepmind/alphagenome_research>

因此，AlphaGenome 的成功并不是“只在一个 ChromBPNet-like 模型尾部插了个轻量 decoder 就变强了”，而是更接近下面这组组合拳：

1. 超长上下文。
2. 全模型 encoder-decoder。
3. 多模态多任务监督。
4. 明确的 teacher/distillation 训练策略。
5. 与任务解析度对齐的输出头设计。

## 5. AlphaGenome 与本项目当前负结果，哪里同、哪里不同

### 5.1 相同点

- 都在乎 base-resolution 或近 base-resolution 的 track 形成。
- 都认为从压缩表示回到高分辨率输出的机制很关键。
- 都会用 skip / multiscale / decoder 结构处理“局部细节 vs 更大上下文”这个张力。

### 5.2 根本不同点

但两者至少在以下关键维度上并不等价：

- 结构放置位置不同
  - 我们测的是 readout-only / profile-branch upgrade。
  - AlphaGenome 用的是整模型级 encoder-decoder 主干。
- 上下文规模不同
  - 我们主任务是 `2114 bp -> 中心 1000 bp`。
  - AlphaGenome 是 `1 Mb`。
- 监督规模不同
  - 我们当前主线更接近单任务或窄任务的 ATAC profile/count。
  - AlphaGenome 是大规模 multimodal track prediction。
- 训练范式不同
  - 我们当前 negative family 主要是 supervised direct training。
  - AlphaGenome 明确使用 teacher/distillation。
- 目标语义不同
  - 我们把 `bias-safe debiased profile + additive fusion + 独立 count head` 当成硬约束。
  - AlphaGenome 并不是围绕 ChromBPNet-style bias decomposition contract 设计的。

## 6. 因而当前最合理的判断边界

综合上面证据，我建议把结论写成下面这句，不再宽也不再窄：

> 当前仓库结果已经足够否掉“在 `corrected B` 上做 readout-only 的 U-Net-like / multiscale decoder 升级，会自然带来 clean gain”这一近端假设；但它还不足以否掉“在更接近 AlphaGenome 的长上下文、全模型 encoder-decoder、teacher/distillation、多模态监督条件下，U-Net-style 结构是有效的”这一更大假设。

这里后半句是推断，不是仓内已验证事实。

## 7. 本次复核新增发现：canonical archive 与 runtime 不完全一致

本次还发现一个复盘严谨性问题：

- `6002` runtime 上的 `src/transchrombp/models/profile_decoder.py` / `transchrombp.py` 包含 `UNetLiteProfileDecoder` 与 `unet_lite_v1` mode。
- 但冻结 worktree `dual-track-20260409` 下对应 `vendor/.../profile_decoder.py` / `transchrombp.py` 没把这段实现完整归档，只保留了 config 名和文档口径。

这不改变 `r3/r4` 的负向结果本身，但会影响后续事后审计：

- runtime 事实是真实跑过 `unet_lite_v1`；
- archival snapshot 目前不足以单独重建这条 runtime 代码路径。

若未来真要重开这条主题，建议先把这部分 runtime delta 以 canonical 方式补回归档，避免再次出现“文档名字在，精确实现不在”的情况。

## 8. 对下一步的建议

当前我不建议做的事：

1. 继续给 `unet_lite_v1` 追加同配方 cheap rerun。
2. 把 AlphaGenome 的正结果口头外推成“我们必须重新押 U-Net”。
3. 把现有负结果写成“U-Net 在基因组任务里不行”。

如果你想因为 AlphaGenome 重新认真看这条方向，我建议只接受一种重开方式：

- 把它定义成显式新 hypothesis，而不是当前 family 的续跑。

这个新 hypothesis 至少应满足下面几条中的 `2-3` 条，而不是只做 readout patch：

1. 更长上下文，而不是继续停留在 `2114 bp`。
2. 更接近全模型 encoder-decoder，而不是只换末端 profile readout。
3. 明确引入 teacher/distillation 或更强的 sequence-to-track supervision。
4. 单独定义新的 success criteria，不能复用 `unet_lite_v1` 的 cheap-screen 口径。

在当前仓库优先级下，更稳妥的动作不是立刻再开跑，而是先决定：

- 我们是否真的要测试“更接近 AlphaGenome 范式”的一个最小切片；
- 还是承认当前项目主线仍不该为 readout-family 再继续烧预算。
