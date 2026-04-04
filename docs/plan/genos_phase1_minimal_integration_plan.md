# Genos Phase 1 Minimal Integration Plan

更新日期：2026-03-21（优先级更新；corrected smoke + G1 batch sweep 已完成）

本文档给出一条基于当前 `vendor/transchrombp/transchrombp` 主线的最小 Genos 接入路线。目标不是一次性把 Genos 线做大，而是在当前有效的 Phase 0 结果基础上，用最小改动回答一个更窄的问题：

`frozen Genos + gated fusion` 在不破坏现有 V2 / bias-safe 语义的前提下，是否能给当前主线带来哪怕很小但可信的净增益？

优先级说明：

- 由于外部研究判断明确建议优先验证“预训练基因组大模型是否有收益”，Genos 线的实验优先级现已提高。
- 这只改变执行顺序，不改变技术门槛：Phase 0 仍然是黄灯，因此 Phase 1 依旧必须是最小、frozen、matched-baseline 的 pilot。

## 1. 当前前提

- Phase 0 已经按新口径重跑完成，旧的 in-sample `AUC=1.0` 结果作废。
- 当前有效结果来自 6000 `genos-1.2b` 环境：
  - `50+50` smoke：最佳层 `layer_3`，`center_cv_auc_mean=0.6580`
  - `500+500` full：最佳层 `layer_6`，`center_cv_auc_mean=0.6719`，`center_holdout_auc=0.6540`
- 结论是 `黄灯`，不是 `绿灯`。
- 资源侧当前可复用结论：
  - Phase 0 独立推理粗估：`batch=8` 双向约 `21.8 GB`，`batch=16` 双向约 `41.0 GB`
  - 但当前 Phase 1 实际训练链路的显存占用明显低于该粗估
- Phase 0 脚本已改成流式统计，host RSS 峰值约 `7.64 GB`；后续不需要再为 Phase 0 的 CPU 内存问题预留额外技术债。
- 工程侧当前已完成：
  - `G0/G1/G2` corrected smoke 全部通过
  - `G1` quick batch sweep（`bs=8/12/16/20`, `dry_run_steps=50`）全部通过
  - sweep 结果：
    - `bs=8`: `35s`, `5432 MiB`
    - `bs=12`: `47s`, `6350 MiB`
    - `bs=16`: `58s`, `7216 MiB`
    - `bs=20`: `71s`, `8090 MiB`
  - 因此当前正式 pilot 的统一 batch 已锁定为 `20`

因此，Phase 1 只能按“最小 frozen pilot”推进，不能按“明显有收益，直接扩展实现范围”的口径推进。

## 2. 接入原则

- 基线只能是当前 `vendor/transchrombp/transchrombp` 主线，不得继续沿用 `.codex_remote_edit/TransChromBP` 旧快照。
- Genos 主体保持：
  - `frozen`
  - `eval()`
  - `torch.no_grad()`
  - 不进入 optimizer
  - 不进入 DDP 包装
  - 不写入 checkpoint
- bias-safe 约束保持不变：
  - 保留 `local_tower`
  - 保留 `profile_bias_stop_gradient: true`
  - 保留 `bias_profile_pool_factor`
  - 保留当前 validation 指标口径
- 第一版只做：
  - `G0 matched baseline`
  - `G1 genos_gate`
  - `G2 genos_mean`
- 第一版默认调度：
  - 不用两卡 DDP
  - `G0` 先单卡跑 matched baseline
  - `G1/G2` 再分配到两张 A6000 上各跑一条单卡 run
- 当前正式 pilot 固定使用 `batch_size_per_gpu=20`
  - 这是 2026-03-21 `G1` quick batch sweep 之后选定的统一值
  - 不再默认使用两卡 DDP 去改变单个 run 的 global batch
- 不做：
  - LoRA
  - cross-attention
  - early exit
  - 长窗口

## 3. 最小代码改动面

### 3.1 新增文件

`vendor/transchrombp/transchrombp/models/genos_adapter.py`

职责拆两块：

- `GenosFeatureExtractor`
  - trainer 侧外部运行的冻结提取器
  - 负责 one-hot -> DNA string -> tokenizer -> fwd + RC -> layer 选择
  - 输出统一为 `[B, L, 1024]`
  - 必须显式使用 `attn_implementation="flash_attention_2"`
- `GenosGatedAdapter`
  - 模型内的可训练小模块
  - 输入：
    - `local_tokens: [B, L, 256]`
    - `genos_feat: [B, L, 1024]`
  - 输出：
    - `fused_tokens: [B, L, 256]`
  - 默认结构：
    - `Linear(1024 -> 256)`
    - `LayerNorm`
    - `gate = sigmoid(Linear(local_tokens) + bias)`
    - `fused = local_tokens + gate * proj(genos_feat)`
  - `gate_bias_init` 建议 `-2.0`

### 3.2 修改文件

`vendor/transchrombp/transchrombp/models/transchrombp.py`

最小改法：

- `__init__` 增加可选 `genos_adapter` 配置解析
- `forward()` 增加可选参数：
  - `genos_feat: Optional[Tensor] = None`
- 融合位置固定在 `local_tower` 之后、Transformer 之前：

```python
conv_feat = self.conv_stem(x)
local_feat = self.local_tower(conv_feat)
tokens = local_feat.transpose(1, 2)

if genos_feat is not None:
    tokens = self.genos_adapter(tokens, genos_feat)

encoded = self.transformer(tokens)
```

约束：

- 没有 `genos_feat` 时，行为必须与当前主线完全一致
- `genos_feat` 绝不能进入 bias branch

`vendor/transchrombp/transchrombp/models/__init__.py`

- 导出新增模块

`vendor/transchrombp/transchrombp/training/train_ddp.py`

最小改法：

- `build_model()` 仍只返回当前主模型
- 新增一层 trainer 侧 runtime 构建逻辑，例如：
  - `build_genos_runtime(model_cfg, device, rank)`
- 在 `run_validation()` 和训练主循环中：
  - 先从 `batch["seq"]` 构建 `genos_feat`
  - 再调用 `outputs = model(seq, genos_feat=genos_feat)`

推荐伪代码：

```python
genos_runtime = build_genos_runtime(model_cfg, dist_env.device, dist_env.rank)

for batch in train_loader:
    seq = batch["seq"].to(dist_env.device, non_blocking=True)
    genos_feat = None
    if genos_runtime is not None:
        with torch.no_grad():
            genos_feat = genos_runtime.extract(seq)
    outputs = model(seq, genos_feat=genos_feat)
```

约束：

- `genos_runtime` 不挂到 DDP model 上
- 每个 rank 各自持有一份 frozen extractor
- 训练/验证都必须走同一条 Genos 特征路径

## 4. one-hot 到 Genos 输入的最小实现

当前数据集返回的 `seq` 来自 `vendor/transchrombp/transchrombp/data/real_data.py`，形状为 `[L, 4]` 的 one-hot，非标准碱基位置是全 0。

因此最小实现可以直接在 `GenosFeatureExtractor` 里做 deterministic 解码：

- `A=[1,0,0,0]`
- `C=[0,1,0,0]`
- `G=[0,0,1,0]`
- `T=[0,0,0,1]`
- 全 0 或异常模式映射到 `N`

不需要先改 dataset，也不需要新增离线 token 缓存。

## 5. 配置设计

### 5.1 model config

建议新增两个配置文件：

- `vendor/transchrombp/transchrombp/configs/model/v2fix_genos_gate.yaml`
- `vendor/transchrombp/transchrombp/configs/model/v2fix_genos_mean.yaml`

二者都从当前 `v2fix_baseline.yaml` 语义出发，新增一个 `genos_branch` 段：

```yaml
genos_branch:
  enabled: true
  model_path: /data1/zhoujiazhen/bylw_atac/foundation_models/genos/Genos-1.2B
  attn_implementation: flash_attention_2
  layer: 6
  hidden_size: 1024
  fusion_mode: gated_add
  pool_mode: none   # G1
  gate_bias_init: -2.0
  freeze: true
```

`G2` 只把 `pool_mode` 改成 `mean`：

- 先对 `[B, L, 1024]` 做全序列均值
- 再广播回 `[B, L, 1024]`

### 5.2 train config

建议新增：

- `vendor/transchrombp/transchrombp/configs/train/train_v2fix_genos_profile_select.yaml`

最小口径：

- `best_metric: peak.profile_target_jsd_full_mean`
- `best_metric_mode: min`
- `early_stop_patience: 8`
- `batch_size_per_gpu: 20`（当前已锁定的正式 pilot 值；如需回退，应显式覆盖）
- `max_epochs: 20`
- 默认是单卡配置；Phase 1 不用两卡 DDP 去放大全局 batch
- 若后续还要改 batch，应先重做 quick sweep，再把选中的 batch 固定用于 `G0/G1/G2`
- 其余尽量继承当前 `train_tutorial_teacher_v2_main_profile_select.yaml` 的稳定设置

## 6. 对照矩阵

- `G0 baseline_v2`
  - 不开 Genos
  - 训练配置与 G1/G2 完全一致
- `G1 genos_gate`
  - `layer=6`
  - `pool_mode=none`
  - 主实验
- `G2 genos_mean`
  - `layer=6`
  - `pool_mode=mean`
  - 判断收益是否只来自粗粒度 global signal

判读优先级：

- `G1 > G0` 且 `G1 > G2`
  - 才能说明逐位置 Genos 特征可能真的有用
- `G1 ≈ G2 > G0`
  - 说明只是多了一个粗粒度 peak indicator
- `G1 ≈ G2 ≈ G0`
  - 说明这条线应尽快收口

## 7. 第一版不该做的事

- 不要先把 Genos 特征接进 bias branch
- 不要先改 checkpoint 结构去保存 Genos 主体
- 不要因为 `layer_6` 是当前最佳层，就把它写成硬编码且不可改
- 不要把“读取中间层 hidden state”写成“early exit”
- 不要把 `G1` 直接和历史 40 epoch fully-trained baseline 硬比

## 8. 实施顺序

1. 先完成 `genos_adapter.py`
2. 再改 `transchrombp.py`，保证 `genos_feat=None` 时行为不变
3. 再改 `train_ddp.py`，只加最小 runtime
4. 配置文件最后加
5. `G0/G1/G2` 的 corrected smoke 已完成
6. `G1` quick batch sweep 已完成，并把正式 pilot 的统一 batch 锁到 `20`
7. 下一步跑 `G0` 单卡 matched baseline
8. 最后把 `G1` / `G2` 分别放到两张 A6000 上并行跑单 seed 20 epoch pilot

## 9. 当前建议（优先级提升后）

在当前有效证据下，最合理的推进方式不是“立即全面并入”，而是“把最小 pilot 提到下一轮最高优先级”：

- 以 `layer_6` 作为当前唯一候选层
- 只做 matched-baseline 的最小 frozen pilot
- 优先验证 `G1 vs G2`，而不是先验证“Genos 能否替代 Transformer”
- 正式 pilot 的 batch 已由 2026-03-21 quick sweep 锁定为 `20`
- 当前已经启动的 `F_s42/G_s42` 不主动中断；它们结束后，A6000 应优先切到 Genos `G0/G1/G2`
- 在 Genos 第一批单 seed pilot 有结论之前：
  - 暂不扩 `A/B` 多 seed
  - 暂不扩 `F/G` 多 seed
  - 暂不启动 LoRA / cross-attention / 长窗口
  - 暂不把“数据扩容 + 模型扩容”作为下一优先级

推荐排程：

1. 当前 run 收口
   - `F_s42` 在 6000 跑完
   - `G_s42` 在 6002 跑完
2. Genos 工程接入与 smoke
   - `genos_adapter.py`
   - `transchrombp.py` 的 `forward(genos_feat=...)`
   - trainer runtime
   - `G0/G1/G2` corrected smoke
   - `G1` quick batch sweep（已完成；正式 batch=`20`）
3. Genos pilot
   - `G0 baseline_v2`：20 epoch，单卡 A6000，耗时仍待以当前 launcher 重新实测
   - `G1 genos_gate`：20 epoch，单卡 A6000，当前 `bs=20` 的 50-step sweep 为 `71s`
   - `G2 genos_mean`：20 epoch，单卡 A6000，预计与 `G1` 同量级
4. 资源排法
   - 先收掉 `G0`
   - 再把两张 A6000 并行给 `G1` / `G2`
   - 不再沿用旧的 `20–24h/run` 估计；应按当前 quick sweep 的 step time 重新预算 wall-clock

当前正式起跑口径：

```bash
cd /data1/zhoujiazhen/bylw_atac/TransChromBP

# 1. 先跑 G0
nohup bash scripts/run_genos_pilot.sh --groups G0 --gpu 0 \
  > /data1/zhoujiazhen/bylw_atac/logs/genos_20260321_G0.launch.log 2>&1 &

# 2. G0 收口后，再并行起 G1 / G2
nohup bash scripts/run_genos_pilot.sh --groups G1 --gpu 0 \
  > /data1/zhoujiazhen/bylw_atac/logs/genos_20260321_G1.launch.log 2>&1 &
nohup bash scripts/run_genos_pilot.sh --groups G2 --gpu 1 \
  > /data1/zhoujiazhen/bylw_atac/logs/genos_20260321_G2.launch.log 2>&1 &
```

- `launcher` 必须显式用 `bash scripts/run_genos_pilot.sh ...` 调用，不要用 `sh`。
- 2026-03-21 已在 6000 上通过 `ssh ... 'bash scripts/run_genos_pilot.sh --groups G0'` 的 `DRY_RUN=1` 验证。

运行时风险说明：

- 当前 launcher 已把正式 pilot 的默认 `BATCH_SIZE_PER_GPU` 改成 `20`，因此上面命令不需要额外传 batch 环境变量。
- `G1` 的 `bs=20` quick sweep 在 stage log 中对应 `step 40/11809 = 50.1s`，按稳定段粗算约 `1.25s/step`。
- 这意味着单个 `G1` run 的训练主体大致在 `~82–93h` 区间；再加上 `validate_every_epochs=5` 带来的 4 次验证开销（当前粗估约 `~1.5–2h`），正式 `20 epoch` run 的总 wall-clock 应按 `~84–95h` 预算。
- 因此，当前真正的阻塞已不再是显存，而是 wall-clock。正式长跑前应明确接受这个量级；如果不可接受，应先下调 `max_epochs`，而不是盲目回退 batch 或改成 DDP。

并行资源说明：

- Genos 执行期间，不要求 `6002` 为 A6000 保留成“评估从属机”。
- 更合理的做法是：`6000` 专心回答 Genos 是否有净增益，`6002` 继续承接 Genos 无关、且不依赖 A6000 checkpoint 同步的独立 TransChromBP 实验。
- 当前最合适的 6002 后续任务是：
  - `GM12878-only` smoke
  - smoke 通过后再跑 `GM12878-only` 单 seed baseline
  - `K562-only` 先做 smoke，不急着直接上完整长跑

## 10. 什么时候降级 Genos

如果出现以下任一情况，Genos 优先级应在首轮 pilot 后立即下降：

- `G1 ≈ G2 ≈ G0`
- `G1` 的 peak 指标没有任何稳定改善，同时 `full/debiased gap` 变大
- gate 长期塌到接近 0，说明模型自动忽略 Genos
- `nonpeak JSD` 明显恶化

这时再回头继续：

- V2fix 多 seed
- 数据扩容 / 迁移实验
- 更大结构搜索
