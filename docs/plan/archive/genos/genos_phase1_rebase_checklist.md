# Genos Phase 1 Rebase Checklist

更新日期：2026-03-21

本文档用于约束 Genos Phase 1 接入的代码基线，避免继续沿用旧版 `.codex_remote_edit/TransChromBP` 快照，把已经确认有效的 V2 / bias-safe 语义一起回退掉。

## 1. 代码基线

- 必须基于当前 `vendor/transchrombp/transchrombp` 主线接入 Genos，不得继续在 `.codex_remote_edit/TransChromBP` 上叠改。
- 接入前先确认以下文件与 6000 当前主线一致：
  - `vendor/transchrombp/transchrombp/models/transchrombp.py`
  - `vendor/transchrombp/transchrombp/models/bias_branch.py`
  - `vendor/transchrombp/transchrombp/training/train_ddp.py`
  - `vendor/transchrombp/transchrombp/configs/model/`
  - `vendor/transchrombp/transchrombp/configs/train/`

## 2. 必须保留的现有语义

- 保留 `ConvStem + local_tower + Transformer` 的现有主干，不得退回旧版“无 local tower”结构。
- 保留 `profile_bias_stop_gradient: true` 的默认安全设置。
- 保留当前 `profile_pool_factor` / bias pooling 的实现与配置解析。
- 保留现有 validation 指标口径：
  - `peak.profile_target_jsd_full_mean`
  - `peak.profile_target_jsd_debiased_mean`
  - `peak.count_pearson_full`
  - `peak.count_pearson_debiased`
  - `peak.count_mae_full`
  - `peak.count_mae_debiased`

## 3. matched baseline 要求

- Genos pilot 的 baseline 必须沿用当前 `bias_pretrain -> main` 流程，不允许把 `bias_branch.pretrained_path` 留空后重新训一个 fresh bias branch。
- `best_metric` 必须与当前主线目标一致，优先使用 `peak.profile_target_jsd_full_mean`，不能回退到仅按 `peak.loss_total` 选 best checkpoint。
- Pilot 对照组至少包括：
  - `G0 baseline_v2`：当前主线同口径 baseline
  - `G1 genos_gate`：逐位置 Genos 注入
  - `G2 genos_mean`：全局均值广播对照
- Phase 1 默认调度是：
  - `G0` 先单卡起跑
  - `G1/G2` 再分别占用两张 A6000 并行
  - 不默认启用 `nproc_per_node=2` 的两卡 DDP
- 2026-03-21 的 corrected smoke 已通过；随后 `G1` quick batch sweep（`bs=8/12/16/20`, `dry_run_steps=50`）也已完成：
  - 峰值显存约 `5.4 / 6.4 / 7.2 / 8.1 GiB`
  - 50-step walltime 约 `35 / 47 / 58 / 71 s`
- 当前正式 `G0/G1/G2` pilot 的统一 batch 已锁定为 `batch_size_per_gpu=20`：
  - 不要再让 launcher / train config 静默回退到旧的 `8`
  - 若后续单独评估 DDP，需要按目标 global batch 反推每卡 batch，而不是直接沿用单卡数值
  - 若模型结构、精度或数据口径发生明显变化，再重新做 quick sweep

## 4. Genos 接入本身的最低要求

- `GenosFeatureExtractor` 必须显式设置 `attn_implementation="flash_attention_2"`，与已验证的 Phase 0/最小前向口径保持一致。
- 当前代码不能声称支持 “early exit 节省 75% 计算”，除非真的实现按层截断前向；仅取中间层 hidden state 不等于 early exit。
- 第一阶段默认保持 `frozen + no_grad + eval()`，不把 Genos 主体放进 optimizer / DDP / checkpoint。
- Genos 特征只能进入 signal branch，不进入 bias branch。

## 5. 进入实现前的检查项

- 先完成新的 Phase 0 重跑，并用 hold-out / CV 结果重新确定：
  - 是否有 go/no-go 信号
  - 候选层位
  - batch / 显存预算
- 如果 Phase 0 结果只是边缘可用，不要直接进入 LoRA 或 cross-attention，先收缩为最小 frozen pilot。

## 6. 合入前验证

- `python -m py_compile` 覆盖新增/修改的 Python 文件。
- 至少完成一次小规模 smoke run，确认：
  - Genos 特征提取链路可跑通
  - matched baseline / genos_gate / genos_mean 三组配置都能被解析
  - validation metrics 正常写出
  - gate 统计和 full/debiased gap 能被记录
