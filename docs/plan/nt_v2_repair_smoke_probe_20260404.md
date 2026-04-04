# NT v2 Repair + Smoke + Probe 执行计划（2026-04-04）

## 目标

在不直接开 tutorial A/B 训练的前提下，用一轮最小但高信息量的执行回答两个问题：

1. `InstaDeepAI/nucleotide-transformer-v2-500m-multi-species` 能否在 6000 上被稳定、本地、可重复地加载并消费 `2114bp` tutorial 窗口？
2. 一个 `bidirectional + masked-LM + token-level hidden state` 的基因组模型，是否比 `Genos` 更像对当前 ATAC 任务有用的外部信号？

本轮不是论文包装动作，而是研究路线筛选动作。

## 已知前提

- `Genos cached P0/P1 restart` 已于 `2026-04-04 05:07 CST` 正常结束，并确认 `P1` 相对新 `P0` 为负结果：
  - `P0`: peak `JSD=0.3348129319`, `count_r=0.8075212978`
  - `P1`: peak `JSD=0.3352303161`, `count_r=0.8031642888`
- `Caduceus-PS tutorial matched A/B` 已于 `2026-04-03 16:31 CST` 收口为 near-null / marginal positive，不足以直接扩到 `GM12878`。
- `2026-04-04 20:52 CST` 复查时，6000 两张 `A6000` 都空闲。
- NT v2 目录当前已包含：
  - `model.safetensors`
  - `pytorch_model.bin`
  - `jax_model/pytree_ckpt.joblib`
- 当前 blocker 已从“缺大文件”转为“HF/PyTorch 加载链路 + 最小 smoke + 实际 probe”。

## 执行范围

### 1. 环境与运行方式

- 机器：`6000`
- GPU：默认 `GPU0` 单卡
- 新环境：`/data1/zhoujiazhen/bylw_atac/.mamba/envs/nucleotide-transformer-py311`
- 基线环境：`/data1/zhoujiazhen/bylw_atac/.mamba/envs/transchrombp`
- 本轮不占满双卡，不启动 DDP，不动 `6002`

### 2. 新增脚本

- `scripts/nt_v2_smoke.py`
  - 负责：
    - 本地加载 tokenizer/model
    - 两条短 DNA 序列 smoke
    - 一条真实 `2114bp` tutorial window smoke
    - 失败时打印 config / 权重 shape 线索
- `scripts/nt_v2_probe.py`
  - 子命令 `extract`
    - 在 `py311 + NT` 环境中运行
    - 输入 tutorial canonical `valid` 的 peak/nonpeak 区域
    - 抽样提取 NT hidden-state summary
    - 同时落盘序列、标签、真值 total count / logcount
  - 子命令 `analyze`
    - 在 `transchrombp` 环境中运行
    - 加载 `corrected-B / center-pool / seed42` baseline
    - 计算 `encoded_only` / `nt_only` / `concat(encoded, nt)` 指标

### 3. Probe 数据与对照

- 数据范围：tutorial canonical `valid`
- 默认采样：`500 peak + 500 nonpeak`
- NT 特征默认覆盖：
  - `global_mean`
  - `bins4_mean`
- 层位默认取 `auto_quartiles`：
  - `1/4`
  - `1/2`
  - `3/4`
  - `last`
- baseline 固定：
  - checkpoint: `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/checkpoints/v2fix_20260320_cpool_s42/best.pt`
  - model config: `transchrombp_teacher_v2_center_pool`
  - data config: `data_tutorial_canonical_v1`

## 判断门槛

### 技术门槛

- `AutoTokenizer` / `AutoModelForMaskedLM` 本地加载通过
- 两条短序列与真实 `2114bp` window 都能产出 hidden states
- 真实窗口 smoke 记录：
  - token 数
  - hidden shape
  - 显存
  - 耗时

### 研究门槛

只有同时满足“可加载 + probe 不弱于 Genos 参考”才进入下一轮 training integration 设计。

优先看三个条件：

1. 至少一个 `layer × feature` 组合的 `nt_only` AUC 明显高于 Genos 参考 `0.6997`
   - 默认希望 `>= 0.73`
2. `concat(encoded, nt)` 不再像 Genos 那样显著拖累 `encoded_only`
   - 默认要求 `concat_auc >= encoded_only_auc - 0.005`
3. count 相关 probe 至少有一项优于 Genos 参考
   - `R² > 0.0167` 或 `Pearson r > 0.1987`

若三条都达不到，则本轮直接按“暂不值得进入 tutorial A/B training integration”收口。

## 输出要求

- `TRACKING.md`：
  - 开始前写 live 状态
  - 启动后台任务后写 run name / GPU / log / ETA
  - 有新结论后再次回写
- `reports/nt_v2_probe_20260404.md`
  - 记录：
    - 加载是否通过
    - smoke 结果
    - probe 指标
    - 是否过 gate
- 若后台任务启动：
  - 日志统一写到 `/data1/zhoujiazhen/bylw_atac/logs/`
  - 结果目录统一写到 `outputs/nt_v2_probe/`

## 明确不做的事

- 不直接做 NT v2 tutorial A/B 训练
- 不在本轮里尝试 base-resolution 直接 token 对齐注入
- 不把 JAX 路径作为本轮 blocker
- 不把 AlphaGenome pilot 和 NT probe 混在一轮
