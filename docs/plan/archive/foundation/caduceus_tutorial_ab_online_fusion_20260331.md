# Caduceus-PS Tutorial Matched A/B 实施计划

## 目标

用一次最小但高信息量的 matched A/B，回答：

> `Genos` 失败，到底是 foundation model 对 ATAC 没用，还是 `causal + coarse summary + 不自然注入` 这条具体路线不行？

首轮只跑 tutorial canonical，不直接跳 GM12878。

## 实验定义

- 数据线：`configs/data/data_tutorial_canonical_v1.yaml`
- A 臂：当前 corrected-B 原始 recipe
  - model: `configs/model/transchrombp_teacher_v2_center_pool.yaml`
  - train: `configs/train/train_tutorial_corrected_b_strict_compare_6000.yaml`
- B 臂：在 `conv stem / local tower` 后、现有 transformer 前接 `Caduceus-PS` token hidden states
  - online extraction
  - `Caduceus-PS` 冻结
  - 只训练 projection/fusion + 原主模型参数
  - 不做 cache
  - 不把 foundation feature 直接打进 count head

## 必须保持不变的比较语义

- `2-rank DDP`
- `batch_size_per_gpu=16`
- global batch `32`
- `peak_max_jitter=500`
- `train_revcomp=true`
- `debiased_profile_weight=2.0`
- `count_pool_mode=center`
- checkpoint selector: `peak.profile_target_jsd_full_mean`

## 代码改动范围

- `vendor/transchrombp/transchrombp/models/caduceus_adapter.py`
  - 新增 `CaduceusFeatureExtractor`
  - 新增 `CaduceusTokenAdapter`
- `vendor/transchrombp/transchrombp/models/transchrombp.py`
  - 新增 `caduceus_feat` 输入
  - 新增 `caduceus_branch` 配置与 pre-transformer token fusion
- `vendor/transchrombp/transchrombp/training/train_ddp.py`
  - 新增 Caduceus online runtime builder
  - 训练 / validation / 日志支持 `caduceus_feat`
  - 保持 `genos_cached` 语义不变
- `vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py`
  - 评估时支持 online Caduceus runtime
- `vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_caduceus_ps.yaml`
  - 以 corrected-B model config 为基线，只加 `caduceus_branch`
- `vendor/transchrombp/transchrombp/scripts/run_caduceus_tutorial_ab.sh`
  - 负责生成 runtime config、注入 bias ckpt、启动 A/B

## 实现约束

- `Caduceus-PS` 是 RC-equivariant；不做 `forward + RC` 平均。
- adapter 采用 residual/gated fusion，零初始化 projection，保证训练起点接近 baseline。
- 若 `genos_branch.enabled` 与 `caduceus_branch.enabled` 同时为 `true`，直接报错。
- 若 `genos_cached.enabled=true` 且任一 online foundation runtime 启用，直接报错。

## Gate

- 安全 gate：B 的 held-out peak `count_r` 相对 A 不允许明显塌陷；若下降超过 `0.02`，直接判失败。
- 推进 gate：在满足安全 gate 的前提下，B 的 held-out peak `mean_jsd` 至少优于 A `0.002`，否则不扩 GM12878。

## 并行资源准备

- 第一优先级：
  - 下载 `Caduceus-PS` 仓库与权重
  - 对 `NT-v2-500M` 做最小加载验证
- 第二优先级：
  - `SegmentNT`
  - `gReLU` model zoo 元数据与 `borzoi-model` / `human-atac-catlas-model`
  - `AlphaGenome` access
  - `HyenaDNA` 中等体量 checkpoint

## 产出要求

- `TRACKING.md` 记录 live 状态
- 训练启动后回写 run name / GPU / 日志路径 / 预计结束时间
- 结果出来后写 `reports/` 报告，不只在聊天里总结

## 2026-03-31 实施进展

- 本地实现已完成，并已同步到 6000 `TransChromBP`：
  - `src/transchrombp/models/caduceus_adapter.py`
  - `src/transchrombp/models/transchrombp.py`
  - `src/transchrombp/models/__init__.py`
  - `src/transchrombp/training/train_ddp.py`
  - `src/transchrombp/evaluation/evaluate_checkpoint.py`
  - `configs/model/transchrombp_teacher_v2_center_pool_caduceus_ps.yaml`
  - `scripts/run_caduceus_tutorial_ab.sh`
  - `scripts/prepare_caduceus_ps_assets.sh`
- 远端已通过：
  - `python3 -m py_compile`（上述 Python 文件）
  - `bash -n scripts/run_caduceus_tutorial_ab.sh`
  - `bash -n scripts/prepare_caduceus_ps_assets.sh`
- `run_caduceus_tutorial_ab.sh` 已修成适配 6000 的真实执行环境：
  - 使用 `VENV_DIR/bin/python` 与 `VENV_DIR/bin/torchrun`
  - `PYTHONPATH` 指向 `ROOT_DIR/src`
- Caduceus 后台准备任务已启动：
  - run name: `caduceus_assets_prep_20260331_205729`
  - log: `/data1/zhoujiazhen/bylw_atac/logs/caduceus_assets_prep_20260331_205729.log`
  - 目标：创建独立 venv、拉官方 repo、下载 `Caduceus-PS`、做最小加载 smoke
- `2026-03-31 21:34 CST` 后补充结论：
  - 6000 侧 `caduceus-ps` 独立环境已经建好，`torchrun`、`causal-conv1d`、`mamba-ssm` 都已安装完成
  - 官方仓库 `caduceus_repo` 已成功 clone 到 `/data1/zhoujiazhen/bylw_atac/foundation_models/caduceus/caduceus_repo`
  - 当前真正 blocker 不是依赖，而是 6000 无法访问 `huggingface.co`；`snapshot_download(...)` 直接报 `Network is unreachable`
  - 因此下一步不再是“重装环境”，而是从可访问 Hugging Face 的机器下载 model snapshot，再同步到 6000 指定目录
- `2026-03-31 21:40 CST` 后补充结论：
  - 用户给出的“本机解析 URL、远端直接下载”路线可行
  - 已在本机解析出 `model.safetensors` 的签名直链，目标 host 为 `cas-bridge.xethub.hf.co`
  - 已在 6000 对该签名 URL 做 `curl -I` 验证，返回 `HTTP/1.1 200 OK`
  - 因此大文件不需要中继；剩余缺口只是把小 metadata（`config.json` / `tokenizer_config.json` / `special_tokens_map.json`）与官方 repo 中的 Python 文件一起组装到模型目录
  - 已在本地把 `prepare_caduceus_ps_assets.sh` 扩展为支持 `CADUCEUS_SAFETENSORS_URL`：当提供签名 URL 时，远端会直接下载 `model.safetensors`，并用官方 repo + 内嵌 metadata 组装出可加载的模型目录
  - 当前未完全落地的唯一原因是 6000 的 SSH 会话间歇性 `Connection reset/timed out during banner exchange`，导致新脚本写回与执行需要在会话稳定时重试
- `2026-03-31 22:29 CST` 后补充进展：
  - 已绕过“先同步脚本再执行”的脆弱路径，直接通过 SSH here-doc 在 6000 上生成并启动临时脚本 `/tmp/caduceus_assets_direct_20260331_222941.sh`
  - 新后台任务：
    - run name: `caduceus_assets_direct_20260331_222941`
    - log: `/data1/zhoujiazhen/bylw_atac/logs/caduceus_assets_direct_20260331_222941.log`
    - 路线：`cas-bridge.xethub.hf.co` 签名 URL 远端直下 `model.safetensors` + 官方 repo Python 文件复制 + 内嵌 metadata 组装 + smoke load
  - 当前日志已确认：
    - `model.safetensors.part` 正在持续增长
    - `MODEL_DIR=/data1/zhoujiazhen/bylw_atac/foundation_models/caduceus/caduceus-ps_seqlen-131k_d_model-256_n_layer-16` 中已经落下 `config.json`、`configuration_caduceus.py`、`modeling_caduceus.py`、`modeling_rcps.py`、`tokenization_caduceus.py`、`special_tokens_map.json`、`tokenizer_config.json`
  - 当前下一步：
    - 等待同一任务完成下载并执行 smoke
    - 若 smoke 成功，立刻转 `bash scripts/run_caduceus_tutorial_ab.sh --smoke`
- `2026-03-31 22:44 CST` 后补充进展：
  - 直链下载任务 `caduceus_assets_direct_20260331_222941` 已完成大文件落盘，`MODEL_DIR` 中的 `model.safetensors` 已存在
  - 原任务自带 smoke 失败，但根因已经定位清楚：
    - 远端 `caduceus-ps` venv 初始装到的是 `numpy 2.2.6`，与 `torch 2.2.0+cu121` 组合会持续触发兼容性警告
    - 更关键的是原 smoke 在 CPU 上做 forward，而 `mamba_ssm` Triton layernorm 不支持 CPU，直接报 `ValueError: Cannot find backend for cpu`
  - 已在 6000 直接执行修复并验证：
    - 将 `caduceus-ps` venv 的 `numpy` 降到 `<2`
    - 用 `CUDA_VISIBLE_DEVICES=0` 实测 `AutoTokenizer.from_pretrained(...)` + `AutoModelForMaskedLM.from_pretrained(...).cuda()`，forward 已通过
  - 同一轮发现新的实现约束：
    - `Caduceus-PS` 暴露出的 token hidden states 实际宽度是 `512`，不是 `config.json` 里直观看到的 `d_model=256`
    - 原因是 RCPS 主干把 forward / reverse 两路状态拼接到暴露给下游的 token states 里
  - 因此已同步修补本地 staging 与 6000 `TransChromBP`：
    - `caduceus_adapter.py` 现在按 `rcps=true => effective hidden_size = 2 * d_model` 推断 extractor 输出维度
    - `train_ddp.py` 新增 `caduceus_branch.hidden_size` 与 extractor 实测 hidden size 的一致性校验，避免静默 shape mismatch
    - `transchrombp_teacher_v2_center_pool_caduceus_ps.yaml` 已改成 `hidden_size: 512`
    - `prepare_caduceus_ps_assets.sh` 已改成固定 `numpy<2`，且 smoke 优先走 CUDA；只有无 CUDA 时才退化为 load-only smoke
  - 当前下一步：
    - 立即起 `bash scripts/run_caduceus_tutorial_ab.sh --smoke --run-label smoke_20260331_2244`
    - 若 A/B smoke 通过，再转完整 matched A/B
- `2026-03-31 22:53 CST` 后补充进展：
  - smoke 共实际跑了三轮，失败链和修复链现在已经闭环：
    - `smoke_20260331_2245`
      - launcher log: `/data1/zhoujiazhen/bylw_atac/logs/caduceus_tutorial_ab_smoke_20260331_2245.launch.log`
      - 结果：A 臂在 import data runtime 阶段失败，报 `ModuleNotFoundError: No module named 'pyBigWig'`
      - 修复：确认 `run_caduceus_tutorial_ab.sh` 是用 `caduceus-ps` venv 跑主训练进程，因此这个 venv 不能只装 foundation 依赖，还要装最小训练依赖；随后补装 `pyBigWig` / `pyfaidx`，并把 `prepare_caduceus_ps_assets.sh` 一并更新为安装这些包
    - `smoke_20260331_2252`
      - launcher log: `/data1/zhoujiazhen/bylw_atac/logs/caduceus_tutorial_ab_smoke_20260331_2252.launch.log`
      - 结果：A 臂已经进入训练初始化，但在 `GradScaler` 构建处失败，报 `AttributeError: module 'torch.amp' has no attribute 'GradScaler'`
      - 修复：确认 `caduceus-ps` venv 是 `torch 2.2.0+cu121`，该版本仍应使用 `torch.cuda.amp.GradScaler`；因此已在 `train_ddp.py` 加入兼容 helper：优先用新 API，若不存在则回退到 `torch.cuda.amp.GradScaler`
    - `smoke_20260331_2256`
      - launcher log: `/data1/zhoujiazhen/bylw_atac/logs/caduceus_tutorial_ab_smoke_20260331_2256.launch.log`
      - A log: `/data1/zhoujiazhen/bylw_atac/logs/caduceus_tutorial_A_corrected_b_s42_smoke_20260331_2256.log`
      - B log: `/data1/zhoujiazhen/bylw_atac/logs/caduceus_tutorial_B_caduceus_ps_s42_smoke_20260331_2256.log`
      - 结果：A/B 两臂都已真实跑完 2-step DDP smoke，checkpoint 与 `run_meta.json` 已写出
  - 成功 smoke 的关键证据：
    - A 臂已正常完成训练、validation、checkpoint 写出
    - B 臂日志已明确打印：
      - `Loaded Caduceus extractor: layer=-1 hidden_size=512 rcps=True rc_equivariant=true`
      - `t_fm=0.9147s`，说明 online foundation extraction 已真实进入 step timing
      - `caduceus_gate_mean=0.1311`、`caduceus_delta_rms=3.218e-06`，说明 gated residual 融合分支已经接通且仍保持接近 baseline 的零初始化状态
  - 基于 smoke 通过，正式 matched A/B 已启动：
    - run label: `matched_20260331_2253`
    - launcher log: `/data1/zhoujiazhen/bylw_atac/logs/caduceus_tutorial_ab_matched_20260331_2253.launch.log`
    - A log: `/data1/zhoujiazhen/bylw_atac/logs/caduceus_tutorial_A_corrected_b_s42_matched_20260331_2253.log`
    - B log: `/data1/zhoujiazhen/bylw_atac/logs/caduceus_tutorial_B_caduceus_ps_s42_matched_20260331_2253.log`
    - 机器 / GPU: `6000`, `GPU 0 + 1`, `2-rank DDP`
    - 预计结束时间：`2026-04-02 12:00 CST`
  - 当前下一步：
    - 先监控 A 臂前几个 epoch 是否稳定推进
    - 待 A/B 全部结束后，按 `peak.profile_target_jsd_full_mean` 选 ckpt，再进入 held-out `test-full`
- `NT-v2-500M` 最小加载验证已完成，但当前不能进入候选对照：
  - `AutoModel.from_pretrained(...)` 直接报 `Linear weight [8192, 1024]` vs config `intermediate_size=4096` 的 size mismatch
  - 说明当前目录不是“文件齐了就可直接用”，后续要么重拉官方目录，要么逐项校正 config/权重组合

- `2026-04-03 16:31 CST` held-out `test-full` 收口结果：
  - A held-out：`checkpoint_epoch=44`，peak `mean_jsd=0.3132558744`、`median_jsd=0.3152509630`、`count_r=0.8388160470`
  - B held-out：`checkpoint_epoch=43`，peak `mean_jsd=0.3126602104`、`median_jsd=0.3145868480`、`count_r=0.8466657043`
  - overall 分类（logcount）侧：A `AUROC/AUPRC/F1=0.8867019489/0.8895893005/0.8215747358`，B `0.8868361809/0.8852656056/0.8235101164`
  - 对照结论：B 没有触发安全 gate，`count_r` 还提升了约 `+0.00785`；但 held-out peak `mean_jsd` 仅比 A 改善约 `0.00060`，明显低于预设的 `0.002` 推进门槛
  - 因此这条 matched A/B 的最终判断应写成：`frozen Caduceus-PS online token fusion` 在 tutorial canonical 单 seed 上是 near-null / marginal positive，而不是足以支持扩线的明确正结果
  - 后续动作：不直接扩到 GM12878；若要继续 foundation-model 主线，应优先改变注入方式 / 模型族 / 训练目标，而不是简单复用当前 recipe 放大数据规模
