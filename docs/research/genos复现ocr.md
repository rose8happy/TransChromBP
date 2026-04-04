# Genos OCR 官方基准复现研究摘要

更新日期：2026-03-23

## 1. 本文用途

这份文档只回答一个问题：

`Genos` 官方公开材料里，和“染色质可及性”最接近、最适合作为仓库内旁线 sanity check 的任务到底是什么？它和当前 `TransChromBP` 主线是什么关系？

结论先写在前面：

- 官方最直接的“可及性”任务不是 ATAC 轨迹回归，而是 `human_ocr_ensembl` 的 OCR 二分类基准。
- 官方 benchmark code 的标准流程是：
  `DNA sequence -> hidden_states[layer] -> attention-mask mean pooling -> MLP/RF/XGB`
- 这条 OCR benchmark 旁线适合做：
  - `Genos-1.2B` 环境、权重、tokenizer、hidden-state 提取链路的 sanity check
  - `layer` 选择是否大体合理的旁证
- 这条 OCR benchmark 旁线不适合做：
  - 当前 `cached-fusion` 主线的 go/no-go gate
  - 对 `profile/count` 轨迹预测收益的直接判定
  - 对 “online full-seq gated fusion” 失败原因的替代解释

## 2. 已核对的官方事实

### 2.1 官方 OCR 任务是什么

官方 benchmark 配置和数据集说明里，`human_ocr_ensembl` 被明确列为支持的数据集，任务类型为二分类：

- 数据集：`human_ocr_ensembl`
- 样本数：`174,756`
- 长度范围：`71-593 bp`
- split：`train/test`
- 任务类型：`classification`

从任务定义看，它表示“开放染色质区域识别”，更准确地说是 OCR 分类，不是 ATAC-seq 连续信号预测。

### 2.2 官方 benchmark code 的真实流程

我核对了官方 `Technical_Notes/benchmarks-code` 的实现，关键链路如下：

1. 输入 JSONL，字段最小为 `{"seq": "...", "label": 0}`
2. `AutoTokenizer.from_pretrained(..., trust_remote_code=True)`
3. `AutoModel.from_pretrained(..., device_map="auto", torch_dtype=torch.bfloat16, trust_remote_code=True)`
4. 前向时打开 `output_hidden_states=True`
5. 取 `outputs.hidden_states[layer]`
6. 用 `attention_mask` 做均值池化
7. 把 pooled embedding 落成 `.pt`
8. 用 `MLP`、`RF` 或 `XGB` 做下游分类

其中默认 `MLP` 是三层隐藏层的轻量分类器，官方默认超参是：

- `mlp_dropout: 0.2`
- `mlp_lr: 1e-4`
- `mlp_epochs: 100`

### 2.3 论文里和 OCR 最相关的公开指标

Genos 论文 Table 2 里，`human_ocr_ensembl` 的 AUC 是：

- `Genos-1.2B: 0.7569`
- `Genos-10B: 0.7623`

这里可以作为“量级对账”的参考值，但不能把“必须 ±0.02 内对齐”当硬门槛，因为论文主文没有把 benchmark code 的所有实现细节逐项锁死。

### 2.4 官方 README 里的 RNA 覆盖度案例是什么

官方 README 另有一个 `RNA Coverage Track Prediction` 应用案例。它的特点是：

- 输入窗口：`32 kb`
- 任务：单碱基 RNA-seq 连续信号回归
- 训练方式：`Genos-1.2B` 全参数 fine-tuning
- 头部：三层 1D conv，末端 `Softplus`
- 优化：`Adafactor + cosine + warmup`

这个案例可以证明官方确实给过“序列到连续轨迹”的配方示例，但它是 RNA-seq，不是官方 ATAC benchmark。

## 3. 对本仓库意味着什么

### 3.1 可以回答的问题

把 `human_ocr_ensembl` 跑通之后，最多可以回答：

- 当前本地 `Genos-1.2B` 权重、tokenizer、`trust_remote_code`、bf16 推理链路是否正常
- `layer 6` 或 `layer 12` 这类候选层，在官方 OCR 分类任务上是否处于合理区间
- `mask mean pooling` 的 summary 表征，是否至少能在官方 OCR benchmark 上给出正常量级的分类性能

### 3.2 不能回答的问题

它不能直接回答：

- `TransChromBP` 当前 `profile/count` 目标上，`cached summary` 是否有真实增益
- 为什么 `online G1/G2` 在 held-out test 上塌陷
- `Genos` 的逐位置特征是否能被当前 `TransChromBP` 主干稳定利用

原因很直接：这两个问题的监督目标、数据语义、模型头、评测指标都不一样。

### 3.3 与当前主线的关系

截至 2026-03-23，仓库里已经有更直接的本地证据：

- `online G1/G2` 的 held-out test 明显落后于 `G0 baseline`
- 当前主线已经切到 `cached-fusion`
- 后续判断应该继续围绕 `summary cache -> probe -> P0/P1/P2`

因此，OCR benchmark 在本仓库里的定位应当是：

- `旁线 sanity check`
- `不改变 cached-fusion 主线优先级`
- `不作为主线 gate`

## 4. 建议的落地边界

### 4.1 当前建议做的

只做一条轻量旁线：

- 模型只用 `Genos-1.2B`
- 先跑 `human_ocr_ensembl`
- 先看 `layer 6` 和 `layer 12`
- 默认分类器先只跑 `MLP`
- 先做 mini smoke，再做全量

### 4.2 当前建议不做的

暂时不做：

- `Genos-10B`
- 上来就全层扫描 `0..12`
- 重新从零搭一套和现有环境平行的新主环境
- 把 OCR benchmark 结果写成 `cached-fusion` 的硬 gate
- 把 RNA-seq 的 full fine-tuning recipe 直接升级成新的 “ATAC 主线”

## 5. 仓库内最小复现链路

建议执行顺序：

1. 复用当前已经验证过的本地 `Genos-1.2B`
2. 拉官方 `benchmarks-code`
3. 用 `genomic-benchmarks` 下载 `human_ocr_ensembl`
4. 用仓库内脚本把目录结构转成 benchmark code 需要的 JSONL
5. 先跑 mini smoke
6. 再跑 `layer 6/12` 的全量
7. 只有前面都干净，才考虑补全层扫描

对应脚本见：

- `scripts/prepare_genos_ocr_jsonl.py`

## 6. 成功与停止条件

### 6.1 可以认为“跑通了”的条件

满足以下条件即可认为这条旁线有效：

- benchmark code 能稳定加载本地 `Genos-1.2B`
- `human_ocr_ensembl` 的 mini smoke 能完整产出 `.pt` 和 `.tsv`
- 全量 `layer 6` 或 `layer 12` 的 `roc_auc` 处于和官方值同量级的区间
  - 实用判断：大致落在 `0.73-0.78`
- 结果没有出现明显反常
  - 例如 `roc_auc < 0.65`
  - 或 train/test 标签映射不一致
  - 或 `trust_remote_code` / bf16 / deterministic 造成流程不稳定

### 6.2 需要立刻停止排查的情况

出现以下情况先停，不要继续补大实验：

- mini smoke 都过不去
- 全量 `layer 6/12` 都明显低于合理区间
- 结果强依赖某个临时手工修补，而不是稳定配置
- 旁线开始挤占 `cached-fusion` 主线的 GPU/工程窗口

## 7. 对执行计划的直接要求

后续执行计划必须满足：

- 默认复用已验证的本地 `Genos-1.2B`
- 数据转换不能再留 `pass`
- 明确把 OCR benchmark 写成 `sidecar sanity check`
- 明确写出“不作为 cached-fusion gate”
- 如果要记录标签语义，必须显式落盘 `label_map.json`

## 8. 官方来源

- Genos 官方 benchmark readme  
  `https://github.com/BGI-HangzhouAI/Genos/blob/main/Technical_Notes/benchmarks-code/readme.md`
- Genos 官方 benchmark config  
  `https://github.com/BGI-HangzhouAI/Genos/blob/main/Technical_Notes/benchmarks-code/config.yaml`
- Genos 官方 benchmark 数据集定义  
  `https://github.com/BGI-HangzhouAI/Genos/blob/main/Technical_Notes/benchmarks-code/benchmarks/datasets_info.yaml`
- Genos 论文  
  `https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giaf132/8296738`
- Genos 官方 README 中的 RNA coverage case  
  `https://github.com/BGI-HangzhouAI/Genos/blob/main/README_en.md`

## 9. 本仓库内相关文件

- 当前 Genos 主线状态：`TRACKING.md`
- 当前 online / cached-fusion 判断：`reports/genos_and_6002_run_analysis_20260323.md`
- 当前 cached-fusion 主方案：`docs/plan/genos_cached_fusion_final_execution_plan_20260323.md`
