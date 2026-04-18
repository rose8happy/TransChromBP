# 6000 环境分工与 Genos 接入记录

本文档记录 6000 机器上与 `chrombpnet`、`TransChromBP`、`Genos-1.2B` 相关的实际可用环境，避免环境信息只散落在 `TRACKING.md`、脚本头部和临时日志里。

更新日期：2026-03-21

## 1. 快速结论

- 当前不要把 `chrombpnet`、`transchrombp`、`genos-1.2b` 三套依赖强行合并成一个公共环境。
- `chrombpnet` 环境用于官方 ChromBPNet 工具链与 TensorFlow 侧流程。
- `transchrombp` 环境用于当前 PyTorch 训练、benchmark、评估脚本。
- `genos-1.2b` 环境用于 `Genos-1.2B` 的加载、推理、embedding 提取和可行性验证。
- 若要做 TransChromBP + Genos 融合实验，建议先走“冻结 Genos 特征提取器”路线，而不是先改造现有 `transchrombp` 主环境。

## 2. 环境总览

| 环境 | 路径 | 核心版本 | 主要用途 |
|---|---|---|---|
| `chrombpnet` | `/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet` | `python 3.9.23`, `tensorflow 2.8.0` | 官方 ChromBPNet CLI、MoDISco、`bedGraphToBigWig` 等旧流程工具链 |
| `transchrombp` | `/data1/zhoujiazhen/bylw_atac/.mamba/envs/transchrombp` | `python 3.10.14`, `torch 2.5.1` | 当前 TransChromBP 训练、benchmark、评估 |
| `genos-1.2b` | `/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b` | `python 3.10.12`, `torch 2.7.1+cu126`, `transformers 4.57.6`, `huggingface_hub 0.36.2`, `flash_attn 2.8.3` | Genos-1.2B 本地加载、GPU 推理、特征提取 |

6000 GPU 基础信息：

- `nvidia-smi`: `Driver 550.67`
- GPU: `NVIDIA RTX A6000 49140 MiB x2`

## 3. 推荐使用方式

### 3.1 `chrombpnet` 环境

用于官方 ChromBPNet 流程，不要拿它跑 Genos 或当前 TransChromBP 训练。

```bash
export CHROMBPNET_ENV=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet
export PATH="$CHROMBPNET_ENV/bin:$PATH"
export LD_LIBRARY_PATH="$CHROMBPNET_ENV/lib:$LD_LIBRARY_PATH"
export PYTHONPATH=/data1/zhoujiazhen/bylw_atac/chromBPNet:$PYTHONPATH
```

### 3.2 `transchrombp` 环境

用于当前 PyTorch 主训练线。

```bash
export TRANSCHROMBP_ENV=/data1/zhoujiazhen/bylw_atac/.mamba/envs/transchrombp
export PATH="$TRANSCHROMBP_ENV/bin:$PATH"
export LD_LIBRARY_PATH="$TRANSCHROMBP_ENV/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH=/data1/zhoujiazhen/bylw_atac/TransChromBP/src:${PYTHONPATH:-}
```

仓库内已有脚本默认按这个环境约定运行，例如：

- `scripts/benchmark/run_training_benchmark.sh`
- `scripts/benchmark/run_profiler.sh`

### 3.3 `genos-1.2b` 环境

用于 Genos 模型推理。若只是日常推理和特征提取，激活虚拟环境即可；只有在重新编译 `flash-attn` 时才需要显式补 `CUDA_HOME`。

```bash
source /data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b/bin/activate
export CUDA_HOME=/usr/local/cuda-12.3
export PATH=/usr/local/cuda-12.3/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-12.3/lib64:${LD_LIBRARY_PATH:-}
```

当前已在这个环境下通过真实验证的能力包括：

- `AutoConfig/AutoTokenizer/AutoModelForCausalLM.from_pretrained(..., local_files_only=True)` 本地加载
- `attn_implementation="flash_attention_2"` 的 GPU 最小前向
- embedding 提取
- next-token 预测
- greedy `generate`

## 4. 关键路径与资源位置

| 项目 | 路径 |
|---|---|
| Genos 官方代码仓库 | `/data1/zhoujiazhen/bylw_atac/foundation_models/genos/Genos` |
| Genos-1.2B 本地模型目录 | `/data1/zhoujiazhen/bylw_atac/foundation_models/genos/Genos-1.2B` |
| TransChromBP 远端代码根目录 | `/data1/zhoujiazhen/bylw_atac/TransChromBP` |
| Genos 可行性脚本 | `scripts/genos_feasibility.py` |
| Genos 可行性输出目录 | `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/genos_feasibility` |

`Genos-1.2B` 本地目录已确认不是“只有单个权重文件”的不完整下载，关键文件已齐全：

- `model.safetensors`
- `config.json`
- `tokenizer.json`
- `tokenizer_config.json`
- `special_tokens_map.json`

## 5. 已验证状态

### 5.1 Genos 接入

截至 2026-03-20，以下事项已完成：

- 官方仓库已克隆到 6000
- 官方 `requirements.txt` 对应依赖已装入隔离环境
- `pip check` 已通过
- 本地模型目录已通过加载验证
- 已在 `cuda:0` 上跑通 `flash_attention_2`
- 已跑通真实 usage 验证

真实 usage 验证摘要：

- 输入：`ACGT` 重复 16 次
- mean pooled hidden state：`(1, 1024)`，范数约 `28.5469`
- 下一位 top-1：`'A'`，概率约 `0.999512`
- greedy `generate`：可继续生成 `ACGT` 周期序列
- 峰值显存：约 `2391.1 MB`

### 5.2 当前未收口项

- 6000 到 `huggingface.co` 的网络可达性在上次检查时不稳定，因此 `hf download` 未完成一次新的远端同步
- 官方 Docker 兜底镜像 `bgigenos/mega:v1` 还没有在本地真正拉下
- 许可证口径存在冲突：GitHub 仓库 `LICENSE` 为 Apache-2.0，但 Hugging Face 页面与论文正文写 MIT；若后续对外分发或纳入正式发布流程，应先人工确认

上述问题不影响当前主机环境已经可以本地加载和使用 `Genos-1.2B`。

## 6. 已知坑与处理办法

### 6.1 官方 `requirements.txt` 的 `torchaudio` 版本约束有误

上游写法：

```text
torchaudio>=0.22.1,<0.23.0
```

按原样安装会报 `No matching distribution found`。当前实际安装时使用的是：

```text
torchaudio==2.7.1
```

后续若上游修正，应优先回归官方约束。

### 6.2 `flash-attn` 不能直接走默认 build isolation

按 `pip install -r requirements.txt` 直接安装时，`flash-attn` 会在构建阶段报：

```text
ModuleNotFoundError: No module named 'torch'
```

当前已验证可行的安装方式是：

1. 先装好 `torch`
2. 补 `ninja`
3. 设定 `CUDA_HOME=/usr/local/cuda-12.3`
4. 单独执行 `flash-attn` 的 `--no-build-isolation` 安装

### 6.3 `generate()` 前要去掉 `token_type_ids`

`Genos-1.2B` 在当前 HF 调用方式下，若把 tokenizer 输出原样传给 `generate()`，会因为未使用的 `model_kwargs` 报错。实际生成时应先移除 `token_type_ids`。

## 7. 如果要纳入 TransChromBP 做实验

当前推荐顺序：

1. 保持现有 `transchrombp` 环境不变
2. 用 `genos-1.2b` 环境先做特征提取或可行性分析
3. 优先尝试“冻结 Genos + 额外特征输入”的 pilot
4. 只有在 pilot 证明有净增益后，再考虑统一环境或在线联合训练

不推荐的做法：

- 直接往现有 `/data1/zhoujiazhen/bylw_atac/.mamba/envs/transchrombp` 里硬装 `transformers`、`flash-attn`、Genos 全套依赖
- 在没有确认收益前就把 Genos 作为可训练主干并入 DDP/optimizer/checkpoint

## 8. 相关文档

- `TRACKING.md`: 当前状态、日志路径、未收口项
- `AGENTS.md`: 6000 上 `chrombpnet` 环境激活约定
- `scripts/genos_feasibility.py`: Genos 可执行验证脚本
- `reports/transchrombp_genos_experiment_plan_20260320.tex`: Genos 融合实验计划
