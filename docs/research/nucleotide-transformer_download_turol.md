## 目标

> 本次下载与 smoke 过程中沉淀下来的小型脚本和 manifest 已归档到 `docs/research/assets/nt_v2_download_20260324/`。

在 Linux 服务器上完成以下任务：

1. 创建一个干净的 Python 3.11 环境。
2. 下载并跑通 `InstaDeepAI/nucleotide-transformer-v2-500m-multi-species`。
3. 用 Hugging Face / PyTorch 跑一个最小 smoke test：
   - 加载 tokenizer 和 model
   - 输入两条短 DNA 序列
   - 输出 logits / hidden states 的 shape
   - 计算 mean pooled sequence embedding
4. 再克隆官方仓库 `instadeepai/nucleotide-transformer`，安装后补跑一个 JAX smoke test。
5. 输出最终环境信息、下载目录、运行日志、失败点和修复动作。

------

## 背景与执行原则

- 默认优先选择 **`InstaDeepAI/nucleotide-transformer-v2-500m-multi-species`**，因为它是公开可下载的 NT v2 模型，适合先验证整条链路。
- **不要求先申请 Hugging Face 权限**。公开仓库通常可直接下载；只有 private / gated 仓库才必须登录或显式提供 token。
- 如果下载时出现 `401/403`，再尝试 `hf auth login` 或设置 `HF_TOKEN`。
- 本任务不要求复现官方原始预训练，因为官方已说明原始预训练脚本包含 proprietary 元素，未完整公开。
- 官方仓库当前推荐 Python 3.11，并支持通过本地安装的方式使用。

------

## 一、前置检查

先执行：

```
set -euxo pipefail

python3 --version || true
python3.11 --version || true
which python3 || true
nvidia-smi || true
df -h
uname -a
```

记录输出到日志文件。

------

## 二、创建工作目录与 Python 环境

```
set -euxo pipefail

mkdir -p ~/work
cd ~/work

python3.11 -m venv nt-env
source ~/work/nt-env/bin/activate

python -m pip install --upgrade pip setuptools wheel
python -V
pip -V
```

------

## 三、安装 HF / PyTorch 路径所需依赖

先安装稳定版：

```
set -euxo pipefail
source ~/work/nt-env/bin/activate

pip install --upgrade huggingface_hub "transformers>=4.52.4" torch
```

如果后续加载 NT 模型时报 `transformers` 兼容问题，再执行：

```
set -euxo pipefail
source ~/work/nt-env/bin/activate

pip install --upgrade git+https://github.com/huggingface/transformers.git
```

`hf download` 支持把 Hub 仓库直接下载到指定本地目录；这正是本任务推荐的下载方式。

------

## 四、可选：登录 Hugging Face

仅在下载公开模型失败时再做。

```
set -euxo pipefail
source ~/work/nt-env/bin/activate

hf auth login
```

或者：

```
export HF_TOKEN="你的token"
```

------

## 五、下载模型到本地目录

目标模型：

- `InstaDeepAI/nucleotide-transformer-v2-500m-multi-species`

执行：

```
set -euxo pipefail
source ~/work/nt-env/bin/activate

mkdir -p ~/models/nt-v2-500m
hf download InstaDeepAI/nucleotide-transformer-v2-500m-multi-species \
  --local-dir ~/models/nt-v2-500m
```

如果网络不稳，可先设置：

```
export HF_HUB_DOWNLOAD_TIMEOUT=30
```

------

## 六、编写 HF/PyTorch smoke test

创建文件 `~/work/smoke_test_hf.py`，内容如下：

```
import os
import torch
from transformers import AutoTokenizer, AutoModelForMaskedLM

MODEL_DIR = os.path.expanduser("~/models/nt-v2-500m")

tokenizer = AutoTokenizer.from_pretrained(
    MODEL_DIR,
    trust_remote_code=True,
)

model = AutoModelForMaskedLM.from_pretrained(
    MODEL_DIR,
    trust_remote_code=True,
)

model.eval()

sequences = [
    "ATTCCGATTCCGATTCCG",
    "ATTTCTCTCTCTCTCTGAGATCGATCGATCGAT",
]

max_length = min(getattr(tokenizer, "model_max_length", 256), 256)

batch = tokenizer.batch_encode_plus(
    sequences,
    return_tensors="pt",
    padding="max_length",
    truncation=True,
    max_length=max_length,
)

attention_mask = batch["input_ids"] != tokenizer.pad_token_id

with torch.no_grad():
    outputs = model(
        batch["input_ids"],
        attention_mask=attention_mask,
        encoder_attention_mask=attention_mask,
        output_hidden_states=True,
    )

hidden = outputs["hidden_states"][-1]
mask = attention_mask.unsqueeze(-1)
mean_emb = (hidden * mask).sum(dim=1) / mask.sum(dim=1)

print("input_ids shape:", tuple(batch["input_ids"].shape))
print("attention_mask shape:", tuple(attention_mask.shape))
print("last_hidden_state shape:", tuple(hidden.shape))
print("mean_sequence_embedding shape:", tuple(mean_emb.shape))
```

> 备注：这里默认开启 `trust_remote_code=True`，因为 NT v2 模型页当前是按自定义实现方式提供使用示例的。

执行：

```
set -euxo pipefail
source ~/work/nt-env/bin/activate

cd ~/work
python smoke_test_hf.py | tee hf_smoke_test.log
```

------

## 七、验收 HF 路径

满足以下条件即视为成功：

1. 无 import error
2. 无模型下载错误
3. 成功打印以下 shape：
   - `input_ids shape`
   - `attention_mask shape`
   - `last_hidden_state shape`
   - `mean_sequence_embedding shape`

若失败，优先排查：

- `transformers` 版本兼容问题
- 本地磁盘不足
- 下载中断 / 网络代理
- 未正确启用虚拟环境
- CUDA / torch 版本不匹配（若使用 GPU）

------

## 八、克隆官方仓库并安装

官方仓库当前提供标准的本地安装方式，适合做第二条验证链路。

执行：

```
set -euxo pipefail

cd ~/work
git clone https://github.com/instadeepai/nucleotide-transformer.git
cd ~/work/nucleotide-transformer

python3.11 -m venv .venv
source ~/work/nucleotide-transformer/.venv/bin/activate

python -m pip install --upgrade pip setuptools wheel
pip install .
```

------

## 九、安装 JAX

### CPU 版

```
set -euxo pipefail
source ~/work/nucleotide-transformer/.venv/bin/activate

pip install --upgrade jax
```

### NVIDIA GPU 版（Linux）

JAX 官方当前推荐优先使用 pip wheels 安装 GPU 版；如果服务器是 Linux + NVIDIA GPU，就先尝试：

```
set -euxo pipefail
source ~/work/nucleotide-transformer/.venv/bin/activate

pip install --upgrade "jax[cuda13]"
```

如果失败，再退回 CPU 版，或者根据机器 CUDA 实际情况改用官方文档对应方式。

------

## 十、编写官方 repo / JAX smoke test

创建文件 `~/work/nucleotide-transformer/smoke_test_jax.py`：

```
import haiku as hk
import jax
import jax.numpy as jnp
from nucleotide_transformer.pretrained import get_pretrained_model

MODEL_NAME = "500M_human_ref"

parameters, forward_fn, tokenizer, config = get_pretrained_model(
    model_name=MODEL_NAME,
    embeddings_layers_to_save=(24,),
    max_positions=128,
)

forward_fn = hk.transform(forward_fn)

sequences = [
    "ATTCCGATTCCGATTCCG",
    "ATTTCTCTCTCTCTCTGAGATCGATCGATCGAT",
]

tokens_ids = [b[1] for b in tokenizer.batch_tokenize(sequences)]
tokens = jnp.asarray(tokens_ids, dtype=jnp.int32)

rng = jax.random.PRNGKey(0)
outs = forward_fn.apply(parameters, rng, tokens)

print("JAX devices:", jax.devices())
print("tokens shape:", tuple(tokens.shape))

for k, v in outs.items():
    try:
        print(k, tuple(v.shape))
    except Exception:
        print(k, type(v))
```

执行：

```
set -euxo pipefail
source ~/work/nucleotide-transformer/.venv/bin/activate

cd ~/work/nucleotide-transformer
python smoke_test_jax.py | tee jax_smoke_test.log
```

官方文档当前提供 `get_pretrained_model(...)` 的使用方式，并列出了 NT 模型的预训练加载接口。

------

## 十一、最终输出

请在任务结束时输出一个简明总结，至少包含：

### 1. 环境信息

```
python --version
pip freeze | sort
nvidia-smi || true
```

### 2. 目录信息

```
du -sh ~/models/nt-v2-500m || true
du -sh ~/work/nucleotide-transformer || true
ls -lah ~/work
```

### 3. 运行日志位置

- `~/work/hf_smoke_test.log`
- `~/work/nucleotide-transformer/jax_smoke_test.log`

### 4. 最终结论

用以下模板输出：

```
[结论]
- HF/PyTorch 路径：成功 / 失败
- 官方 repo/JAX 路径：成功 / 失败

[成功项]
- 已完成 …
- 已验证 …

[失败项]
- 失败步骤：
- 报错摘要：
- 已尝试修复：
- 下一步建议：
```

------

## 十二、故障修复优先级

按这个顺序排查：

1. Python 版本不对
2. 没激活正确虚拟环境
3. `transformers` 版本兼容问题
4. Hugging Face 下载失败
5. 磁盘空间不足
6. GPU / CUDA / JAX 不兼容
7. JAX GPU 安装失败则退回 CPU 验证

------

## 十三、不要做的事

- 不要尝试“完整复现官方原始预训练流程”。
- 不要默认直接上 2.5B 模型。
- 不要假设 `git clone` 会把大权重一起拉下来。
- 不要在没验证公开下载是否可用前就强制要求 HF token。

------

## 十四、最低成功标准

只要满足以下两条，就算本轮任务完成：

1. `smoke_test_hf.py` 成功运行并打印 embedding shape
2. `smoke_test_jax.py` 至少能成功加载模型并打印设备与输出结构
