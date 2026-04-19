# TransChromBP V2 代码改进消融实验方案

> 创建日期：2026-03-19
> 更新日期：2026-03-20（结合 ablation_tf_20260318 + bias_safe 分析报告调整优先级）
> 背景：基于 ChatGPT 5.4 Pro 对 V2 代码的审查建议 + Claude 对实际代码的交叉验证，
> 筛选出 5 项值得用实验验证的改动。每项改动设计为独立消融，改动前先跑实验看效果。

---

## 〇、结合 ablation_tf 结果的优先级调整（2026-03-20 新增）

分析报告 `reports/paper/transchrombp_ablation_tf_20260318_analysis_20260320.tex` 和
scale 轨迹图揭示了几个关键事实，直接影响实验优先级：

### 发现 1：Profile scale 已自行衰减到接近 0

αp 从 ~0.5 起步，5 epoch 内降到 <0.05，最终 ~0.01。
**结论：bias branch 对 profile 预测几乎没有实际贡献。**

- **实验 C（profile_scale_init 0.1）降级为可选**：
  αp 从 1.0 起步时已在 5 epoch 内自行衰减到 <0.05。
  从 0.1 起步只加速前 2-3 epoch 的衰减，最终终点相同。预期效果可忽略。
- **实验 A（freeze fix）的 profile 路径影响极小**：
  frozen bias 的 profile 输出被乘以 ~0.01 的 scale，即使 BN drift 也几乎不可见。

### 发现 2：Count scale 有意义且 V2-full 更低

αc 先升后降：peak ~2.0-2.4 → 30 epoch 后 full ~1.0、noTF ~1.4。
**结论：Transformer 在逐步接管 count 预测，减少对 bias 的依赖。**

- **实验 B（center pool）仍然有价值**：它改善的恰好是 signal branch 自身的 count 预测质量，
  方向与 Transformer 的自然趋势一致。
- **实验 A（freeze fix）的 count 路径可能有轻微影响**：
  αc ~1.0-1.4，frozen bias 的 count 输出会被乘以这个 scale，BN drift 可能造成轻微波动。

### 发现 3：Full/debiased gap 已经极小（0.00006）

说明当前 bias-safe 框架非常有效，profile shortcut 已被充分抑制。
- 实验 C 和 D 原本针对的"shortcut 回潮"问题在当前数据中不存在。

### 发现 4：30 epoch 偏紧，best epoch 集中在 27-30

六个 run 的 best 全在 27-30，曲线到 30 仍在下降。
- **所有新实验应改为 40 epoch**，否则可能因训练不充分而掩盖真实效果。
- 同时也使 checkpoint soup（实验 E）更有意义：更多末尾 epoch 可平均。

### 发现 5：Bias-Safe 2×2 消融表明 `stop_gradient` 应成为后续默认

`reports/paper/transchrombp_bias_safe_ablation_20260320.tex` 显示：

- 当前 `v2fix` 源配置没有显式设置 `fusion.profile_bias_stop_gradient`，运行时会落回代码默认 `false`
- 在 `2×2` 析因中，`sg=true` 对最终 `peak.profile_target_jsd_*` 几乎没有可见影响，但能明显改善 full/debiased 隔离
- 因此，**后续新 run（A/B 多种子、F/G、以及保留兜底的 C/D）统一改成 `sg=true` 更合理**

边界约束：

- 已经启动的 `v2fix_20260320_cpool_s42` / `v2fix_20260320_freeze_s42` 不重启
- `C/D` 继续保留为低优先级兜底，而不是硬删除；Bias-Safe 结果只是在方向上进一步支持它们降权

### 修订后的优先级

| 原优先级 | 实验 | 修订优先级 | 理由 |
|---|---|---|---|
| P0 | E: Soup | **P0 不变** | 零成本，更多 epoch 数据使 soup 更有效 |
| P0 | A: Freeze fix | **P1 降级** | full-debiased gap 已 0.00006，BN drift 实际影响极小；仍值得做因为是正确性修复 |
| P1 | B: Count pool | **P0 升级** | 直接改善 count 路径，与 Transformer 趋势一致，预期回报最高 |
| P1 | C: pscale=0.1 | **P3 降级/可选** | αp 已自行衰减到 ~0.01，改初始值无实际效果 |
| P2 | D: BN→GN | **P3 降级/可选** | 前提假设（BN drift 有害）被数据否定 |
| — | 40 epoch 复试 | **P0 新增** | 报告第一建议；与 B 合并可以一石二鸟 |

### 修订后的执行计划

**第零批**：E（soup 现有 checkpoint）
**第一批**：B（center pool）+ 40 epoch 复试（合并为一个实验）+ A（freeze fix，也用 40 epoch）
**统一约束**：除已在跑的 `cpool_s42` / `freeze_s42` 外，后续所有新启动的 v2fix run 默认在 `fusion` 段加 `profile_bias_stop_gradient: true`
**第二批（可选）**：C、D — 仅在第一批出现意外结果时补做

---

## 一、总体设计原则

- **单变量原则**：每个消融只改一个因素，其余与 V2 baseline 完全一致
- **共享 bias pretrain**：实验 A/B 可复用当前 `ablation_tf` 的 shared bias checkpoint（bias 架构不变）
- **3 seeds**：42, 1234, 2024（与现有 ablation 一致）
- **best_metric**：`peak.profile_target_jsd_debiased_mean` (min)
- **max_epochs**：**40**（从 30 提高，基于报告建议，给模型充足收敛空间）
- **评估指标**（按重要性排序）：
  1. `peak.profile_target_jsd_full_mean` — profile 质量（论文主指标）
  2. `peak.profile_target_jsd_debiased_mean` — signal branch 自身能力
  3. `count_pearson_full` / `count_mae_full` — count 预测质量（**升级为第三**）
  4. `profile_full_debiased_jsd` gap — shortcut 回潮指标（已确认极小，降为参考）
  5. `nonpeak.profile_target_jsd_full_mean` — 背景区质量
  6. seed 间标准差 — 稳定性指标

---

## 二、命名规范

### 实验 tag

统一前缀 `v2fix_YYYYMMDD`，后缀为实验代号：

| 实验 | run_name 模板 | 示例 |
|---|---|---|
| A: Freeze 修复 | `v2fix_{date}_freeze_s{seed}` | `v2fix_20260321_freeze_s42` |
| B: Count pooling | `v2fix_{date}_cpool_s{seed}` | `v2fix_20260321_cpool_s42` |
| C: profile_scale | `v2fix_{date}_pscale01_s{seed}` | `v2fix_20260321_pscale01_s42` |
| D: BN→GN | `v2fix_{date}_gn_s{seed}` | `v2fix_20260321_gn_s42` |
| D bias pretrain | `v2fix_{date}_gn_shared_bias` | — |
| E: Soup | 无需训练，后处理 | — |

### 目录结构

```
TransChromBP/outputs/
├── checkpoints/v2fix_{date}_freeze_s42/     # checkpoint 目录
│   ├── epoch_001.pt ... epoch_030.pt
│   └── best.pt
├── logs/v2fix_{date}_freeze_s42/
│   ├── epoch_metrics.jsonl                   # 每 epoch 验证指标
│   └── run_meta.json                         # 训练完成标记
└── runtime/v2fix_{date}/
    ├── configs/                              # 运行时生成的 config 副本
    └── meta/                                 # 实验元信息
```

---

## 三、Baseline 定义

Baseline 来自 `ablation_tf_20260318` 的 **V2-full** 组（3 seeds），已全部完成。
数据来源：`reports/assets/ablation_tf_20260318/summary_table.csv`，held-out test 评估。

| run_name | best_epoch | peak JSD full | peak JSD debiased | full-debiased gap | peak count r | nonpeak JSD |
|---|---|---|---|---|---|---|
| `full_s42` | 29 | 0.31425 | 0.31430 | 0.00006 | 0.83536 | 0.55939 |
| `full_s1234` | 29 | 0.31403 | 0.31408 | 0.00005 | 0.83107 | 0.55381 |
| `full_s2024` | 28 | 0.31430 | 0.31437 | 0.00007 | 0.84385 | 0.55684 |
| **均值±std** | 28.7 | **0.31419±0.00012** | **0.31425±0.00012** | **0.00006±0.00001** | **0.83676±0.00531** | **0.55668±0.00228** |

noTF 对照组（参考）：

| run_name | best_epoch | peak JSD full | peak JSD debiased | full-debiased gap | peak count r | nonpeak JSD |
|---|---|---|---|---|---|---|
| `noTF_s42` | 30 | 0.32391 | 0.32401 | 0.00010 | 0.82863 | 0.56268 |
| `noTF_s1234` | 29 | 0.32368 | 0.32378 | 0.00010 | 0.81812 | 0.55685 |
| `noTF_s2024` | 27 | 0.32419 | 0.32430 | 0.00011 | 0.82835 | 0.55988 |
| **均值±std** | 28.7 | **0.32393±0.00021** | **0.32403±0.00021** | **0.00010±0.00000** | **0.82503±0.00489** | **0.55980±0.00238** |

**关键观察**（来自 scale 轨迹图）：
- αp（profile scale）：从 ~0.5 衰减到 ~0.01（epoch 30），bias 对 profile 几乎无贡献
- αc（count scale）：full ~1.0，noTF ~1.4（epoch 30），bias 对 count 仍有贡献但在减少
- 所有 best epoch 在 27-30，训练远未收敛

关键路径：
- bias pretrain checkpoint: `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/checkpoints/ablation_tf_20260318_shared_bias/best.pt`
- epoch metrics: `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/logs/ablation_tf_20260318_full_s{seed}/epoch_metrics.jsonl`
- 分析报告: `reports/paper/transchrombp_ablation_tf_20260318_analysis_20260320.tex`

---

## 四、各实验详细设计

### 实验 A：freeze_bias_core eval 语义修复

#### 动机

当前 `freeze_core()` 只设 `requires_grad=False`，但 bias branch 里的 BatchNorm1d 和 Dropout
在 `model.train()` 下仍会更新 running stats / 随机 drop。这意味着冻结的 bias branch
每次 forward 的输出不完全一致，破坏了"预训练 → 冻结"的语义。

#### 需要修改的文件与精确位置

**文件 1**：`src/transchrombp/models/transchrombp.py`

改动 1 — `__init__` 里新增状态变量（约 line 144 附近，`freeze_bias_core: bool = False` 之后）：
```python
# 在 super().__init__() 之后新增
self._bias_branch_frozen = False
```

改动 2 — `freeze_bias_core()` 方法（line 350-354）：
```python
def freeze_bias_core(self) -> None:
    """Freeze bias branch parameters and set to eval mode."""
    if self.bias_branch is None:
        return
    self.bias_branch.freeze_core()
    self._bias_branch_frozen = True
    self.bias_branch.eval()  # 关键：关掉 dropout / 冻结 BN 统计
```

改动 3 — `unfreeze_bias_core()` 方法（line 356-360）：
```python
def unfreeze_bias_core(self) -> None:
    """Unfreeze bias branch parameters."""
    if self.bias_branch is None:
        return
    self.bias_branch.unfreeze_core()
    self._bias_branch_frozen = False
```

改动 4 — 新增 `train()` override（在 `unfreeze_bias_core` 之后）：
```python
def train(self, mode: bool = True):
    """Override train() to keep frozen bias branch in eval mode."""
    super().train(mode)
    if self.bias_branch is not None and self._bias_branch_frozen:
        self.bias_branch.eval()
    return self
```

改动 5 — `forward()` 里的 Case 2（line 305-309）：
```python
elif self.use_bias_decomposition:
    if self.bias_branch is None:
        raise RuntimeError("Bias decomposition is enabled but bias_branch is not initialized")
    if self._bias_branch_frozen:
        with torch.no_grad():
            profile_bias, count_bias = self.bias_branch(seq_onehot)
    else:
        profile_bias, count_bias = self.bias_branch(seq_onehot)
```

#### 需要创建的 config

不需要新的 model config 或 train config（与 baseline 完全一致），
只需要部署修改后的 `transchrombp.py` 到远端。

Train config 直接复用 `train_ablation_v2_main.yaml`（seed 和 run_name 由 orchestrator 脚本覆盖）。

#### 对照设计

| 组别 | 描述 | bias branch 行为 |
|---|---|---|
| baseline | `ablation_tf_20260318_full_s*` | freeze 时 BN/Dropout 仍 train mode |
| fix_freeze | `v2fix_*_freeze_s*` | freeze 时强制 eval + no_grad |

#### 预期

- 如果 fix 有效：debiased 指标的 seed 间方差减小，full-debiased gap 可能缩小
- 如果差异不大：说明当前 bias branch 的 BN/Dropout 影响较小（4 层 dilated conv + BN 在 eval vs train 差异有限）

#### 验证修复正确性的快速 smoke test

部署代码后，先做一个不训练的 forward 一致性检查：
```python
model.freeze_bias_core()
model.train()  # 模拟 trainer 行为
x = torch.randn(2, 2114, 4)
y1 = model(x)
y2 = model(x)
assert torch.allclose(y1.profile_bias, y2.profile_bias), "frozen bias should be deterministic"
assert torch.allclose(y1.count_bias, y2.count_bias), "frozen bias should be deterministic"
```

---

### 实验 B：Count Head Pooling 改法

#### 动机

当前 count 分支用 `encoded.mean(dim=1)` 对全长 2114 token 做 mean pooling，
但监督的 count 对应中心 1000bp 窗口。全长 pooling 会把两端各 557 token 的上下文一起平均，
稀释中心信号。Profile 分支已经用了 `center_crop_1d(output_len=1000)`，count 应对齐。

#### 需要修改的文件与精确位置

**文件**：`src/transchrombp/models/transchrombp.py`

改动（line 296，forward 方法内）：
```python
# 原始
pooled = encoded.mean(dim=1)

# 改为（center_crop_1d 操作在 last dim，encoded 是 [B, L, D]，需转换）
start = (encoded.size(1) - self.output_len) // 2
pooled = encoded[:, start : start + self.output_len, :].mean(dim=1)
```

**注意**：这里直接用 slice 而不是调用 `center_crop_1d`，因为 `center_crop_1d` 对 last dim 操作，
而 `encoded` 的 seq 维度在 dim=1。用 slice 更直观、无需 transpose。

#### 需要创建的 config

不需要新 config。Model config 和 train config 均复用 baseline。

#### 对照设计

| 组别 | pooling 范围 | 对应代码 |
|---|---|---|
| baseline | 全长 2114 token | `encoded.mean(dim=1)` |
| center_pool | 中心 1000 token | `encoded[:, start:end, :].mean(dim=1)` |

#### 预期

- `count_pearson_full` 提升（中心信号更聚焦）
- `count_mae_full` 下降
- profile 指标不应受影响（profile head 不变）
- 如果 count 指标反而变差：说明全长上下文对 count 预测有正面贡献（transformer 已学会利用上下文）

---

### 实验 C：profile_scale_init 调低

#### 动机

当前 `profile_scale_init=1.0`，训练开始时 profile bias 以等权参与融合。
V2 的实验已经表明 signal branch 自己就能学好 profile（debiased 接近 full）。
调低 profile_scale_init = 用初始化编码先验：bias 主要校准 count，不主导 profile。

#### 需要修改的文件

无代码修改。只需要创建一份新的 model config。

#### 需要创建的 config

**文件**：`configs/model/transchrombp_teacher_v2_pscale01.yaml`

基于 `transchrombp_teacher_v2.yaml`，只改一处：
```yaml
fusion:
  profile_fusion: add
  count_fusion: logsumexp
  profile_scale_init: 0.1   # ← 唯一改动，从 1.0 改为 0.1
  count_scale_init: 1.0
  learnable_scales: true
  positive_scales: true
```

Train config 复用 `train_ablation_v2_main.yaml`。

#### 对照设计

| 组别 | profile_scale_init | count_scale_init |
|---|---|---|
| baseline | 1.0 | 1.0 |
| pscale01 | 0.1 | 1.0 |

#### 预期

- `profile_full_debiased_jsd` gap 缩小（bias 对 profile 贡献更小）
- `peak.profile_target_jsd_full_mean` 不掉或略好（signal branch 被迫独立学好 profile）
- count 指标不受影响（count_scale 未改）

---

### 实验 D：Bias Branch BatchNorm1d → GroupNorm

#### 动机

`ResidualDilatedConvBlock` 用 `BatchNorm1d`，在小 batch + peak/nonpeak 混合比例波动时
统计量不稳定。`GroupNorm` 不依赖 batch 统计，训练和推理行为完全一致，对冻结语义更友好。

#### 需要修改的文件与精确位置

**文件**：`src/transchrombp/models/bias_branch.py`

改动（line 60，`ResidualDilatedConvBlock.__init__`）：
```python
# 原始
self.norm = nn.BatchNorm1d(channels)

# 改为
self.norm = nn.GroupNorm(num_groups=min(8, channels // 4), num_channels=channels)
# 对于 hidden_channels=128: GroupNorm(8, 128)，即每组 16 个 channel
```

#### 需要创建的 config

**文件**：`configs/model/transchrombp_teacher_v2_gn.yaml`

基于 `transchrombp_teacher_v2.yaml`，在 bias_branch 段新增一个参数：
```yaml
bias_branch:
  type: chrombpnet_bias
  enabled: true
  hidden_channels: 128
  kernel_size: 21
  n_dil_layers: 4
  norm_type: groupnorm    # ← 新增参数，需要在代码中支持
```

> **注意**：当前代码中 `ResidualDilatedConvBlock` 没有 `norm_type` 参数，
> 需要同时改 `bias_branch.py` 支持 config 传入的 norm 类型选择。
> 或者更简单的做法：直接硬编码为 GroupNorm，不加 config 参数。
> **建议采用简单做法**：直接在代码里改 `nn.BatchNorm1d` → `nn.GroupNorm`，
> 因为这个实验本身就是为了验证 GN 是否更好。

#### 特殊注意

- **需要重新训练 bias pretrain**（bias 架构变了，state_dict key 不同，不能复用现有 bias checkpoint）
- 建议同时应用实验 A 的 freeze 修复（避免混淆变量）
- bias pretrain 用的 train config 和 data config 与 baseline 完全一致

#### 对照设计

| 组别 | bias branch norm | freeze 语义 |
|---|---|---|
| baseline | BatchNorm1d | 旧版（无 eval 强制） |
| fix_freeze (实验 A) | BatchNorm1d | 新版（eval 强制） |
| gn_bias (实验 D) | GroupNorm(8, 128) | 新版（eval 强制） |

> 实验 D 的准确对照应该是"实验 A 的结果"，而不是原始 baseline。
> 这样才能隔离 norm 类型这个单一变量。

#### 预期

- bias pretrain 阶段：训练曲线更稳（seed 间方差更小）
- 主训练阶段：冻结 bias 的输出在 eval 和 train 之间完全一致（GN 没有 running stats）
- 整体指标可能持平或略好

---

### 实验 E：Checkpoint Soup（后处理，不需要新训练）

#### 动机

多个 seed 的 best epoch 都贴近训练末尾（如 `noTF_s1234` best at epoch 29/30），
说明模型仍在持续改善区间。对末尾 K 个 epoch 做 weight average（soup）
通常可以在不增加训练成本的前提下提升泛化。

#### 实现

需要写一个独立脚本 `scripts/checkpoint_soup.py`：

```python
"""Checkpoint soup: average weights from multiple checkpoints and evaluate."""
import argparse
import torch
from pathlib import Path

def load_model_state(ckpt_path):
    ckpt = torch.load(ckpt_path, map_location="cpu")
    return ckpt["model_state"]

def average_states(state_dicts):
    avg = {}
    for key in state_dicts[0]:
        stacked = torch.stack([sd[key].float() for sd in state_dicts])
        avg[key] = stacked.mean(dim=0)
    return avg

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ckpt-dir", required=True, help="checkpoints/run_name/")
    parser.add_argument("--mode", choices=["last_k", "top_k_jsd"], default="last_k")
    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--metrics-jsonl", help="epoch_metrics.jsonl for top_k selection")
    parser.add_argument("--output", required=True, help="output soup checkpoint path")
    args = parser.parse_args()

    ckpt_dir = Path(args.ckpt_dir)
    all_ckpts = sorted(ckpt_dir.glob("epoch_*.pt"))

    if args.mode == "last_k":
        selected = all_ckpts[-args.k:]
    elif args.mode == "top_k_jsd":
        # 从 metrics jsonl 读取每个 epoch 的 JSD，选最好的 K 个
        import json
        metrics = []
        with open(args.metrics_jsonl) as f:
            for line in f:
                metrics.append(json.loads(line))
        # 按 peak JSD debiased 排序（越小越好）
        metrics.sort(key=lambda m: m["val"]["peak"]["profile_target_jsd_debiased_mean"])
        top_epochs = [m["epoch"] for m in metrics[:args.k]]
        selected = [ckpt_dir / f"epoch_{e:03d}.pt" for e in top_epochs]
        selected = [p for p in selected if p.exists()]

    print(f"Averaging {len(selected)} checkpoints: {[p.name for p in selected]}")
    states = [load_model_state(p) for p in selected]
    avg_state = average_states(states)

    # 保存为与原始 checkpoint 相同格式
    base_ckpt = torch.load(selected[-1], map_location="cpu")
    base_ckpt["model_state"] = avg_state
    base_ckpt["soup_sources"] = [str(p) for p in selected]
    torch.save(base_ckpt, args.output)
    print(f"Saved soup to {args.output}")

if __name__ == "__main__":
    main()
```

#### 对照设计

对每个已完成的 baseline run（`ablation_tf_20260318_full_s*`）分别做：

| 组别 | 描述 | checkpoint |
|---|---|---|
| best.pt | 当前 best checkpoint（单点） | `best.pt` |
| soup_last5 | 最后 5 个 epoch 的 weight average | `epoch_{26..30}.pt` 平均 |
| soup_top3_jsd | validation JSD 最好的 3 个 epoch | 按 `peak.profile_target_jsd_debiased_mean` 选 |

#### 评估方式

用训练中相同的 validation pipeline 评估 soup checkpoint：
```bash
# 伪命令 — 需要一个 evaluate-only 入口
python -m transchrombp.training.evaluate \
  --model-config configs/model/transchrombp_teacher_v2.yaml \
  --checkpoint outputs/checkpoints/run_name/soup_last5.pt \
  --data-config configs/data/data_tutorial_canonical_v1.yaml \
  --split test
```

> **注意**：当前 `train_ddp.py` 没有独立的 evaluate-only 入口。
> 需要写一个轻量的 `scripts/evaluate_checkpoint.py` 复用 `run_validation()` 逻辑。
> 或者直接在 soup 脚本里内嵌评估。

#### 预期

- soup 通常比单点 best 提升 0.5–2%
- 如果训练末尾已经很平坦，soup 效果更明显
- 如果 soup 反而变差：说明模型末尾在震荡，weight space 不平坦

---

### 实验 F：Count Head Attention Pooling（2026-03-20 新增）

> 来源：ChatGPT Pro 第 2 层建议 — "让模型自己学习哪些位置对 count 最重要"

#### 动机

实验 B 的 center pooling 通过缩小 pooling 范围来聚焦中心信号，但方式比较"硬"：
固定取中心 1000 token。Attention pooling 更灵活——模型自己学习每个位置的权重，
理论上可以适应不同 region 的信息分布（某些 region 的关键信号可能不完全在中心）。

如果 transformer 确实让某些位置学到了跨区域上下文信息，attention pooling 能直接捕获。

#### 代码改动

**文件**：`src/transchrombp/models/transchrombp.py`（已实现）

- `__init__` 中 `count_pool_mode` 新增 `"attention"` 选项
- 新增 `self.count_pool_proj = nn.Linear(d_model, 1)` 投影层
- `_pool_for_count()` 新增 attention 分支：

```python
if self.count_pool_mode == "attention":
    attn_logits = self.count_pool_proj(encoded).squeeze(-1)          # [B, L]
    attn_weights = torch.softmax(attn_logits, dim=1)                 # [B, L]
    return torch.sum(encoded * attn_weights.unsqueeze(-1), dim=1)    # [B, D]
```

- `build_transchrombp_from_config` 已支持通过 `heads.count_pool_mode: "attention"` 传入

#### Config

不需要新 model config 文件。只需在运行时覆盖 count_pool_mode：

```yaml
heads:
  count_pool_mode: "attention"  # ← 唯一改动
```

#### 对照设计

| 组别 | count_pool_mode | 备注 |
|---|---|---|
| baseline | full | 全长 mean pooling (2114 tokens) |
| center_pool (实验 B) | center | 中心 1000 tokens mean |
| attn_pool (实验 F) | attention | 学习的 attention pooling |

#### 预期

- `count_pearson_full` 提升（attention 权重自动学会聚焦信息密集位置）
- 如果 attention pooling > center pooling：说明固定中心不够灵活
- 如果 attention pooling ≈ center pooling：说明信号主要集中在中心，两种方式等价
- 如果 attention pooling < center pooling：说明自由度太大导致过拟合，center 的硬约束更好

#### 诊断工具

训练完成后，可以可视化 attention weights 的分布：
```python
# 检查 attention 是否集中在中心
with torch.no_grad():
    attn_logits = model.count_pool_proj(encoded).squeeze(-1)
    attn_weights = torch.softmax(attn_logits, dim=1)
    # 画 attn_weights 的热力图，看是否聚焦中心
```

---

### 实验 G：Profile Head 轻量 Conv Refinement（2026-03-20 新增）

> 来源：ChatGPT Pro 第 2 层建议 — "transformer 负责看远处，conv 负责把局部峰形磨尖"

#### 动机

当前 profile head 是 `Linear(d_model, 1)`，直接将 transformer 的每个 token 映射到标量。
这对全局上下文友好，但对 profile 的高分辨率峰形（sharp peaks）有点"太平"。

更好的设计是在 transformer 输出后加一个轻量的 depthwise-separable conv，
让全局上下文先被 transformer 建模，再由 conv 精修局部峰形。
这样 long-range information 更容易真正落到 JSD 指标上。

#### 代码改动

**文件**：`src/transchrombp/models/transchrombp.py`（已实现）

- `__init__` 中新增 `use_profile_refine: bool = False` 参数
- 当启用时，用 depthwise conv + pointwise conv + 残差连接替代单层 Linear：

```python
if self.use_profile_refine:
    self.profile_refine = nn.Sequential(
        nn.Conv1d(d_model, d_model, kernel_size=5, padding=2, groups=d_model),  # depthwise
        nn.GELU(),
        nn.Conv1d(d_model, d_model, kernel_size=1),  # pointwise
        nn.GELU(),
    )
    self.profile_signal_head = nn.Conv1d(d_model, 1, kernel_size=1)
else:
    self.profile_signal_head = nn.Linear(d_model, 1)
```

- forward 中的 profile 分支用残差连接：

```python
if self.use_profile_refine:
    feat = encoded.transpose(1, 2)                                  # [B, D, L]
    feat = feat + self.profile_refine(feat)                         # 残差连接
    profile_signal = self.profile_signal_head(feat).squeeze(1)      # [B, L]
else:
    profile_signal = self.profile_signal_head(encoded).squeeze(-1)  # [B, L]
```

- `build_transchrombp_from_config` 已支持通过 `heads.use_profile_refine: true` 传入

#### 设计选择说明

- **Depthwise separable conv** 而非普通 conv：参数量极小（~d_model×5 + d_model² ≈ 67K for d_model=256），不会引入过拟合风险
- **残差连接**：保证模型至少不比 baseline 差——refinement 可以学为零
- **kernel_size=5**：覆盖 ±2bp 的局部上下文，足够精修峰形但不会太大
- **Conv1d(d, 1) 替代 Linear(d, 1) 作为最终 head**：在 use_profile_refine 模式下保持 channel-first 一致性

#### Config

```yaml
heads:
  use_profile_refine: true  # ← 唯一改动
```

#### 对照设计

| 组别 | profile head | 新增参数 |
|---|---|---|
| baseline | Linear(d_model, 1) | 0 |
| profile_refine (实验 G) | depthwise conv + pointwise conv + Conv1d(d,1) + 残差 | ~67K |

#### 预期

- `peak.profile_target_jsd` 改善（峰形更尖锐，JSD 更低）
- 参数量增加极小，不影响训练速度
- 如果效果不大：说明 Linear head 已经足够，或者瓶颈不在局部峰形

---

## 五、统一评估协议

### 评估入口

所有实验用统一的评估脚本，避免指标计算口径不一致。

**方案 A**（推荐）：写一个 `scripts/evaluate_checkpoint.py`，复用 `train_ddp.py` 中的
`run_validation()` + `compute_selection_metric_sums()` + `finalize_selection_metrics()` 逻辑。

**方案 B**：直接从训练日志的 `epoch_metrics.jsonl` 提取 best epoch 的验证指标
（训练时已经计算了所有需要的指标）。

> 对于 A/B/C/D/F/G 实验，方案 B 已经足够——训练时的 validation 已经算好了。
> 对于实验 E（soup），必须用方案 A——soup checkpoint 没有对应的训练日志。

### 命名规范（补充）

| 实验 | run_name 模板 | 示例 |
|---|---|---|
| F: Attention pooling | `v2fix_{date}_attnpool_s{seed}` | `v2fix_20260321_attnpool_s42` |
| G: Profile refine | `v2fix_{date}_profref_s{seed}` | `v2fix_20260321_profref_s42` |

### 结果汇总表

所有实验完成后，填充以下汇总表：

| 实验 | seed | peak JSD full ↓ | peak JSD debiased ↓ | full-debiased gap ↓ | count_pearson ↑ | count_mae ↓ | nonpeak JSD ↓ | best_epoch |
|---|---|---|---|---|---|---|---|---|
| baseline_s42 | 42 | | | | | | | |
| baseline_s1234 | 1234 | | | | | | | |
| baseline_s2024 | 2024 | | | | | | | |
| **baseline 均值** | — | | | | | | | |
| **baseline std** | — | | | | | | | |
| fix_freeze_s42 | 42 | | | | | | | |
| attn_pool_s42 | 42 | | | | | | | |
| prof_refine_s42 | 42 | | | | | | | |
| ... | | | | | | | | |

### 合并决策规则

每个实验跑完 3 seeds 后，按以下规则决定是否合入代码：

1. **必要条件**：3 seeds 的 `peak.profile_target_jsd_debiased_mean` 均值不劣于 baseline
   - 判定方法：paired t-test（同 seed 对比），p < 0.1 为显著劣化；否则视为"不劣于"
2. **充分条件**（满足任一即合入）：
   - peak JSD debiased 均值显著更好（p < 0.1）
   - `count_pearson_full` 显著提升 **且** peak JSD 不劣
   - `full-debiased gap` 缩小 > 10% **且** peak JSD 不劣
   - seed 间标准差显著减小（稳定性提升）
3. **组合叠加**：如果多个实验都通过，按 A → B/F → G → C 顺序依次叠加，每叠加一个后确认组合效果不冲突

---

## 六、执行时间线与依赖（2026-03-20 修订 v2）

```
当前 (2026-03-20)
│
├── [已完成] ablation_tf_20260318 (6000)
│     6 runs 全部完成，held-out test 评估已出
│     分析报告已写
│
├── [进行中] bias_safe_ablation (6002)
│     B3/B4 仍在跑
│
▼ ─────────── 第零批（立即可做）──────────
│
├── [E] Checkpoint Soup (零训练成本)
│     脚本已写好：checkpoint_soup.py + evaluate_checkpoint.py
│     耗时：~2h（部署 + 6 个 run 评估）
│
▼ ─────────── 第一批（核心实验）──────────
│
├── [B+40ep] Count Center Pooling + 40 epoch 复试 (6000, ~7h)
│     合并原实验 B 与"训练轮数复试"
│     config: heads.count_pool_mode=center + max_epochs=40
│     3 seeds 串行 × ~2.1h/run ≈ 6.4h
│
├── [A+40ep] Freeze 修复 + 40 epoch (6000 or 6002, ~7h/~58h)
│     正确性修复 + 更长训练作为新 baseline
│
▼ ─────────── 第一批结果分析 ──────────────
│
├── 对比 4 组（baseline 30ep / A 40ep / B 40ep / soup）
│     如果 B 的 count 改善明确 → 在 B 基础上叠加 F 和 G
│     如果 B 改善不大 → 直接跑 F 替代 B
│
▼ ─────────── 第二批（读出头改进）──────────
│
├── [F] Attention Pooling + 40ep (~7h)
│     若 B 不够 → 替代方案
│     若 B 好 → 与 B 对比哪个更优
│
├── [G] Profile Refinement + 40ep (~7h)
│     与 B 或 F 的最佳方案正交，可叠加
│
▼ ─────────── 第二批结果分析 ──────────────
│
├── 汇总 best 组合
│     如果收益仍小 → 考虑第 4 层改动（降采样 TF、gated residual）
│     如果收益明确 → 合入代码
│
▼ ─────────── 第三批（可选）──────────
│
├── [C] pscale=0.1 — 仅在 αp 未自行衰减时补做
├── [D] BN→GN — 仅在 A 修复后稳定性仍不足时
├── [E+] EMA — 如果 soup 有效但想进一步提升
└── 第 4 层改动 — 仅 1-3 层全部无效时考虑
```

### 资源分配建议

| 服务器 | GPU | 可分配实验 | 估计耗时（40 epoch） |
|---|---|---|---|
| 6000 | 2×A6000 | B+40ep → A+40ep → F+40ep → G+40ep | ~28h 串行 |
| 6002 | 1×RTX3080 | 可并行分担部分实验 | ~19h/run（串行 3 seeds） |

> 推荐方案：6000 串行跑 B → A → F → G，总计 ~28h。
> 或 6000 跑 B+A（~14h），6002 并行跑 F+G。

### 为什么把 C 和 D 降为可选

**实验 C（pscale=0.1）**：scale 轨迹图明确显示 αp 从 1.0 起步已在 5 epoch 内自行衰减到 <0.05。
将初始值从 1.0 改到 0.1 只改变前 2-3 epoch 的过渡行为，对最终性能没有实质影响。
除非出现 αp 不衰减的新情况（当前数据中从未出现），否则不值得花 3 seeds 的 GPU 时间。

**实验 D（BN→GN）**：ChatGPT 建议的前提是"BN 在冻结后仍会漂移导致不稳定"。
但 full-debiased gap 已经是 0.00006（本质上为零），说明即使有漂移也没有造成实际问题。
实验 A 修复 freeze 语义后，BN 的 running stats 被彻底冻结在 eval 值，
GN 的唯一额外好处（消除 train/eval 行为差异）就完全没有必要了。

---

## 七、风险与应急

| 风险 | 影响 | 应对 |
|---|---|---|
| 40 epoch 仍不够（best 仍贴尾） | 仍未摸到上限 | 补 50 epoch 单 seed 探针 |
| 实验 B count pooling 反而变差 | Transformer 利用了全长上下文做 count | 改用实验 F attention pooling |
| 实验 F attention pooling 过拟合 | 自由度太大 | 加 dropout 或约束 attention 分布，或退回 center pool |
| 实验 G profile refine 无效果 | 瓶颈不在局部峰形 | kernel_size 从 5 试到 9，或加第二层 conv |
| 实验 A freeze 修复无差异 | BN drift 影响确实可忽略（数据已暗示） | 仍然合入代码（正确性修复），不影响其它结论 |
| A+40ep 与 baseline 30ep 不可直接对比（epoch 不同） | 无法隔离 freeze fix 的效果 | 可从 A 40ep 的 epoch 28-30 处截取对比 30ep baseline |
| F 和 G 叠加时相互干扰 | 无法判断各自贡献 | 先单独验证 F 和 G，再做叠加实验 |
