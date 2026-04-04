# Strict Compare 后续分析评估（2026-03-27）

## 总结

对这轮新增分析，我的结论是：

- `A`：**成立，而且是当前除 selector 泄漏外最值得优先处理的问题。**
- `B`：**成立，但它更像 limitation / 归因边界，不是实现 bug。**
- `C`：**大方向成立。原来的验证方法不够好，应该改成 runtime 级别验证。**
- `D`：**“存在语义歧义”不成立；代码已经明确说明 warmup_ratio 会生效。**
- `E`：**基本成立，属于文档记录项，不是 blocker。**

---

## A. 训练 region 集合不一致

### 结论

**成立。**

而且这不是“理论上可能不同”，而是我在 6000 上已经直接核对出差异。

### 证据

官方 controlled 臂代码：

- `ChromBPNet` 训练直接吃原始 peaks `overlap.bed.gz`，并现场 `prep nonpeaks`：
  - `scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh:205-213`
  - `scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh:233-246`

自研 controlled 臂配置：

- `TransChromBP` tutorial data config 直接指向自家预处理产物：
  - `vendor/transchrombp/transchrombp/configs/data/data_tutorial_canonical_v1.yaml:6-7`

6000 上的实测结果：

- peaks：
  - 官方原始 peaks：`269800`
  - `TransChromBP filtered.peaks.bed`：`267175`
  - 以 `(chr,start,end,summit)` 为 key 比较，`filtered.peaks` 是原始 peaks 的**严格子集**，少了 `2625` 条
- nonpeaks：
  - 官方 candidate nonpeaks：`508429`
  - `TransChromBP filtered.nonpeaks.bed`：`267175`
  - `filtered.nonpeaks` 是官方 candidate nonpeaks 的**严格子集**，少了 `241254` 条

所以“same frozen dataset”如果指的是：

- 同一原始 BAM / genome / fold：**成立**
- 同一实际 candidate peaks/nonpeaks 集合：**不成立**

### 影响

这会削弱“严格控制变量”的说服力，尤其是当我们想把结果写成“只差 architecture / recipe”时。

### 建议

在 tutorial / GM12878 的公平主表前，必须统一 region 集合。更稳的两种做法是：

1. 让 `TransChromBP` 直接吃 official arm 用的 peaks + nonpeaks
2. 或先冻结一份共享 `filtered.peaks + filtered.nonpeaks`，再让两边都用这份

---

## B. 优化器 / schedule / bf16 差异

### 结论

**成立，但性质是“归因边界”而不是“工程 bug”。**

### 证据

- `ChromBPNet` 模型 compile 用 `Adam(learning_rate=args.learning_rate)`：
  - `chrombpnet/training/models/chrombpnet_with_bias_model.py:134`
- `TransChromBP` strict-compare config 用：
  - `precision: bf16`
  - `optimizer: adamw`
  - `learning_rate: 5e-4`
  - `weight_decay: 0.01`
  - `schedule: cosine + warmup_ratio 0.06`
  - 文件：`vendor/transchrombp/transchrombp/configs/train/train_tutorial_corrected_b_strict_compare_6000.yaml`

### 判断

如果当前目标是比较 **full systems under matched data/hardware/batch budget**，那这不算 bug。  
但如果要把优势写成“Transformer 架构本身更强”，这就是未控制混杂变量。

### 建议

- 主文必须明确这是 **system-vs-system** 比较
- Table / limitation 里显式列出这些未锁定变量
- 若后面上 `noTF-controlled`，必须沿用与 `corrected-B-controlled` 完全相同的 optimizer / schedule / precision

---

## C. MirroredStrategy global batch 验证方法

### 结论

**成立。**

你的修正比之前“看 `steps_per_epoch`”更准确。

### 证据

- generator 的 `__len__()` 只返回 `ceil(N / batch_size)`：
  - `chrombpnet/training/data_generators/batchgen_generator.py:53-55`
- `__getitem__()` 也是按这个 batch size 切片：
  - `chrombpnet/training/data_generators/batchgen_generator.py:99-107`
- 多卡只是外层包 `MirroredStrategy()`：
  - `chrombpnet/training/train.py:95-117`

因此仅看 `steps_per_epoch`，确实无法区分“global 32”还是“per-replica 32”。

### 建议

把验证改成 runtime probe：

- 最直接的是打印 replica 看到的 `inputs.shape`
- 或者插一个 callback / model-side debug print，只跑极短 smoke

“global batch 大概率是 32”我同意，但这仍然是**推断**，不是已闭环证据。

---

## D. `warmup_steps=0 + warmup_ratio=0.06` 是否有歧义

### 结论

**不成立。这里没有语义歧义，代码已经明确。**

### 证据

`train_ddp.py` 的 scheduler builder：

- `warmup_steps_cfg = int(... warmup_steps ...)`
- `warmup_ratio = float(... warmup_ratio ...)`
- `warmup_steps = warmup_steps_cfg if warmup_steps_cfg > 0 else int(total_steps * warmup_ratio)`
  - `vendor/transchrombp/transchrombp/training/train_ddp.py:1338-1341`

也就是说，当前 strict-compare config 里：

- `warmup_steps: 0`
- `warmup_ratio: 0.06`

实际语义就是：

- **启用 6% 的 warmup**

这不是“可能有 warmup”，而是**确定有 warmup**。

### 正确表述

这里应改写为：

- “TransChromBP 明确使用 warmup，而 ChromBPNet 没有”
- 这属于 `B` 中的优化差异，不需要再单列成“代码语义待确认”

---

## E. JSD 计算微差

### 结论

**基本成立，影响低。**

### 证据

- `ChromBPNet` 用 `scipy.spatial.distance.jensenshannon(...)`
  - `chrombpnet/training/metrics.py:70`
- `TransChromBP` 手写 KL 后再 `sqrt`
  - `vendor/transchrombp/transchrombp/training/train_ddp.py:239-246`
- 两边最终都落在 `sqrt(JSD)` 这个量上

### 判断

对 peak 区域选择指标，这种差别更接近“实现细节”，不是主风险。  
我同意它应该被记录，但不认为它会改变当前 strict-compare 的 go/no-go 决策。

---

## 更新后的建议优先级

1. 先修 / 先控：
   - official selector 泄漏
   - region 集合统一
2. 再验证：
   - official multi-GPU global batch runtime smoke
3. 再写口径：
   - system-vs-system 的 optimizer / schedule / bf16 limitation
4. 不再单独追：
   - split mismatch（已证伪）
   - `patience=0` 语义（已证伪）
   - warmup 歧义（已证伪，实际就是有 warmup）
