# Strict Compare 代码审计（2026-03-27）

## 结论摘要

- `P0`：**成立**。ChromBPNet 官方 selector 当前确实在 `test` split 上逐 epoch 选 best，存在 test 泄漏。
- `P1`：**部分成立**。MirroredStrategy 下的实际 global batch 还没有做 runtime 实证；这是验证缺口，但目前还不是已证实 bug。
- `P1`：**不成立（已核实对齐）**。两侧 tutorial split 并不存在实际不一致；6000 上 `fold_0.json` 与 `chrombpnet_tutorial/data/folds.json` 已程序化比对为完全一致。
- `P2`：**成立**。`run_tutorial_strict_compare_official.sh` 当前仍会重复传 `--max-parallel`。
- `P3`：**成立**。`metrics_only` 模式下 `counts_metrics()` 仍会写出散点图 PNG。
- `P3`：**不成立**。TransChromBP 的 `early_stop_patience: 0` 在当前 trainer 里确实等价于禁用 early stop。

## 2026-03-27 13:22 CST 后续状态

这轮审计后，以下修复已经在本地主仓库落地，等待继续同步到 6000 并重启 official 臂：

- `ChromBPNet selector` 已新增 `--split`，默认改为 `valid`
- `predict.py` 已支持显式 `--split`
- `metrics_only` 模式下不再额外写 count 散点图
- `run_tutorial_strict_compare_official.sh` 已去掉重复 `--max-parallel`

尚未完成的只剩一项高优先验证：

- 对 official arm 做一次 multi-GPU runtime smoke，确认当前 TF/Keras 组合下 generator `batch_size=32` 的实际 global batch 语义

补充澄清：

- `tutorial_official_controlled_s42_20260327_013246` 的第一次失败，确实是一个会导致 official controlled arm **直接秒退** 的真实 wrapper bug，而不是单纯“代码看起来不优雅”。
- 更准确地说，致命问题是：wrapper 当时向下游同时传了 `--multi-gpu-train` 和不兼容的 `--max-parallel=2`，而下游脚本本身明确禁止这个组合。
- 与此不同，`--max-parallel` 的“重复传递”属于单独的代码卫生问题；它本身未必一定触发失败，但会让参数覆盖关系变脆弱、难以推断。
- 2026-03-27 后续 code review 已确认：Claude 的修复把下游从 `hard error` 改成了 `warn + auto-coerce to 1`，同时 tutorial wrapper 的调试输出补上了 `max_parallel`；两份脚本均已通过本地 `bash -n`。这属于合理的鲁棒性修复，没有引入新的明显行为风险。
- 2026-03-27 15:53 CST 已完成更快的 synthetic runtime smoke：新增 [debug_mirrored_sequence_batch.py](/home/zhengwei/project/python/TransChromBP/scripts/paper_aligned_repro/debug_mirrored_sequence_batch.py)，在 6000 的 `CUDA_VISIBLE_DEVICES=0,1` 下直接验证 `Keras Sequence + MirroredStrategy`。输出显示：
  - `Using MirroredStrategy with 2 replicas`
  - `[batch-debug] replica 0 input_shape [16 2114 4]`
  - `[batch-debug] replica 1 input_shape [16 2114 4]`
  因此，当前官方 multi-GPU 链路里 `batch_size=32` 的语义已经有 runtime 证据支持为 **global batch 32**，而不是每卡 32。

## 逐项判断

### 1. ChromBPNet selector 在 test 集上选 best

**成立，且是当前最严重的问题。**

证据：

- 官方 selector 只把 `fold_json` 传给 `predict.py`，没有任何 `split` 参数：[select_best_epoch.py](/home/zhengwei/project/python/TransChromBP/scripts/paper_aligned_repro/select_best_epoch.py#L66)
- `predict.py` 在入口里把评估 split 硬编码为 `mode="test"`：[predict.py](/home/zhengwei/project/python/TransChromBP/chrombpnet/training/predict.py#L106)

因此当前 official selector 的行为是：

1. 对每个 `chrombpnet.epoch_*.h5`
2. 在 `fold_json["test"]` 上跑评估
3. 用 test 指标挑 best epoch

这会直接破坏“统一外部 validation 选 best”的核心声称。

### 2. `--max-parallel` 重复传递

**成立。**

证据：

- `CMD` 初始数组已经写入 `--max-parallel "${MAX_PARALLEL}"`：[run_tutorial_strict_compare_official.sh](/home/zhengwei/project/python/TransChromBP/scripts/paper_aligned_repro/run_tutorial_strict_compare_official.sh#L90)
- 当 `MULTI_GPU_TRAIN=1` 时又追加一次 `--max-parallel 1`：[run_tutorial_strict_compare_official.sh](/home/zhengwei/project/python/TransChromBP/scripts/paper_aligned_repro/run_tutorial_strict_compare_official.sh#L113)

当前行为依赖“后面的值覆盖前面的值”，能跑通，但确实脆弱。

需要额外说明的是，第一次 official controlled run 的真正秒退原因并不是“重复传参”本身，而是 wrapper 在旧版本里把 `--multi-gpu-train` 与 `--max-parallel=2` 这一对不兼容参数一起传给了下游。失败日志中可以直接看到：

```text
ERROR: --multi-gpu-train requires --max-parallel=1
```

而下游约束就在 [run_paper_aligned_fast_1seed.sh](/home/zhengwei/project/python/TransChromBP/scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh#L143)。

### 3. MirroredStrategy 下 global batch 需要实证验证

**部分成立。**

我没有发现“已经确定 batch 被放大到 64”的证据，但你指出的验证缺口确实存在。

现有代码证据：

- ChromBPNet generator 直接用 `args.batch_size` 构造 batch：[initializers.py](/home/zhengwei/project/python/TransChromBP/chrombpnet/training/data_generators/initializers.py#L80)
- `Sequence.__len__` 和 `__getitem__` 都按这个 batch size 切数据：[batchgen_generator.py](/home/zhengwei/project/python/TransChromBP/chrombpnet/training/data_generators/batchgen_generator.py#L53), [batchgen_generator.py](/home/zhengwei/project/python/TransChromBP/chrombpnet/training/data_generators/batchgen_generator.py#L99)
- 多卡训练仅是在 `model.fit(...)` 外包了一层 `tf.distribute.MirroredStrategy()`：[train.py](/home/zhengwei/project/python/TransChromBP/chrombpnet/training/train.py#L95)

这说明 `32` 是 generator 产出的 batch 大小，但**还没有 runtime 证据**证明 TF 在当前环境下把它当作 global batch，而不是 per-replica batch。

2026-03-27 15:53 CST 的 synthetic smoke 已经把这件事从“推断”推进成了“runtime 证据”：

- `Visible GPUs: 2`
- `Using MirroredStrategy with 2 replicas`
- `[batch-debug] replica 0 input_shape [16 2114 4]`
- `[batch-debug] replica 1 input_shape [16 2114 4]`

所以更准确的当前结论是：

- “已经有 global batch=64 bug”目前**不成立**
- “当前 TF/Keras + Keras Sequence 组合下，official arm 的 `batch_size=32` 会被分发成每卡 `16`”已经**有实证支持**

### 4. 两侧数据 split 对齐未显式验证

**作为实际 mismatch 指控，不成立。**

我已在 6000 上直接比较：

- `/data1/zhoujiazhen/bylw_atac/chrombpnet_strict_compare/folds/fold_0.json`
- `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/folds.json`

结果：

- `equal True`
- `train 454 == 454`
- `valid 1 == 1`
- `test 1 == 1`

同时，自研 strict-compare config 明确指向 tutorial 的 `folds.json`：[data_tutorial_canonical_v1.yaml](/home/zhengwei/project/python/TransChromBP/vendor/transchrombp/transchrombp/configs/data/data_tutorial_canonical_v1.yaml#L1)

所以：

- “之前缺少显式 diff 这一步”这个批评是合理的
- 但“当前两侧 split 可能不一致”这件事，经核对后**不成立**

### 5. `counts_metrics()` 在 metrics_only 模式下仍可能写文件

**成立。**

证据：

- `predict.py` 在 `metrics_only=True` 时仍无条件调用 `counts_metrics()`：[predict.py](/home/zhengwei/project/python/TransChromBP/chrombpnet/training/predict.py#L119), [predict.py](/home/zhengwei/project/python/TransChromBP/chrombpnet/training/predict.py#L138), [predict.py](/home/zhengwei/project/python/TransChromBP/chrombpnet/training/predict.py#L156)
- `counts_metrics()` 内部无条件 `plt.savefig(...)`：[metrics.py](/home/zhengwei/project/python/TransChromBP/chrombpnet/training/metrics.py#L19)

这不会改变指标数值，但会在逐 epoch selector 里额外产出大量 PNG。

### 6. `early_stop_patience: 0` 的语义

**不成立。当前实现里 `0` 的确表示禁用 early stop。**

证据：

- strict-compare config 把 `early_stop_patience` 设为 `0`：[train_tutorial_corrected_b_strict_compare_6000.yaml](/home/zhengwei/project/python/TransChromBP/vendor/transchrombp/transchrombp/configs/train/train_tutorial_corrected_b_strict_compare_6000.yaml#L19)
- trainer 里只有当 `early_stop_patience > 0` 时才可能触发早停：[train_ddp.py](/home/zhengwei/project/python/TransChromBP/vendor/transchrombp/transchrombp/training/train_ddp.py#L1927)

因此这不是风险点。

## 建议优先级

1. 先修 `P0`：给官方 selector/predict 增加显式 `split=valid` 支持；在这之前不要使用 official selector 输出任何 best-epoch 结果。
2. 既然 global batch 语义已闭环，就可以恢复 official controlled arm。
3. 顺手修 `--max-parallel` 重复传参和 `metrics_only` 额外写 PNG。
4. split 对齐和 `patience=0` 可以从“怀疑项”移出。
