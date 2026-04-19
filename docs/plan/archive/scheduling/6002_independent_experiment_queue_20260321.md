# 6002 独立实验线执行清单（Genos 提优后的非 Genos 实验）

> 日期：2026-03-21  
> 适用机器：`6002 / RTX 3080 12GB`  
> 目标：在 `A6000` 优先承接 `Genos Phase 1` 的同时，让 `6002` 持续推进 **Genos 无关、且不依赖 A6000 checkpoint 同步** 的独立 `TransChromBP` 实验。

---

## 一、边界与原则

### 1.1 这条线负责什么

- 承接 `Genos` 无关的独立 `TransChromBP` 实验。
- 优先回答“换独立数据集后，当前主线结构是否仍然稳定有效”。
- 实验必须能在 `RTX 3080 12GB` 单卡上独立完成。

### 1.2 这条线不负责什么

- 不承担 `Genos` 训练或推理。
- 不承担 `A6000` 训练权重的默认从属评估。
- 不为了评估而频繁同步大 checkpoint。
- 不在 `3080` 上一次性串行排多个 `40ep` 长训练。

### 1.3 当前前提

- `6002` 当前仍在跑：`v2fix_20260321_profref_s42`
- `6002` 当前仍在跑低优先级 CPU 预处理：`GM12878/K562`
- 因此这份清单的默认执行时点是：
  - `G_s42` 收口后
  - 且目标数据集的 `prep_v1` 关键产物已齐全后

---

## 二、当前实地状态

### 2.1 现有可复用训练入口

- 启动脚本：
  - `scripts/launch_teacher_v2_6002_single.sh`
- 单卡 train config：
  - `configs/train/train_tutorial_teacher_v2_main_6002_single.yaml`
- 现有 profile-select train config：
  - `configs/train/train_ablation_v2_main_profile_select.yaml`
- 当前数据 config 模板：
  - `configs/data/data_tutorial_canonical_v1_6002.yaml`

### 2.2 当前数据落盘状态（2026-03-21 实查）

`GM12878`
- 已有：`rep*.bam`、`overlap.bed.gz`、`merged.bam`、`merged.bam.bai`
- 未确认完成：`prep_v1/bigwig/merged_unstranded.bw`
- 未确认完成：`prep_v1/background/tutorial_folds/nonpeaks_tutorial_folds.bed`

`K562`
- 已有：`rep*.bam`、`overlap.bed.gz`
- 尚未完成：`merged.bam`
- 尚未完成：`prep_v1/bigwig/merged_unstranded.bw`
- 尚未完成：`prep_v1/background/tutorial_folds/nonpeaks_tutorial_folds.bed`

### 2.3 直接结论

- `GM12878` 是 6002 收口后最适合接的第一条独立实验线。
- `K562` 目前只适合排在 `GM12878` 之后做 smoke，不适合立即上完整长跑。

---

## 三、推荐实验顺序

### 3.1 第一优先级：GM12878-only smoke

目标：
- 验证 `GM12878` 数据在 `6002` 上可以完整走通当前主线 `TransChromBP` 训练链路。
- 在不拉长 wall-clock 的前提下，先回答“独立数据集 + 单卡 3080”是否工程可用。

建议配置：
- backbone：默认先用 `v2fix_baseline.yaml`
- 如果 `F/G` 在启动前已经明确分出胜负：
  - 可把 smoke 的 model config 改成当前非 Genos 最优 readout 版本
- 训练长度：
  - `5-10 epoch`
- batch：
  - 从 `batch_size_per_gpu=8` 起步
- workers：
  - `num_workers=2`

预计耗时：
- 代码/配置落地：`0.5 天`
- smoke 运行：`2-5 小时`

go/no-go 条件：
- 训练可稳定进入 validation
- `best.pt`、`last.pt`、metrics JSON 正常写出
- 不出现明显的 OOM / data loader 卡死 / region 为空等链路错误
- peak 指标至少达到“可比较”水平，不要求一次性超过 tutorial 主线

### 3.2 第二优先级：GM12878-only 单 seed baseline

目标：
- 在 `GM12878-only` 上拿到一条和 tutorial 口径接近的单 seed 基线。

建议配置：
- model：
  - 默认 `v2fix_baseline.yaml`
  - 若 `F/G` 收口后有明确更优非 Genos 版本，可切换为赢家
- train：
  - 基于 `train_ablation_v2_main_profile_select.yaml`
  - 单卡改为 `batch_size_per_gpu=8`
  - `config_path` 改到新的 `GM12878` data config
  - `logging.output_dir` 保持 `6002` 本地路径

预计耗时：
- 正式 `40ep`：`19-22 小时`
- 若先做 `20ep` pilot：`10-12 小时`

建议判定：
- 先做 `10ep` 或 `20ep` 小 pilot 更稳
- 只有 smoke 指标和吞吐都正常，再决定是否上完整 `40ep`

### 3.3 第三优先级：K562-only smoke

目标：
- 只验证 `K562` 是否具备进入正式训练的最小工程条件。

预计耗时：
- 预处理收口后：`1-3 小时`

建议：
- 不直接上 `40ep`
- 先等 `GM12878` 跑出一条有效结果，再判断是否值得扩 `K562`

---

## 四、执行前置条件清单

### 4.1 不打断现有任务

- 保持 `v2fix_20260321_profref_s42` 自然收口
- 保持 `chrombpnet_dataset_prep_6002_20260321_012824` 自然收口

### 4.2 GM12878 启动前必须检查的文件

必须存在：

- `/home/zhengwei/bylw_atac/chrombpnet_datasets/GM12878/merged.bam`
- `/home/zhengwei/bylw_atac/chrombpnet_datasets/GM12878/merged.bam.bai`
- `/home/zhengwei/bylw_atac/chrombpnet_datasets/GM12878/prep_v1/inputs/overlap.bed`
- `/home/zhengwei/bylw_atac/chrombpnet_datasets/GM12878/prep_v1/bigwig/merged_unstranded.bw`
- `/home/zhengwei/bylw_atac/chrombpnet_datasets/GM12878/prep_v1/background/tutorial_folds/nonpeaks_tutorial_folds.bed`

若缺任一项：
- 不启动 `GM12878-only` 训练
- 先查 `chrombpnet_dataset_prep_6002_20260321_012824.log`

### 4.3 K562 启动前必须检查的文件

必须存在：

- `/home/zhengwei/bylw_atac/chrombpnet_datasets/K562/merged.bam`
- `/home/zhengwei/bylw_atac/chrombpnet_datasets/K562/merged.bam.bai`
- `/home/zhengwei/bylw_atac/chrombpnet_datasets/K562/prep_v1/inputs/overlap.bed`
- `/home/zhengwei/bylw_atac/chrombpnet_datasets/K562/prep_v1/bigwig/merged_unstranded.bw`
- `/home/zhengwei/bylw_atac/chrombpnet_datasets/K562/prep_v1/background/tutorial_folds/nonpeaks_tutorial_folds.bed`

---

## 五、需要补的配置与脚本

### 5.1 新增 data config

优先新增：

- `configs/data/data_gm12878_canonical_v1_6002.yaml`
- `configs/data/data_k562_canonical_v1_6002.yaml`

格式直接参照：

- `configs/data/data_tutorial_canonical_v1_6002.yaml`

其中只替换：

- `input.peaks_bed`
- `input.nonpeaks_bed`
- `input.bigwig`

其余窗口参数先保持：

- `input_len=2114`
- `output_len=1000`
- `max_jitter=500`

### 5.2 新增单卡 train config

建议新增：

- `configs/train/train_gm12878_v2_6002_smoke.yaml`
- `configs/train/train_gm12878_v2_6002_single.yaml`
- `configs/train/train_k562_v2_6002_smoke.yaml`

基线模板：

- `configs/train/train_ablation_v2_main_profile_select.yaml`

需要修改的字段：

- `seed`
- `max_epochs`
- `data.config_path`
- `data.batch_size_per_gpu=8`
- `data.num_workers=2`
- `logging.output_dir=/home/zhengwei/bylw_atac/TransChromBP/outputs`
- `logging.run_name`

### 5.3 新增启动脚本

建议新增一个通用单卡 launcher，例如：

- `scripts/launch_v2fix_6002_single_dataset.sh`

职责：

- 激活 `6002` 环境
- 设置 `PYTHONPATH=/home/zhengwei/bylw_atac/TransChromBP/src`
- 接收 `MODEL_CONFIG`、`TRAIN_CONFIG`、`RUN_NAME`
- 用 `torchrun --nproc_per_node=1` 或 `python -m transchrombp.training.train_ddp` 启动

备注：
- 现有 `launch_teacher_v2_6002_single.sh` 偏向 tutorial 的 bias->main 串联流程
- 对 `GM12878/K562-only` 更合适的是直接做单段 main train launcher

---

## 六、建议执行顺序（可直接照单推进）

### 6.1 阶段 A：等待当前运行收口

1. 等 `G_s42` 跑完
2. 等 `GM12878/K562` CPU 预处理至少完成 `GM12878`
3. 检查 `GM12878` 五个关键文件是否齐全

预计耗时：
- 取决于当前后台任务剩余时间

### 6.2 阶段 B：配置落地

1. 新建 `data_gm12878_canonical_v1_6002.yaml`
2. 新建 `train_gm12878_v2_6002_smoke.yaml`
3. 新建 `train_gm12878_v2_6002_single.yaml`
4. 新建 `launch_v2fix_6002_single_dataset.sh`
5. `py_compile` / `bash -n` 检查

预计耗时：
- `1-2 小时`

### 6.3 阶段 C：GM12878 smoke

1. 运行 `5-10 epoch` smoke
2. 记录：
  - 吞吐
  - 显存
  - validation 指标
  - 日志路径
3. 判断是否进入正式单 seed

预计耗时：
- `2-5 小时`

### 6.4 阶段 D：GM12878 单 seed

1. 若 smoke 通过，启动 `20ep` 或 `40ep`
2. 训练结束后立即做 held-out 评估
3. 把结果写回 `TRACKING.md`

预计耗时：
- `20ep`：`10-12 小时`
- `40ep`：`19-22 小时`

### 6.5 阶段 E：K562 smoke

1. 确认 `K562` 关键文件齐全
2. 复制 `GM12878` 的 smoke 配置改成 `K562`
3. 跑 `5-10 epoch` smoke

预计耗时：
- `1-3 小时`

---

## 七、推荐的具体排队顺序

```text
当前：
  6002 继续跑 v2fix_20260321_profref_s42
  6002 继续跑 chrombpnet_dataset_prep_6002_20260321_012824

下一步：
  1. GM12878-only smoke
  2. GM12878-only 单 seed baseline
  3. K562-only smoke

暂缓：
  4. K562-only 完整长跑
  5. 任何依赖 A6000 checkpoint 同步的评估流程
  6. 3080 上的 Genos 相关实验
```

---

## 八、与 A6000 主线的协同关系

- `6000 / A6000` 负责回答：`Genos` 是否对当前主线有净增益。
- `6002 / 3080` 负责回答：在 `Genos` 之外，独立数据集上的 `TransChromBP` 链路是否能稳定复用。
- 两条线并行，但互不依赖大 checkpoint 同步。
- 只有当 `Genos` 给出明确负结论，或 `GM12878` 结果显著优于 tutorial 主线时，才重新讨论两条线的资源再平衡。

---

## 九、相关文档

- `TRACKING.md`
- `docs/plan/archive/scheduling/experiment_arrangement_20260321.md`
- `docs/plan/archive/genos/genos_phase1_minimal_integration_plan.md`
