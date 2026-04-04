# Tutorial L3 Launch 复核（Codex，2026-03-29）

## Findings

### P1: 当前已启动的任务不是完整的 L3 shared-region compare，只是 TransChromBP 单侧 shared-region 重训

按执行稿定义，L3 需要在 **官方 controlled** 和 **TransChromBP controlled** 两侧同时锁死同一份 `shared_filtered_peaks.bed` / `shared_filtered_nonpeaks.bed` / `shared_unstranded.bw`。[strict_chrombpnet_official_comparison_execution_20260327.md](/home/zhengwei/project/python/chromBPNet/docs/plan/strict_chrombpnet_official_comparison_execution_20260327.md#L103)

但当前仓库里的官方 tutorial wrapper 仍然只接受：

- `--mode`
- `--work-root`
- `--seed`
- `--folds`
- `--gpus`
- `--max-parallel`
- `--name`

它**没有**实现执行稿要求的：

- `--shared-peaks`
- `--shared-nonpeaks`
- `--run-suffix`

见 [run_tutorial_strict_compare_official.sh](/home/zhengwei/project/python/chromBPNet/scripts/paper_aligned_repro/run_tutorial_strict_compare_official.sh#L36) 和执行稿的 [L3 启动命令](/home/zhengwei/project/python/chromBPNet/docs/plan/strict_chrombpnet_official_comparison_execution_20260327.md#L380)。

进一步检查 `scripts/paper_aligned_repro/` 当前代码也没有任何 `shared-peaks` / `shared-nonpeaks` / `run-suffix` 接口痕迹；这说明问题不只是 wrapper 没转发，而是官方 L3 启动链路本身还未实现。

更关键的是，当前底层脚本 [run_paper_aligned_fast_1seed.sh](/home/zhengwei/project/python/chromBPNet/scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh#L195) 会在每个 fold 内固定执行 `chrombpnet prep nonpeaks`，然后把新生成的 `nonpeaks_negatives.bed` 继续传给后续 bias/chrombpnet train 与 predict。[run_paper_aligned_fast_1seed.sh](/home/zhengwei/project/python/chromBPNet/scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh#L205) [run_paper_aligned_fast_1seed.sh](/home/zhengwei/project/python/chromBPNet/scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh#L233) [run_paper_aligned_fast_1seed.sh](/home/zhengwei/project/python/chromBPNet/scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh#L279)

这意味着就算现在给 wrapper 补了参数，官方臂也仍然不能真正“直接吃外部 frozen shared nonpeaks”；要把 L3 做成定义中的 shared-region compare，至少还需要改到底层训练脚本，让它支持跳过 `prep nonpeaks` 并显式使用外部给定的 shared peaks/nonpeaks。

因此，Claude 当前写的“官方侧不需重训，L2 结果即 L3 baseline”并不符合仓库中原本定义的 L3 protocol。当前实际启动的是：

- `TransChromBP` 吃 shared-region 的新训练
- 对照仍然是 **official L2**

这可以作为有价值的 side experiment，但不能直接记成完整 `L3 shared-region compare`。

### P2: shared nonpeaks 的文档口径前后不一致，但当前 run 的实际输入来源已经可以核实

执行稿在 L3 方案 A 的原则描述里写的是：

- 从 official `auxiliary/` 取 `filtered.peaks.bed`
- 以及 `candidate_nonpeaks.bed`

但同一份执行稿在后面的实际命令又改成复制：

- `filtered.peaks.bed`
- `filtered.nonpeaks.bed`

见 [strict_chrombpnet_official_comparison_execution_20260327.md](/home/zhengwei/project/python/chromBPNet/docs/plan/strict_chrombpnet_official_comparison_execution_20260327.md#L116) 与 [strict_chrombpnet_official_comparison_execution_20260327.md](/home/zhengwei/project/python/chromBPNet/docs/plan/strict_chrombpnet_official_comparison_execution_20260327.md#L361)。

而前面的 follow-up 审计又明确记录过：

- 官方 `candidate nonpeaks = 508,429`
- 自研 `filtered.nonpeaks = 267,175`

见 [strict_compare_followup_assessment_20260327.md](/home/zhengwei/project/python/chromBPNet/reports/strict_compare_followup_assessment_20260327.md#L36)。

但本轮进一步到 6000 实地核查后，当前 shared-region 文件已经可以确认：

- `shared_filtered_peaks.bed` 与 official controlled auxiliary 下的 `filtered.peaks.bed` **完全一致**
- `shared_filtered_nonpeaks.bed` 与 official controlled auxiliary 下的 `filtered.nonpeaks.bed` **完全一致**
- `shared_unstranded.bw` 与 official controlled auxiliary 下的 `data_unstranded.bw` **完全一致**

并且三者的行数/文件规模也对得上：

- peaks：`269,773`
- nonpeaks：`269,773`

所以，当前这条已经启动的 TransChromBP run 的 region 定义本身是清楚的；真正的问题变成了：

1. 执行稿前文仍保留了 `candidate_nonpeaks` 的旧写法，和当前实际落地不一致
2. 官方侧尚未具备按同样 shared files 重跑的能力

这会影响“能否称为完整 L3”，但**不再构成必须立刻停掉当前 TransChromBP run 的理由**。

### P2: 配置和汇报名称已经出现漂移，增加了操作失误风险

当前仓库里的实际 data config 名称是：

- `configs/data/data_tutorial_shared_region_L3.yaml`

训练 config 也确实指向了它。[train_tutorial_corrected_b_strict_compare_L3_6000.yaml](/home/zhengwei/project/python/chromBPNet/vendor/transchrombp/transchrombp/configs/train/train_tutorial_corrected_b_strict_compare_L3_6000.yaml#L40) [data_tutorial_shared_region_L3.yaml](/home/zhengwei/project/python/chromBPNet/vendor/transchrombp/transchrombp/configs/data/data_tutorial_shared_region_L3.yaml#L1)

此外，这份训练 config 里的默认 `run_name` 实际是 `tutorial_corrected_b_strict_compare_L3`，并不带 Claude 汇报里的 `_s42` 后缀；`train_ddp.py` 只有在命令行显式传 `--run-name` 时才会覆盖这个值。不过本轮到 6000 核查进程后，已确认实际 launcher 的确显式传了 `--run-name tutorial_corrected_b_strict_compare_L3_s42`，因此当前正在运行的目录和日志口径是一致的。日志路径实际为 `/data1/zhoujiazhen/bylw_atac/logs/tutorial_corrected_b_strict_compare_L3_s42.log`。

但前面方法学文档和 Claude 汇报里多次写成：

- `data_tutorial_shared_region_v1.yaml`

这类命名漂移本身不会让当前 run 立刻失败，但在远端部署、补 selector、复盘日志时容易导致“文档说的不是实际跑的那份 config”。

## 判断

- **任务方向合理**：在 locked follow-up 之后，优先推进 shared-region 方向是合理的。
- **当前执行形态不够严谨**：它更准确地说是“已正确定义输入的 TransChromBP 单侧 shared-region 重训”，而不是已经成立的完整 L3 protocol。
- **现阶段最大漏洞不是训练代码本身，而是 protocol 落地不完整**。

## 建议

1. 当前这条 run 可以继续跑，不需要因为上述问题立刻停掉。
2. 但结果出来后不要直接写成“L3 baseline 对比完成”。
3. 2026-03-29 晚些时候，Codex 已在本地补齐官方 shared-region 启动链路：
   - `run_tutorial_strict_compare_official.sh` 现已支持 `--run-suffix` / `--shared-peaks` / `--shared-nonpeaks` / `--shared-bigwig`
   - `run_paper_aligned_fast_1seed.sh` 现已支持跳过内部 `prep nonpeaks`，并在评估阶段复用输入 peaks 与共享 bigWig
4. 因此当前真正剩余的 blocker 已经从“代码缺口”收缩为“是否实际重跑 official L3 臂”。
5. 如果短期内不想重跑官方侧，那当前 run 最多只能写成：
   - `ours-on-shared-region vs official-L2 reference`
   - 不能写成完整 `shared-region system compare`
