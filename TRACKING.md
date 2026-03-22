# 项目追踪清单（工作清单 + 下载资源清单）

本文件用于减少支线过多导致的信息遗漏。每次开始新任务或汇报前，先查看并按最新状态更新本文件。

> 已完成 / 已合并事项已迁至 [TRACKING_archive.md](TRACKING_archive.md)。
> 第一节原则上只保留 `进行中` / `待处理` / `待验证` 条目；阶段性完成或已并入其它主线的事项转入 `TRACKING_archive.md`。

状态约定：
- `已完成`：可直接使用
- `进行中`：任务仍在跑或尚未收尾
- `待处理`：明确要做但还没开始
- `待验证`：文件/目录已在，但需要实际加载或运行验证

## 一、在做事情清单（实时更新）

| 事项 | 当前状态 | 当前结论 / 进度 | 下一步 | 关键路径 / 日志 |
|---|---|---|---|---|
| 本地 Git 主档案与远端运行目录分工收口 | 进行中 | 2026-03-21 已定稿：本地仓库 `/home/zhengwei/project/python/chromBPNet` 是唯一主档案，GitHub 管长期历史；6000/6002 仅作训练或部署副本；规则已写入 `DEVELOPMENT.md` | 继续把远端仍有效的 launcher/config/文档回收到本地，并收紧本地与 6000 的脏工作树 | `DEVELOPMENT.md`、`TRACKING.md`、`scripts/sync_project.sh` |
| Git 工作树清理与实验归档方案 | 进行中 | 已完成清理分流：低风险代码/文档/摘要优先版本化，`tmp_remote_edit/` 暂保留，重文件与噪音已纳入忽略；当前进入批次 A | 按 `docs/plan/git_cleanup_and_archive_plan_20260321.md` 先收口低风险文件，再单独处理 `tmp_remote_edit/` 与根目录杂项 | `docs/plan/git_cleanup_and_archive_plan_20260321.md`、`DEVELOPMENT.md`、`.gitignore` |
| TransChromBP validation 指标对齐改造（JSD/count\_r/MAE） | 进行中 | `tmp_remote_edit/transchrombp` 已与 6000 当前源码/配置重新同步并通过关键文件 `py_compile`；新口径训练器草案已补 `JSD/count_r/MAE` 统计，以及 `best_metric=peak.profile_target_jsd_full_mean` | 待 6000 训练窗口释放后原子写回 `train_ddp.py` 与新 config，并做 `py_compile + dry-run` | `tmp_remote_edit/transchrombp/training/train_ddp.py`、`tmp_remote_edit/transchrombp/configs/train/train_tutorial_teacher_v2_main_profile_select.yaml` |
| 6002 上 A6000 历史数据集缓存补齐 | 待验证 | 2026-03-21 13:06 CST 已确认 6002 上 `GM12878/K562` 的 `merged.bam(.bai)`、`merged_unstranded.bw` 与 `nonpeaks_tutorial_folds.bed` 均已落盘；`chrombpnet_dataset_prep_6002_20260321_012824.log` 已在 `02:12:52` 收口，`HeLa-S3` / `ATAC fastq` 仍未推进 | 抽样核对关键文件可读性并在独立实验里实际用一次；`HeLa-S3` / `ATAC fastq` 继续排后 | `scripts/start_6002_dataset_cache_downloads.sh`、`/home/zhengwei/bylw_atac/logs/dataset_cache_downloads_20260318_022755.log` |
| 6000/6002 ChromBPNet 数据集 CPU 预处理（GM12878/K562） | 待验证 | 2026-03-21 13:06 CST 已确认两边后台作业均已结束：6000 日志在 `02:41:06`、6002 日志在 `02:12:52` 写出 `all requested datasets prepared`；两边 `GM12878/K562` 的 `merged.bam(.bai)`、`merged_unstranded.bw`、`nonpeaks_tutorial_folds.bed` 均已生成 | 用实际 loader / smoke run 抽样验证一轮后归档，再决定是否继续补 `HeLa-S3` / `ATAC fastq` | `scripts/run_remote_chrombpnet_dataset_prep.sh`、`/data1/zhoujiazhen/bylw_atac/logs/chrombpnet_dataset_prep_6000_20260321_012824.log`、`/home/zhengwei/bylw_atac/logs/chrombpnet_dataset_prep_6002_20260321_012824.log` |
| 6002 Genos 无关独立实验线（GM12878/K562） | 进行中 | `G_s42=profref` 已收口；`GM12878` 五个关键文件、data/train config 与 launcher 均已齐全。`gm12878_v2_6002_smoke` 的 `Invalid interval bounds` 根因已定位为 nonpeak BED 含 13 条 `bigWig` 缺失染色体记录（`chrM`×1、`chrUn_KI270302v1`×1、`chrUn_GL000226v1`×4、`chrY_KI270740v1_random`×7），并已在本地主档案与 6002 运行目录的 `real_data.py` 加入按 `bigWig` 染色体集合的过滤/防护；`GM12878/K562` 的 6002 data config fallback `sampling.nonpeak_ratio` 也已对齐到 `0.1`。修复后远端已通过 `py_compile + dataset 构造 + 多样本 __getitem__` 验证。`gm12878_v2_6002_smoke_retry_20260322_003645` 已于 2026-03-22 03:56 CST 在 `epoch 7` early-stop 收口，`best epoch=3`、`best peak.profile_target_jsd_full_mean=0.4351`，最后一次验证为 `val:peak profile_target_jsd_full_mean=0.4371`、`count_pearson_full=0.7970`。随后 `k562_v2_6002_smoke_20260322_161602` 已于 17:55 CST 按 `max_epochs=10` 正常收口；`best epoch=9`、`best peak.profile_target_jsd_full_mean=0.6331`，最后一次验证为 `val:peak profile_target_jsd_full_mean=0.6354`、`count_pearson_full=0.8411`。当前 pilot20 已按既定队列在 2026-03-22 19:17 CST 于 6002 后台启动；截至 2026-03-22 23:19 CST，主训练进程 `PID 33613` 仍在跑，`nvidia-smi pmon` 显示 GPU0 `sm≈97%`、显存约 `1.9 GiB`，标准输出已推进到 `epoch 9 step 15220/31482`。`epoch_metrics.jsonl` 已写到 `epoch 8`，当前仍是 `best epoch=3`、`best peak.profile_target_jsd_full_mean=0.4352`；最近一次验证（`epoch 8`）为 `val:peak profile_target_jsd_full_mean=0.4404`、`count_pearson_full=0.7771`，`best.pt` 最近更新时间仍停在 `20:43 CST`。为把 6002 总收口时间压到 2026-03-23 白天，6002 上后台挂起的 `queue_after_gm12878_pilot20_to_k562_single_20260322.sh`（PID `38088`）仍在守候；其日志已持续写到 `2026-03-22 23:18 CST` 的 `current run still active`，会在当前 run 结束且不早于 `2026-03-23 06:00 CST` 时自动启动 `train_k562_v2_6002_single.yaml` 对应的 `K562` 单 seed 长跑。启动期日志里的 `RuntimeWarning` 仅表示按 `bigWig` 染色体集合过滤了 `13` 条 `GM12878` nonpeak 不匹配记录（`chrM`、`chrUn_GL000226v1`、`chrUn_KI270302v1`、`chrY_KI270740v1_random`），训练链路未被中断 | 继续盯 `gm12878_v2_6002_pilot20_20260322_191746` 是否在 `epoch 10+` 或收尾前刷新 `best.pt`；若夜里提前结束，明早确认 queue 是否在 `06:00 CST` 后按预期接棒启动 `K562`，并核对其 launch/train log | `docs/plan/6002_independent_experiment_queue_20260321.md`、`tmp_remote_edit/TransChromBP/scripts/launch_v2fix_6002_single_dataset.sh`、`tmp_remote_edit/TransChromBP/scripts/queue_after_gm12878_pilot20_to_k562_single_20260322.sh`、`tmp_remote_edit/transchrombp/data/real_data.py`、`tmp_remote_edit/transchrombp/configs/data/data_gm12878_canonical_v1_6002.yaml`、`tmp_remote_edit/transchrombp/configs/data/data_k562_canonical_v1_6002.yaml`、`tmp_remote_edit/transchrombp/configs/train/train_gm12878_v2_6002_pilot20.yaml`、`tmp_remote_edit/transchrombp/configs/train/train_k562_v2_6002_single.yaml`、`/home/zhengwei/bylw_atac/logs/queue_after_gm12878_pilot20_to_k562_single_20260322.sh`、`/home/zhengwei/bylw_atac/logs/queue_after_gm12878_pilot20_to_k562_single_20260322.log`、`/home/zhengwei/bylw_atac/logs/gm12878_v2_6002_smoke_retry_20260322_003645.log`、`/home/zhengwei/bylw_atac/logs/gm12878_v2_6002_smoke_retry_20260322_003645.launch.log`、`/home/zhengwei/bylw_atac/logs/k562_v2_6002_smoke_20260322_161602.log`、`/home/zhengwei/bylw_atac/logs/k562_v2_6002_smoke_20260322_161602.launch.log`、`/home/zhengwei/bylw_atac/logs/gm12878_v2_6002_pilot20_20260322_191746.log`、`/home/zhengwei/bylw_atac/logs/gm12878_v2_6002_pilot20_20260322_191746.launch.log` |
| 6000 A6000 vs 6002 RTX3080 实时训练速度粗对比 | 进行中 | 历史 A6000 链路评估与分析已收口，结论保持：`V2-full` 相比 `V2-noTF` 在 held-out peak 上稳定更优，且未见 bias leakage 回潮；本条目前只保留为跨机器对照入口 | 需要时把 6000 `ablation_tf_20260318` 与 6002 `Bias-Safe B1-B4` 汇成统一对照表 | `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/reports/ablation_tf_20260318/summary_table.csv`、`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/reports/ablation_tf_20260318/ablation_summary.txt` |
| v2fix 读出头实验设计与待部署代码审查 | 待处理 | 2026-03-21 13:06 CST：6000 `F_s42=attnpool` 已在 `epoch 40/40` 收口，最终 `val:peak JSD=0.3318`、`count_r=0.7991`；6002 `G_s42=profref` 已于 `11:16 CST` 在 `epoch 25` early-stop 收口，最后一次验证为 `val:peak JSD=0.4502`、`count_r=0.8028`，RTX 3080 已释放 | 整理 `F/G` 最终指标与 checkpoint，对照后决定是否还需要追加 readout 头方案或把结论并入 `40 epoch` 复盘 | `/data1/zhoujiazhen/bylw_atac/logs/v2fix_20260321_attnpool_s42.log`、`/home/zhengwei/bylw_atac/logs/v2fix_20260321_profref_s42_launch.log` |
| V2 代码改进消融实验（修订版） | 待处理 | 第零批与第一批单 seed 已收口到 test 口径；`freeze_s42` 与 `cpool_s42` 近乎平手。现 `F_s42` 已完成、`G_s42` 也已于 2026-03-21 11:16 CST early-stop 收口，因此 V2fix 运行阶段已结束，当前进入结果汇总与是否继续扩 seed 的判断 | 汇总 `F_s42/G_s42` 的 epoch/test 指标，再决定是否扩多 seed；并将 6000/6002 资源转给 Genos 与独立数据集实验 | `docs/plan/v2_code_improvement_ablations.md`、`docs/plan/experiment_arrangement_20260321.md` |
| 40 epoch 尾部收益复盘（F_s42/G_s42 收口后） | 待处理 | 已约定不默认继续加 epoch；当前 `F_s42` 已在 `epoch 40` 收口且 `best epoch=40`，`G_s42` 也已在 `epoch 25` early-stop 收口，因此现在可以统一复盘 `30 -> 40` 的 held-out 收益、末段验证曲线，以及继续拉长训练对 Genos 窗口的机会成本 | 汇总 `epoch_metrics`、best checkpoint 与 held-out test，输出是否补 `50 epoch` 的结论 | `docs/plan/v2_code_improvement_ablations.md`、`docs/plan/v2fix_readout_head_experiment_design.md` |
| TransChromBP tutorial 真实数据基线重训（bs16 + peak-best） | 待处理 | 2026-03-15 run 曾在 `epoch 1 step 10120/13491` 被 SIGTERM；当前已被新的 `bias_pretrain -> main` 主线替代 | 仅在需要“无预训练 bias”对照时再恢复 | `logs/tutorial_real_bs16_peakbest_20260315_021233.log` |
| GM12878 自定义 ChromBPNet 复现（全流程） | 待处理 | 误启动产物已清理；该数据集当前仅作为 TransChromBP 输入，不再单独扩 ChromBPNet 复现 | 如后续需求变化再恢复独立复现线 | `/data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/GM12878/` |
| Genos-1.2B 本地接入准备 | 进行中 | 本地 `Genos-1.2B` 已跑通加载、最小前向、embedding、next-token 和 generate；当前剩余阻塞只剩 HF 连通性、官方 Docker 兜底镜像，以及许可证口径待确认 | 网络恢复后重跑 `hf download` 并继续拉 Docker 镜像；现有 host 环境和使用记录可直接复用 | `/data1/zhoujiazhen/bylw_atac/foundation_models/genos/Genos-1.2B`、`docs/env/transchrombp_genos_env.md` |
| Genos 融合实验设计与 Phase 0 快速验证 | 进行中 | Phase 0 结论稳定：唯一候选层仍是 `layer_6`，Phase 1 只能做最小 frozen pilot。2026-03-21 已完成起跑前脚本/配置修正并在 6000 重启主线；`G1` quick batch sweep 也已收口，`bs=8/12/16/20` 全部成功，峰值显存约 `5.4/6.4/7.2/8.1 GiB`，50-step walltime 为 `35/47/58/71 s`。截至 2026-03-23 00:06 CST，`G0(batch=20)` 已于 2026-03-21 18:52 CST 在 `epoch 20/20` 收口，最后一次 `val:peak profile_target_jsd_full_mean=0.3337`、`count_pearson_full=0.7998`，`best epoch=20`。`G1(batch=20)` 仍在 GPU1 上运行，标准输出已推进到 `epoch 9 step 5160/11809`、`loss=16.97586`、`lr=3.546e-04`；进程 `PID 1825216` 已运行约 `34.8 h`，`ps` 显示 `%CPU≈112`，`nvidia-smi pmon` 观察到 GPU1 `sm≈92%`、显存约 `8.1 GiB`。`G1` 的 `epoch_metrics.jsonl` 仍只写到 `epoch 8`，最近一次 validation 还是 `epoch 5`：`val:peak profile_target_jsd_full_mean=0.3420`、`count_pearson_full=0.7898`，`best.pt` 仍停在 `2026-03-22 09:05 CST`。`G2(batch=20)` 也在 GPU0 上正常推进，标准输出已到 `epoch 6 step 5660/11809`、`loss=17.02354`、`lr=4.449e-04`；进程 `PID 1889714` 已运行约 `23.3 h`，`ps` 显示 `%CPU≈112`，`nvidia-smi pmon` 观察到 GPU0 `sm≈94%`、显存约 `8.1 GiB`。`G2` 的 `epoch_metrics.jsonl` 仍只写到 `epoch 5`，当前最好验证依旧是 `val:peak profile_target_jsd_full_mean=0.3414`、`count_pearson_full=0.7955`，`best.pt` 更新时间仍为 `2026-03-22 22:04 CST`。当前两张 A6000 只见 `G1/G2` 自身占用，且到午夜为止仍无新 validation 节点落盘。cached-fusion 实现已做一轮工程 review：已修 data-config 解析与 `manifest/input_len/split/features` 校验、probe 改为 `5-fold OOF`、launcher 改为相对路径+自动发现 `G0 best.pt`，并在 `train_ddp.py` 加入 cached mode 缺 cache 时的硬报错；本轮又补上了训练期 `step_time_total/data_to_device/genos_fetch/forward_backward/optimizer` 与 batch 级 `FiLM/count injection` 使用统计，还补齐了缺失的 `configs/data/data_tutorial_canonical_v1.yaml` 与 cached pilot 脚本的导入路径/前置校验。`tmp_remote_edit/transchrombp` 相关关键文件现已形成本地 Git 专项快照，版本管理收口到可追溯状态 | 在远端按最终执行版先跑 `build_genos_summary_cache.py --dry_run`，确认 `train_regions/train_epoch_regions/valid_regions` 与 cache 预算；随后再跑 probe 与 `P0/P1/P2` | `docs/plan/genos_phase1_rebase_checklist.md`、`docs/plan/genos_phase1_minimal_integration_plan.md`、`docs/plan/experiment_arrangement_20260321.md`、`docs/plan/genos_roi_pivot_experiment_design_20260323.md`、`docs/plan/genos_cached_fusion_experiment_design_20260323.md`、`docs/plan/genos_cached_fusion_refined_experiment_design_20260323.md`、`docs/plan/genos_cached_fusion_final_execution_plan_20260323.md`、`tmp_remote_edit/transchrombp/scripts/run_genos_pilot.sh`、`tmp_remote_edit/transchrombp/scripts/run_genos_cached_pilot.sh`、`tmp_remote_edit/transchrombp/scripts/build_genos_summary_cache.py`、`tmp_remote_edit/transchrombp/scripts/run_genos_summary_probe.py`、`tmp_remote_edit/transchrombp/training/train_ddp.py`、`tmp_remote_edit/transchrombp/data/real_data.py`、`tmp_remote_edit/transchrombp/models/transchrombp.py`、`tmp_remote_edit/transchrombp/models/genos_adapter.py`、`tmp_remote_edit/transchrombp/configs/data/data_tutorial_canonical_v1.yaml`、`/data1/zhoujiazhen/bylw_atac/logs/genos_G0_20260321_104549.launch.log`、`/data1/zhoujiazhen/bylw_atac/logs/genos_20260321_baseline_s42.log`、`/data1/zhoujiazhen/bylw_atac/logs/genos_G1_20260321_launch.log`、`/data1/zhoujiazhen/bylw_atac/logs/genos_20260321_gate_s42.log`、`/data1/zhoujiazhen/bylw_atac/logs/genos_G2_20260322_004704.launch.log`、`/data1/zhoujiazhen/bylw_atac/logs/genos_20260322_mean_s42.log`、`/data1/zhoujiazhen/bylw_atac/logs/genos_batch_sweep_20260321_100046.log`、`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/logs/genos_batch_sweep_20260321_100046.tsv` |
| foundation_models 下载完整性验证 | 进行中 | `Genos-1.2B` 已完成本地完整性与真实使用验证；当前只剩 HF 可达性恢复后的同步复查，以及其它 foundation model 的同类核验 | 网络恢复后重跑 `hf download`，再继续复核其他模型 | `/data1/zhoujiazhen/bylw_atac/foundation_models/logs/` |
| nanochat / modded-nanogpt 学习环境准备 | 待处理 | 两个仓库已克隆，`modded-nanogpt` 仍是浅克隆 | 确认是否需要独立环境并安装依赖 | `/data1/zhoujiazhen/bylw_atac/nanochat`、`/data1/zhoujiazhen/bylw_atac/modded-nanogpt` |

## 二、下载资源清单（模型与数据）

说明：
- 这里记录"下载/准备了什么、从哪来、放在哪、用于什么"。
- 复现产物（如训练好的 `chrombpnet.h5`）虽不是"下载"，也记录在清单里，便于后续直接复用。

### A. 参考基因组与教程数据（ChromBPNet）

| 名称 | 服务器路径 | 来源（在哪里下载/获得） | 用途 | 当前状态 |
|---|---|---|---|---|
| hg38 参考基因组 | `/data1/zhoujiazhen/bylw_atac/hg38.fa` | 本地工作目录已有（教程中复用） | ChromBPNet 训练/解释、后续自定义模型数据处理 | 已完成 |
| hg38.fa 索引 | `/data1/zhoujiazhen/bylw_atac/hg38.fa.fai` | 本地工作目录已有 | 支持按染色体随机访问 | 已完成 |
| hg38 chrom sizes | `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/hg38.chrom.sizes` | 教程 step1 从 ENCODE 下载 | `chrombpnet_makebigwig`、`prep nonpeaks`、`prep splits` 的必需输入 | 已完成 |
| ENCODE blacklist | `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/blacklist.bed.gz` | 教程 step1 从 ENCODE 下载 `ENCFF356LFX` | 生成 `nonpeaks` 时作为排除区 | 已完成 |
| 教程数据目录 | `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial` | `chromBPNet` tutorial workflow | 官方教程复现基线 | 已完成 |
| 教程 folds 配置 | `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/data/folds.json` | 教程流程生成/准备 | 染色体划分（train/valid/test） | 已完成 |
| 预生成 genomewide GC bins | `/data1/zhoujiazhen/bylw_atac/chrombpnet_refs/genomewide_gc_hg38_stride_1000_inputlen_2114.bed` | 本地预生成 | 频繁制备 nonpeaks 时节省时间 | 待处理 |

### A2. 6002 直连官方缓存（TransChromBP 准备中）

| 名称 | 服务器路径 | 来源 | 用途 | 当前状态 |
|---|---|---|---|---|
| 6002 参考基因组与基础 refs | `/home/zhengwei/bylw_atac/hg38.fa` 等 | ENCODE 公开下载 | 6002 上共用参考输入 | 进行中（hg38.fa 已重建完成） |
| 6002 教程原始 peak 输入 | `/home/zhengwei/bylw_atac/chrombpnet_tutorial/data/overlap.bed.gz` | ENCODE `ENCFF333TAT` | tutorial preprocessing 输入 | 已完成 |
| 6002 示例 BAM/peaks/nonpeaks | `/home/zhengwei/bylw_atac/chrombpnet_tutorial/data/ENCSR868FGK_merged.bam` | `chrombpnet_data` 公共示例 | 快速 smoke test | 已完成（文件在，但 BAM 不完整，不用于 tutorial 正式预处理） |
| 6002 tutorial raw BAM（rep1/2/3 + merged） | `/home/zhengwei/bylw_atac/chrombpnet_tutorial/data/rep{1,2,3}.bam`、`merged.bam` | ENCODE tutorial step1 原始 BAM | 复现 A6000 tutorial preprocessing / 生成 `merged_unstranded.bw` | 已完成 |
| 6002 tutorial bigWig | `/home/zhengwei/bylw_atac/TransChromBP/outputs/preprocessing/tutorial_canonical_v1/step2_bigwig/merged_unstranded.bw` | 6002 本地由 `rep1/2/3 -> merged.bam -> bedGraphToBigWig` 生成 | 单卡 tutorial 训练直接输入 | 已完成 |
| 6002 folds / bias model 公共缓存 | `/home/zhengwei/bylw_atac/chrombpnet_tutorial/data/folds.json`、`bias.h5` 等 | Zenodo 小文件改由 6000 补齐 | 快速复用 split 与预训练 bias model | 已完成 |

### B. ChromBPNet 复现产物（可直接复用）

| 名称 | 服务器路径 | 来源（如何得到） | 用途 | 当前状态 |
|---|---|---|---|---|
| 自定义训练可复用 bias 模型 | `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/bias/models/bias.h5` | 教程 bias 训练（step4） | `chrombpnet pipeline -b` 输入候选 | 已完成 |
| ChromBPNet 主模型 | `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/chrombpnet/chrombpnet.h5` | 教程训练（step6b 双卡） | 教程评估、解释 | 已完成 |
| ChromBPNet 无 bias 模型 | `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/chrombpnet/chrombpnet_nobias.h5` | 教程训练自动导出 | 解释和对比分析 | 已完成 |
| 缩放后 bias 模型 | `/data1/zhoujiazhen/bylw_atac/chrombpnet_tutorial/outputs/chrombpnet/models/bias_model_scaled.h5` | 教程 bias 训练与预处理 | ChromBPNet 因子化训练输入 | 已完成 |
| SHAP（counts/profile） | `.../evaluation/corrected.{counts,profile}_scores.h5` | 双卡并行解释脚本 | TF-MoDISco / motif 分析 | 已完成 |
| MoDISco 结果+报告（counts/profile） | `.../evaluation/modisco_{results,reports}_{counts,profile}/` | MoDISco 运行 | motif 聚类与可视化 | 已完成 |

### C. 自定义 ChromBPNet 训练候选数据

| 名称 | 服务器路径 | 来源（在哪里下载） | 用途 | 当前状态 |
|---|---|---|---|---|
| `chrombpnet_datasets/GM12878` | `/data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/GM12878` | ENCODE/BAM + overlap peaks | 后续自定义训练 | 已完成（rep1/2/3.bam + overlap.bed.gz） |
| `chrombpnet_datasets/K562` | `/data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/K562` | ENCODE/BAM + overlap peaks | 后续自定义训练 | 已完成（rep1/2.bam + overlap.bed.gz） |
| `chrombpnet_datasets/HeLa-S3` | `/data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/HeLa-S3` | SRA `prefetch` | 需自行对齐/建 bigwig/出 peaks | 部分完成（.sra/.sralite 在，未见 BAM） |
| `TransChromBP` | `/data1/zhoujiazhen/bylw_atac/TransChromBP` | 手动创建 | 自研模型统一输出目录 | 已完成 |
| `6002/chrombpnet_datasets/GM12878` | `/home/zhengwei/bylw_atac/chrombpnet_datasets/GM12878` | ENCODE 公开下载 | 6002 缓存 | 进行中 |
| `6002/chrombpnet_datasets/K562` | `/home/zhengwei/bylw_atac/chrombpnet_datasets/K562` | ENCODE 公开下载 | 6002 缓存 | 进行中 |

### D. 通用 ATAC 数据（自建数据池）

| 名称 | 服务器路径 | 来源 | 用途 | 当前状态 |
|---|---|---|---|---|
| `ATACseq/peaks` | `/data1/zhoujiazhen/bylw_atac/ATACseq/peaks` | `download_atac_datasets.sh` | 候选 peak 集合 | 已完成（~7 个峰文件） |
| `ATACseq/fastq` | `/data1/zhoujiazhen/bylw_atac/ATACseq/fastq` | SRA (`prefetch`/`fasterq-dump`) | 自主对齐 | 已完成（~18 个 FASTQ 文件） |
| `6002/ATACseq/peaks` | `/home/zhengwei/bylw_atac/ATACseq/peaks` | 同上 | 6002 缓存 | 待处理 |

### E. 基因组基础模型（foundation_models）

| 名称 | 服务器路径 | 来源 | 用途 | 当前状态 |
|---|---|---|---|---|
| AlphaGenome 代码+权重 | `/data1/zhoujiazhen/bylw_atac/foundation_models/alphagenome/` | GitHub + HuggingFace | AlphaGenome 推理/对比 | 代码已完成；权重待验证 |
| Genos (1.2B/10B/10B-v2) | `/data1/zhoujiazhen/bylw_atac/foundation_models/genos/` | HuggingFace + GitHub `BGI-HangzhouAI/Genos` | 基因组基础模型对比 | 进行中（1.2B 已完成关键文件核对、隔离环境安装、`pip check`、本地加载、GPU `flash_attention_2` 前向与真实 usage 验证；剩余主要是 HF 可达性与 Docker 兜底镜像） |
| Gengram-10B | `/data1/zhoujiazhen/bylw_atac/foundation_models/gengram/` | HuggingFace | 基因组基础模型对比 | 待验证 |
| AlphaGenome 专用环境 | `/data1/zhoujiazhen/bylw_atac/.mamba/envs/alphagenome` | 6000 本地创建 | Python>=3.10 环境 | 已完成 |

### F. Transformer 学习代码仓库

| 名称 | 服务器路径 | 来源 | 当前状态 |
|---|---|---|---|
| nanochat | `/data1/zhoujiazhen/bylw_atac/nanochat` | GitHub `karpathy/nanochat` | 已完成（仓库已克隆） |
| modded-nanogpt | `/data1/zhoujiazhen/bylw_atac/modded-nanogpt` | GitHub `KellerJordan/modded-nanogpt` | 已完成（浅克隆，依赖未装） |

## 三、已知问题与判断规则（避免重复踩坑）

- 自定义训练不要为 `bigWig` 单独安排"下载任务"；`chrombpnet pipeline` 会从 `bam/fragment/tagAlign` 现场生成。真正值得提前缓存的是 `hg38.fa(+.fai)`、`hg38.chrom.sizes`、`blacklist.bed.gz`、`folds.json`、可复用 `bias.h5`。
- AlphaGenome 适合做外部 black-box 对比，不适合把 API 输出拿来训练学生模型；优先做 `1000s` 级别小规模对照。
- 6000 上 `chrombpnet` 环境是 `Python 3.9.23`，与 AlphaGenome SDK 的 `Python >= 3.10` 不兼容；AlphaGenome 使用独立环境 `/data1/zhoujiazhen/bylw_atac/.mamba/envs/alphagenome`。
- `chrombpnet_deepshap`、`chrombpnet_convert_html_to_pdf` 在当前环境里可能不存在；优先使用 Python 模块入口。
- `bias_qc` / `chrombpnet_qc` 历史失败常见原因是环境未激活（`modisco: not found`）；必须先执行 `CLAUDE.md` 中的环境激活命令。
- `foundation_models` 的 HF 下载目录里的 `.lock` 文件不等于下载失败；是否可用以"最小加载验证"为准。
- Genos 官方 `requirements.txt` 当前把 `torchaudio` 写成了 `>=0.22.1,<0.23.0`；按原样安装会直接报 `No matching distribution found`，当前临时修正为 `torchaudio==2.7.1` 继续装，后续如官方更新应回归上游文件为准。
- Genos 安装 `flash-attn` 时，直接按 `pip install -r requirements.txt` 走 build isolation 会在构建前置阶段报 `ModuleNotFoundError: No module named 'torch'`；更稳的处理路径应是先装好 `torch`，再对 `flash-attn` 单独使用 `--no-build-isolation`。

## 四、2026-03-17: 开发环境与版本管理配置

### 完成事项
1. Git 仓库初始化：本地与远程服务器 (6000) 均已关联 GitHub 仓库 `rose8happy/TransChromBP`；`.gitignore` 排除 `*.h5, *.pt`；统一 Git 用户身份 `yangmeisuan <345687960@qq.com>`。
2. 自动化同步脚本：`scripts/sync_project.sh` 支持 `deploy` / `download_results`。
3. 创建了 `DEVELOPMENT.md` 详细说明开发工作流。
