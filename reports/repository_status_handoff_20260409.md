# 仓库现况总览（2026-04-09 21:30 CST）

## 一句话结论

截至 `2026-04-09 21:30-21:31 CST`，本地档案仓的真实状态已经不再是“`msdec_v1_s2 + skipprobe_wide` 双机并跑”，而是：

- 本地 canonical git 档案仓仍是 `/home/zhengwei/project/python/chromBPNet`
- 当前活跃实验只剩 `6000 / A6000x2` 上的 `NT v2 teacher-distill short10`
- `6002 / RTX 3080` 当前空闲，`U-Net-lite` 方案只完成到本地 worktree 级别，还未在远端实际起跑

---

## 1. 角色边界

### 1.1 本地仓库

- 本地仓库 `master` 是当前唯一有完整 git 历史的 canonical 档案仓。
- 当前状态：`upstream/master` 之上 **ahead 20 commits**，且工作树仍是明显脏状态：
  - tracked 修改 `23` 个文件
  - untracked 文件 `31` 个
- 已提交历史的主轴不是最近双机实验本身，而是：
  - `official chrombpnet` 外置化
  - 仓库档案结构收口
  - vendor snapshot / evaluation tooling

### 1.2 6000 运行仓

- 路径：`/data1/zhoujiazhen/bylw_atac/TransChromBP`
- 现状：目录下存在 `.git`，但 `git status` 直接返回 `No commits yet on main`
- 结论：它更像长期运行目录，而不是一个有可追溯提交历史的正式 git 仓库

### 1.3 6002 运行仓

- 路径：`/home/zhengwei/bylw_atac/TransChromBP`
- 现状：目录下 **没有 `.git`**
- 结论：`6002` 当前只是运行副本，不具备本地仓库级别的分支/提交可追溯性

---

## 2. 本地分支与 worktree 盘点

| 分支 / worktree | 基线 | 当前状态 | 主要职责 |
|---|---|---|---|
| `master` | `4c62096` | 脏工作树；docs / paper / foundation / multiscale / externalization 混在一起 | 主档案仓，当前也是默认汇总入口 |
| `dual-track-20260409` | `4c62096` | tracked 修改 `7`、untracked `26` | 最新实验 pivot：`teacher-distill + U-Net-lite` |
| `multiscale-decoder-probe-20260407` | `4c62096` | clean | `msdec/skipprobe` 旧 readout probe 的干净切片 |
| `autonomy/20260406-chrombpnet-externalization` | `909ebef` | untracked `1` | `official chrombpnet` 外置化收尾文档 |
| `autonomy/20260406-structure` | `12c7c8c` | tracked 修改 `4`、untracked `3` | `foundation cache contract` 抽象层重构 |

### 2.1 `master` 当前混合了哪些改动簇

- docs / tracking / paper：
  - `TRACKING.md`
  - `docs/plan/project_roadmap_20260330.md`
  - `reports/paper_claim_evidence_matrix_20260326.md`
  - `reports/transchrombp_paper_{cn,draft}_v1.tex`
- official ChromBPNet 外置化：
  - `tests/test_chrombpnet_official_externalization.sh`
  - `workflows/tutorial/step3_get_background_regions.sh`
  - `scripts/start_6000_chrombpnet_dataset_prep.sh`
- foundation / training contract：
  - `vendor/transchrombp/transchrombp/training/train_ddp.py`
  - `vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py`
  - `vendor/transchrombp/transchrombp/scripts/build_foundation_cache.py`
  - `vendor/transchrombp/transchrombp/models/foundation_adapter.py`
- readout probe：
  - `vendor/transchrombp/transchrombp/models/transchrombp.py`
  - `vendor/transchrombp/transchrombp/models/profile_decoder.py`
  - `vendor/transchrombp/transchrombp/scripts/run_multiscale_decoder_probe.sh`
  - `vendor/transchrombp/transchrombp/scripts/run_short10_no_foundation_control.sh`

### 2.2 `dual-track-20260409` 持有的新增代码

这条 worktree 才是当前活跃实验路线的真正代码入口，核心内容包括：

- `NT v2 teacher-distill`
  - `teacher_cache_export.py`
  - `transchrombp_teacher_v2_center_pool_ntv2_distill.yaml`
  - `train_tutorial_teacher_v2_ntv2_distill_{short10,full}.yaml`
  - `run_ntv2_teacher_distill.sh`
- `U-Net-lite readout`
  - `profile_decoder.py` 内 `UNetLiteProfileDecoder`
  - `transchrombp_teacher_v2_center_pool_unet_lite_v1.yaml`
  - `run_unet_lite_decoder_probe.sh`
- 对应验证：
  - `tests/test_teacher_cache_export.py`
  - `tests/test_ntv2_teacher_distill_launcher.py`
  - `tests/test_unet_lite_launcher.py`

### 2.3 `autonomy/20260406-structure` 的定位

这条 worktree 不是当前 live experiment 主线，而是基础设施整理线。它的主要价值在于把 foundation cache 相关的重复合同层抽成共享 helper：

- `vendor/transchrombp/transchrombp/utils/foundation_contract.py`
- `tests/test_foundation_contract.py`

当前这条线还没有并入 `master`。

---

## 3. 实验现况

### 3.1 已收口的旧线

#### 6000：`msdec_v1_s2`

- run：`teacher_v2_center_pool_msdec_v1_s2_s42_20260409_full1`
- 日志：`/data1/zhoujiazhen/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_msdec_v1_s2_s42_20260409_full1.log`
- 截至 `2026-04-09 21:30 CST` 的最新可见结果：
  - peak valid `profile_target_jsd_full_mean=0.3314`
  - peak valid `count_pearson_full=0.8022`
  - overall valid `profile_target_jsd_full_mean=0.3541`
  - overall valid `count_pearson_full=0.8085`
- 结论：
  - run 已完成，不再活跃
  - 相比 paper-facing `corrected B` 仍偏弱，因此 `msdec_v1` 不再作为默认扩线 family

#### 6002：`skipprobe_v1_wide`

- run：`teacher_v2_center_pool_skipprobe_v1_wide_short10_s42_6002_r1`
- 日志：`/home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_skipprobe_v1_wide_short10_s42_6002_r1.log`
- 截至 `2026-04-09 21:30 CST` 的最新可见结果：
  - peak valid `profile_target_jsd_full_mean=0.4494`
  - peak valid `count_pearson_full=0.8075`
  - overall valid `profile_target_jsd_full_mean=0.4692`
  - overall valid `count_pearson_full=0.8135`
- 结论：
  - run 已完成，不再活跃
  - 这条 `skipprobe` 高带宽 cheap-screen 没有形成值得保留的正向信号

### 3.2 当前活跃线

#### 6000：`NT v2 teacher-distill short10`

- run：`ntv2_teacher_distill_short10_s42_6000_20260409_r2`
- 机器：`6000 / A6000 x2`
- 证据：
  - `nvidia-smi` 显示两张 A6000 均有活跃 python 进程
  - 进程 PID：`862231`、`862232`
  - 命令行直接指向 `transchrombp.training.train_ddp` + `transchrombp_teacher_v2_center_pool_ntv2_distill.yaml`
  - 日志 `/data1/zhoujiazhen/bylw_atac/TransChromBP/logs/ntv2_teacher_distill_short10_s42_6000_20260409_r2.log` 在 `2026-04-09 21:30 CST` 已推进到 `epoch=1 step=3380/11809`
- 配置事实：
  - 本地 `dual-track-20260409` worktree 的 `train_tutorial_teacher_v2_ntv2_distill_short10.yaml` 明确写的是 `max_epochs: 10`
- 粗略估计：
  - 以当前 step 速度估算，训练部分大致还有 `~2` 小时量级
  - 若不早停，预计结束时间窗约为 `2026-04-09 23:45` 到 `2026-04-10 00:15 CST`

### 3.3 尚未真正起跑的线

#### 6002：`U-Net-lite readout`

- 代码、配置、launcher 已存在于 `dual-track-20260409` worktree
- 但截至 `2026-04-09 21:31 CST`：
  - `6002` 的 3080 已空闲
  - `ps -ef | grep transchrombp.training.train_ddp` 没有活跃训练进程
- 结论：
  - 这条线目前仍停留在“本地已实现、远端未起跑”的状态

---

## 4. 当前仍需注意的不一致

- `TRACKING.md` 已在本轮改成指向本报告，但旧的 `reports/session_handoff_multiscale_and_next_tasks_20260409.md` 仍保留着“`msdec_v1_s2` / `skipprobe_wide` 是 active run”的旧快照
- 真实现场已经切到：
  - `6000` 正在跑 `teacher-distill`
  - `6002` 当前空闲
- `dual-track-20260409` 的计划与 pivot 记录还没有完整回流到主档案仓，因此“最新实验路线”实际分散在：
  - 远端日志
  - worktree 内部文件
  - 主仓库旧 handoff 文档

---

## 5. 当前最值得做的事

1. 继续把 `TRACKING.md` 视为 live 入口，而把本报告当成当前仓库 / worktree / 双机状态的统一快照。
2. 等 `ntv2_teacher_distill_short10_s42_6000_20260409_r2` 收口后，立即刷新本报告与 `TRACKING.md`，完成一次 closeout：
   - short10 gate 是否通过
   - 是否值得升 full
   - `6002` 是否真的要起 `U-Net-lite`
3. 在实验结论稳定后，决定 `dual-track-20260409` 是否要回流到主档案仓；否则主仓库会继续处于“旧实验已归档、新实验只存在于 worktree”的半分叉状态。

---

## 6. 建议的阅读顺序

1. 本文件：回答“现在真实在发生什么”
2. `reports/project_plan_code_review_20260405.md`：回答“最近代码为什么会改成这样”
3. `reports/session_handoff_multiscale_and_next_tasks_20260409.md`：保留为旧双机排程快照
4. `dual-track-20260409` worktree 下两份文件：
   - `docs/plan/dual_track_ntv2_distill_unet_lite_20260409.md`
   - `reports/dual_track_experiment_pivot_20260409.md`
   用来回答“为什么从 `msdec/skipprobe` 切到 `teacher-distill + U-Net-lite`”
