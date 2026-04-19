# 仓库现况总览（2026-04-10 10:10 CST）

> 新 charter、`TRACKING.md` 和 `docs/experiments/registry.md` / `runs.csv` 已分别承担规则源、live 入口和 canonical 索引；本文只负责快照，不再承担排程裁决角色。

## 一句话结论

截至 `2026-04-10 10:10 CST`，本地档案仓的真实状态已经切换成终态快照，而不是“`msdec_v1_s2 + skipprobe_wide` 双机并跑”那套旧口径：

- `6000` 的 A6000 formal gate `teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1` 已收口并判为 `fail`
- `6000` 的 AlphaGenome v2 sidecar `alphagenome_matched_raw_track_slice_v2_20260410` 已收口并判为 `pass`
- `6002` 的 `U-Net-lite v1` 已完成 `r4` 确认性复跑并收成 `no-go / stop`
- 当前没有 active 的 `6000` / `6002` training run；本报告只记录现场快照，不再给出“双机必须串行推进”的默认动作

---

## 1. 角色边界

### 1.1 本地仓库

- 本地仓库 `master` 是当前唯一有完整 git 历史的 canonical 档案仓。
- 当前状态：`master` 已重新回到 clean 状态，并保持为相对 `upstream/master` 的本地领先档案分支。
- 本轮已经用连续两条 docs 提交把原先混在主仓里的 docs / handoff / archive 内容收口进 canonical 档案仓。
- 当前实验 family / run 的 canonical 索引已经单独落到 `docs/experiments/registry.md` 与 `docs/experiments/runs.csv`，不再只散落在 `TRACKING.md`、handoff 与远端日志里。
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

## 2. 本地分支与 archive snapshot 盘点

> 下表保留的是 `2026-04-10` 当天本地可恢复的 branch/snapshot 现场，不等于当前仍挂载的 live worktree 集合。当前 mounted worktree 以 `TRACKING.md` 与 `docs/experiments/registry.md` 为准。

| 分支 / snapshot | 基线 | 当前状态 | 主要职责 |
|---|---|---|---|
| `master` | 本轮 docs closeout 后的最新状态 | clean | 主档案仓，承载 canonical docs / tracking / paper / reports |
| `dual-track-20260409` | `4ffc6e7` | clean | 近期 dual-track pivot 快照：`teacher-distill + U-Net-lite` |
| `multiscale-decoder-probe-20260407` | `4c62096` | clean | `msdec/skipprobe` 旧 readout probe 的干净切片 |
| `autonomy/20260406-chrombpnet-externalization` | `4aadffa` | clean | `official chrombpnet` 外置化收尾与工具链 bridge |
| `autonomy/20260406-structure` | `764b8c0` | clean | `foundation cache contract` 抽象层重构 |

### 2.1 `master` 当前保留什么

- 现在 `master` 已不再混着实验性代码与工具链 WIP。
- 本轮保留在 `master` 的内容只有 canonical docs / archive：
  - `TRACKING.md`
  - `docs/experiments/*`
  - `reports/repository_status_handoff_20260409.md`
  - `docs/plan/*.md`
  - `reports/*.md`
  - `reports/*.tex`
  - `docs/superpowers/` 下的历史 agent 设计 / 执行归档（当前 live 计划以 `docs/plan/` 为准）

### 2.2 `dual-track-20260409` 的历史 pivot 代码入口

这条 snapshot branch 已固定在 clean 提交 `4ffc6e7 (wip: snapshot dual track pivot)`，是近期 dual-track pivot 的历史代码入口。当前 master-side 的 canonical 说明见 `reports/dual_track_pivot_snapshot_20260409.md`，核心内容包括：

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

这条 branch 不是当前 live experiment 主线，而是基础设施整理线。本轮整理后，它已经形成两条本地提交：

- `24649aa (wip: snapshot foundation cache contract)`
- `764b8c0 (wip: add foundation cache alignment regression)`

它的主要价值在于把 foundation cache 相关的重复合同层抽成共享 helper：

- `vendor/transchrombp/transchrombp/utils/foundation_contract.py`
- `tests/test_foundation_contract.py`

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

### 3.2 已收口的 6000 / 6002 运行线

#### 6000：A6000 formal gate

- run：`teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1`
- 完成时间：`2026-04-10 10:10:10 CST` 左右
- 结论：`fail`
- best epoch：`22`
- best peak `profile_target_jsd_full_mean=0.3337811803218466`
- best peak `count_pearson_full=0.7948227405497477`
- 与历史 `corrected B` comparator 相比，JSD / count 都明显落后，因此不做 promotion

#### 6000：AlphaGenome v2 sidecar

- run：`alphagenome_matched_raw_track_slice_v2_20260410`
- 完成时间：`2026-04-10 03:34:42 CST`
- 结论：`pass`
- `16` 个 loci 全部完成，且每个位点都保留 `1` 条可用 `ATAC` track
- 这是已完成的 technical / external-coordinate sidecar，不占 active slot，也不改写 A6000 formal gate

#### 6002：U-Net-lite v1

- run：`teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r4`
- 完成时间：`2026-04-10 03:34:42 CST` 左右
- 结论：`no-go / stop`
- `r1/r2` 仍是启动失败，`r3/r4` 是有效 run，但 `r4` 仍没有把这条 family 推上 shortlist
- 当前没有 active 6002 training run

---

## 4. 当前仍需注意的不一致

- `TRACKING.md` 已把规则源指向 `docs/plan/2026-04-09_dual_machine_experiment_charter.md`，并把它作为 live 入口，而不是让 handoff 自己承担优先级裁决
- 旧的 `reports/session_handoff_multiscale_and_next_tasks_20260409.md` 仍保留着 `msdec_v1_s2 / skipprobe_wide` 的历史快照
- 真实现场已经切到：
  - `6000` 的 A6000 formal gate 已收口并判 `fail`
  - `6000` 的 AlphaGenome v2 sidecar 已收口并判 `pass`
  - `6002` 的 `U-Net-lite v1` 已收口并判 `no-go / stop`
- `dual-track-20260409` 的 pivot 计划与切换记录现已在 `reports/dual_track_pivot_snapshot_20260409.md` 里给出 canonical 恢复入口；不再要求依赖 live mounted worktree

---

## 5. 按当前 live 文档读取

1. `6000` 当前没有 active training run；`A6000 formal gate` 已 `fail`，AlphaGenome v1/v2 sidecar 都已完成 closeout。
2. `6002` 当前没有 active training run；`U-Net-lite v1` 已收成 `no-go / stop`。
3. 若未来要改默认下一步，以 charter / `TRACKING.md` 为准，本报告只记录终态快照，不再裁决。

---

## 6. 建议的阅读顺序

1. `docs/plan/2026-04-09_dual_machine_experiment_charter.md`：回答“当前什么实验合法、什么优先、什么时候才允许切回论文”
2. `TRACKING.md`：回答“当前 live 进度和下一步是什么”
3. 本文件：回答“现在真实在发生什么”
4. `reports/project_plan_code_review_20260405.md`：回答“最近代码为什么会改成这样”
5. `reports/session_handoff_multiscale_and_next_tasks_20260409.md`：保留为旧双机排程快照
6. `reports/dual_track_pivot_snapshot_20260409.md`
   用来回答“为什么从 `msdec/skipprobe` 切到 `teacher-distill + U-Net-lite`，以及现在如何从 snapshot 恢复那条 branch”
