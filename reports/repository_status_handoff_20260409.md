# 仓库现况总览（2026-04-09 22:46 CST）

> 双机实验规则源统一看 `docs/plan/2026-04-09_dual_machine_experiment_charter.md`；本文只负责快照，不负责裁决默认下一步。

## 一句话结论

截至 `2026-04-09 22:46 CST`，本地档案仓的真实状态已经不再是“`msdec_v1_s2 + skipprobe_wide` 双机并跑”，而是：

- 本地 canonical git 档案仓仍是 `/home/zhengwei/project/python/chromBPNet`，且本轮已完成 worktree cleanup
- `6000 / A6000x2` 上的 `NT v2 teacher-distill short10` 已完成，并已正式判负
- `6002 / RTX 3080` 已重新起跑 `U-Net-lite v1` 的确认性 cheap rerun `r4`，用来确认 `r3` 的 `provisional no-go` 是否稳定

---

## 1. 角色边界

### 1.1 本地仓库

- 本地仓库 `master` 是当前唯一有完整 git 历史的 canonical 档案仓。
- 当前状态：`master` 已重新回到 clean 状态，并保持为相对 `upstream/master` 的本地领先档案分支。
- 本轮已经用连续两条 docs 提交把原先混在主仓里的 docs / handoff / archive 内容收口进 canonical 档案仓。
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
| `master` | 本轮 docs closeout 后的最新状态 | clean | 主档案仓，承载 canonical docs / tracking / paper / reports |
| `dual-track-20260409` | `4ffc6e7` | clean | 最新实验 pivot：`teacher-distill + U-Net-lite` |
| `multiscale-decoder-probe-20260407` | `4c62096` | clean | `msdec/skipprobe` 旧 readout probe 的干净切片 |
| `autonomy/20260406-chrombpnet-externalization` | `4aadffa` | clean | `official chrombpnet` 外置化收尾与工具链 bridge |
| `autonomy/20260406-structure` | `764b8c0` | clean | `foundation cache contract` 抽象层重构 |

### 2.1 `master` 当前保留什么

- 现在 `master` 已不再混着实验性代码与工具链 WIP。
- 本轮保留在 `master` 的内容只有 canonical docs / archive：
  - `TRACKING.md`
  - `reports/repository_status_handoff_20260409.md`
  - `docs/plan/*.md`
  - `reports/*.md`
  - `reports/*.tex`
  - `docs/superpowers/` 下的设计 / 执行文档

### 2.2 `dual-track-20260409` 持有的新增代码

这条 worktree 已在本轮整理后落成 clean 提交 `4ffc6e7 (wip: snapshot dual track pivot)`，是当前活跃实验路线的真正代码入口，核心内容包括：

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

这条 worktree 不是当前 live experiment 主线，而是基础设施整理线。本轮整理后，它已经形成两条本地提交：

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

### 3.2 当前活跃线

#### 6000：`NT v2 teacher-distill short10`

- run：`ntv2_teacher_distill_short10_s42_6000_20260409_r2`
- 机器：`6000 / A6000 x2`
- 截至 `2026-04-09 22:43 CST` 的证据：
  - `nvidia-smi` 显示两张 A6000 都已回到 `13 MiB / 0% util`
  - `ps` 中已无该 run 的活跃训练进程
  - 日志 `/data1/zhoujiazhen/bylw_atac/TransChromBP/logs/ntv2_teacher_distill_short10_s42_6000_20260409_r2.log` 尾部明确写出 `[early-stop] epoch=5` 与 `NT v2 teacher distill completed.`
  - `run_meta.json` 已写出：
    - `best_epoch=2`
    - `best_metric_value=0.345324...`
    - `stopped_early=true`
  - log 中 `epoch=2` 对应的 `val:peak` 为：
    - `profile_target_jsd_full_mean=0.3453`
    - `count_pearson_full=0.5781`
- gate 判读：
  - matched no-foundation control 的同口径 `valid peak` best 为 `0.336346 / 0.8006`
  - distill `r2` 的 best `valid peak` 为 `0.345324 / 0.5781`
  - 因此既没有达到 `peak JSD` 改善 `>= 0.002`，也没有达到 `count_r` 不下降超过 `0.005`
- 结论：
  - run 已收口
  - formal gate 已判负
  - tutorial `teacher-distill` 线当前不升 full、不扩 seed，直接停表

### 3.3 已完成 cheap-screen、待正式判读的线

#### 6002：`U-Net-lite readout`

- 实际日志：
  - `/home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r1.log`
  - `/home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r2.log`
  - `/home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r3.log`
- 截至 `2026-04-09 22:41 CST`：
  - `6002` 的 3080 已重新被 `r4` 占用
  - 活跃训练进程已指向 `teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r4`
  - `r1` 是 loader/config contract 失败，`r2` 是 bigWig 路径 remap 失败
  - 当前真正有效的只有 `r3`
  - `r3` 日志里的 best 出现在 `epoch=2`，峰区 `profile_target_jsd_full_mean=0.45525`
  - 日志末尾已打印 `U-Net-lite decoder probe completed.`
- 结论：
  - 这条线不是“尚未起跑”，而是“已有一个有效 cheap-screen run，当前正在用 `r4` 做确认性补强”

---

## 4. 当前仍需注意的不一致

- `TRACKING.md` 现在应把规则源改指向 `docs/plan/2026-04-09_dual_machine_experiment_charter.md`，而不是让 handoff 自己承担优先级裁决
- 旧的 `reports/session_handoff_multiscale_and_next_tasks_20260409.md` 仍保留着“`msdec_v1_s2` / `skipprobe_wide` 是 active run”的旧快照
- 真实现场已经切到：
  - `6000 teacher-distill r2` 已结束，并已正式判负
  - `6002` 正在跑 `U-Net-lite r4` 确认性 cheap rerun
- `dual-track-20260409` 的计划与 pivot 记录已经以 `wip` 提交固定在对应 worktree，但尚未合回 `master`

---

## 5. 当前最值得做的事

1. 先把 `docs/plan/2026-04-09_dual_machine_experiment_charter.md` 视为双机实验规则源，把 `TRACKING.md` 视为 live 入口，而把本报告当成仓库 / worktree / 双机状态快照。
2. 等 `teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r4` 收口后补齐第二个 closeout：
   - `U-Net-lite v1` 在 `r4` 后能否正式停表
3. 先按 [reports/unet_lite_v1_rigor_review_20260409.md](reports/unet_lite_v1_rigor_review_20260409.md) 的口径收紧 `6002 U-Net-lite`：当前已从 `single-run provisional no-go` 进入确认性 cheap rerun 阶段，不能再把 `r1/r2/r3` 直接当成三次负向重复。

---

## 6. 建议的阅读顺序

1. `docs/plan/2026-04-09_dual_machine_experiment_charter.md`：回答“当前什么实验合法、什么优先、什么时候才允许切回论文”
2. `TRACKING.md`：回答“当前 live 进度和下一步是什么”
3. 本文件：回答“现在真实在发生什么”
4. `reports/project_plan_code_review_20260405.md`：回答“最近代码为什么会改成这样”
5. `reports/session_handoff_multiscale_and_next_tasks_20260409.md`：保留为旧双机排程快照
6. `dual-track-20260409` worktree 下两份文件：
   - `docs/plan/dual_track_ntv2_distill_unet_lite_20260409.md`
   - `reports/dual_track_experiment_pivot_20260409.md`
   用来回答“为什么从 `msdec/skipprobe` 切到 `teacher-distill + U-Net-lite`”
