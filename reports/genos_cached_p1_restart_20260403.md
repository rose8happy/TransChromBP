# Genos Cached `P0/P1` Restart 执行记录（2026-04-03）

## 背景

- `2026-04-03 16:35 CST` 实查显示：`6000` 两张 `RTX A6000` 都空闲，`6002` 的 `RTX 3080` 被 `ComfyUI` 常驻进程占住约 `8.2 GiB` 显存且本轮不动。
- `Caduceus-PS` tutorial matched A/B 已在同日下午按 `near-null / marginal positive` 收口，不适合继续直接扩到 `GM12878`。
- `SegmentNT / gReLU / HyenaDNA` 在 `6000` 上没有本地可用资产；`AlphaGenome` 更适合小规模 black-box pilot，不适合作为下一条占满 `A6000` 的主线。
- 因此在“方法探索优先”的前提下，当前唯一能当天启动、且能在 `1` 个晚上内给出 go/no-go 的路线，是重启一条最小 `Genos cached summary` pilot。

## 更正说明

- `16:59 CST` 第一版 launcher `genos_cached_restart_p0p1_20260403_165748` 已经启动，但 `step0 build_genos_summary_cache.py` 被我错误地写成单进程单卡。
- 用户指出仓库约定后，我已在 `17:14 CST` 主动停掉该 launcher，并把未完成 cache 目录归档为 `tutorial_global_20260403_165748_aborted_singlegpu_*`。
- 当前有效 launcher 已更正为：`genos_cached_restart_p0p1_dualcache_20260403_171304`。

## 当前有效执行内容

### 运行环境

- 机器：`6000`
- GPU：`GPU 0 + 1`
- 拓扑：`step0` 为两路单卡并行，`step2/step3` 为 `2-rank DDP`
- venv：`/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b`
- 根目录：`/data1/zhoujiazhen/bylw_atac/TransChromBP`

### 启动信息

- launcher：`genos_cached_restart_p0p1_dualcache_20260403_171304`
- 启动时间：`2026-04-03 17:14:33 CST`
- PID：`3659968`
- launcher log：
  - `/data1/zhoujiazhen/bylw_atac/logs/genos_cached_restart_p0p1_dualcache_20260403_171304.launch.log`
- cache 子日志：
  - `/data1/zhoujiazhen/bylw_atac/logs/genos_cached_restart_p0p1_dualcache_20260403_171304.cache_train.log`
  - `/data1/zhoujiazhen/bylw_atac/logs/genos_cached_restart_p0p1_dualcache_20260403_171304.cache_valid.log`
- 临时脚本：
  - `/tmp/genos_cached_restart_dualcache_20260403_171304.sh`
- 预计结束时间：
  - `step0 cache` 已于 `2026-04-04 00:13 CST` 实际完成；当前已切入 `P0` 双卡训练。按 `00:24 CST` 时的 `epoch=1 step=3100/11809` 与当前 step-time 粗算，`P0` 整体量级约还需 `2-2.5` 小时，随后 launcher 才会进入 `P1`。

### 当前双卡占用

- `17:14 CST` 启动后，launcher 已实际拉起两路 cache：
  - `train cache pid=3659976`
  - `valid cache pid=3659977`
- 这两路分别对应：
  - `GPU0 -> train`
  - `GPU1 -> valid`
- `2026-04-03 18:19 CST` 复查：
  - `valid` 已完成并写出 `manifest_valid.json`（时间戳 `17:37 CST`）
  - `train cache` 仍在 `GPU0` 上运行
  - 因此 `P0/P1` 还没有开始，当前阶段仍是 cache build

### 新产物路径

- cache 目录：
  - `/data1/zhoujiazhen/bylw_atac/TransChromBP/genos_cache/tutorial_global_20260403_171304`
- probe 输出：
  - `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/genos_cached_probes/20260403_171304_p0p1_restart`
- 运行时 train config：
  - `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/runtime/20260403_171304_genos_cached_restart/train_genos_cached_short10_runtime.yaml`

### 2026-04-03 晚间调度更新

- `18:38-18:40 CST` 已在 `GPU1` 完成 `GM12878_6000 best test-full` sidecar，结论是 6000 侧 held-out 偏弱可复现。
- `20:38-20:40 CST` 已在 `GPU1` 完成 `K562_6000 epoch_019 test-full` sidecar，结论是后期 checkpoint 只是在继续放大 profile/count tradeoff，不改写 `best.pt` 的 profile-first 选择。
- 因为这两个高信息量、几分钟级的本地 sidecar 都已收口，而 `Genos train cache` 仍在推进，所以从 `2026-04-03 20:47 CST` 起不再继续给 `GPU1` 排新任务，改为显式预留给 `P0/P1` 的双卡切换。
- 当前最新可见 cache 进度为 `216000/426545 records`，最近一次 ETA `12403s`，约 `3.4` 小时；但仍应把它视为 rough ETA，而不是硬结束时间。

### 2026-04-04 凌晨进度更新

- `00:13 CST`：`manifest_train.json` 写出，说明 `train cache` 已完整结束。
- 随后 `step1 summary probe` 正常完成，probe 结果继续落盘到 `outputs/genos_cached_probes/20260403_171304_p0p1_restart/probe_results.json`。
- `00:20 CST` 左右：launcher 正式进入 `[step2]`，以双卡 `torchrun` 启动 `genos_cached_P0_baseline_rerun_20260403_171304`。
- `00:24 CST` 复查：两张 `A6000` 都在工作（约 `94-96% util`、各 `3.0 GiB` 显存），launcher 日志已推进到 `epoch=1 step=3100/11809`。当前仍未进入 `P1`。

### 子 run

- `P0`：
  - `genos_cached_P0_baseline_rerun_20260403_171304`
- `P1`：
  - `genos_cached_P1_late_film_restart_20260403_171304`

## 当前 launcher 实际做什么

1. 把 `global_mean cache` 拆成两路：`train` 在 `GPU0`，`valid` 在 `GPU1`，并行构建。
2. 基于历史 `G0 best.pt` 重跑一次 `genos summary probe`，确保 cache 与当前代码路径兼容。
3. 跑一条 fresh `P0 baseline_short10`，避免直接拿 `2026-03-23` 的旧 `P0` 当对照。
4. 跑一条 `P1 global_late_film`。

本轮没有包含：

- `P2 global_count_only`
- `online Genos G1/G2`
- `GM12878`
- `6002`
- 任何新的 foundation model 下载

## 判断门槛

- 先看 `epoch 5-8` 的同 epoch 对比。
- 继续推进：
  - `P1` 相对本轮新 `P0` 的 `peak.profile_target_jsd_full_mean` 改善 `>= 0.002`
  - 或 `JSD >= 0.0015` 且 `peak.count_pearson_full >= +0.005`
- 直接收口：
  - 到 `epoch 4` 仍基本追不平 `P0`
  - 或 step-time overhead 明显偏高
  - 或 `P1` 的增益只停留在噪声级别

## 本轮完成后的下一步

- 若 `P1` 过 gate：补 `held-out test-full`，再决定是否继续 `offset-aware global summary`。
- 若 `P1` 不过 gate：把 `Genos` foundation 主线正式降级，不再继续占用 `A6000` 主训练窗口。

## 进度检查命令

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1   'tail -f /data1/zhoujiazhen/bylw_atac/logs/genos_cached_restart_p0p1_dualcache_20260403_171304.launch.log'
```

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1   'tail -f /data1/zhoujiazhen/bylw_atac/logs/genos_cached_restart_p0p1_dualcache_20260403_171304.cache_train.log'
```

```bash
ssh -p 6000 zhoujiazhen@127.0.0.1   'tail -f /data1/zhoujiazhen/bylw_atac/logs/genos_cached_restart_p0p1_dualcache_20260403_171304.cache_valid.log'
```
