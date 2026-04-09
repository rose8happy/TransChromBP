# 2026-04-09 双线实验切换记录

## 结论

- `msdec_v1` 与 `skipprobe` 两条 readout 小变体链在 `2026-04-09` 都已经拿到了足够明确的负向证据，因此不再继续往这条 family 上堆窄调参。
- 新的占卡方案正式切成两条互不重叠的线：
  - `6000 / A6000x2`：`NT v2 -> teacher distill`
  - `6002 / RTX 3080`：`U-Net-lite readout`
- 本地代码、配置、launcher 和文档都已补齐；接下来的关键动作只剩远端同步、起跑和 gate 判断。

## 已收口 run

### 1. 6000: `msdec_v1_s2`

- run name：`teacher_v2_center_pool_msdec_v1_s2_s42_20260409_full1`
- 日志：`/data1/zhoujiazhen/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_msdec_v1_s2_s42_20260409_full1.log`
- 关键结果：
  - peak valid `profile_target_jsd_full_mean=0.3314`
  - peak valid `count_pearson_full=0.8022`
  - overall valid `profile_target_jsd_full_mean=0.3541`
  - overall valid `count_pearson_full=0.8085`
- 判断：
  - 仍明显落后于 `corrected B` 双 seed `0.3146 / 0.8496`
  - 因此 `msdec_v1` 不具备“替代默认 readout”的价值

### 2. 6002: `skipprobe_v1_wide`

- run name：`teacher_v2_center_pool_skipprobe_v1_wide_short10_s42_6002_r1`
- 日志：`/home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_skipprobe_v1_wide_short10_s42_6002_r1.log`
- 关键结果：
  - peak valid `profile_target_jsd_full_mean=0.4494`
  - peak valid `count_pearson_full=0.8075`
  - overall valid `profile_target_jsd_full_mean=0.4692`
  - overall valid `count_pearson_full=0.8135`
- 判断：
  - 作为 3080 cheap-screen 也明显差于 matched `short10_nofoundation_control` 的 peak `0.3193 / 0.8298`
  - 因此 `skipprobe` 不再保留为 3080 默认候选

## 为什么改成“蒸馏 + 新 readout family”

- 之前的 `NT v2 residual` 结论只能否掉：
  - summary/coarse residual
  - bins16 residual
- 它并没有否掉“foundation 只做 teacher”的路线。
- 与此同时，`msdec/skipprobe` 这类小 readout 变体已经给出足够多的负向证据，再继续扩线只会重复低价值搜索。
- 所以新的划分更干净：
  - A6000 用来回答“teacher supervision 是否能带来 clean gain”
  - 3080 用来回答“真正不同的 readout family 是否能在低成本筛选里赢下来”

## 本轮新增实现

### 1. NT v2 teacher-distill 基础设施

- `teacher_cache_export.py`
  - 用 dataset-aligned cache 导出 `profile16/logcount`
  - 同步写出 `teacher_manifest_<split>.json`
- 新配置：
  - `transchrombp_teacher_v2_center_pool_ntv2_distill.yaml`
  - `train_tutorial_teacher_v2_ntv2_distill_short10.yaml`
  - `train_tutorial_teacher_v2_ntv2_distill_full.yaml`
- 新 launcher：
  - `run_ntv2_teacher_distill.sh`

### 2. U-Net-lite readout family

- `UNetLiteProfileDecoder`
  - encoder-decoder + skip fusion
  - 只负责 debiased profile readout
- 新配置：
  - `transchrombp_teacher_v2_center_pool_unet_lite_v1.yaml`
  - `train_tutorial_teacher_v2_readout_short10.yaml`
- 新 launcher：
  - `run_unet_lite_decoder_probe.sh`

### 3. 验证

- 本地回归：
  - `13 passed in 3.25s`
- 额外静态检查：
  - 两个 launcher 的 `bash -n`
  - `teacher_cache_export.py` 的 `py_compile`

## 当前 gate

### 6000 蒸馏线

- 先跑 short10
- 只有满足以下条件才升 full：
  - `peak JSD` 改善 `>= 0.002`
  - `peak count_r` 不下降超过 `0.005`
  - gap 不失控，默认 `<= 0.010`

### 6002 U-Net-lite 线

- 只跑 short10
- 只有满足以下条件才保留：
  - `JSD` 改善 `>= 0.002`
  - 或 `JSD` 在 `0.001` 内持平且 `count_r` 提升 `>= 0.01`

## 当前状态

- `2026-04-09` 实时复核时：
  - 6000 两张 A6000：`13 MiB / 0% util`
  - 6002 的 3080：`27 MiB / 0% util`
- 因此当前没有资源阻塞；真正的 blocker 只剩远端代码同步与 launcher 起跑。
