# 2026-04-09 双线实验执行计划：NT v2 蒸馏 + U-Net-lite

## 1. 触发原因

- `6000 / A6000x2` 上的 `teacher_v2_center_pool_msdec_v1_s2_s42_20260409_full1` 已收口：
  - 日志：`/data1/zhoujiazhen/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_msdec_v1_s2_s42_20260409_full1.log`
  - peak valid：`profile_target_jsd_full_mean=0.3314`、`count_pearson_full=0.8022`
  - 结论：仍落后于 paper-facing `corrected B` 双 seed 口径 `0.3146 / 0.8496`，不值得继续给 `msdec_v1` family 扩 full-budget 变体。
- `6002 / RTX 3080` 上的 `teacher_v2_center_pool_skipprobe_v1_wide_short10_s42_6002_r1` 也已收口：
  - 日志：`/home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_skipprobe_v1_wide_short10_s42_6002_r1.log`
  - peak valid：`profile_target_jsd_full_mean=0.4494`、`count_pearson_full=0.8075`
  - 结论：作为 cheap-screen 也没有逼近 matched `short10_nofoundation_control` 的 peak `0.3193 / 0.8298`，因此 `skipprobe` 不再作为 3080 主线。
- `2026-04-09` 再次实时复核后，6000 两张 A6000 与 6002 的 3080 均为空闲，因此可以把新的双线方案直接排上队列。

## 2. 已完成本地实现

### 2.1 NT v2 teacher-distill 线

- 新增 `teacher cache` 导出模块：
  - `src/transchrombp/evaluation/teacher_cache_export.py`
  - 从 dataset-aligned NT v2 cache 生成：
    - `teacher_manifest_<split>.json`
    - `<split>_profile16.f32.npy`
    - `<split>_logcount.f32.npy`
- 新增蒸馏模型与训练配置：
  - `configs/model/transchrombp_teacher_v2_center_pool_ntv2_distill.yaml`
  - `configs/train/train_tutorial_teacher_v2_ntv2_distill_short10.yaml`
  - `configs/train/train_tutorial_teacher_v2_ntv2_distill_full.yaml`
- 新增 launcher：
  - `scripts/run_ntv2_teacher_distill.sh`
  - 串行阶段固定为：
    1. `build_foundation_cache.py`
    2. `teacher_cache_export`
    3. `train_ddp`

### 2.2 U-Net-lite readout 线

- 在 `src/transchrombp/models/profile_decoder.py` 新增 `UNetLiteProfileDecoder`
- 在 `src/transchrombp/models/transchrombp.py` 接入 `profile_readout_mode=unet_lite_v1`
- 新增模型/训练配置：
  - `configs/model/transchrombp_teacher_v2_center_pool_unet_lite_v1.yaml`
  - `configs/train/train_tutorial_teacher_v2_readout_short10.yaml`
- 新增 launcher：
  - `scripts/run_unet_lite_decoder_probe.sh`

### 2.3 本地验证

- 单测：
  - `tests/test_teacher_cache_export.py`
  - `tests/test_teacher_distill_configs.py`
  - `tests/test_ntv2_teacher_distill_launcher.py`
  - `tests/test_multiscale_decoder_probe.py`
  - `tests/test_ntv2_bins16_residual_head.py`
  - `tests/test_unet_lite_launcher.py`
- 当前结果：
  - `13 passed in 3.25s`
  - `bash -n scripts/run_ntv2_teacher_distill.sh`
  - `bash -n scripts/run_unet_lite_decoder_probe.sh`
  - `python -m py_compile src/transchrombp/evaluation/teacher_cache_export.py`

## 3. 6000 主线：NT v2 teacher-distill

### 3.1 固定边界

- 学生模型保持 `corrected B` 语义：
  - bias branch 保留
  - `full/debiased` additive contract 不变
  - count head 仍是 `center-pool`
- foundation 只作为 offline teacher：
  - `foundation_model.mode=distill_only`
  - 不再启用 `residual_head` / `cross_attention`
- 第一轮 teacher targets 固定为：
  - `profile16`
  - `logcount`
- 第一轮不启用 `rank`

### 3.2 启动顺序

1. 先跑 short10 gate：
   - launcher：`bash scripts/run_ntv2_teacher_distill.sh`
   - train config：`train_tutorial_teacher_v2_ntv2_distill_short10.yaml`
2. 只有 short10 通过 gate 才升 full：
   - train config：`train_tutorial_teacher_v2_ntv2_distill_full.yaml`

### 3.3 gate

- 相对 matched no-foundation short10 control，必须同时满足：
  - `peak JSD` 改善 `>= 0.002`
  - `peak count_r` 不下降超过 `0.005`
  - profile/count 的优势不能靠单侧明显塌陷换来，默认要求 `gap <= 0.010`
- 不过 gate：
  - 不扩 seed
  - 不升 full-budget
  - 直接停表

### 3.4 预计耗时

- 如果 cache 需要重建，short10 全流程预计 `4-6` 小时：
  - cache build：约 `1-2` 小时
  - teacher cache export：约 `10-30` 分钟
  - short10 训练：约 `3-4` 小时
- full-budget 只有在 gate 通过后才排，预计额外 `8-12` 小时。

## 4. 6002 主线：U-Net-lite readout

### 4.1 固定边界

- 只替换 debiased profile readout
- backbone、bias branch、count head 全部保持 `corrected B`
- skip 只来自 local/encoded feature，不把 bias 特征回灌进 decoder
- full profile 仍通过现有 additive fusion 形成

### 4.2 启动顺序

- 只跑单卡 short10 cheap-screen：
  - launcher：`bash scripts/run_unet_lite_decoder_probe.sh`
  - train config：`train_tutorial_teacher_v2_readout_short10.yaml`
- 这条线不自动晋级 6000，也不自动扩 seed。

### 4.3 gate

- 相对当前 `skipprobe/msdec` shortlist，满足任一条件即可保留：
  - `JSD` 改善 `>= 0.002`
  - `JSD` 在 `0.001` 内持平且 `count_r` 提升 `>= 0.01`
- 如果达不到：
  - 3080 readout 新 family 直接停表

### 4.4 预计耗时

- 单卡 short10 预计 `4.5-5.5` 小时，按当前 `skipprobe_wide` 日志估算与之同量级。

## 5. 执行与同步注意事项

- 本地归档仓使用 `vendor/transchrombp/transchrombp/...`；6000/6002 运行仓使用 `src/transchrombp/...`，同步时必须按路径映射写回。
- `TRACKING.md` 必须在三次时点同步：
  - 本轮本地实现与 pivot 结论形成后
  - 远端后台任务刚启动后
  - short10 gate 收口后
- 如果 launch 成功，必须立刻补：
  - run name
  - 机器 / GPU
  - 日志路径
  - 预计结束时间窗
  - 下一步 gate 动作
