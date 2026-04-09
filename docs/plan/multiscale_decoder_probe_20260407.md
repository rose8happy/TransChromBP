# `corrected B + multiscale decoder probe` 实施计划（2026-04-07）

## 目标

把“下一代路线”的第一优先级压缩成一条可执行的结构实验树：

- 固定 `corrected B` 作为唯一主比较基线
- 保留 `bias branch`、`full/debiased` 语义、`center pool` count head
- 只替换 profile readout
- 让 `6000 / A6000x2` 与 `6002 / RTX 3080 12G` 都做结构探索，但彼此不再互相等待

## 已落地代码

1. 新增 [profile_decoder.py](/home/zhengwei/project/python/chromBPNet/vendor/transchrombp/transchrombp/models/profile_decoder.py)
   - `MultiScaleProfileDecoder`
   - `LocalSkipProfileProbe`

2. 扩展 [transchrombp.py](/home/zhengwei/project/python/chromBPNet/vendor/transchrombp/transchrombp/models/transchrombp.py)
   - `profile_readout_mode` 现支持：
     - `linear`
     - `multiscale_decoder_v1`
     - `skip_probe_v1`

3. 新增模型配置
   - [transchrombp_teacher_v2_center_pool_msdec_v1.yaml](/home/zhengwei/project/python/chromBPNet/vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_msdec_v1.yaml)
   - [transchrombp_teacher_v2_center_pool_msdec_v1_s2.yaml](/home/zhengwei/project/python/chromBPNet/vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_msdec_v1_s2.yaml)
   - [transchrombp_teacher_v2_center_pool_msdec_v1_narrow.yaml](/home/zhengwei/project/python/chromBPNet/vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_msdec_v1_narrow.yaml)
   - [transchrombp_teacher_v2_center_pool_skipprobe_v1.yaml](/home/zhengwei/project/python/chromBPNet/vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool_skipprobe_v1.yaml)

4. 复用 launcher [run_multiscale_decoder_probe.sh](/home/zhengwei/project/python/chromBPNet/vendor/transchrombp/transchrombp/scripts/run_multiscale_decoder_probe.sh)
   - `6000` 用它跑 full-budget 主线
   - `6002` 用它跑单卡 `short10` 结构筛选树

5. 扩展测试 [test_multiscale_decoder_probe.py](/home/zhengwei/project/python/chromBPNet/tests/test_multiscale_decoder_probe.py)
   - multiscale shape contract
   - `full/debiased` additive contract
   - `skip_probe_v1` shape contract
   - `narrow / skipprobe` YAML 构建路径

## 当前验证

本地已通过：

```bash
./.venv/bin/python -m pytest tests/test_multiscale_decoder_probe.py -q
```

远端已同步并通过：

```bash
/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b/bin/python -m py_compile ...
/home/zhengwei/bylw_atac/.mamba/envs/transchrombp/bin/python -m py_compile ...
```

## 已得到的结构证据

`6000` 第一条结构 smoke 已完成：

- run name: `teacher_v2_center_pool_msdec_v1_short10_smoke_s42_20260407`
- best epoch: `9`
- peak `profile_target_jsd_full_mean=0.33617`
- peak `count_pearson_full=0.79695`
- 结论：`msdec_v1` 至少 non-collapse，可升级到 tutorial canonical full-budget

`6002` 旧对照 `s2` 已不再作为当前决策 gate：

- 已于 `2026-04-07 17:33 CST` 停止
- 原因：它把 3080 绑定成“给 6000 做从属对照”的角色，不符合当前双机解耦目标

## 新的双机分工

### 硬规则

- `6000` 和 `6002` 各自维护独立结构探索队列。
- 两边可以研究相邻 family，但**不互相等待**，也不把对方结果当作新任务启动前置。
- `6000` 只按自己的 `A6000 full-budget` 价值排序；`6002` 只按自己的 `short10 / 单卡` 信息增益排序。

### 6000 / A6000 x2：正式 verdict 队列

只承担高价值、双卡才值得跑的正式 verdict；默认优先完整 decoder family，不用它去追单卡 cheap gate。

当前 active：

1. `teacher_v2_center_pool_msdec_v1_narrow_s42_20260408_full1`
   - 日志：`/data1/zhoujiazhen/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_msdec_v1_narrow_s42_20260408_full1.log`
   - 配置：
     - `TRAIN_CONFIG=configs/train/train_tutorial_teacher_v2_main.yaml`
     - `DATA_CONFIG=configs/data/data_tutorial_canonical_v1.yaml`
     - `MODEL_CONFIG=configs/model/transchrombp_teacher_v2_center_pool_msdec_v1_narrow.yaml`
     - `TRAIN_SEED=42`
     - `TRAIN_GPU_IDS=0,1`
     - `BATCH_SIZE_PER_GPU=12`

当前 launch-ready backlog：

2. `msdec_v1_s2 full-budget`
   - 结构含义：保持 `multiscale_decoder_v1` family 不变，只把 decoder scales 从 `3 -> 2`
   - 当前价值：这是一个真正不同于 `narrow` 的 decoder-depth ablation，适合放在 `6000` 自己的完整 decoder 队列里
   - 启动条件：只取决于 `6000` 侧的队列优先级，不等待 `6002`

3. `msdec_v2_wide/deep`（待代码/配置落地）
   - 结构含义：仍留在 bias-safe multiscale decoder family，但往“更强 decoder”方向走，而不是回到 adapter family
   - 当前状态：还不是 launch-ready；等 `narrow / s2` 给出方向后再落配置

`6000` 当前明确不默认排入：

- `skipprobe full-budget` 的镜像复现
- 任何 `ntv2_*` adapter family 变体
- `GM12878 / K562` 新数据线长训练

### 6002 / RTX 3080 12G：单卡 cheap-screen 队列

只承担同机口径下的廉价结构排序；它的目标是快速排掉不值得继续在 3080 深挖的 readout 变体，不再承担 A6000 的 gate 职责。

当前 active：

1. `teacher_v2_center_pool_skipprobe_v1_short10_s42_6002_r1`
   - 日志：`/home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_skipprobe_v1_short10_s42_6002_r1.log`
   - 配置：
     - `TRAIN_CONFIG=configs/train/train_tutorial_foundation_short10.yaml`
     - `DATA_CONFIG=configs/data/data_tutorial_canonical_v1_6002.yaml`
     - `MODEL_CONFIG=configs/model/transchrombp_teacher_v2_center_pool_skipprobe_v1.yaml`
     - `TRAIN_SEED=42`
     - `TRAIN_GPU_IDS=0`
     - `BATCH_SIZE_PER_GPU=16`

当前 launch-ready backlog：

2. `msdec_v1_s2 short10`
   - 结构含义：保持 multiscale decoder family，但把 decoder scales 从 `3 -> 2`
   - 当前价值：它是 6002 上现成、尚未跑过、且不同于 `skipprobe` 的 cheap ablation
   - 启动条件：只取决于 `6002` 自己的队列，不等待 `6000`

3. `6002-only cheap readout ablation`
   - `2026-04-09` 已优先落地为 `skipprobe_wide`
   - 当前 run：`teacher_v2_center_pool_skipprobe_v1_wide_short10_s42_6002_r1`
   - 配置：`configs/model/transchrombp_teacher_v2_center_pool_skipprobe_v1_wide.yaml`
   - 结构含义：保持 `skip_probe_v1` family 不变，只把 readout hidden width 从 `128 -> 256`
   - 目的：继续让 3080 做“便宜但结构上真的不同”的筛选，而不是重复 `6000` 正式 verdict

`6002` 当前明确不默认排入：

- 为 `6000` 做镜像 full-budget 对照
- 任何 `ntv2_*` adapter family 变体
- 跨到 `GM12878 / K562` 新数据线

## Gate 规则

### 6000 gate

- 只和历史 `corrected B` 比，不和 `6002` 比
- 相对历史 `corrected B`，若 peak `profile_target_jsd_full_mean` 恶化超过 `0.003`，判不过
- 若 peak `count_pearson_full` 恶化超过 `0.01`，判不过
- 若 `full/debiased gap` 明显放大，判不过
- `6002` 的结果不能阻塞或触发 `6000` 的启动

### 6002 gate

- 只和 `6002` 自己的 anchor / 已完成单卡基线比较
- `skipprobe` / `s2` 只服务于 `6002` 自己的 cheap-screen 排序
- 若 JSD 改善至少 `0.002`，或 JSD 基本持平且 count 更稳，可进入 `6002` 自己的 shortlist
- `6002` 的 shortlist 不直接等于 `6000` 的下一条任务

## `2026-04-08`：6000 `full1` gate 结果

`teacher_v2_center_pool_msdec_v1_s42_20260407_full1` 已完成正式判定。

### 证据

- `6000 full1` best epoch：`22`
- best peak `profile_target_jsd_full_mean=0.3345583154971703`
- best peak `count_pearson_full=0.8095856875678048`
- best peak `full/debiased gap ~= -6.4e-06`（安全性本身没有异常放大）
- 历史 `corrected B` held-out 两 seed：
  - `B_s42`：`0.3146706970465267 / 0.8503032921152159`
  - `B_s1234`：`0.3145302712082182 / 0.8488473843266611`
  - 两 seed mean 约：`0.3146 / 0.8496`

### 判读

- 相对 `corrected B` 两 seed mean，`msdec_v1 full1` 恶化约：
  - peak JSD `+0.01996`
  - peak count `-0.04001`
- 这不仅超过了本计划里“JSD 恶化 `>0.003` / count 恶化 `>0.01`”的 6000 gate，
  也直接触发了既有 paper-facing 硬阈值：
  - `peak.profile_target_jsd_full_mean > 0.3197`
  - `peak.count_pearson_full < 0.8400`

### 结论

- `teacher_v2_center_pool_msdec_v1_s42_20260407_full1` **不过 gate**
- 不再开 `teacher_v2_center_pool_msdec_v1_s1234`
- `6000` 的旧 `msdec_v1` 主线到此结束；当前 A6000 队列已经切换到 `narrow full-budget`

## 当前建议顺序

1. `6000` 继续跑完 `teacher_v2_center_pool_msdec_v1_narrow_s42_20260408_full1`
2. `6002` 继续跑完 `teacher_v2_center_pool_skipprobe_v1_short10_s42_6002_r1`
3. 两边收口后，各自只按自己的 gate 口径判读，不互相卡下一条
4. `6000` 的下一条默认从自己的 backlog 里选 `msdec_v1_s2 full-budget` 或新的 decoder 配置
5. `6002` 的下一条默认从自己的 backlog 里选 `msdec_v1_s2 short10` 或新的 6002-only cheap ablation
