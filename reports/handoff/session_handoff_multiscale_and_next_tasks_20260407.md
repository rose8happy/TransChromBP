# 当前工作接续说明：Multiscale Decoder 双机解耦版（2026-04-07）

## 目的

这份 handoff 只回答三个问题：

1. 现在两台机器各自在跑什么？
2. 跑完后各自下一步是什么？
3. 哪些事明确不要再让两台机器互相等待？

它替代旧的 “`6000 s3` 等 `6002 s2` 再决策” 口径。

从 `2026-04-08 11:3x CST` 起，再补一条更硬的调度规则：

> `6000` 和 `6002` 各自维护独立结构探索队列；两边可以主题相关，但不互相等待，也不把对方结果当作新任务启动前置。

## 一、当前状态（`2026-04-08 22:53 CST` 复核）

### 1. 6000 `narrow full-budget` 已收口

- 机器：`6000 / A6000 x2`
- run name: `teacher_v2_center_pool_msdec_v1_narrow_s42_20260408_full1`
- 日志：`/data1/zhoujiazhen/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_msdec_v1_narrow_s42_20260408_full1.log`
- 模型：`msdec_v1_narrow`
- 配置：
  - `train_tutorial_teacher_v2_main.yaml`
  - `data_tutorial_canonical_v1.yaml`
  - `transchrombp_teacher_v2_center_pool_msdec_v1_narrow.yaml`
  - `TRAIN_SEED=42`
  - `BATCH_SIZE_PER_GPU=12`
- 启动时间：`2026-04-08 11:32 CST`
- 收口证据：
  - `2026-04-08 22:53 CST` 复核时两张 A6000 均空闲，未见训练进程
  - `run_meta.json` / log 最后写入时间约 `2026-04-08 19:13:56-19:13:57 CST`
  - 日志末尾已写出 `epoch=30`、`[best]` 与 `Multiscale decoder probe completed.`
- 最终读数：
  - best epoch `30`
  - best peak `loss_total=16.92284`
  - best peak `profile_target_jsd_full_mean=0.33934`
  - best peak `count_pearson_full=0.79843`
  - best overall `profile_target_jsd_full_mean=0.36188`
  - best overall `count_pearson_full=0.80052`
- 当前判读：
  - 相对历史 `corrected B` held-out 两 seed mean（约 `0.3146 / 0.8496`）仍明显偏弱
  - 当前应按既定 6000 gate 收口，倾向延续 `msdec_v1 full1` 的负 verdict，不再在这条小变体链上盲扩

### 2. 6002 `s2 short10` 已收口

- 机器：`6002 / RTX 3080 12G`
- run name: `teacher_v2_center_pool_msdec_v1_s2_short10_s42_6002_r1`
- 日志：`/home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_msdec_v1_s2_short10_s42_6002_r1.log`
- 模型：`msdec_v1_s2`
- 配置：
  - `train_tutorial_foundation_short10.yaml`
  - `data_tutorial_canonical_v1_6002.yaml`
  - `transchrombp_teacher_v2_center_pool_msdec_v1_s2.yaml`
  - `TRAIN_SEED=42`
  - `BATCH_SIZE_PER_GPU=16`
- 启动时间：`2026-04-08 14:59 CST`
- 收口证据：
  - `2026-04-08 22:53 CST` 复核时 3080 已空闲，未见训练进程
  - `run_meta.json` / log 最后写入时间约 `2026-04-08 20:18:52 CST`
  - 日志末尾已写出 `epoch=10`、`[best]` 与 `Multiscale decoder probe completed.`
- 最终读数：
  - best epoch `10`
  - best peak `profile_target_jsd_full_mean=0.44732`
  - best peak `count_pearson_full=0.80733`
  - best overall `profile_target_jsd_full_mean=0.46716`
  - best overall `count_pearson_full=0.81335`
- 当前判读：
  - 相对 `6002 anchor / narrow` 是明确正向改进
  - 相对 `skipprobe` 是“JSD 基本持平但 count 更稳”的 tradeoff，符合进入 `6002` shortlist 的 cheap-screen 口径

## 二、已经结束或被主动取消的任务

### 1. 6000 上一条 `msdec_v1 full1` 已完成并判负

- run name: `teacher_v2_center_pool_msdec_v1_s42_20260407_full1`
- 收口证据：`2026-04-08 11:24 CST` 复核时无训练进程，log 已完整写出 `epoch=28`、`early-stop` 与完成标记
- best epoch: `22`
- best peak `profile_target_jsd_full_mean=0.33456`
- best peak `count_pearson_full=0.80959`
- best overall `profile_target_jsd_full_mean=0.35714`
- best overall `count_pearson_full=0.81544`
- 结论：
  - 相对历史 `corrected B` held-out 两 seed mean（约 `0.3146/0.8496`）恶化约 `+0.01996 JSD / -0.04001 count`
  - 触发 paper-facing 硬阈值 `peak JSD > 0.3197` 与 `peak count_r < 0.8400`
  - 不再开 `teacher_v2_center_pool_msdec_v1_s1234`

### 2. 6000 结构 smoke 已完成

- run name: `teacher_v2_center_pool_msdec_v1_short10_smoke_s42_20260407`
- 关键结论：
  - best epoch `9`
  - peak `profile_target_jsd_full_mean=0.33617`
  - peak `count_pearson_full=0.79695`
  - 可以授权进入 A6000 full-budget

### 3. 6002 旧 `s2` 对照已停止

- run name: `teacher_v2_center_pool_msdec_v1_s2_short10_smoke_s42_6002_20260407_fix1`
- 停止时间：`2026-04-07 17:33 CST`
- 停止原因：
  - 它把 `6002` 绑定成 `6000` 的从属 gate
  - 当前目标改成“双机都做结构探索，但不互相等待”

### 4. 6002 anchor 已完成

- run name: `teacher_v2_center_pool_msdec_v1_short10_s42_6002_anchor_20260407`
- 收口证据：`2026-04-07 20:20 CST` 左右 log 与 `run_meta.json` 停止写入，当前无训练进程
- best epoch: `2`
- best peak `profile_target_jsd_full_mean=0.45232`
- best peak `count_pearson_full=0.76995`
- best overall `profile_target_jsd_full_mean=0.47167`
- best overall `count_pearson_full=0.77610`
- 最终状态：
  - `epoch=5` 触发 early-stop
  - 这条 run 已经固定了 6002 同机比较基线

### 5. 6002 narrow 已完成

- run name: `teacher_v2_center_pool_msdec_v1_narrow_short10_s42_6002_r1`
- 收口证据：`2026-04-08 09:43 CST` 复核时 3080 空闲，log 已写出 `epoch=8` early-stop 与完成标记
- best epoch: `5`
- best peak `profile_target_jsd_full_mean=0.44927`
- best peak `count_pearson_full=0.79079`
- best overall `profile_target_jsd_full_mean=0.46895`
- best overall `count_pearson_full=0.79626`
- 结论：
  - 按既定 6002 gate，这是相对 anchor 的正向改进
  - 可以进入 6002 shortlist，并作为 `6002` 自己后续单卡结构探索的参考

### 6. 6002 skipprobe 已完成

- run name: `teacher_v2_center_pool_skipprobe_v1_short10_s42_6002_r1`
- 启动时间：`2026-04-08 09:54 CST`
- 收口时间：`2026-04-08 14:53 CST` 左右
- 日志：`/home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_skipprobe_v1_short10_s42_6002_r1.log`
- 最终结论：
  - `epoch=9` early-stop
  - best epoch `6`
  - best peak `profile_target_jsd_full_mean=0.44535`
  - best peak `count_pearson_full=0.79626`
  - best overall `profile_target_jsd_full_mean=0.46523`
  - best overall `count_pearson_full=0.80196`
  - 相对 `6002 narrow`，它是当前更优的单卡 cheap-screen 结果

### 7. 6002 `s2 short10` 已完成

- run name: `teacher_v2_center_pool_msdec_v1_s2_short10_s42_6002_r1`
- 启动时间：`2026-04-08 14:59 CST`
- 日志：`/home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_msdec_v1_s2_short10_s42_6002_r1.log`
- 收口时间：`2026-04-08 20:18 CST` 左右
- 最终结论：
  - best epoch `10`
  - best peak `profile_target_jsd_full_mean=0.44732`
  - best peak `count_pearson_full=0.80733`
  - best overall `profile_target_jsd_full_mean=0.46716`
  - best overall `count_pearson_full=0.81335`
  - 相对 `6002 anchor / narrow` 为正向改进
  - 相对 `skipprobe` 则是 JSD 近乎持平、count 更稳，符合进入 `6002` shortlist 的 cheap-screen 口径

## 三、当前研究计划

### 6000 的职责

- `msdec_v1 full1` 已完成职责，并给出负 verdict：当前不值得成为正式下一代 readout 候选
- `narrow full-budget` 也已跑完；当前只剩相对 `corrected B` 的正式 gate 与收口记录
- 从这条之后开始，`6000` 的后续候选由 `6000` 自己的 A6000 队列决定，不再要求来自 `6002` 的 shortlist 才能启动
- 它不再为 `msdec_v1` 打开 `s1234`

### 6002 的职责

- 只回答：在单卡、短预算、同机口径下，哪种 readout 结构最值得进入下一轮 shortlist？
- 它不等待 `6000` 的 full-budget verdict 才继续跑自己的队列
- 它的结果从当前轮开始也不再反向决定 `6000` 下一条开谁；当前 `skipprobe` 与 `s2` 都只服务于 `6002` 自己的 shortlist 排序

一句话说：

> `6000` 负责自己的正式 verdict 队列，`6002` 负责自己的便宜排序队列；两边可以研究相邻 family，但不再共享同一个 stop/go 节点，也不互相卡启动。

## 四、跑完后各自下一步

### 6000 当前下一步

1. 记录 `teacher_v2_center_pool_msdec_v1_narrow_s42_20260408_full1` 已于 `2026-04-08 19:13 CST` 左右收口
2. 只与历史 `corrected B` 做正式 gate
3. 如果 `narrow full-budget` 仍不过 gate，则不要继续在这条小变体链上盲扩，而要先切到新的 readout family
4. 按独立队列口径整理 `6000` 自己的下一组 A6000 候选，不等待 `6002`

### 6002 当前下一步

不需要等 6000；当前 `s2 short10` 已经收口，继续维护 `6002` 自己的单卡筛选队列。

固定顺序直接转：

1. 记录 `skipprobe` 已优于 `6002 narrow`
2. 记录 `s2 short10` 已收口，并与 `skipprobe` 形成“JSD 近乎持平、count 更稳”的 tradeoff
3. 明确 `6002` 内部 shortlist 排序后，继续下一条单卡 cheap-screen；是否进入 `6000` 另由 `6000` 自己的独立队列决定

## 五、复核命令

### 6000 复核 `narrow full-budget`

```bash
ssh zhoujiazhen@127.0.0.1 -p 6000 'tail -n 80 /data1/zhoujiazhen/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_msdec_v1_narrow_s42_20260408_full1.log'
ssh zhoujiazhen@127.0.0.1 -p 6000 'grep -i "Multiscale decoder probe completed\\|early-stop\\|\\[val\\]" /data1/zhoujiazhen/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_msdec_v1_narrow_s42_20260408_full1.log | tail -n 20'
```

### 6002 复核 `s2 short10`

```bash
ssh -i /home/zhengwei/.ssh/codex_6002_ed25519 zhengwei@127.0.0.1 -p 6002 'tail -n 80 /home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_msdec_v1_s2_short10_s42_6002_r1.log'
ssh -i /home/zhengwei/.ssh/codex_6002_ed25519 zhengwei@127.0.0.1 -p 6002 'grep -i "Multiscale decoder probe completed\\|early-stop\\|\\[val\\]" /home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_msdec_v1_s2_short10_s42_6002_r1.log | tail -n 20'
```

## 六、当前明确不要做的事

1. 不要再把 `6002` 当成“先跑完再决定 6000 能不能继续”的 gate 机
2. 不要现在切到 `GM12878 / K562` 数据线
3. 不要回到旧的 foundation adapter family
4. 不要因为 `skipprobe / s2` 都已收口，就继续在同一条 `msdec` 小变体链上无门槛盲扩；先完成 `corrected B` gate 与 `6002` shortlist 排序

## 一句话总结

**截至 `2026-04-08 22:53 CST`，`6000 narrow full-budget` 与 `6002 s2 short10` 都已实际收口；从这一刻起，两边都进入各自的收口与下一轮候选排序，不再互相等待。**
