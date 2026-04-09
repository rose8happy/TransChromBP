# 双机实验执行宪章（2026-04-09）

> 生效时间：`2026-04-09 22:11 CST`
> 本文件只负责回答一件事：当前 `6000` / `6002` 上什么实验合法、什么实验优先、什么情况下才允许把默认动作切回“先写论文”。

## 1. 角色与权威顺序

当前双机实验决策按以下顺序裁决：

1. 最新实验事实与日志
2. 本文件
3. `TRACKING.md` 的 live 状态
4. `reports/repository_status_handoff_20260409.md` 等快照 / handoff / 历史计划文档

三份常用文档的职责固定如下：

- `docs/plan/2026-04-09_dual_machine_experiment_charter.md`：规则源
- `TRACKING.md`：live 状态与下一步
- `reports/repository_status_handoff_20260409.md`：仓库 / worktree / 双机现场快照

除非本文件被新的 charter 明确替换，否则任何旧文档都不能单独把默认动作改成“先写论文”。

## 2. 实验优先门

只要当前 shortlist 中仍存在以下任一情况，就默认继续实验，不允许把“先写论文”替代为默认下一步：

1. 活跃 run 尚未收口
2. run 已收口但 gate 尚未正式判读
3. 白名单 family 中存在已经定义好的下一条实验

论文写作、supplementary 收口、图表整理仍可并行，但它们只占“实验空档”，不抢当前白名单队列的默认下一步。

## 3. 机器职责固定门

### 3.1 `6000 / A6000 x2`

- 默认承担 full-budget、正式 gate、长任务，或明显受益于更高吞吐的实验
- foundation 与 non-foundation 都可以在这里跑，是否 foundation 不是硬门
- 当前默认动作仍是不等 `6002`，而是维护自己的独立高价值队列
- `AlphaGenome matched raw-track slice v1` 已在 `2026-04-10 00:22:55 CST` 完成 technical/alignment gate closeout 并通过
- 若继续 AlphaGenome 这条线，唯一允许下一步是扩到约 `12-20 loci` 的稍大 matched panel；它仍是窄 external coordinate，不是大 benchmark

### 3.2 `6002 / RTX 3080`

- 默认承担 short10 cheap-screen、确认性 rerun、轻量结构探索、工具/数据小验证
- foundation 与 non-foundation 都可以在这里跑，是否 foundation 不是硬门
- 结果只能作为建议，不能自动变成 `6000` 发车命令

## 4. 当前允许与禁止的 family

### 4.1 白名单

当前只有以下三类实验允许自然推进：

1. 只有 genuinely new teacher / init hypothesis 的 foundation-assisted family 才允许重新进入；当前 tutorial `teacher-distill` 近邻变体不在自然推进范围内
2. `multiscale / dense decoder / U-Net-lite readout family`
3. `AlphaGenome matched raw-track slice` 这类窄外部对照

### 4.2 黑名单

当前默认禁止自然扩线的 family：

1. `summary/token/coarse residual/bins16 residual` 整条 measured adapter family
2. `msdec/skipprobe` 旧 readout 小变体链
3. 没有 genuinely new hypothesis 的 foundation bolt-on 变体

## 5. 当前双机队列定义

### 5.1 `6000 / A6000 x2`

截至 `2026-04-10 02:48:42 CST`，`teacher-distill tutorial` 线已停表，而 `6000` 已切到一条独立的双轨执行面：一条 A6000 formal gate + 一条 AlphaGenome sidecar。当前事实如下：

- `teacher-distill` 历史 formal gate verdict 仍为 `fail`
- `teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1`
- launch time：`2026-04-10 02:48:42 CST`
- machine / GPU：`6000 / A6000 x2`
- purpose：在保持 `corrected B` 语义不变的前提下，对 `multiscale/local-skip decoder v2` 做正式 gate
- hard compare target：只和历史 `corrected B` 比；不和 `6002`、不和旧 `msdec/skipprobe` 小变体链比
- log：`/data1/zhoujiazhen/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1.log`
- monitoring：`tail -f /data1/zhoujiazhen/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1.log`
- completion check：`grep -i "Finished" /data1/zhoujiazhen/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1.log`
- initial ETA window：`2026-04-10 07:48:42 CST` 到 `2026-04-10 10:48:42 CST`
- fresh launch evidence：发车后日志已进入真实 training step，`nvidia-smi` 已看到双 A6000 `96% util / 3664 MiB`
- `alphagenome_matched_raw_track_slice_v2_20260410`
- launch time：`2026-04-10 02:48:42 CST`
- machine / env：`6000 / alphagenome env`
- purpose：把已通过 v1 technical/alignment gate 的 AlphaGenome matched slice 扩到约 `12-20 loci` 的稍大 matched panel
- role：仍是窄 external coordinate / pilot，不是大 benchmark，也不作为 A6000 formal gate 的等待前置
- log：`/data1/zhoujiazhen/bylw_atac/TransChromBP/logs/alphagenome_matched_raw_track_slice_v2_20260410.log`
- monitoring：`tail -f /data1/zhoujiazhen/bylw_atac/TransChromBP/logs/alphagenome_matched_raw_track_slice_v2_20260410.log`
- completion check：`grep -i "Finished" /data1/zhoujiazhen/bylw_atac/TransChromBP/logs/alphagenome_matched_raw_track_slice_v2_20260410.log`
- fresh launch evidence：后台 PID 与 `run_alphagenome_pilot.py` 子进程已存在，输出目录已开始写入 `profiles/*.npz` 与 `region_metadata.jsonl`
- `alphagenome_matched_raw_track_slice_v1_20260410`
- start time：`2026-04-10 00:22:40 CST`
- end time：`2026-04-10 00:22:55 CST`
- hard verdict：`pass`
- 4 个固定 loci 全部成功，且每个位点都保留 `1` 条可比较的 `ATAC` track
- 已成功产出 `summary.csv`、`region_metadata.jsonl`、`profiles/*.npz`、`run_meta.json`、`merged_locus_totals.csv`
- 详细证据与边界结论见 `reports/alphagenome_matched_raw_track_slice_v1_closeout_20260410.md`

规则如下：

- `teacher-distill tutorial` 线已停表，不再为这条线分配后续 `6000` run
- `6000` 继续独立维护自己的高价值 backlog，不等待 `6002`
- 当前 `6000` 的 active GPU run 是 `teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1`
- 当前 `6000` 的 active sidecar 是 `alphagenome_matched_raw_track_slice_v2_20260410`
- 两条线互不等待；AlphaGenome sidecar 不构成 A6000 formal gate 的共享前置物
- A6000 formal gate 只和历史 `corrected B` 比，不和 `6002` 或旧 `msdec` 小变体链比较
- 明确禁止为了同步去复制 `6002` 的确认性 rerun 或镜像 `6002` 的 cheap-screen
- 后续若补充 `6000` 候选，必须是独立高价值任务，而不是旧 tutorial 线的续跑
- `AlphaGenome matched raw-track slice v1` 已经完成并作为已通过的 technical/alignment gate 基线保留；当前 active 的是 `v2` 扩面 sidecar
- 若继续 AlphaGenome 线，唯一允许下一步是当前这次 `12-20 loci` 的小幅扩面；禁止把这次 `pass` 或 `v2` 扩面改写成模型质量胜利或直接扩成大 benchmark

### 5.2 `6002 / RTX 3080`

截至 `2026-04-09 23:16:50 CST`，当前已确认的事实是：

- `U-Net-lite` 已在 `6002` 实际跑过 `r1/r2/r3`
- 经严谨性复核后，`r1/r2` 被确认为启动失败，当前唯一有效历史 run 是 `r3`
- 用户已批准补一条同配方确认性 cheap rerun：`teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r4`
- `6002` 当前仍在跑 `r4`
- 整条 `r4` 收口 ETA 窗口约为 `2026-04-10 00:50 CST` 到 `2026-04-10 01:30 CST`

规则如下：

- `6002` 的队列只保留本次 `r4` 收口与最终判读，不在 `r4` 结束前追加新 `6002` run
- `r1/r2/r3/r4` 的 rigor 背景保留，但不把 `r4` 结束前的任何外部 backlog 当成已发车项
- `r4` 结束后再基于最终判读决定下一步；在此之前不追加新的 `6002` run
- 若 `r4` 与 `r3` 仍落在同一区间，则 `U-Net-lite v1` 正式收成 `no-go` 并停表
- 若 `r4` 明显改善，则唯一允许下一步是一次确认性 cheap rerun；确认性 rerun 通过后，才可把它保留为新的 `3080 shortlist` 候选

### 5.3 双机协同门

- 两台机器各自维护 backlog，不等待对方 run 收口才决定自己的下一条
- 一台机器的结果只能提供建议，不能成为另一台机器发车的必需前提
- 只有共享前置物缺失时，才允许出现等待
- 明确禁止为了同步而做镜像 run 与共用唯一下一步

## 6. 晋级与停表规则

每条白名单实验线都必须写死：

1. 过 gate 后唯一允许的下一步
2. 不过 gate 后必须停到哪里
3. 明确禁止“没过也再试一个小改动”

当前两条 active family 的 stop-rule 已经在第 5 节固定，不允许再回到“先小修一版试试看”的松散模式。

`AlphaGenome matched raw-track slice` 的最小规则如下：

- 这是一条 `6000` 独立 backlog 的窄 external coordinate / pilot，不是大 benchmark
- 过 gate 的最小条件：在固定 matched loci 面板上成功产出可用的 AlphaGenome raw-track 输出；“可用”指所有选定位点都生成 summary / metadata / profile 结果，ontology / filter 后仍保留至少一条可比较的 ATAC track，且没有因为 API / track matching 失败让大部分位点失效
- 过 gate 后唯一下一步：扩大到约 `12-20 loci` 的稍大 matched panel
- 不过 gate 后：这条 pilot 停表，不扩成大 benchmark，`6000` 转去下一个 genuinely new high-value family

`6002 U-Net-lite v1` 的最小规则如下：

- 当前只保留 `r4` 收口与最终判读
- 若 `r4` 与 `r3` 同区间，则 `U-Net-lite v1` 正式收成 `no-go` 并停表
- 若 `r4` 明显改善，则唯一允许下一步是一次确认性 cheap rerun

## 7. 论文让路条件

只有在双机各自 backlog 都没有合法下一条时，论文类事项才允许成为全局默认下一步：

1. 双机各自都没有合法下一条实验
2. 当前白名单 family 已全部触发停表
3. 双机都没有共享前置物之外的可发车项

某台机器正在跑，不等于另一台机器的 backlog 自动清空。

论文工作只能占各自机器的空档，不再通过“另一台机器在跑”来顶掉 backlog。

## 8. 文档同步门

每次实验收口或 stop/go 变化后，必须同轮更新三处：

1. 本文件里的当前队列状态
2. `TRACKING.md` 顶部 live 行
3. 对应 handoff / plan / report

只更新口头回复、不更新文档，不算完成；但同步动作是为了保持事实一致，不是要求两台机器必须串行等对方收口。

## 9. 新实验准入门

任何想插入当前白名单之外的新实验，必须先在本文件里补齐以下 5 个字段：

1. hypothesis
2. 机器 / GPU
3. 为什么不属于当前黑名单 family
4. 过 gate 后唯一下一步
5. 不过 gate 后停在哪里

未补齐前，不允许发车。

## 10. 当前默认读法

当前阶段的默认读法固定为：

1. `6000` 当前 active GPU run 是 `teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1`；它是双卡 A6000 formal gate，只和历史 `corrected B` 比，不等待 `6002`
2. `6000` 当前 active sidecar 是 `alphagenome_matched_raw_track_slice_v2_20260410`；它是基于 `v1 pass` 的 `12-20 loci` 小幅扩面，仍只算窄 external coordinate，不升级成大 benchmark，也不和 A6000 gate 互相等待
3. `6002` 继续按自己的 `U-Net-lite r4` rigor closeout 线收口；一台机器的正负结果都不能自动改写另一台机器当前已发车项的默认下一步
