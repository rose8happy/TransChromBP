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
- 当前没有 active training run；`teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1` 已在 `2026-04-10 10:10:10 CST` 左右完成 formal gate closeout，并判定 `fail`
- `alphagenome_matched_raw_track_slice_v1_20260410` 已在 `2026-04-10 00:22:55 CST` 完成 technical/alignment gate closeout 并通过
- `alphagenome_matched_raw_track_slice_v2_20260410` 已在 `2026-04-10 03:34:42 CST` 完成 technical / external-coordinate sidecar closeout 并通过
- AlphaGenome v1/v2 都已经收口；它们现在只算窄 external coordinate，不再占 active slot，也不自动升级成 benchmark

### 3.2 `6002 / RTX 3080`

- 默认承担 short10 cheap-screen、确认性 rerun、轻量结构探索、工具/数据小验证
- foundation 与 non-foundation 都可以在这里跑，是否 foundation 不是硬门
- 当前没有 active training run；`teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r4` 已在 `2026-04-10 03:34:42 CST` 左右收口，`r1/r2` 仍是无效启动，`r3/r4` 是有效 run，但最终 verdict 仍是 `no-go / stop`
- 结果只能作为建议，不能自动变成 `6000` 发车命令

## 4. 当前允许与禁止的 family

### 4.1 白名单

当前没有 active family。若未来要重新推进任何新 run，必须先通过第 9 节新增准入。历史上曾允许自然推进、但现在都已收口的 family 仅作为参考保留：

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

截至 `2026-04-10 10:10:10 CST`，`6000` 没有 active training run。当前事实如下：

- `teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1`
  - start time：`2026-04-10 02:48:42 CST`
  - end time：`2026-04-10 10:10:10 CST` 左右
  - hard verdict：`fail`
  - best epoch：`22`
  - best peak `profile_target_jsd_full_mean=0.3337811803218466`
  - best peak `count_pearson_full=0.7948227405497477`
  - formal gate 只和 corrected B comparator 比，不和 `6002` 或旧 `msdec/skipprobe` 链比较
- `alphagenome_matched_raw_track_slice_v1_20260410`
  - hard verdict：`pass`
  - 4 个固定 loci 全部成功，且每个位点都保留 `1` 条可比较的 `ATAC` track
  - 已成功产出 `summary.csv`、`region_metadata.jsonl`、`profiles/*.npz`、`run_meta.json`、`merged_locus_totals.csv`
- `alphagenome_matched_raw_track_slice_v2_20260410`
  - hard verdict：`pass`
  - `16` 个 loci 全部成功落盘，且每个位点都保留 `1` 条可用 `ATAC` track
  - 这是已完成的 external-coordinate sidecar，不占 active slot

规则如下：

- `6000` 现在没有 active run，不再等待 `6002`
- 明确禁止为了同步去复制 `6002` 的确认性 rerun 或镜像 `6002` 的 cheap-screen
- 后续若补充 `6000` 候选，必须是显式新 hypothesis，而不是当前已关闭 run 的续跑
- `AlphaGenome matched raw-track slice v1` / `v2` 都已经完成，不再按“待启动 / active run”口径描述
- AlphaGenome 的两个 closeout 都只算 external coordinate，不自动改写 A6000 formal gate 语义

### 5.2 `6002 / RTX 3080`

截至 `2026-04-10 03:34:42 CST`，`6002` 没有 active training run。当前已确认的事实是：

- `U-Net-lite` 已在 `6002` 实际跑过 `r1/r2/r3/r4`
- 经严谨性复核后，`r1/r2` 被确认为启动失败，`r3/r4` 是有效 run
- `r4` 已完成确认性复跑，但仍不足以把 `U-Net-lite v1` 推上 shortlist
- 当前 family verdict 保持 `no-go / stop`

规则如下：

- `6002` 现在没有 active run，不再追加同配方 `U-Net-lite v1` run
- `r1/r2/r3/r4` 的 rigor 背景保留，但不把任何已关闭 run 当成当前 active backlog
- 若未来要再动 `6002`，必须是显式新 hypothesis，而不是当前 `U-Net-lite v1` 的续跑
- 当前没有任何理由把这条 family 写成正向 shortlist

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
- 当前 v1/v2 都已完成；后续若再推进，必须作为新的显式 hypothesis 重新准入，而不是自动续跑
- 不过 gate 后：这条 pilot 停表，不扩成大 benchmark，也不再占 active slot

`6002 U-Net-lite v1` 的最小规则如下：

- 当前只保留已完成的 `r1/r2/r3/r4` closeout 结论
- `U-Net-lite v1` 已正式收成 `no-go / stop`
- 若未来要继续，只能是一个新的显式 hypothesis，而不是同配方补跑

## 7. 论文让路条件

只有在双机各自 backlog 都没有合法下一条时，论文类事项才允许成为全局默认下一步：

1. 双机各自都没有合法下一条实验
2. 当前白名单 family 已全部触发停表
3. 双机都没有共享前置物之外的可发车项

现在两台机器都没有 active training run，因此论文工作、归档整理和新 hypothesis 准入可以并行推进；任何未来新 run 都必须单独经过本 charter 的准入门。

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

1. `6000` 没有 active training run；`teacher_v2_center_pool_msdls_v2_30ep_s42_6000_20260410_r1` 已正式 `fail`，`alphagenome_matched_raw_track_slice_v1_20260410` / `v2_20260410` 都已完成 closeout
2. `6002` 没有 active training run；`U-Net-lite v1` 已正式收成 `no-go / stop`
3. 现在默认动作不是“等对方收口”，而是按各自已关闭的事实继续做文档收口、归档和新 hypothesis 准入
