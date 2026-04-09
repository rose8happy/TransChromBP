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

- 只负责高价值正式 gate
- 当前优先承担 `NT v2 teacher-distill` 这类需要双卡正式 verdict 的线
- 不承担 readout 小变体 cheap-screen

### 3.2 `6002 / RTX 3080`

- 只负责 cheap-screen、readout probe、窄外部对照
- 不承担 `6000` 的 full-budget verdict 职责
- 结果只能作为晋级候选，不能自动变成 `6000` 发车命令

## 4. 当前允许与禁止的 family

### 4.1 白名单

当前只有以下三类实验允许自然推进：

1. `NT v2 teacher-distill`
2. `multiscale / dense decoder / U-Net-lite readout family`
3. `AlphaGenome matched raw-track slice` 这类窄外部对照

### 4.2 黑名单

当前默认禁止自然扩线的 family：

1. `summary/token/coarse residual/bins16 residual` 整条 measured adapter family
2. `msdec/skipprobe` 旧 readout 小变体链
3. 没有 genuinely new hypothesis 的 foundation bolt-on 变体

## 5. 当前双机队列定义

### 5.1 `6000 / A6000 x2`

截至 `2026-04-09 22:11 CST`，当前唯一合法 active run 是：

- `ntv2_teacher_distill_short10_s42_6000_20260409_r2`

规则如下：

- 若 short10 过 gate：唯一允许的下一步是同配方升 `full-budget distill`
- 若 short10 不过 gate：`teacher-distill tutorial` 线直接停表
- 明确禁止因为旧路线惯性而跳去 `msdec/skipprobe/residual/cross-attention`
- 即便 short10 判负，也不能仅因为“论文更该收口”就跳过 `6002` 的白名单队列

### 5.2 `6002 / RTX 3080`

截至 `2026-04-09 22:11 CST`，当前已确认的事实是：

- `U-Net-lite` 已在 `6002` 实际跑过 `r1/r2/r3`
- `6002` 当前 GPU 空闲
- 当前最先要做的不是“先写论文”，而是把 `U-Net-lite v1` 收成正式 verdict

规则如下：

- 若 `U-Net-lite v1` 过 gate：唯一允许的下一步是同 family 的一次确认性 cheap rerun
- 若 `U-Net-lite v1` 不过 gate：下一条仍必须来自 `dense decoder / readout family` 的 genuinely new 变体
- 只有在没有 ready 的下一个 readout 变体时，才允许把 `6002` 转去做 `AlphaGenome matched raw-track slice`

### 5.3 双机协同门

- `6000` 负责正式 gate，`6002` 负责 cheap-screen，不互相镜像跑同一条 full-budget
- `6002` 的正结果只能提供晋级候选，不能自动触发 `6000`
- `6000` 的负结果也不能直接把默认动作改成“开始写论文”；只要 `6002` 白名单里还有下一条，就继续实验

## 6. 晋级与停表规则

每条白名单实验线都必须写死：

1. 过 gate 后唯一允许的下一步
2. 不过 gate 后必须停到哪里
3. 明确禁止“没过也再试一个小改动”

当前两条 active family 的 stop-rule 已经在第 5 节固定，不允许再回到“先小修一版试试看”的松散模式。

## 7. 论文让路条件

只有在以下任一条件满足时，论文类事项才允许成为默认下一步：

1. 双机都没有合法下一条实验
2. 当前白名单 family 已全部触发停表
3. 当前 run 正在跑且没有合理并行 cheap-screen 可开

换句话说，论文工作默认是填空档，不是拿来顶掉仍然合法的实验队列。

## 8. 文档同步门

每次实验收口或 stop/go 变化后，必须同轮更新三处：

1. 本文件里的当前队列状态
2. `TRACKING.md` 顶部 live 行
3. 对应 handoff / plan / report

只更新口头回复、不更新文档，不算完成。

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

1. 先等 `6000 teacher-distill short10` 收口并正式判 gate
2. 先把 `6002 U-Net-lite r1/r2/r3` 收成 verdict，而不是把它继续写成“尚未起跑”
3. 即便 `6000` 判负，只要 `6002` 白名单队列还在，就继续实验，不默认切回“先写论文”
