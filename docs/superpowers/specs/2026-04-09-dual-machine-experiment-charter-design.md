# Dual-Machine Experiment Charter Design

## Goal

把当前双机实验决策从“实验与论文并行、容易临时切回写作”改成“实验队列强约束优先”，使 `/home/zhengwei/project/python/TransChromBP` 后续在 `6000` 与 `6002` 上的执行默认遵循同一套白名单、黑名单、stop-rule 与插队规则。

## Why This Exists

当前仓库里已经有三层足够明确的外部路线共识：

1. 旧的 measured foundation adapter family 默认停表；
2. 下一代优先级应放在 `multiscale / dense decoder / U-Net-like readout` 与 `foundation as teacher/distill`；
3. 论文主线应收口，但不应继续反过来吞掉尚未跑完的实验 shortlist。

现在缺的不是更多口头共识，而是一份会压过“先写论文”默认分叉的硬规则文件。

## Scope

本设计只覆盖三类内容：

1. 当前双机实验的全局硬规则；
2. `6000` 与 `6002` 的职责分工和允许队列；
3. 文档同步与新实验准入门。

## Non-Goals

本设计不做以下事情：

1. 不重写旧实验结论；
2. 不直接替代 `TRACKING.md` 的 live 状态职责；
3. 不引入复杂打分系统、队列评分器或自动调度器；
4. 不把所有未来实验都预先写死，只约束当前阶段的合法动作边界。

## Deliverables

实施时只改两处：

1. 新增 `docs/plan/2026-04-09_dual_machine_experiment_charter.md`
2. 更新 `TRACKING.md` 顶部引用与 live 行，使其把“规则源”显式指向该 charter

## Authority Model

新 charter 建成后，实验优先级判断遵循以下权威顺序：

1. 最新实验事实与日志
2. `docs/plan/2026-04-09_dual_machine_experiment_charter.md`
3. `TRACKING.md` live 状态
4. 其他 handoff / plan / report

这意味着：

- `TRACKING.md` 不再自行定义“默认下一步”
- `repository_status_handoff_20260409.md` 与 `dual-track` 计划文档只能引用 charter
- 若某份旧文档与 charter 冲突，以 charter 为准

## Hard Rules

### 1. Experiment-First Gate

只要当前 shortlist 中仍存在以下任一情况，就默认继续实验，不允许把“先写论文”替代为默认下一步：

1. 活跃 run 尚未收口
2. run 已收口但 gate 尚未正式判读
3. 白名单 family 中存在已经定义好的下一条实验

### 2. Fixed Machine Roles

`6000 / A6000x2`：

- 只负责高价值正式 gate
- 当前优先承担 `NT v2 teacher-distill` 这类需要双卡结论的线
- 不负责 readout 小变体 cheap-screen

`6002 / RTX 3080`：

- 只负责 cheap-screen、readout probe、窄外部对照
- 不承担 `6000` 的 full-budget verdict 职责
- 结果只能作为候选晋级信号，不能自动升级为 `6000` 发车命令

### 3. Allowed Families

当前唯一允许自然推进的 family 只有三类：

1. `NT v2 teacher-distill`
2. `multiscale / dense decoder / U-Net-lite readout family`
3. `AlphaGenome matched raw-track slice` 这类窄外部对照

不在白名单中的实验，默认不合法。

### 4. Blocked Families

当前默认禁止自然扩线的 family：

1. `summary/token/coarse residual/bins16 residual` 整条 measured adapter family
2. `msdec/skipprobe` 旧 readout 小变体链
3. 没有 genuinely new hypothesis 的 foundation bolt-on 变体

### 5. Promotion And Stop Rules

每条白名单实验线都必须明确：

1. 过 gate 后唯一允许的下一步
2. 不过 gate 后必须停到哪里
3. 明确禁止“没过也再试一个小改动”

### 6. Writing-Yields-To-Experiment Rule

只有在以下任一条件满足时，论文类事项才允许成为默认下一步：

1. 双机都没有合法下一条实验
2. 当前白名单 family 已全部触发停表
3. 当前 run 正在跑且没有合理并行 cheap-screen 可开

## Current Queue Definition

### `6000 / A6000x2`

当前唯一合法 active run：

- `ntv2_teacher_distill_short10_s42_6000_20260409_r2`

后续规则：

- 若 short10 过 gate：唯一下一步是同配方升 `full-budget distill`
- 若 short10 不过 gate：`teacher-distill tutorial` 线直接停表
- 明确禁止因为论文写作或旧 family 惯性而跳去 `msdec/skipprobe/residual/cross-attention`

### `6002 / RTX 3080`

当前已知事实：

- `U-Net-lite` 已经至少实际跑过 `r1/r2/r3`
- 当前 GPU 空闲
- 下一步不是“先写论文”，而是把 `U-Net-lite v1` 线收成正式 verdict

后续规则：

- 若 `U-Net-lite v1` 过 gate：唯一下一步是同 family 的一次确认性 cheap rerun
- 若 `U-Net-lite v1` 不过 gate：下一条仍必须来自 `dense decoder / readout family` 的 genuinely new 变体
- 只有在没有 ready 的下一个 readout 变体时，才允许把 `6002` 转去做 `AlphaGenome matched raw-track slice`

## Synchronization Rules

每次实验收口或 stop/go 变化后，必须同轮更新三处：

1. charter 里的当前队列状态
2. `TRACKING.md` 顶部 live 行
3. 对应 handoff / plan 文档

若只更新口头回复、不更新文档，则不算完成。

## New Experiment Admission Gate

任何想插入当前白名单之外的新实验，必须先在 charter 中补齐以下 5 个字段：

1. hypothesis
2. 机器 / GPU
3. 为什么不属于当前黑名单 family
4. 过 gate 后唯一下一步
5. 不过 gate 后停在哪里

未补齐前，不允许启动。

## Success Criteria

实施完成后，应满足：

1. `TRACKING.md` 不再把“先写论文”写成默认下一步
2. 后续讨论双机实验时能直接引用 charter，而不是重新解释优先级
3. 即便 `6000 teacher-distill` 判负，也不会自动切回写论文，而是继续执行 `6002` 白名单队列
4. 新实验若不在白名单内，必须先补准入字段才可发车

## Risks And Tradeoffs

### 风险

1. 规则太弱，会再次被旧 handoff 文档冲淡
2. 规则太强，如果白名单写得过窄，可能阻塞 genuinely new hypothesis

### 取舍

当前优先选择“规则偏强”，因为现在真正的问题不是实验太多，而是实验队列经常被论文写作默认打断。

## Implementation Preference

本设计的推荐实现是：

1. 新建 `docs/plan/2026-04-09_dual_machine_experiment_charter.md` 作为单一权威规则源
2. 在 `TRACKING.md` 顶部显式引用它
3. 保持 `TRACKING.md` 只写 live 状态，不再承担路线裁决职责
