# 2026-04-11 Factor Ladder 无人值守串行矩阵设计

## 1. 设计目标

本设计服务于一个明确约束：

- 接下来一段时间，用户不会持续盯盘；
- 但 `alphagenome_factor_ladder` 这条实验线仍希望继续推进，不希望每跑完一段就停下来等待人工确认。

因此，本设计的目标不是“让所有实验立刻并发开满”，而是：

1. 在**无人值守**前提下，把 factor-ladder 的最小高价值矩阵按固定顺序继续跑完。
2. 不因中间指标好坏停表，只因**技术故障**停表。
3. 让用户回来时，能够用最少上下文成本看懂：
   - 当前跑到哪一段；
   - 哪一段完成了；
   - 如果停了，是停在哪、为什么停。

本设计的成功标准是：

- 形成一条可执行的单机串行队列，能从当前正在运行的 `E2 short10` 开始，继续接完后续矩阵；
- 过程中所有 stage 都有清晰的启动条件、技术成功条件、技术熔断条件和产物约定；
- 队列状态被同步写回 canonical 文档，而不是只存在于远端 shell 历史里。

## 2. 当前上下文

当前 `alphagenome_factor_ladder` 的 live 状态如下：

1. `6000` remote isolated repo：`/data1/zhoujiazhen/bylw_atac/.codex_jobs/alphagenome_factor_ladder_20260411`
2. remote runtime commit：`5877077`
3. canonical archival commit：`079b603`
4. 当前正在运行的 formal gate：
   - `teacher_v2_hierdec4096_short10_s42_6000_20260411_r1`
   - dual A6000
   - formal `E2 short10`

这意味着：

1. 队列头已经存在，不能通过“重启整个矩阵”来实现无人值守；
2. 新设计必须支持**接管后续 stage**，而不是重启当前 `E2`；
3. `E3` 依赖 `E2 teacher30` 的 teacher checkpoint，因此矩阵天然不是完全独立的全交叉组合。

## 3. 方案选择

本轮考虑 3 类调度方案。

### 3.1 方案 A：6000 单机串行队列

- 当前 `E2 short10` 跑完后，依次接：
  - `E1 short10`
  - `E2 teacher30`
  - `teacher cache export`
  - `E3 distill short10`

优点：

- 最稳，最适合无人值守；
- 日志、产物、状态文件都集中在 `6000` 单一运行域；
- 不会引入 `6000/6002` 环境漂移和双机文档分叉。

缺点：

- 总墙钟时间较长；
- 不能利用 `6002` 并行加速。

### 3.2 方案 B：双机拆分并行

- `6000`：`E2 -> teacher30 -> E3`
- `6002`：独立跑 `E1`

优点：

- 完成更快。

缺点：

- 无人值守时更容易出现环境差异、路径不同步、文档回写分叉；
- `6002` 不是当前 factor-ladder 主运行域，恢复和排障成本更高。

### 3.3 方案 C：更大的全交叉矩阵

- 在上述基础上再加多 seed、补长 budget comparator、或者额外结构分支。

优点：

- 覆盖面最大。

缺点：

- 明显超出“最小高价值矩阵”的范围；
- 信息密度开始下降，不适合当前无人值守约束。

### 3.4 选型结论

本设计采用 **方案 A：6000 单机串行队列**。

原因：

1. 它最符合“用户暂时不盯盘”的现实约束；
2. 它保留了 factor-ladder 的解释性；
3. 它把复杂度压在一个运行域内，更利于状态沉淀和后续 closeout。

## 4. 固定范围与非目标

### 4.1 本设计覆盖的矩阵范围

固定串行队列只包含 4 段：

1. 当前已在跑的 `E2 short10`
2. `E1 short10`
3. `E2 teacher30`
4. `E3 distill short10`

其中 `teacher cache export` 不单独算作一个实验 family run，但算作一个必须显式建模的中间 stage。

### 4.2 本设计明确不做的事

1. 不把 `6002` 纳入同一条无人值守主队列
2. 不增加额外 seed
3. 不扩成 `8192/16384` 第二层矩阵
4. 不因为中间结果“看起来不好”就自动停表
5. 不因为“看起来很好”就自动扩线

## 5. 队列顺序

队列的固定顺序如下：

1. `S0`: 等待当前 `E2 short10` 完成
   - run id: `teacher_v2_hierdec4096_short10_s42_6000_20260411_r1`

2. `S1`: 启动 `E1 short10`
   - 目标：测 `longctx4096 + corrected-B readout`

3. `S2`: 启动 `E2 teacher30`
   - 目标：生成可用于 distill 的 teacher checkpoint

4. `S3`: 导出 teacher cache
   - 目标：生成 `train/valid` manifests 与对应 teacher targets

5. `S4`: 启动 `E3 distill short10`
   - 目标：测 distillation 的边际价值

这里故意把 `E1` 放在当前 `E2` 之后，而不是强行前置，原因只有一个：

- 当前 `E2` 已经在跑，强行改顺序没有收益；更合理的做法是接管其后续阶段。

## 6. “无条件全跑”的精确定义

本设计中的“无条件全跑”定义为：

- **科学上不停表**：
  - 不因为某个 run 的指标是 `positive / flat / negative` 而停；
  - 不因为某个阶段“看起来已经足够说明问题”而停。

- **技术上允许熔断**：
  - 只有在技术故障导致下一阶段无法安全继续时，队列才停。

这一定义是本设计的核心。它允许无人值守，同时避免把“无人值守”误解成“无论发生什么都盲跑到底”。

## 7. 每段 stage 的触发条件与技术成功条件

### 7.1 `S0 -> S1`：当前 `E2 short10` 结束后启动 `E1`

触发条件：

1. 当前 `E2 short10` 日志出现完成标记；
2. 对应 metrics / meta 文件已落盘；
3. 至少有一个可读 checkpoint 存在。

技术成功条件：

1. run 正常退出；
2. 日志和产物完整；
3. 不出现 `OOM / NaN / DDP crash`。

### 7.2 `S1 -> S2`：`E1 short10` 完成后启动 `E2 teacher30`

触发条件：

1. `E1 short10` 正常完成；
2. 只看技术成功，不看指标判优。

技术成功条件：

1. `E1` run 正常退出；
2. `epoch_metrics.jsonl` 和 `run_meta.json` 存在；
3. checkpoint 存在。

### 7.3 `S2 -> S3`：`E2 teacher30` 完成后导出 teacher cache

触发条件：

1. `E2 teacher30` 正常完成；
2. 至少存在一个 teacher checkpoint。

checkpoint 选择策略固定为：

1. 优先 `best.pt`
2. 如果没有 `best.pt`，退到最新 `epoch_*.pt`
3. 两者都没有则停表

### 7.4 `S3 -> S4`：teacher cache 导出成功后启动 `E3`

触发条件：

1. `train` / `valid` 两份 teacher manifest 都存在；
2. `n_records` 校验通过；
3. `teacher_target_names` 与 distill config 对齐。

技术成功条件：

1. cache 文件完整；
2. manifest 与数据配置一致；
3. distill student 能正常启动。

## 8. 技术熔断条件

以下情况允许并要求队列停表：

### 8.1 训练类熔断

1. `OOM`
2. `loss=nan` 或明显数值发散
3. `DDP` 初始化失败
4. rank 异常退出
5. 数据路径不存在

### 8.2 checkpoint / export 类熔断

1. `E2 teacher30` 结束后没有任何可读 checkpoint
2. teacher cache 缺失 `train` 或 `valid` manifest
3. `n_records` mismatch
4. teacher targets 与 distill config 不匹配

### 8.3 自动恢复的唯一例外

只允许一个窄场景自动重试：

- `master_port` 冲突

恢复策略：

1. 自动换一个预留端口再试一次
2. 若再次失败，则停表

除此之外，不做无限自动重试。

## 9. 队列执行器设计

### 9.1 形态

本设计不采用一条超长 `bash && bash && bash` 链，而采用一个薄的 queue sidecar。

理由：

1. 无人值守时，最重要的是“停在哪里能看懂”，而不是“命令最短”；
2. sidecar 更适合显式记录 stage 状态和失败原因；
3. 它能接管当前已经在跑的 `E2`，而不是强制重新启动矩阵。

### 9.2 组成

1. 一个远端队列控制脚本
   - 建议路径：`scripts/run_factor_ladder_unattended_matrix.sh`
2. 一份显式 queue config
   - 建议路径：`configs/queues/factor_ladder_unattended_20260411.yaml`
3. 一组队列状态文件
   - 建议目录：`outputs/queue/factor_ladder_unattended_20260411/`

### 9.3 控制脚本职责

控制脚本只负责：

1. 监视当前 stage 是否完成
2. 判断是否满足下一阶段启动条件
3. 启动下一阶段
4. 记录状态事件
5. 在技术故障时停表

它不负责重新实现训练逻辑；训练仍通过现有 launcher / `torchrun` 路径执行。

## 10. 队列状态文件设计

建议至少保留以下文件：

1. `queue_state.json`
   - 当前 stage
   - 整体状态：`waiting / running / halted / completed`
   - 最近更新时间
   - 当前关注的 log 路径

2. `events.jsonl`
   - 追加写入每个关键事件：
     - `stage_started`
     - `stage_passed`
     - `stage_failed`
     - `queue_halted`
     - `queue_completed`

3. `summary.md`
   - 给人看，不给程序看
   - 用户回来后优先打开这一份

4. `queue.log`
   - orchestrator 自己的短日志
   - 只记录“何时启动了谁、为何停表”

## 11. run 与日志组织

每个正式 run 继续保留自己的独立日志，不混成一个总日志。

建议 run 命名如下：

1. 当前 `E2 short10`
   - `teacher_v2_hierdec4096_short10_s42_6000_20260411_r1`
2. `E1 short10`
   - `teacher_v2_center_pool_longctx4096_short10_s42_6000_20260411_r1`
3. `E2 teacher30`
   - `teacher_v2_hierdec4096_teacher30_s42_6000_20260411_r1`
4. `E3 distill short10`
   - `teacher_v2_hierdec4096_distill_short10_s42_6000_20260411_r1`

这样做的目的：

1. run manifest 容易登记
2. 单条日志容易排障
3. closeout 时不需要从一个巨型 orchestrator log 里反推训练细节

## 12. 启动前一次性检查

真正启动队列前，必须完成以下预检：

1. 当前 `E2 short10` 的完成探针可靠
2. 后续 stage 的 launcher 都通过 `bash -n`
3. teacher checkpoint 选择规则写死
4. teacher cache 成功判据写死
5. 每段启动前都做 GPU 资源检查
6. 预留一组 `master_port` fallback
7. 输出目录和磁盘空间预检通过
8. 队列锁文件存在，防止重复启动

其中第 1 条最关键：

- 队列的第一段不是“重新启动 `E2`”，而是“等待现有 `E2` 正常收口”。

## 13. 文档同步规则

无人值守不意味着文档可以滞后。

本设计要求：

1. 在 canonical 文档里登记这条 unattended queue
2. 在 `runs.csv` 中把后续 planned runs 先注册为 `queued`
3. stage 真正启动后，再把对应 row 改成 `running`
4. 若技术熔断，则立即回写：
   - 停在第几段
   - 最后一段 log 路径
   - 失败类型
5. 整条队列跑完后，再做 formal gate / closeout 判读

## 14. 预计时间窗

按当前 `6000` 双卡吞吐与既有 `short10` 历史，先记录以下保守时间估计：

1. 当前 `E2 short10`
   - 预计结束：`2026-04-11 23:20-23:45 CST`
2. `E1 short10`
   - 约 `1.5-2.5h`
3. `E2 teacher30`
   - 约 `4.5-6.5h`
4. `teacher cache export`
   - 约 `0.5-1.5h`
5. `E3 distill short10`
   - 约 `1.5-2.5h`

整条无人值守矩阵的总完成时间窗先记为：

- `2026-04-12 08:00-12:00 CST`

这只是排程估计，不是科学结论。

## 15. 设计结论

本设计把 factor-ladder 的下一步推进定义成：

- 一个以 `6000` 为唯一执行域、
- 从当前 active `E2 short10` 开始接管、
- 科学上不停表、技术上才熔断、
- 且带有显式状态文件和 canonical 文档同步的无人值守串行矩阵。

这套设计的价值不在于“最省命令”，而在于：

1. 用户暂时离开后，实验不会失去推进能力；
2. 任何中断都能被明确定位；
3. 用户回来时，可以低成本恢复上下文并继续做 formal 判读。
