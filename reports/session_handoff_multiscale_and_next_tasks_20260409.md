# 当前工作接续说明：双机实验排程快照（2026-04-09）

> **Superseded / 历史快照**
>
> 这份文档记录的是 `2026-04-09` 较早时点的一次双机排程快照，**不能再作为当前调度入口**。
> 当前规则源和 live 入口请以以下文件为准：
> - [docs/plan/2026-04-09_dual_machine_experiment_charter.md](../docs/plan/2026-04-09_dual_machine_experiment_charter.md)
> - [TRACKING.md](../TRACKING.md)
> - [reports/repository_status_handoff_20260409.md](./repository_status_handoff_20260409.md)

## 目的

这份快照用于替代仍以 `20260407` 命名、但内容已被多次追加更新的旧 handoff 文档，给 `2026-04-09` 当前时点一个更直观的双机排程入口。

它只回答四个问题：

1. `6000` 和 `6002` 现在是否仍有任务在跑？
2. 两台机器当前各自承担什么职责？
3. 哪些 family 已经收口，不应再按“默认下一条”盲扩？
4. 如果今天继续推进，双机各自最合理的下一步是什么？

若本文件与更新的实机日志、`TRACKING.md` 或新实验事实冲突，以更新事实为准。

---

## 一、当前实机状态（`2026-04-09 02:40 CST` 复核）

| 机器 | GPU 状态 | 最新收口 run | 当前结论 |
|---|---|---|---|
| `6000 / A6000 x2` | 双卡活跃（约 `3630 MiB / 49140 MiB`, `94-95% util`） | `teacher_v2_center_pool_msdec_v1_s2_s42_20260409_full1` | 已按既有 backlog 实际启动 `msdec_v1_s2 full-budget` 双卡 run，用于把 `msdec_v1` family 的 decoder-depth ablation 正式收口；预计完成时间窗 `2026-04-09 09:30-10:30 CST` |
| `6002 / RTX 3080 x1` | 单卡活跃（约 `4123 MiB / 12288 MiB`, `98% util`） | `teacher_v2_center_pool_skipprobe_v1_wide_short10_s42_6002_r1` | 已补齐并启动 `skipprobe_wide short10`，用于给 `6002` 增加一条不同于 `s2` 的 `skip_probe_v1` 高带宽 cheap ablation；预计完成时间窗 `2026-04-09 07:30-08:30 CST` |

补充证据：

- `6000` 当前 active run：`teacher_v2_center_pool_msdec_v1_s2_s42_20260409_full1`
- `6000` 日志：`/data1/zhoujiazhen/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_msdec_v1_s2_s42_20260409_full1.log`
- `6000` launcher pid：`723040`
- `6000` 日志已写出 launcher header、`torchrun` 警告、dataset setup 与 `epoch=1` 训练行，说明不是秒挂
- `6002` 当前 active run：`teacher_v2_center_pool_skipprobe_v1_wide_short10_s42_6002_r1`
- `6002` 日志：`/home/zhengwei/bylw_atac/TransChromBP/logs/teacher_v2_center_pool_skipprobe_v1_wide_short10_s42_6002_r1.log`
- `6002` wrapper pid：`1729594`
- `6002` 日志已写出 launcher header、dataset setup 与 single-card setup 行，说明新配置已被实际消费且不是秒挂
- 当前双机状态已切换为：`6000` 跑双卡 full-budget verdict，`6002` 跑单卡 cheap-screen，重新回到真正的双机并行

---

## 二、当前有效的双机分工

### 1. `6000 / A6000 x2`

职责固定为：

- 承担 `full-budget / 正式 verdict` 级别的高价值实验
- 只按 `6000` 自己的 paper-facing gate 和 A6000 信息增益排序
- 不等待 `6002` 的单卡结果才决定自己下一条

当前不应再默认做的事：

- 不把 `msdec_v1_s2 full-budget` 当作无条件自然下一条
- 不继续在 `msdec_v1 full1 -> narrow -> s2` 这条小变体链上盲扩
- 不回到任何 `ntv2_*` / foundation adapter family

### 2. `6002 / RTX 3080 x1`

职责固定为：

- 承担同机口径下的 `cheap-screen / 单卡结构排序`
- 只回答“哪条 readout 变体值得留在 6002 shortlist”
- 结果可为 `6000` 提供旁证，但不能直接触发或阻塞 `6000` 的启动

当前不应再默认做的事：

- 不做 `6000` 的镜像 full-budget 对照
- 不回头跑 `ntv2_*` / foundation adapter family
- 不因为 `skipprobe` 与 `s2` 都收口了，就继续在同一条 `msdec` 小变体链上无门槛加 run

---

## 三、已经稳定收口的判断

### 1. `foundation adapter family`

当前稳定读法仍是：

- `summary fusion`
- `token fusion`
- `coarse residual`
- `bins16 center-aligned residual`

这条 family 默认停表。若未来重开，必须提出一个不属于当前 adapter family 的新 hypothesis，而不是“再试一个小变体”。

### 2. `6000` 的当前 `msdec_v1` 主线

已知事实：

- `teacher_v2_center_pool_msdec_v1_s42_20260407_full1` 已正式不过 gate
- `teacher_v2_center_pool_msdec_v1_narrow_s42_20260408_full1` 也已收口，当前读数相对历史 `corrected B` 仍明显偏弱

因此当前默认动作不是继续顺着 `msdec_v1` 小变体链发新 run，而是：

1. 先把 `narrow` 相对 `corrected B` 的正式 gate 与 closeout 写死
2. 若 verdict 仍为负，则切到新的 `A6000 readout family` 讨论与排程

### 3. `6002` 的当前 shortlist

当前已完成并可复用的单卡结果：

- `anchor`
- `narrow`
- `skipprobe`
- `s2`

当前稳定结论是：

- `skipprobe` 优于 `6002 narrow`
- `s2` 相对 `skipprobe` 是“JSD 近乎持平、count 更稳”的 tradeoff
- 因此 `skipprobe` 与 `s2` 都可留在 `6002` shortlist

---

## 四、今天继续推进时的默认下一步

### `6000`

1. 等待 `teacher_v2_center_pool_msdec_v1_s2_s42_20260409_full1` 收口
2. 用同一口径对历史 `corrected B` 做正式 gate
3. 若 `s2 full-budget` 仍不过 gate，则停止把 `msdec_v1` 小变体链当默认 backlog，并切到新的 `A6000 readout family`

### `6002`

1. 等待 `teacher_v2_center_pool_skipprobe_v1_wide_short10_s42_6002_r1` 收口
2. 用同机口径把它与 `skipprobe` / `s2` 做 shortlist 排序
3. 是否值得升到 `6000`，由 `6000` 自己重新判定，不自动继承

### 双机共同规则

1. 当前 `6000` 已有活跃后台训练；在它收口前，不再并发启动新的 `6000` 训练
2. 两台机器仍按“各自维护队列、不互相等待”执行
3. 若没有新的高信息量候选 family，优先做收口、归档与 paper-facing 文档同步，而不是为了占 GPU 继续发车

---

## 五、与旧文档的关系

- 本文件是 `2026-04-09` 时点的双机排程快照
- `20260407` 命名的两份文档仍保留作为历史上下文与结构细节来源：
  - `docs/plan/multiscale_decoder_probe_20260407.md`
  - `reports/session_handoff_multiscale_and_next_tasks_20260407.md`
- `TRACKING.md` 现在应把本文件视为历史快照，不再把它当作双机实验 live 入口
