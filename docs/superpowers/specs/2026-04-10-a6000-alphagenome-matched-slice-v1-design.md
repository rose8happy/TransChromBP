# 2026-04-10 A6000 AlphaGenome Matched Slice v1 设计

## 1. 背景与问题

当前双机调度已经改成“吞吐优先 + 完全解耦”：

- `6000 / A6000 x2` 不再等待 `6002` 的 cheap-screen 结果；
- `6002 / RTX 3080` 继续独立完成 `U-Net-lite v1 r4` 的 rigor closeout；
- `6000` 需要自己的独立 backlog，避免 A6000 空转。

当前已知事实是：

- tutorial `teacher-distill` 线已经正式停表；
- `6002` 仍在跑 `teacher_v2_center_pool_unet_lite_v1_short10_s42_6002_20260409_r4`；
- 现有仓库已经具备 AlphaGenome pilot 脚手架：
  - [run_alphagenome_pilot.py](/home/zhengwei/project/python/chromBPNet/scripts/alphagenome_pilot/run_alphagenome_pilot.py)
  - [regions_k562_tutorial_selected_loci.csv](/home/zhengwei/project/python/chromBPNet/scripts/alphagenome_pilot/regions_k562_tutorial_selected_loci.csv)
- `6000` 上的 AlphaGenome SDK、认证 smoke、metadata 查询都已经通过，当前 blocker 不再是环境安装，而是把第一轮 matched slice 定义清楚。

因此，本设计要解决的问题不是“要不要继续写论文”，也不是“要不要立刻上新的训练 family”，而是：

> 如何把 `6000` 的下一条任务收成一个独立、可快速收口、不会重新膨胀成 benchmark 的 `AlphaGenome matched raw-track slice v1`。

---

## 2. 目标与非目标

### 2.1 目标

本设计的目标有四个：

1. 给 `6000` 安排一条不依赖 `6002` 的独立任务。
2. 用最小闭环回答：AlphaGenome 的 matched raw-track slice 能否稳定跑通，并形成有解释力的外部坐标。
3. 复用现有 pilot 脚手架，不为 `v1` 额外引入大工程。
4. 把 pass/fail gate、输出落点、closeout 方式一次写死，避免“先跑点输出以后再说”。

### 2.2 非目标

本设计明确不做以下事情：

1. 不把 `v1` 扩成大样本 benchmark。
2. 不做 bias-only stress test。
3. 不把 AlphaGenome 输出拿来训练任何学生模型。
4. 不把这条任务包装成新的双卡长训练 family。
5. 不在本设计中决定更大规模的 AlphaGenome 第二阶段实现。

---

## 3. 方案比较

### 3.1 方案 A：`corrected B + local-skip dense decoder v2`

优点：

- 更贴近“直接提升模型增益”的主问题；
- 与外部报告里“multiscale / dense decoder”主假设更一致；
- 若成功，后续更容易进入正式训练晋级。

缺点：

- 需要先补新实现，不是立即可发车；
- 工程风险和设计不确定性都更高；
- 会把当前这轮任务重新拉回“先写新模型再看”。

### 3.2 方案 B：新的 `teacher / init` 高价值 family

优点：

- 潜在上限高；
- 若形成 genuinely new hypothesis，可能带来新的 foundation 方向。

缺点：

- 当前仓里没有 ready-to-run 的具体实现；
- 概念风险和工程风险都更大；
- 不适合作为 A6000 当下立刻接上的下一条独立任务。

### 3.3 方案 C：`AlphaGenome matched raw-track slice v1`

优点：

- 现有前置最完整，最快可闭环；
- 与 `6002` 当前的 `U-Net-lite r4` 完全解耦；
- 能先拿到一个窄而准的 external coordinate；
- 不会回到已判负的 adapter family。

缺点：

- 它回答的是“外部坐标是否能干净建立”，不是“模型训练是否直接提升”；
- 它主要占用的是 `6000` 的实验窗口，而不是长期占满双 A6000 的训练吞吐；
- 若后续要扩成更大规模比较，还需要单独二阶段设计。

### 3.4 结论

本次采用 **方案 C：`AlphaGenome matched raw-track slice v1`**。

这里的含义不是说 AlphaGenome 已经成为项目主线，而是：

- 当前需要给 `6000` 一个独立、可立即执行、信息密度足够高的任务；
- `v1` 先回答“这条窄 slice 能否稳定形成外部坐标”；
- 更重的训练 family 仍可在后续单独设计，不与本轮 pilot 混写。

---

## 4. `v1` 的任务边界

### 4.1 任务性质

`AlphaGenome matched raw-track slice v1` 被定义为：

- `6000` 上独立推进的 `external-coordinate pilot`
- 不等待 `6002 r4`
- 不依赖 `U-Net-lite v1` 的 stop/go 结论
- 基于 API/推理脚手架的小规模对照

它不是：

- 新的训练 family
- 双卡长训练
- 大 benchmark
- bias stress test

### 4.2 当前限定范围

`v1` 只做以下限定任务：

1. 物种/细胞类型：`K562`
2. 输出类型：`ATAC`
3. 任务形式：`matched raw-track compare`
4. 数据规模：固定小面板，不扩 held-out 大样本

### 4.3 允许复用的现有资产

`v1` 默认复用：

1. [run_alphagenome_pilot.py](/home/zhengwei/project/python/chromBPNet/scripts/alphagenome_pilot/run_alphagenome_pilot.py)
2. [merge_locus_totals.py](/home/zhengwei/project/python/chromBPNet/scripts/alphagenome_pilot/merge_locus_totals.py)
3. [regions_k562_tutorial_selected_loci.csv](/home/zhengwei/project/python/chromBPNet/scripts/alphagenome_pilot/regions_k562_tutorial_selected_loci.csv)
4. `6000` 上现有 `alphagenome` 环境与 API key 读取方式

---

## 5. `matched loci` 面板

`v1` 固定使用现成的 4 个 tutorial 代表位点：

1. `peak_high`
2. `peak_mid`
3. `nonpeak_high`
4. `nonpeak_mid`

设计理由：

- 这 4 个位点已经在现有脚手架里固化；
- 它们覆盖 peak / nonpeak 与高 / 中信号两类代表情形；
- 对 `v1` 来说，关键是闭环完整，不是样本量大；
- 若连这 4 个位点都无法稳定形成对照，直接说明这条 pilot 当前不该扩。

`v1` 明确不做：

- 新建 held-out 大 panel
- 临时扩到几十或几百个位点
- 混入新的细胞类型或新的输出类型

---

## 6. Gate 与 stop-rule

### 6.1 `v1` 的 gate 类型

`v1` 的 gate 是 **技术与对齐 gate**，不是性能胜负 gate。

它回答的是：

- AlphaGenome 这条 matched slice 能否稳定跑通；
- 能否形成统一、可读、可复用的 raw-track 对照产物。

它不回答：

- AlphaGenome 是否在大样本上赢过当前模型；
- 是否应该立刻把 AlphaGenome 变成主 benchmark。

### 6.2 `pass` 条件

`v1` 只有在以下条件全部满足时才算 `pass`：

1. 4 个位点全部成功返回 AlphaGenome 输出。
2. 每个位点都成功生成：
   - `summary.csv`
   - `region_metadata.jsonl`
   - `profiles/*.npz`
3. ontology/filter 后，每个位点仍保留至少 1 条可比较的 `ATAC` track。
4. 可以与本地已有 locus totals/predictions 成功合并，产出一张统一的 merged 表。
5. 输出不是明显坏结果，例如：
   - 全零
   - 全部缺 track
   - 大部分位点 API / track matching 失败

### 6.3 `fail` 条件

只要出现以下任一情况，就记为 `fail` 并停表：

1. 4 个位点里有明显比例无法稳定返回结果。
2. ontology/filter 后经常没有可比较 track。
3. 输出格式不稳定，无法形成统一 merged 表。
4. 需要额外补大量新工程，才能勉强继续。

### 6.4 `pass` 后唯一允许的下一步

若 `v1 pass`，唯一允许的下一步是：

- 扩成一个稍大的 matched panel，规模约 `12-20 loci`
- 仍然定位为窄 external coordinate
- 仍然不升级成大 benchmark

### 6.5 `fail` 后的动作

若 `v1 fail`，则：

- 这条 pilot 直接停表
- 不伪装成“先修一版再无限补跑”
- `6000` 转去下一个 genuinely new high-value task

---

## 7. 输出与 closeout

### 7.1 远端原始输出

原始输出仍放在 `6000` 远端的 `outputs/alphagenome_pilot/...` 下，至少包括：

1. `summary.csv`
2. `region_metadata.jsonl`
3. `profiles/*.npz`
4. `run_meta.json`
5. `merged_locus_totals.csv`

### 7.2 本仓库必须保留的产物

本仓库不保存大体积中间物，但必须保留下列可复用结论：

1. 一份 `reports/` closeout
2. `TRACKING.md` 的 live 状态更新
3. `charter` 中对 `6000` 当前 backlog 的状态更新

### 7.3 closeout 必须回答的问题

closeout 至少要明确写出：

1. 4 个位点是否全部成功
2. ontology/filter 后保留了哪些 track
3. merged 表是否成功形成
4. `v1` 最终判定为 `pass` 还是 `fail`
5. 若 `pass`，下一步只能扩到 `12-20 loci`
6. 若 `fail`，停表原因是什么

---

## 8. 执行顺序与耗时预期

### 8.1 执行顺序

`v1` 按以下顺序执行：

1. 在 `6000` 核验 `alphagenome` 环境、API key、脚手架可运行
2. 跑 4 loci 的 `v1` pilot
3. 立即生成 merged 表
4. 同轮写 closeout 与文档同步

### 8.2 耗时预期

因为这是 API/推理 pilot，不是双卡训练，所以预期应按短闭环管理：

1. 启动与环境核验：`10-20 min`
2. 4 loci 实际调用与落盘：`10-40 min`
3. 合并表与 closeout：`20-40 min`

整体预期：

- 若前置正常，`1-2 h` 内应收口；
- 若 API / track matching 出问题，也应在同一时间窗内尽快判成 `fail`；
- 不允许把这条 pilot 拖成长任务。

---

## 9. 文档同步与失败处理

### 9.1 文档同步门

`v1` 一旦发车或收口，必须同轮更新：

1. [TRACKING.md](/home/zhengwei/project/python/chromBPNet/TRACKING.md)
2. [2026-04-09_dual_machine_experiment_charter.md](/home/zhengwei/project/python/chromBPNet/docs/plan/2026-04-09_dual_machine_experiment_charter.md)
3. 一份新的 `reports/` closeout

### 9.2 启动前失败

若环境 / API key / 基础脚手架在启动前就失败，这条任务应记为 `launch-blocked`，而不是伪装成实验结果。

### 9.3 运行后失败

若 pilot 已启动，但发现 track/filter 口径无法形成稳定可比输出，则按 `fail` 停表。

### 9.4 唯一允许的工程性修正

若只是小型路径或输出整理问题，且不改变 hypothesis，本轮允许一次工程性修正。

明确禁止：

- 因为结果不好而无限补跑
- 借工程修正名义偷偷改任务定义

---

## 10. 行动化结论

当前 `6000` 的下一条任务应收成一句话：

> 立刻在 `6000` 上发起一个 `1-2 h` 内可收口的 `AlphaGenome matched raw-track slice v1` pilot，用 4 个固定 tutorial loci 验证 AlphaGenome raw-track 对照能否干净跑通并形成 merged closeout；若通过，只扩到 `12-20 loci`，若失败，直接停表并切回下一个 genuinely new high-value task。
