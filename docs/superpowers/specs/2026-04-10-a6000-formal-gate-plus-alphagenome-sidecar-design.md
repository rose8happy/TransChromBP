# 2026-04-10 A6000 正式训练 Gate + AlphaGenome Sidecar 设计

## 1. 背景与问题

`AlphaGenome matched raw-track slice v1` 已经在 `6000` 上完成 closeout，并以 **technical/alignment gate** `pass` 收口。这个结果只说明：

- `AlphaGenome` 这条窄 external coordinate 路径已经可以稳定跑通；
- 可以形成可复用的 `summary / metadata / profiles / merged` 产物；
- 若继续推进，唯一允许的下一步是扩到 `12-20 loci` 的稍大 matched panel。

它**不**说明：

- `AlphaGenome` 在模型质量上优于当前方案；
- `AlphaGenome` 应该成为 `6000` 的主 GPU 实验；
- A6000 下一步应该继续空着等这条 sidecar 扩面。

当前新的现实约束是：

1. `AlphaGenome` 这条线继续扩面，大概率仍然走 SDK/API + CPU/网络路径，不会真正吃掉 A6000 训练算力。
2. `6000 / A6000 x2` 当前空闲，应该承担一条真正使用 GPU 的正式训练实验。
3. 外部共识仍然稳定指向：
   - `corrected B` 作为 bias-safe 基座
   - 下一代最值得押注的方向是 `multiscale / local-skip / coarse-to-fine dense decoder`
   - 若继续 foundation，只允许 genuinely new `teacher / init` 范式，不回到旧 adapter family。

因此，本设计要解决的核心问题是：

> 如何让 `AlphaGenome 12-20 loci` 扩面作为 sidecar 独立推进，同时给 `6000` 设计一条真正吃 GPU 的正式训练 gate，而且两条线互不等待、互不污染解释口径。

---

## 2. 目标与非目标

### 2.1 目标

本设计的目标有四个：

1. 把下一步明确拆成两条并行轨道，而不是继续共享一个“唯一下一步”。
2. 让 `AlphaGenome` 扩面继续作为 sidecar，保持 external coordinate 角色。
3. 让 A6000 真正承担一条正式、吃 GPU 的训练 gate。
4. 把这条正式训练 gate 的边界、budget 和 stop-rule 一次写死，避免跑完后再临场解释。

### 2.2 非目标

本设计明确不做以下事情：

1. 不把 `AlphaGenome` 扩面升级成大 benchmark。
2. 不把 `AlphaGenome` 结果当作 A6000 训练 gate 的前置。
3. 不重开旧的 `teacher-distill tutorial` 近邻变体。
4. 不复活旧 `msdec / skipprobe / U-Net-lite` 的命名链作为正式主线。
5. 不在本设计中直接实现新 decoder；这里只定义下一条正式实验该是什么。

---

## 3. 方案比较

### 3.1 方案 A：只继续 `AlphaGenome 12-20 loci` 扩面

优点：

- 几乎不需要新的训练实现；
- 可以快速继续积累外部坐标。

缺点：

- 仍然主要消耗 CPU / API / I/O，而不是 A6000 训练算力；
- 无法回答“当前 ceiling 是否卡在 profile readout / multiscale 结构”这个主问题；
- 会让 A6000 继续空着。

### 3.2 方案 B：`AlphaGenome` sidecar + 新的 `teacher / init` GPU family

优点：

- 也能让 A6000 真正承担正式训练；
- 若 hypothesis 成立，潜在上限较高。

缺点：

- 当前工程与概念风险都更大；
- 不如 multiscale decoder 那样贴近当前外部共识的第一优先级；
- 当前仓库中也没有一个 ready-to-run 的新 family 实现锚点。

### 3.3 方案 C：`AlphaGenome` sidecar + 新的 `multiscale/local-skip decoder` 正式 gate

优点：

- 最符合外部共识的第一优先级；
- 直接回答当前主问题：`corrected B` 的 ceiling 是否被 profile readout / multiscale 结构限制；
- 可以把 A6000 用在真正吃 GPU 的正式训练上；
- `AlphaGenome` sidecar 与其在角色上天然解耦。

缺点：

- 需要新一轮 clean family 设计；
- 不能简单拿旧 `msdec / skipprobe` 名字复用，否则会把旧历史包袱带回来。

### 3.4 结论

本次采用 **方案 C：`AlphaGenome sidecar` + `A6000 multiscale/local-skip decoder formal gate`**。

---

## 4. 双轨并行结构

### 4.1 轨道 A：`AlphaGenome matched raw-track slice v2`

这条轨道被定义为：

- `AlphaGenome matched raw-track slice v1` 的唯一允许续步
- 规模从 `4 loci` 扩到约 `12-20 loci`
- 仍是窄 `external coordinate`
- 仍不升级成大 benchmark
- 仍不用于训练学生模型
- 不是 `6000` 的主 GPU 实验

### 4.2 轨道 B：`corrected B + multiscale/local-skip decoder formal gate`

这条轨道被定义为：

- `6000 / A6000 x2` 的正式 GPU 训练实验
- 作用是回答“当前 ceiling 是否真由 profile readout / multiscale 结构限制”
- 它是一个新的、干净的 `v2 family`
- 它不是：
  - `U-Net-lite v1` 的镜像 rerun
  - 旧 `msdec / skipprobe` 命名链的自然续跑
  - 旧 adapter family 的 disguised revival

### 4.3 解耦硬规则

两条轨道必须满足：

1. 互不等待。
2. `AlphaGenome` sidecar 不能成为 A6000 正式训练的启动前置。
3. A6000 正式训练的结果也不能改写 `AlphaGenome` sidecar 的技术/对齐结论。
4. 两条线分别 closeout，分别判读，分别写结论。

---

## 5. A6000 正式训练 Gate 的最小实现边界

### 5.1 基座固定

正式训练 gate 的基座固定为：

- [transchrombp_teacher_v2_center_pool.yaml](/home/zhengwei/project/python/TransChromBP/vendor/transchrombp/transchrombp/configs/model/transchrombp_teacher_v2_center_pool.yaml)

这就是当前 `corrected B` 基座，必须保留：

1. `sequence_encoder`
2. `conv_stem`
3. `local_tower`
4. `bias_branch`
5. `fusion`
6. `count_head`
7. `full/debiased` 语义

### 5.2 唯一允许变化的部分

这条实验只允许改 `profile` 支路。

具体边界：

1. 引入 `2-3` 级 coarse-to-fine decoder
2. 至少一条 local-skip 路径
3. 输出仍回到 `1000 bp`
4. `bias fusion` 仍发生在 debiased profile 之后
5. 不允许 skip / decoder 把 bias 语义偷渡回主路径

### 5.3 明确禁止混入的变量

在这条 gate 里明确禁止再同时混入：

1. foundation token / summary / cross-attention
2. 新数据线
3. 新 count loss
4. 新训练范式
5. 多个 decoder family 并行混跑

换句话说，这条 A6000 实验只回答一个问题：

> 在 `corrected B` 其余语义完全不动的前提下，单独升级 `profile readout` 为 `multiscale/local-skip decoder`，能不能给出 clean gain？

---

## 6. 30 Epoch Formal Gate 规则

### 6.1 budget 设计

这条新 family 的第一次正式判定，直接定义为：

- `30 epoch`
- `single-seed`
- `seed=42`
- `6000 / A6000 x2`
- tutorial canonical 数据线

它不是：

- `short10 gate`
- full campaign
- multi-seed expansion

### 6.2 为什么直接用 30 epoch

用户已明确建议直接跑 `30 epoch`，而不是再先做 `short10`。

在本设计里，这个选择的含义是：

- 第一次正式 run 就要给出足够强的 verdict；
- 但仍然通过 `single-seed formal gate` 控制预算，不直接跳到更长、更重的 full-budget。

### 6.3 过 gate 后唯一允许的下一步

若 `30 epoch formal gate` 通过，唯一允许的下一步才是：

- 第二个 seed，或
- 更完整 budget 的升级

### 6.4 不过 gate 后的动作

若 `30 epoch formal gate` 不通过，则：

- 这条 `v2 family` 直接停表
- 不补一串 `v2a / v2b / v2c` 小变体
- 不把“不太差”解释成继续扩线理由

---

## 7. 判读阈值

### 7.1 主比较对象

这条 A6000 formal gate 只和当前 `corrected B` 比，不和以下对象比：

1. 旧 `U-Net-lite`
2. 旧 `msdec`
3. `6002` 的单卡结果
4. `AlphaGenome` sidecar

### 7.2 `pass` 条件

推荐的 `pass` 规则分两档。

**强通过条件：**

1. `peak profile JSD` 改善 `>= 0.003`
2. `peak count_r` 不下降超过 `0.005`
3. `full/debiased gap` 不明显放大

**弱通过条件：**

1. `peak profile JSD` 基本持平，定义为不差于 `0.001`
2. `peak count_r` 提升 `>= 0.01`
3. `full/debiased gap` 仍稳定

只有满足其中一档，才可视为 `pass`。

### 7.3 直接 `fail`

只要满足任一条，就记为 `fail` 并停表：

1. `peak profile JSD` 恶化超过 `0.003`
2. `peak count_r` 恶化超过 `0.01`
3. `full/debiased gap` 明显恶化
4. 训练出现非偶然不稳定 / collapse

---

## 8. `AlphaGenome` sidecar 的并行规则

`AlphaGenome matched raw-track slice v2` 可以和 A6000 正式训练并行推进，但必须满足：

1. 它只是一条 sidecar
2. 它不占 A6000 主实验位
3. 它不能在结果解释上盖过 A6000 正式训练
4. 它也不能被 A6000 正式训练的结果自动取消

如果继续 `AlphaGenome` 线，唯一允许的下一步仍然是：

- `12-20 loci`
- 小幅扩面
- 窄 external coordinate

不允许：

- 直接扩成大 benchmark
- 把这次 `pass` 漂移成质量结论

---

## 9. 行动化结论

当前应把下一步固定为两条并行轨道：

1. `AlphaGenome matched raw-track slice v2`
   - 作为 sidecar
   - 扩到约 `12-20 loci`
   - 保持窄 external coordinate 定位

2. `corrected B + multiscale/local-skip decoder v2`
   - 作为 `6000 / A6000 x2` 的正式 GPU 训练 gate
   - 直接跑 `30 epoch single-seed formal gate`
   - 只比较 `corrected B`
   - 过 gate 才允许二次扩展

一句话压缩：

> `AlphaGenome` 继续做 sidecar 坐标扩面，A6000 同时承担一条真正吃 GPU 的 `corrected B + multiscale/local-skip decoder` 正式训练 gate；两条线互不等待、互不借口径、各自 closeout。
