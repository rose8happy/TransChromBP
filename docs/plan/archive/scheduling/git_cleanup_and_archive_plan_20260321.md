# Git 工作树清理与实验归档方案（2026-03-21）

> 目标：把当前本地仓库收紧成真正的“唯一主档案”，同时不破坏 6000/6002 正在运行的实验。

---

## 一、当前判断

截至 2026-03-21，本地仓库的脏工作树主要由四类内容构成：

1. **应正式版本化的代码与文档**
2. **本地/远端运行所需，但尚未整理提交的实验资产**
3. **明显属于本地噪音或生成物的文件**
4. **命名或位置还不理想，需迁移后再提交的中间材料**

清理目标不是“把 `git status` 立刻清零”，而是先把这些内容分流清楚，再按批次收口。

---

## 二、当前脏工作树的结构化分类

### 2.1 可直接纳入 Git 的内容

这部分应视为正式项目资产，建议进入后续提交批次：

- `docs/`
  - `docs/learning/`：根目录旧学习资料已实际迁移到这里
  - `docs/plan/`：实验计划、调度、Genos 方案、6002 独立实验线
  - `docs/research/`：研究笔记与小型可复用附件
  - `docs/env/`：环境说明
- `scripts/`
  - `scripts/benchmark/`
  - `scripts/genos_feasibility.py`
  - `scripts/run_remote_chrombpnet_dataset_prep.sh`
  - `scripts/start_6000_chrombpnet_dataset_prep.sh`
  - `scripts/start_6002_*`
  - `scripts/prepare_6002_transchrombp_single_gpu.sh`
  - `scripts/setup_report_env.sh`
- `requirements-report.txt`
- `requirements-report-notebook.txt`
- `TRACKING.md`
- `TRACKING_archive.md`
- `DEVELOPMENT.md`
- `reports/*.tex`
- `reports/assets/` 下的小型摘要文件
  - `csv`
  - `png`
  - `json`

### 2.2 应保留，但更适合迁移/整理后再提交的内容

这些文件本身有信息价值，但当前位置或命名不理想：

- `TransChromBP_Teacher_V2_Code.txt`
- `TransChromBP_Teacher_V2_Evaluation.txt`
- `startup_copy.md`

建议处理方式：

- 迁移到 `reports/` 或 `docs/research/assets/`
- 重命名成带日期和主题的文件名
- 保留内容，但不要长期悬挂在仓库根目录

### 2.3 当前最需要明确策略的内容

#### `tmp_remote_edit/`

当前事实：

- 体量很小，约 `772K`
- 文件数约 `68`
- 里面包含当前真实在用的 `TransChromBP` 代码、配置和脚本草案

结论：

- 这部分**不能丢**，因为它已经承载了当前主线实验代码改造
- 但目录名 `tmp_remote_edit/` 明显是临时态命名，不适合作为长期结构

建议：

- **短期**：允许其继续留在仓库并进入版本管理，优先保证不丢
- **中期**：迁移到更正式的命名空间，例如 `vendor/transchrombp/`、`worktrees/transchrombp/` 或独立仓库

#### 根目录旧学习资料删除

当前 `git` 里被标记为删除的 8 份中文学习文档：

- `BPNet与ChromBPNet模型对比.md`
- `BPNet模型小白教程.md`
- `ChromBPNet学习资料总索引.md`
- `ChromBPNet完整项目深度解析.md`
- `ChromBPNet数据处理完全指南.md`
- `学习材料索引.md`
- `学习资料说明.md`
- `小白学习指南.md`

这些文件现在已经在 `docs/learning/` 下存在对应新位置。

结论：

- 这不是“意外删除”
- 应当按“根目录迁移到 `docs/learning/`”处理

### 2.4 应忽略或不进入 Git 的内容

这部分不应继续占用主仓库注意力：

- `.claude/`
- `.agents/`
- `server_public_key.txt`
- `*.deb`
- LaTeX 构建垃圾：
  - `*.aux`
  - `*.log`
  - `*.out`
  - `*.toc`
  - `*.fls`
  - `*.fdb_latexmk`
  - `*.xdv`
  - `*.synctex.gz`
- 原始训练日志目录：`logs/`
- checkpoint / runtime 大文件：`outputs/`

### 2.5 关于 PDF 的建议

当前有两类 PDF：

1. 外部参考论文 PDF
2. 本地由 LaTeX 生成的报告 PDF

建议默认策略：

- **报告源码** 版本化：`reports/*.tex`
- **报告 PDF** 默认不作为必须版本化对象
- 若某份 PDF 需要对外发送或固定归档，再单独决定是否提交

这能避免“源码已在、PDF 又重复占状态”的噪音。

---

## 三、建议执行批次

### 批次 A：低风险收口

目标：先把不涉及语义判断的部分收干净。

建议动作：

1. 保留并提交：
  - `docs/`
  - `scripts/`
  - `requirements-report*.txt`
  - `TRACKING*.md`
  - `DEVELOPMENT.md`
2. 把根目录旧学习资料的删除，与 `docs/learning/` 的新增一起作为一次“迁移”提交
3. 提交 `reports/*.tex` 与 `reports/assets/` 下的小型摘要
4. 保持 `outputs/`、`logs/`、大数据文件继续不进 Git

### 批次 B：临时工作树保全

目标：确保当前实验代码草案不丢。

建议动作：

1. 保留 `tmp_remote_edit/` 中真正有用的代码与配置
2. 不纳入下列明显 scratch 文件：
  - `tmp_remote_edit/tracking_copy.md`
  - `tmp_remote_edit/training_analysis_copy.md`
  - `*.checked`
3. 之后再决定是否整体改名迁移

### 批次 C：根目录整理

目标：减少根目录噪音。

建议动作：

1. 把有价值的说明文本迁入 `docs/`
2. 把明显本地材料加入忽略规则
3. 根目录只保留：
  - 仓库元数据
  - 核心说明文件
  - 少量稳定入口文件

---

## 四、我建议的具体处理顺序

### 第一步

先完成“低风险版本化批次”：

- `docs/`
- `scripts/`
- `reports/*.tex`
- `reports/assets/*` 小型摘要
- `requirements-report*.txt`
- `TRACKING*.md`
- `DEVELOPMENT.md`

### 第二步

把根目录 8 份旧学习资料按“迁移到 `docs/learning/`”收口。

### 第三步

单独处理 `tmp_remote_edit/`：

- 保留代码
- 剔除 scratch copy
- 再决定是否改目录名

### 第四步

最后再决定以下内容是否正式归档：

- `TransChromBP_Teacher_V2_Code.txt`
- `TransChromBP_Teacher_V2_Evaluation.txt`
- `startup_copy.md`

---

## 五、当前不建议立刻做的事

- 不建议直接在 6000 上做大规模 `git clean`
- 不建议把 6002 运行目录临时变成 Git 主档案
- 不建议现在就重构 `tmp_remote_edit/` 为新仓库或 submodule
- 不建议把原始训练日志全文纳入 Git

---

## 六、收口标准

当以下条件满足时，才算“本地 Git 主档案”这件事真正落地：

1. 当前活跃实验所需的代码、配置、launcher、计划文档都已在本地 Git 中
2. 远端运行目录不再承担唯一信息源角色
3. 小型实验摘要可以从本地仓库直接追溯
4. `git status` 中剩余未收口内容只剩少量明确知道如何处理的项

---

## 七、下一步建议

下一轮执行建议直接按下面顺序来：

1. 做一次“批次 A”收口
2. 再单独审一次 `tmp_remote_edit/`
3. 最后清理根目录剩余杂项

这样风险最低，也最符合当前实验仍在运行的现实约束。
