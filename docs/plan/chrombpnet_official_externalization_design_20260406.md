# 官方 ChromBPNet 外置化设计（2026-04-06）

## 1. 目标

把当前根仓库从“官方 `chrombpnet` clone 上叠加项目档案与自研代码”收口为：

> `TransChromBP` 主仓 + 项目总档案仓

同时把官方 `ChromBPNet` 的源码查阅与官方复现职责，明确外置到 6000（A6000）上的独立参考仓，不再把官方源码长期留在本仓库内。

本设计解决的不是磁盘容量问题，而是以下三个长期混淆：

1. 根仓当前仍保留官方 `chrombpnet/` 包、`setup.py`、`README.md` 与 `pip install -e .` 语义，导致“这个仓库到底是不是官方仓”不清楚。
2. 当前项目文档与代码主线已经把 `vendor/transchrombp/` 定义为正式 `TransChromBP` snapshot，但根目录结构仍会误导维护者把官方代码当成本地主线。
3. 一部分 official-compare / dataset-prep 脚本还直接依赖本地 `chrombpnet/` 路径，因此本地官方 clone 还没有真正退役。

## 2. 当前状态

当前根仓库的结构同时承载了三类语义：

- 官方 `ChromBPNet` 代码：
  - `chrombpnet/`
  - `chrombpnet.egg-info/`
  - `setup.py`
  - `MANIFEST.in`
  - 根目录 `README.md`
- 自研 `TransChromBP` 正式 snapshot：
  - `vendor/transchrombp/transchrombp/`
  - `vendor/transchrombp/scripts/`
- 项目级文档、论文、实验计划与结果分析：
  - `docs/`
  - `reports/`
  - `TRACKING.md`

与此同时，6000 上已经具备承担“官方复现 + 官方源码查阅”的条件：有 `chrombpnet` 环境、A6000 资源和长期运行目录。因此，官方代码没有必要继续作为本仓库的常驻源码树存在。

## 3. 决策结论

采用“两阶段外置化”，不做一次性硬切。

### 3.1 最终目标态

根仓库只保留以下角色：

- `TransChromBP` 主线代码与版本化 snapshot
- 项目级脚本、文档、论文、计划、追踪
- 与官方复现相关的桥接脚本、patch 台账、执行说明

根仓库不再保留以下内容：

- 本地 `chrombpnet/` 官方源码树
- `chrombpnet.egg-info/`
- “本仓库可直接 `pip install -e .` 并作为官方 ChromBPNet 使用”的包装语义

官方 `ChromBPNet` 的 canonical 参考位置改为：

- 6000：`/data1/zhoujiazhen/bylw_atac/chrombpnet_official`

后续所有“看官方源码”和“跑官方复现”的需求，都默认在该 A6000 参考仓完成，而不是在本仓完成。

### 3.2 选择该方案的原因

相比“直接删本地 `chrombpnet/`”的硬切，这个方案更稳，因为当前仍有一批本地脚本和文档依赖本地官方路径；如果不先桥接这些依赖，硬切会把仓库从“结构混乱”直接推成“链路断裂”。

相比“继续共存，只在文档里声明不是主仓”，这个方案边界更清楚，因为它真正移除了最容易误导人的本地官方源码树和 packaging 语义，而不是仅靠口头说明维持约束。

## 4. 迁移范围

### 4.1 在本轮设计内的内容

1. 盘清并替换本地对 `chrombpnet/` 的显式路径依赖。
2. 固化 6000 官方参考仓的路径契约与使用方式。
3. 把当前依赖官方源码的脚本改成：
   - 显式使用 6000 官方仓；
   - 或显式失败并提示正确远端路径；
   - 不再默认吃本地 `chrombpnet/`。
4. 修改根仓文档与规范，使其不再把 `chrombpnet/` 说成本仓核心代码。
5. 在依赖切断后，删除本地 `chrombpnet/`、`chrombpnet.egg-info/` 和相关官方 packaging 入口。

### 4.2 明确不在本轮处理的内容

1. 不重跑任何官方实验。
2. 不搬迁 6000/6002 既有结果目录。
3. 不把 `vendor/transchrombp/` 再拆成独立仓。
4. 不在本轮解决 `tmp_remote_edit/TransChromBP/scripts/` 的 11 个历史 helper script 问题；那是并行 cleanup 线，不是本设计的主目标。

## 5. 需要切断的依赖类型

### 5.1 运行时脚本依赖

这类脚本目前直接引用本地 `chrombpnet/` 路径，必须先改：

- `scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh`
- `scripts/paper_aligned_repro/select_best_epoch.py`
- `scripts/start_6000_chrombpnet_dataset_prep.sh`
- `scripts/start_6002_chrombpnet_dataset_prep.sh`
- `workflows/tutorial/step3_get_background_regions.sh`

处理原则：

- 本仓脚本不再假设 `REPO_ROOT/chrombpnet/...` 存在。
- 任何官方 compare / official prep 行为，都显式走 6000 官方仓。
- 若脚本在非 6000 环境被调用，应给出清晰错误，而不是静默走错路径。

### 5.2 文档与规范依赖

以下文件当前仍把本地 `chrombpnet/` 视为本仓核心组成，需要改口径：

- `AGENTS.md`
- `DEVELOPMENT.md`
- 根目录 `README.md`
- 与 official compare 相关的当前权威计划/说明文档

处理原则：

- 历史报告不强行全文改写，只修当前权威入口。
- 历史报告里若引用本地 `chrombpnet/...` 路径，保留为“历史上下文”即可。
- 当前 canonical 文档必须明确：
  - 本仓是 `TransChromBP` 主仓；
  - 官方 ChromBPNet 在 6000 外置存在。

### 5.3 Packaging 语义依赖

根仓当前仍有官方 packaging 入口：

- `setup.py`
- `MANIFEST.in`
- `chrombpnet.egg-info/`

这套入口与“本仓不再承载官方源码”直接冲突，因此迁移完成后必须一起退役或改写。

本设计选择：

- 退役根仓的官方 ChromBPNet packaging 语义；
- 不在本轮把 `vendor/transchrombp/` 补成新的 installable package；
- 本地导入 `transchrombp` 继续沿用当前约定：通过 `PYTHONPATH=vendor/transchrombp` 或远端 runtime 工作区使用。

## 6. 目标架构

### 6.1 本仓角色

本仓保留：

- `vendor/transchrombp/`：版本化 `TransChromBP` snapshot
- `scripts/`：项目脚本、桥接脚本、官方复现 wrapper
- `docs/` / `reports/` / `TRACKING.md`：项目文档与证据链
- `tests/`：针对自研代码或桥接脚本的检查

本仓移除：

- `chrombpnet/`
- `chrombpnet.egg-info/`
- 根仓作为官方 ChromBPNet 安装包的入口

### 6.2 6000 官方参考仓角色

6000 上的 `/data1/zhoujiazhen/bylw_atac/chrombpnet_official` 成为唯一官方参考仓，承担：

- 官方源码查阅
- official compare 所需的 patched `predict.py` / `metrics.py` / helper code
- 官方 tutorial / controlled compare / dataset prep 的真实执行

该仓的约束是：

- 本仓不再镜像其完整源码；
- 若我们对官方代码有本地修补需求，只保留“patch 台账 / wrapper / 说明”，不把整份官方源码继续带回主仓。

### 6.3 桥接方式

本仓涉及官方代码的脚本，只允许两种方式之一：

1. 显式 SSH 到 6000，在官方仓路径下执行；
2. 本地先组装命令，再把命令或脚本发到 6000 的官方仓执行。

不再允许：

- 直接引用 `REPO_ROOT/chrombpnet/...`
- 默认把当前根仓当成官方源码根目录

## 7. 迁移步骤

### 阶段 A：依赖盘点与路径合同冻结

先做一次面向实现的硬盘点，把所有仍引用本地 `chrombpnet/` 的脚本、文档、流程列全。

同时固定唯一远端路径合同：

- `CHROMBPNET_OFFICIAL_ROOT=/data1/zhoujiazhen/bylw_atac/chrombpnet_official`

任何桥接脚本都必须围绕这个路径合同工作，不允许各自再发明第二套路径。

### 阶段 B：改桥接脚本

把仍依赖本地 `chrombpnet/` 的 active 脚本改成远端官方仓调用。

重点目标是：

- `scripts/paper_aligned_repro/*`
- `scripts/start_6000_chrombpnet_dataset_prep.sh`
- `scripts/start_6002_chrombpnet_dataset_prep.sh`
- `workflows/tutorial/step3_get_background_regions.sh`

处理后，这些脚本即使仍留在本仓，也不会再依赖本地官方源码树。

### 阶段 C：固化官方 patch 台账

对于我们实际依赖过的官方改动点，要留下“差异说明”，而不是继续保留整份源码树。

最低必须固化的点：

- `chrombpnet/training/predict.py`
- `chrombpnet/training/metrics.py`
- `chrombpnet/helpers/make_gc_matched_negatives/*`

这部分的产物形式应是：

- 当前权威设计/执行文档中的差异说明；
- 必要时补一份 patch 台账文档；
- 而不是继续把官方完整源码留在本仓。

### 阶段 D：修改根仓口径

把当前权威文档与规范同步改成新角色：

- 本仓不是官方 ChromBPNet 仓
- 本仓也不再提供官方 ChromBPNet 本地安装能力
- 官方源码查阅与 official reproduction 默认去 6000 外置仓

### 阶段 E：删除本地官方源码与包装

只有在阶段 A-D 完成后，才允许删除：

- `chrombpnet/`
- `chrombpnet.egg-info/`
- 与官方 package 强绑定的根仓 packaging 入口

这一步做完后，根仓结构才真正完成角色切换。

## 8. 风险与防护

### 风险 1：脚本静默失效

如果只删本地 `chrombpnet/`，很多脚本会在运行时才报“找不到文件”，甚至可能在错误目录下执行。

防护：

- 先改脚本，再删目录；
- 脚本改完后加显式路径检查和清晰报错。

### 风险 2：文档口径落后于仓库现实

如果删了目录但不改 `AGENTS.md` / `DEVELOPMENT.md` / `README.md`，维护者仍会被旧口径误导。

防护：

- 删除本地官方源码的同一轮，必须同步修改权威文档。

### 风险 3：官方 patch 证据丢失

如果本仓不再保留官方源码，而我们又没有把关键 patch 点沉淀成说明，后续会不知道 official compare 为什么和原始 upstream 不同。

防护：

- 在删本地官方源码前，先写清 patched points 和用途。

### 风险 4：根仓 packaging 变成半残状态

如果只删 `chrombpnet/` 但保留 `setup.py` / `MANIFEST.in` 原语义，根仓会变成“看起来可安装、实际不可用”的半残仓。

防护：

- 把 packaging 退役动作明确纳入同一条迁移链，而不是留作顺手处理。

## 9. 验证标准

迁移完成后，至少满足以下条件才算成功：

1. 本仓 `rg` 不再出现 active 脚本对 `REPO_ROOT/chrombpnet/...` 的依赖。
2. 当前权威文档不再把 `chrombpnet/` 描述为本仓核心代码。
3. official compare / official prep 相关脚本能明确解析到 6000 官方仓路径，或在缺路径时给出清晰错误。
4. 本仓根目录不再包含 `chrombpnet/` 与 `chrombpnet.egg-info/`。
5. 根仓不再保留“官方 ChromBPNet 可直接 `pip install -e .`”的文档与包装语义。

## 10. 一句话执行原则

这次迁移的核心不是“把一个目录删掉”，而是：

> 先把本仓对官方源码的活依赖切到 6000 外置仓，再删除本地官方源码与官方 packaging，让根仓真正只代表 `TransChromBP` 主线。
