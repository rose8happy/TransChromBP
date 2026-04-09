# TransChromBP 开发工作流指南 (Development Workflow)

本文档描述了 TransChromBP 项目的开发、同步和版本管理流程。我们采用了 **"Git 版本控制 + Rsync/Scp 快速同步"** 的混合模式，以兼顾代码安全性和开发效率。

## 1. 核心流程图解

*   **本地 (Local)**: 代码编写、文档撰写、结果分析。
*   **服务器 (Remote)**: 模型训练、数据处理 (GPU 环境)。
*   **GitHub**: 代码托管、版本备份。

```mermaid
graph TD
    A[本地开发 Local] -->|Git Push| B(GitHub 仓库)
    B -->|Git Pull| C[服务器 Remote]
    A -->|Deploy (Rsync)| C
    C -->|Download Results (Rsync)| A
```

## 2. 版本管理与官方源码外化

从 2026-04-06 起，本项目把“TransChromBP 主仓”和“官方 ChromBPNet 源码”分开管理：

1. **本仓库** 是 TransChromBP 主仓 + 项目档案，不是官方 ChromBPNet 源码主库。
2. **官方 ChromBPNet 源码**、`setup.py`、`MANIFEST.in` 和 `chrombpnet.egg-info` 正在外化到 6000 上的 `/data1/zhoujiazhen/bylw_atac/chrombpnet_official`。
3. **官方代码查找、官方复现、官方对比** 默认先去 6000 外部仓，不要把当前仓库当作官方来源。
4. **GitHub** 负责长期备份与版本历史，不是实验运行目录。
5. **6000 工作区** 是训练/评估与官方源码对照的远端工作区。
6. **6002 工作区** 当前按部署/运行目录管理，默认不作为 Git 档案节点。

这意味着：

- 代码、配置、实验计划、launcher、结果摘要，应先落本地仓库，再同步到远端。
- 远端允许存在运行期临时改动，但必须回收并整理回本地。
- 模型权重、原始数据、原始训练日志这类大文件，不作为 Git 版本管理对象。
- 需要核对官方 ChromBPNet 行为时，优先在 `/data1/zhoujiazhen/bylw_atac/chrombpnet_official` 里完成代码查找和复现，再把结论带回本仓库。

### 2.1 当前三处工作区的角色

| 位置 | 当前角色 | 是否应视为主档案 |
| :--- | :--- | :---: |
| 本地 `/home/zhengwei/project/python/chromBPNet` | TransChromBP 主开发仓库、文档与分析汇总中心 | ✅ |
| 6000 `/data1/zhoujiazhen/bylw_atac/chrombpnet_official` | 官方 ChromBPNet 外部源码仓与复现基准 | ✅（官方源码） |
| 6000 `/data1/zhoujiazhen/bylw_atac/TransChromBP` | TransChromBP 运行/开发工作区 | ❌ |
| 6002 `/home/zhengwei/bylw_atac/TransChromBP` | 运行/部署目录 | ❌ |

### 2.2 哪些内容必须回收到本地版本管理

以下内容默认应进入本地仓库并参与 Git 版本管理：

- `scripts/`、`docs/`、`tests/` 等代码与文档目录
- `chrombpnet/` 如仍需参考，只视为外化期间的兼容路径，不应作为官方源码主档案
- `vendor/transchrombp/` 中的版本化 TransChromBP snapshot 与辅助脚本
- `configs/`、launcher、runtime config 模板、实验清单、计划文档
- `reports/*.tex`
- `reports/assets/` 下的小型结果摘要，例如 `csv`、`png`、`json`
- `TRACKING.md`、`DEVELOPMENT.md`、实验安排文档

以下内容默认不进入 Git：

- checkpoint、权重、`*.h5`、`*.pt`
- 原始训练数据与派生大文件：`bam/bw/bigWig/fa/fai/bed.gz`
- 远端 `outputs/` 下的大型中间产物
- 原始长训练日志全文

### 2.3 原始日志与摘要的处理规则

- `logs/` 目录下的原始训练日志仍然默认 **不进 Git**。
- 但日志中的关键信息不应只留在远端：
  - 需要整理成 `md/txt/csv/json/png` 摘要后回收到本地仓库。
- 如果一份 `json` 很小、且是结果摘要而不是大规模中间产物，则应允许版本化。

### 2.4 每次实验的最小归档动作

在启动实验前，至少保证以下文件已经存在于本地仓库：

1. 对应的 launcher
2. 对应的 `model/train/data config`
3. 对应的计划文档或实验清单

在实验结束后，至少补齐以下动作：

1. 将关键结果摘要带回本地
2. 更新 `TRACKING.md`
3. 如有必要，补 `reports/` 或 `docs/plan/` 中的结论文档
4. 再决定是否提交 Git

## 3. 常用操作命令

常用的本地-远端同步操作通常通过封装好的脚本 `./scripts/sync_project.sh` 完成；下面也保留了直接调用的工作流与测试命令。

### 场景 A: 常规开发与提交 (Git)
当你完成了一个功能模块或修复了一个 Bug，并希望永久保存记录时：

1.  **提交并推送 (Local)**:
    ```bash
    ./scripts/sync_project.sh push
    ```
    *提示：脚本会检测未提交的更改，并询问是否立即提交。*

2.  **拉取更新 (Remote)**:
    登录服务器后：
    ```bash
    ./scripts/sync_project.sh pull
    ```

### 场景 B: 快速调试 (Fast Deploy)
当你正在频繁修改代码（如调整 `train.py` 中的打印语句或逻辑），不想产生大量琐碎的 Git 提交记录时：

1.  **一键部署 (Local -> Remote)**:
    ```bash
    ./scripts/sync_project.sh deploy
    ```
    *   **原理**: 使用 `rsync` 将本地代码直接覆盖到服务器。
    *   **注意**: 不会覆盖服务器上的 `logs/`, `outputs/`, `data/` 等数据目录。
    *   **默认目标**: 该脚本当前默认写入 `/data1/zhoujiazhen/bylw_atac/chromBPNet`；6000 上的 TransChromBP 运行工作区是 `/data1/zhoujiazhen/bylw_atac/TransChromBP`，两者不要混为一谈。

### 场景 C: 查看训练结果 (Download Results)
当服务器上的训练任务完成后，你希望在本地查看日志曲线或 HTML 报告时：

1.  **回传结果 (Remote -> Local)**:
    ```bash
    ./scripts/sync_project.sh download_results
    ```
    *   **功能**: 将服务器的 `logs/` 和 `reports/` 目录同步回本地对应目录。

## 4. 目录结构与同步规则

为了防止数据丢失和冲突，我们定义了严格的同步规则：

| 目录/文件 | Git 托管? | Deploy 同步? (Local->Remote) | 结果回传? (Remote->Local) | 说明 |
| :--- | :---: | :---: | :---: | :--- |
| `chrombpnet/` | ⚠️ 过渡兼容 | ✅ | ❌ | 仅用于外化期间的兼容路径，不是官方源码主档案 |
| `scripts/` | ✅ | ✅ | ❌ | 脚本工具 |
| `docs/` | ✅ | ✅ | ✅（必要时） | 计划、清单、结果说明 |
| `logs/` | 原始日志 ❌；摘要 ✅ | ❌ (排除) | ✅ (下载) | 原始日志不入 Git，摘要需整理回本地 |
| `reports/` | ✅ (源码与摘要) | ✅ | ✅ (下载) | 保留 `tex/md/assets` 源文件；PDF 与 LaTeX 构建产物默认不入 Git |
| `vendor/transchrombp/` | ✅ | ✅ | ❌ | 版本化的 TransChromBP 本地 snapshot |
| `outputs/` | 大文件 ❌；小摘要视情况 ✅ | ❌ (排除) | ❌ | checkpoint 不进 Git；小型汇总需手动回收 |
| `data/` | ❌ | ❌ (排除) | ❌ | 训练数据 (只在服务器存在) |
| `*.h5`, `*.bw` | ❌ | ❌ (排除) | ❌ | 大文件 |

## 5. 环境配置信息

*   **本地 Git 用户**: `yangmeisuan <345687960@qq.com>`
*   **服务器 Git 用户**: `yangmeisuan <345687960@qq.com>`
*   **服务器地址**: `127.0.0.1` (Port 6000)
*   **脚本默认部署路径**: `/data1/zhoujiazhen/bylw_atac/chromBPNet`
*   **6000 运行工作区**: `/data1/zhoujiazhen/bylw_atac/TransChromBP`

### 6000 上的多环境分工

6000 上 `chrombpnet`、`transchrombp`、`genos-1.2b` 三套环境的用途、版本、激活方式与已知坑，已集中整理到：

- `docs/env/transchrombp_genos_env.md`

### 本地报告/画图独立环境

为避免把 `matplotlib`、`pandas`、`seaborn`、`plotly` 等报告依赖混入训练环境，仓库提供了独立环境脚本：

```bash
bash scripts/setup_report_env.sh
source .venv-report/bin/activate
```

默认安装的是画图与表格处理核心包，适合生成报告、CSV 汇总和可视化脚本。

如果还需要 Jupyter / ipykernel，再执行：

```bash
INSTALL_NOTEBOOK=1 bash scripts/setup_report_env.sh
```

## 6. 故障排查

*   **SSH 连接失败**: 确保你的 SSH 密钥已添加到服务器，且可以通过 `ssh zhoujiazhen@127.0.0.1 -p 6000` 直接登录。
*   **权限错误**: 确保你对服务器上的目标目录有写权限。
*   **Git 冲突**: 如果 `pull` 时出现冲突，请先手动解决冲突 (`git merge` 或 `git rebase`)，然后再继续。
*   **远端代码与本地不一致**: 以本地仓库为准，先把远端有效改动回收成本地文件，再决定 `commit` 或重新部署；不要把 6000/6002 工作目录当成最终档案源。
