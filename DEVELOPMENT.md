# TransChromBP 开发工作流指南 (Development Workflow)

本文档描述了 TransChromBP 项目的开发、同步和版本管理流程。我们采用了 **"Git 版本控制 + Rsync 快速同步"** 的混合模式，以兼顾代码安全性和开发效率。

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

## 2. 常用操作命令

所有操作均通过封装好的脚本 `./scripts/sync_project.sh` 完成。

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

### 场景 C: 查看训练结果 (Download Results)
当服务器上的训练任务完成后，你希望在本地查看日志曲线或 HTML 报告时：

1.  **回传结果 (Remote -> Local)**:
    ```bash
    ./scripts/sync_project.sh download_results
    ```
    *   **功能**: 将服务器的 `logs/` 和 `reports/` 目录同步回本地对应目录。

## 3. 目录结构与同步规则

为了防止数据丢失和冲突，我们定义了严格的同步规则：

| 目录/文件 | Git 托管? | Deploy 同步? (Local->Remote) | 结果回传? (Remote->Local) | 说明 |
| :--- | :---: | :---: | :---: | :--- |
| `chrombpnet/` | ✅ | ✅ | ❌ | 核心代码库 |
| `scripts/` | ✅ | ✅ | ❌ | 脚本工具 |
| `logs/` | ❌ | ❌ (排除) | ✅ (下载) | 训练日志 |
| `reports/` | ✅ (部分) | ✅ | ✅ (下载) | 评估报告 |
| `outputs/` | ❌ | ❌ (排除) | ❌ | 模型检查点 (过大，手动管理) |
| `data/` | ❌ | ❌ (排除) | ❌ | 训练数据 (只在服务器存在) |
| `*.h5`, `*.bw` | ❌ | ❌ (排除) | ❌ | 大文件 |

## 4. 环境配置信息

*   **本地 Git 用户**: `yangmeisuan <345687960@qq.com>`
*   **服务器 Git 用户**: `yangmeisuan <345687960@qq.com>`
*   **服务器地址**: `127.0.0.1` (Port 6000)
*   **部署路径**: `/data1/zhoujiazhen/bylw_atac/chromBPNet`

## 5. 故障排查

*   **SSH 连接失败**: 确保你的 SSH 密钥已添加到服务器，且可以通过 `ssh zhoujiazhen@127.0.0.1 -p 6000` 直接登录。
*   **权限错误**: 确保你对服务器上的目标目录有写权限。
*   **Git 冲突**: 如果 `pull` 时出现冲突，请先手动解决冲突 (`git merge` 或 `git rebase`)，然后再继续。
