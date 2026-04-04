# TransChromBP Snapshot

该目录保存当前版本化的 `TransChromBP` 代码快照与辅助脚本。

约定：

- `transchrombp/`：Python 包、配置、训练与评估入口。
- `scripts/`：与该快照配套的辅助脚本，例如 external selector、ablation launcher、评估脚本。
- 本目录是本仓库内关于 `TransChromBP` 的正式归档位置；不要再把 `tmp_remote_edit/` 当成权威代码源。

使用提醒：

- 若需要本地直接导入 `transchrombp`，应将 `vendor/transchrombp` 加入 `PYTHONPATH`。
- 远端运行目录仍可按需部署到独立工作区，但回收时以本目录为准。
