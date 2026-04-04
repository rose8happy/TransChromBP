# CLAUDE.md

本文件只补充 Claude Code 的代理适配信息。仓库通用规则、目录结构、命令、文档协议、远端路径与长任务规范，一律以 [AGENTS.md](AGENTS.md) 为准。

## 使用顺序

1. 先读 [AGENTS.md](AGENTS.md) 获取仓库主规范
2. 再读 [TRACKING.md](TRACKING.md) 获取当前 live 状态
3. 需要完整证据链时，再看 `reports/` 与 `docs/plan/`

## Claude Code 特有提醒

- 本地工作区主要用于代码编辑、报告生成和同步；训练、评估和大规模数据处理主要发生在远端 6000 / 6002。
- 判断 GPU 或训练状态时，优先 `ssh` 到对应远端；不要只看本地 `nvidia-smi`。
- 通过 Claude 执行远端写回时，继续遵循 `AGENTS.md` 中的“本地临时文件 -> `scp`/上传 -> 远端原子替换 -> 立即验证”流程。

## 远端快速入口

- 6000（A6000 × 2）: `ssh zhoujiazhen@127.0.0.1 -p 6000` → `/data1/zhoujiazhen/bylw_atac`
- 6002（RTX 3080）: `ssh -i /home/zhengwei/.ssh/codex_6002_ed25519 zhengwei@127.0.0.1 -p 6002` → `/home/zhengwei/bylw_atac`

## 回复语言

- 面向用户统一使用中文。
