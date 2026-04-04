# GEMINI.md

本文件只补充 Gemini CLI 的代理适配信息。仓库通用规则、目录结构、命令、文档协议、远端路径与长任务规范，一律以 [AGENTS.md](AGENTS.md) 为准。

## 使用顺序

1. 先读 [AGENTS.md](AGENTS.md) 获取仓库主规范
2. 再读 [TRACKING.md](TRACKING.md) 获取当前 live 状态
3. 需要完整证据链时，再看 `reports/` 与 `docs/plan/`

## Gemini CLI 特有提醒

- 当前项目的大量训练和数据处理任务在远端服务器上执行，本地主要用于代码编辑和报告生成。
- 判断训练/GPU 状态时，优先登录对应远端；不要只在本地运行 `nvidia-smi`。
- 涉及远端写回时，继续遵循 `AGENTS.md` 中“本地临时文件 -> `scp`/上传 -> 远端原子替换 -> 立即验证”的流程。

## 远端快速入口

- 6000（A6000 × 2）: `ssh zhoujiazhen@127.0.0.1 -p 6000` → `/data1/zhoujiazhen/bylw_atac`
- 6002（RTX 3080）: `ssh -i /home/zhengwei/.ssh/codex_6002_ed25519 zhengwei@127.0.0.1 -p 6002` → `/home/zhengwei/bylw_atac`

## 回复语言

- 面向用户统一使用中文。
