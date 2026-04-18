# 仓库治理与同步收口（2026-04-19）

## 目标

把当前仓库从“本地 / 6000 / 6002 混合心智”收成明确的单一真源模型，避免继续把 runtime 镜像、官方 lookup 仓和 canonical trunk 当成同一类节点。

## 本轮落地

1. 新增 `docs/env/repository_governance.md` 作为 repository topology 与 sync contract 的集中入口。
2. 重写 `scripts/sync_project.sh` 的接口：
   - 新增 `publish-runtime-6000`
   - 新增 `publish-runtime-6002`
   - 新增 `pull-results-6000`
   - 新增 `pull-results-6002`
   - 新增 `status-all`
   - 旧的 `deploy` / `download_results` 改为直接报 deprecated 提示
3. 更新 `README.md`、`DEVELOPMENT.md`、`AGENTS.md`，统一为“本地 canonical trunk + 官方 lookup root + 两个 runtime”的口径。
4. 为 `reports/`、`references/`、`references/local-only/` 补齐 README，明确 `md/tex/assets`、`local-only` 和 staging 边界。
5. 把几个不应长期留在仓库根目录的单文件归位：
   - `TRAINING_ANALYSIS.md` → `reports/chrombpnet_tutorial_bias_training_analysis_20260123.md`
   - `REPRODUCTION_NOTES.md` → `docs/research/chrombpnet_reproduction_notes.md`
   - `bpnet_model_annotated.py` → `docs/learning/assets/bpnet_model_annotated.py`
   - `visualize_bpnet.py` → `scripts/visualize_bpnet.py`

## 当前结论

- 本地 `/home/zhengwei/project/python/chromBPNet` 继续作为唯一 canonical trunk。
- 6000 `/data1/zhoujiazhen/bylw_atac/chromBPNet` 只承担官方 ChromBPNet 查阅 / 复现。
- 6000 与 6002 的 `TransChromBP` 目录默认都只是 runtime/mirror，远端热修必须先回收进本地仓。
- 目录物理改名值得做，但不应和本轮同步接口重写绑在一起。

## 后续最值得做的事

1. 对现有 mounted worktree 做一轮 registry 对账，决定哪些 family 已经满足 closeout 条件、可以收起。
2. 如果仍想把本地根目录从 `chromBPNet` 改成更直观的名字，单开一轮 rename，把绝对路径、worktree 挂载和文档链接一起处理。
