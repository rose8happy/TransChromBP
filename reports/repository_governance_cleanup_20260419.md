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

## 后续补收口

1. 从临时保护分支 `master-wip-20260419-pre-governance-merge` 只挑回真正仍有价值的内容，而不是整支回滚：
   - `TRACKING.md`、`docs/experiments/registry.md`、`docs/experiments/runs.csv` 的 factor-ladder / loss-balance 最新状态
   - `reports/chatgpt_bundle_loss_balance_20260416/` 外发专题包
   - `reports/repository_hygiene_cleanup_20260417.md`
   - `reports/assets/chrombpnet_bias_training_20260123/` 小型正式资产
   - `vendor/transchrombp/scripts/summarize_loss_balance_selectors.py`
2. 明确不回流会冲掉治理主线的旧内容，例如旧版 `scripts/sync_project.sh`、回顶层的历史单文件路径、以及删除 `docs/env/repository_governance.md` 的改动。
3. 修正 `.gitignore` 对 `docs/env/` 的误伤：保留 `env/` 目录默认忽略的意图，但显式放开 `docs/env/**`，并把 `docs/env/transchrombp_genos_env.md` 纳入 Git。
4. 完成一轮 mounted worktree 审计并写入 `reports/worktree_registry_audit_20260419.md`：
   - 冗余的 `loss-balance-bundle-20260416` 挂载已删除；
   - dirty 的 `loss-balance-20260417` 先冻结成提交 `c9066b6`，补齐 `run/`、`closeout/`、`snapshot/` tags，再卸载 live worktree；
   - `loss_balance_curriculum` 后续若要复查 exact local dynamic-count 实现，不再依赖 live worktree，而是走 `snapshot/loss-balance-20260417/20260419`。

## 当前结论

- 本地 canonical naming / entrypoint 已统一收口为 `TransChromBP`；当前应把 `/home/zhengwei/project/python/TransChromBP` 视为规范化本地入口。
- 6000 `/data1/zhoujiazhen/bylw_atac/chromBPNet` 只承担官方 ChromBPNet 查阅 / 复现。
- 6000 与 6002 的 `TransChromBP` 目录默认都只是 runtime/mirror，远端热修必须先回收进本地仓。
- 本地目录的物理 rename / symlink 过渡由独立的 `local-repo-rename` workstream 处理，不应和本轮同步接口重写绑在一起。
- 保护分支 `master-wip-20260419-pre-governance-merge` 只剩“把治理主线打回去”的旧内容；在完成上述 selective curation 之后，已于 `2026-04-19` 删除。

## 后续最值得做的事

1. 对现有 mounted worktree 做一轮 registry 对账，决定哪些 family 已经满足 closeout 条件、可以收起。
2. 如果仍想把本地根目录从 `chromBPNet` 改成更直观的名字，单开一轮 rename，把绝对路径、worktree 挂载和文档链接一起处理。
