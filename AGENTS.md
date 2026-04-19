# Repository Guidelines

## Instruction Hierarchy
- `AGENTS.md` is the canonical repository instruction file. Stable, cross-agent rules belong here.
- `CLAUDE.md` and `GEMINI.md` are thin agent-specific wrappers. They should only add tool-specific reminders and point back to `AGENTS.md` for repository rules.
- Repository topology, single-source-of-truth wording, and sync command names are summarized in `docs/env/repository_governance.md`; if that document and `AGENTS.md` diverge, `AGENTS.md` wins.
- Keep live status in `TRACKING.md`, experiment registry in `docs/experiments/registry.md` and `docs/experiments/runs.csv`, reusable full analysis in categorized `reports/`, and medium/large execution plans in `docs/plan/`.
- `docs/superpowers/` is a legacy archive for historical agent-generated design/execution docs. New or still-live plans/specs must live under `docs/plan/`, not under `docs/superpowers/`.
- When repository reality changes, update `AGENTS.md` first and then trim or refresh any wrapper files that depend on it.
- Remove stale guidance when it becomes misleading instead of endlessly appending new caveats below it.

## Project Structure & Module Organization
- Active local canonical trunk is `/home/zhengwei/project/python/TransChromBP`. If the old `/home/zhengwei/project/python/chromBPNet` path still exists, treat it only as a compatibility symlink / rollback path, not as the canonical root.
- This repo is the TransChromBP main repository and project archive; official ChromBPNet source lookup and reproduction currently use the 6000 external repo at `/data1/zhoujiazhen/bylw_atac/chromBPNet`.
- Local official ChromBPNet payload is retired; official ChromBPNet source lookup and reproduction live only in the 6000 external repo at `/data1/zhoujiazhen/bylw_atac/chromBPNet`. Do not assume a `/data1/zhoujiazhen/bylw_atac/chrombpnet_official` alias exists unless `TRACKING.md` explicitly says it has been validated.
- The repository topology is fixed as: local canonical trunk, 6000 official lookup root, 6000 runtime workspace, and 6002 runtime mirror. Do not describe this as three-way peer sync; it is a single-source-of-truth plus one-way publish model.
- `scripts/` holds utility scripts: `sync_project.sh`, `benchmark/`, data prep launchers, `setup_report_env.sh`.
- `workflows/` holds end-to-end bash workflows plus `tutorial/` step scripts.
- `tests/` contains shell-based integration checks.
- `docs/` organizes durable documentation by topic: `plan/`（live canonical 实验计划）, `research/`（稳定研究摘要）, `archive/ops/`（执行日志 / runbook 归档）, `env/`（环境配置）, `learning/`（学习资料）. 按需新增子目录，不预建空目录。
- `docs/plan/` 根目录只保留仍在使用的 canonical 计划；历史计划统一迁入 `docs/plan/archive/<topic>/`。
- `docs/superpowers/` only keeps historical agent-planning artifacts that still matter for auditability; treat it as an archive namespace, not a second live planning root.
- `docs/superpowers/archive_manifest.md` records archive-path to live-successor mapping when historical agent docs are still worth keeping for auditability.
- `docs/experiments/registry.md` is the canonical family/workstream index; `docs/experiments/runs.csv` is the canonical run-level manifest.
- `reports/` stores reusable analysis sources (`.tex`/`.md`) with `assets/` for result summaries (csv/png/json). Root `reports/` should normally only expose `README.md`; reusable content belongs under `analysis/`, `closeout/`, `governance/`, `handoff/`, `external/`, or `paper/`. LaTeX build products and PDF outputs stay local-only unless explicitly needed for delivery.
- `references/` stores lightweight indexes for local-only lookup material; keep the actual heavy files under `references/local-only/`.
- `vendor/transchrombp/` is the versioned local snapshot of the TransChromBP codebase and helper scripts. Local TransChromBP imports should come from this snapshot or from the relevant 6000/6002 runtime workspace, not from the root repo as an official ChromBPNet package.
- 6000 的实际 TransChromBP 运行仓使用 `src/transchrombp/` 布局，不是本地归档仓的 `vendor/transchrombp/` 布局；核对远端代码、运行脚本或环境时不要机械套本地路径。
- `tmp_remote_edit/` is only the staging area for remote file edits and transient copies — not the final archive.
- `images/` and `README.md` provide documentation assets and usage notes.

## Build, Test, and Development Commands
- `pip install -r requirements.txt` installs shared Python deps for the repo.
- If you need local TransChromBP imports, set `PYTHONPATH=vendor/transchrombp:$PYTHONPATH` in the active workspace, or use the relevant remote runtime workspace on 6000/6002.
- Official ChromBPNet CLI, code lookup, and reproduction should be done in `/data1/zhoujiazhen/bylw_atac/chromBPNet` on 6000, or through an explicitly exported `CHROMBPNET_OFFICIAL_ROOT` that points to the current validated official root.
- `./scripts/sync_project.sh publish-runtime-6000|publish-runtime-6002 [--dry-run]` are the only supported publish commands; do not use ambiguous `deploy` wording.
- `./scripts/sync_project.sh pull-results-6000|pull-results-6002 [--dry-run]` are the only supported result-sync commands; do not sync `outputs/` wholesale.
- `bash workflows/tutorial/step1_download_bams_and_peaks.sh /path/to/data` downloads tutorial inputs.
- `bash tests/full_workflow.sh 0` executes the full tutorial workflow on GPU 0 (long, GPU-heavy).

## 环境激活（6000）
- 6000 上至少区分两套环境，不要混用：
  - 官方 ChromBPNet CLI、`modisco`、`bedGraphToBigWig`、官方 data prep / reproduction：用 `chrombpnet` 环境。
  - TransChromBP / foundation 训练、评估、代码自检：优先用 `.venvs/genos-1.2b`，并指向远端运行仓的 `src/` 布局。
- `chrombpnet` 环境适合官方工具链，但当前**不包含** foundation 训练所需的 `torch`；凡是要跑 `TransChromBP` foundation 代码、最小复现实验或 runtime 自检，默认不要用这套环境。

### 6000: 官方 ChromBPNet / 工具链
```
export CHROMBPNET_ENV=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet
export PATH="$CHROMBPNET_ENV/bin:$PATH"
export LD_LIBRARY_PATH="$CHROMBPNET_ENV/lib:$LD_LIBRARY_PATH"
export PYTHONPATH=/data1/zhoujiazhen/bylw_atac/TransChromBP/vendor/transchrombp:$PYTHONPATH
```

### 6000: TransChromBP / foundation runtime
```
export TRANSCHROMBP_ENV=/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b
export PATH="$TRANSCHROMBP_ENV/bin:$PATH"
export LD_LIBRARY_PATH="$TRANSCHROMBP_ENV/lib:$LD_LIBRARY_PATH"
export PYTHONPATH=/data1/zhoujiazhen/bylw_atac/TransChromBP/src:$PYTHONPATH
```

## Coding Style & Naming Conventions
- Python is the primary language; preserve existing indentation (tabs are used in several modules).
- Keep function/module naming in `snake_case`; avoid reformatting files unless required.
- Bash scripts live under `scripts/`, `workflows/`, and `tests/`; prefer descriptive, verb-led script names.

## Testing Guidelines
- Tests are shell scripts in `tests/` and require external downloads and GPU resources.
- `bash tests/genomewide_gc_bin_test.sh` validates genomewide GC binning via md5 checks.
- `bash tests/test_pred_to_bigwig.sh` depends on outputs from the full workflow.

## Commit & Pull Request Guidelines
- Commit messages are short and lowercase, typically imperative (e.g., "update readme").
- Keep commits focused and note data/model changes explicitly.
- PRs should include: a brief summary, commands run, hardware notes, and links to any data sources.

## Long-running Jobs (Downloads/Training)
- Do not wait interactively for long downloads or training runs.
- Run jobs in the background with logs, for example:
  `nohup bash workflows/tutorial/step6_train_chrombpnet_model.sh ... > logs/train.log 2>&1 &`
- 即使已有程序在跑，只要目标 GPU 还有明显余量、没有跑满，就可以继续启动新任务；启动前先用 `nvidia-smi` 确认可用余量，并根据情况显式指定 `CUDA_VISIBLE_DEVICES`、适当下调 `batch_size`，避免把现有任务直接挤爆。
- 在 6000 上启动长跑训练时，如果代码路径已支持 DDP 且两张 A6000 都空闲或仍有明显余量，默认优先双卡；不要无说明地把本应双卡的长跑训练落成单卡。
- 启动长任务前，必须把任务预计耗时纳入是否“现在就启动”的判断：结合当前时间、GPU 余量、串并行队列、以及后续需要人工跟进的时间窗，明确评估现在开跑是否合适；不要只因为 GPU 空闲就直接启动。
- 若一个串行任务包含“单卡预处理/缓存构建 + 双卡训练”两个阶段，必须在 launcher、日志和 `TRACKING.md` 明确写出“当前单卡是预处理，后续训练会切到双卡”，避免误判为训练没有用双卡。
- 启动后台任务后，向用户汇报时必须包含预计结束时间；若只能估区间，也要明确给出预计完成时间窗，不能只报日志路径和监控命令。
- Tell the user how to check progress: `tail -f logs/train.log`, `grep -i "Finished" logs/train.log`.
- Pause for user confirmation that the run finished before proceeding.
- 开始处理新任务前先查看并按最新进度更新 `TRACKING.md`（尤其”在做事情清单”和”下载资源清单”）。
- 文档同步是硬 gate，不是可选收尾：开始实质分析前、启动后台任务后、以及形成新结论后，都必须在同一轮内更新文档。
- 若启动了新的后台任务，必须立即把以下信息写入 `TRACKING.md`：任务目的、run name、机器/ GPU、关键配置差异、日志路径、预计结束时间、下一步。
- 若本轮产出了超过“状态 + 下一步 + 一句话结论”粒度的分析、实验设计、证据链或论文口径调整，必须同步新增或更新 `docs/plan/` / `reports/` 文档；不能只在回复里口头说明。
- 一个任务只有在“实际操作已执行 + 对应文档已更新”后才算完成；如果代码/训练已跑但文档未回写，应明确标记为未完成并继续补文档。
- 维护时遵循三层分工：
  - `TRACKING.md`：live 状态与短结论，只保留 `进行中` / `待处理` / `待验证` 条目
  - `TRACKING_archive.md`：阶段性完成事项归档
  - `docs/experiments/registry.md` / `docs/experiments/runs.csv`：family/workstream 与 run 的 canonical 索引
  - `reports/`：可复用的完整分析报告（`.tex`/`.md`），含 `assets/` 下的结果摘要
- `已完成` / `已合并` 条目堆积时先做归档整理，再继续追加新状态。
- 超出”状态 + 下一步 + 一句话结论”粒度的分析，应改写为 `reports/` 下的独立报告，`TRACKING.md` 只保留高层结论与报告链接。
- 适合单独出报告的内容：多 run 对照、gate/停跑原因、跨实验线比较、路线结论链、后续实验建议。
- 启动或收口任何 run 时，除了更新 `TRACKING.md` / `reports/`，还必须同步更新 `docs/experiments/runs.csv`；若 family/workstream 的状态、canonical branch 或 next action 变化，再同步更新 `docs/experiments/registry.md`。
- `runs.csv` 的 `commit_sha` 默认记录本地 canonical archival commit；当 6000/6002 运行仓缺乏可靠 git 历史时，不要伪造“远端真实 commit”，而要在 `notes` 中注明这是归档快照 commit。
- 关键实验里程碑使用 annotated tag：`run/<machine>/<run_id>`、`closeout/<family_id>/<yyyymmdd>`、`snapshot/<branch_slug>/<yyyymmdd>`。
- 挂载中的 worktree 必须先在 `docs/experiments/registry.md` 里有对应 family/workstream 记录，确认 closeout tag、report 与 run manifest 已齐全后，才允许收起明显过时的 worktree。

## Agent Response Language
- Reply in Chinese for all user-facing responses and explanations.
- 在直接回答用户当前问题后，若上下文允许，继续思考并补 1-2 句“下一步最值得做什么”的建议；如果当前没有合理下一步，也要明确说明暂时不建议推进。

## 工作目录记录
- 6000（A6000 × 2）: `ssh zhoujiazhen@127.0.0.1 -p 6000` → `/data1/zhoujiazhen/bylw_atac`
- 6002（RTX 3080）: `ssh zhengwei@127.0.0.1 -p 6002`（密钥 `/home/zhengwei/.ssh/codex_6002_ed25519`）
- TransChromBP 远端代码：6000 `/data1/zhoujiazhen/bylw_atac/TransChromBP`、6002 `/home/zhengwei/bylw_atac/TransChromBP`
- 自定义数据集：6000 `/data1/zhoujiazhen/bylw_atac/chrombpnet_datasets/`、6002 `/home/zhengwei/bylw_atac/chrombpnet_datasets/`

## 远端文件传输与写回规范
- 从 Windows/PowerShell 向 6000/6002 写远端文件时，首选“先生成本地临时文件，再用 `scp`/重定向上传，再在远端原子替换”的方式；不要把大段正文直接内嵌进 `ssh "..."` 命令字符串。
- 小型 ASCII shell 脚本可使用单引号 here-doc：`ssh ... 'cat > /path/file <<'"'"'EOF'"'"' ... EOF'`；正文中只要包含中文、Markdown、反引号或大量引号，就不要用这种内联方式。
- 对包含中文、Markdown、反引号、`\n`、反斜杠或多层引号的内容，优先使用 `base64` 负载或 `scp` 传输；避免让 PowerShell 参与正文转义。
- 禁止使用 PowerShell 双引号字符串把整段正文直接包进 `ssh` 命令，尤其是正文里含 `` ` ``、`\n`、中文或 YAML/Markdown 时；这类写法容易把反引号解释成转义并污染文件内容。
- 需要在远端改现有文本文件时，优先上传完整新文件后再 `mv` 覆盖，或在远端运行简短的 `python3`/`perl` 脚本处理；不要在本地命令串里手写复杂替换逻辑。
- 每次远端写回后必须立即验证：至少检查 `wc -l`、`sed -n '1,40p'` 或 `tail -n 20`；如果是脚本，还要补 `bash -n`/解释器语法检查，确认没有出现整文件一行、空文件或乱码。
