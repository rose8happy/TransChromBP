# Repository Guidelines

## Project Structure & Module Organization
- `chrombpnet/` contains the Python package: CLI entrypoints, pipelines, helpers, training, and evaluation code.
- `workflows/` holds end-to-end bash workflows plus `tutorial/` step scripts.
- `tests/` contains shell-based integration checks.
- `images/` and `README.md` provide documentation assets and usage notes.

## Build, Test, and Development Commands
- `pip install -r requirements.txt` installs Python deps; use `pip install -e .` for editable dev.
- `chrombpnet pipeline ...` runs bias-factorized ChromBPNet training (see `README.md` for full flags).
- `chrombpnet bias pipeline ...` trains a bias model first.
- `bash workflows/tutorial/step1_download_bams_and_peaks.sh /path/to/data` downloads tutorial inputs.
- `bash tests/full_workflow.sh 0` executes the full tutorial workflow on GPU 0 (long, GPU-heavy).

## 环境激活（6000）
- 在 6000 机器上跑脚本前先加载 chrombpnet 环境，避免找不到 `modisco`、`bedGraphToBigWig` 等工具。
```
export CHROMBPNET_ENV=/data1/zhoujiazhen/bylw_atac/.mamba/envs/chrombpnet
export PATH="$CHROMBPNET_ENV/bin:$PATH"
export LD_LIBRARY_PATH="$CHROMBPNET_ENV/lib:$LD_LIBRARY_PATH"
export PYTHONPATH=/data1/zhoujiazhen/bylw_atac/chromBPNet:$PYTHONPATH
```

## Coding Style & Naming Conventions
- Python is the primary language; preserve existing indentation (tabs are used in several modules).
- Keep function/module naming in `snake_case`; avoid reformatting files unless required.
- Bash scripts live under `workflows/` and `tests/`; prefer descriptive, verb-led script names.

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
- Tell the user how to check progress: `tail -f logs/train.log`, `grep -i "Finished" logs/train.log`.
- Pause for user confirmation that the run finished before proceeding.
- 开始处理新任务前先查看并按最新进度更新仓库根目录的 `TRACKING.md`（尤其“在做事情清单”和“下载资源清单”）。

## Agent Response Language
- Reply in Chinese for all user-facing responses and explanations.

## 工作目录记录
- 使用 `ssh zhoujiazhen@127.0.0.1 -p 6000` 访问 `/data1/zhoujiazhen/bylw_atac`，这是我们的工作目录。
- 可用的另一台机器：`ssh zhengwei@127.0.0.1 -p 6002`。
- 连接 6002 使用的密钥：`/home/zhengwei/.ssh/codex_6002_ed25519`（公钥：`/home/zhengwei/.ssh/codex_6002_ed25519.pub`）。

## 远端文件传输与写回规范
- 从 Windows/PowerShell 向 6000/6002 写远端文件时，首选“先生成本地临时文件，再用 `scp`/重定向上传，再在远端原子替换”的方式；不要把大段正文直接内嵌进 `ssh "..."` 命令字符串。
- 小型 ASCII shell 脚本可使用单引号 here-doc：`ssh ... 'cat > /path/file <<'"'"'EOF'"'"' ... EOF'`；正文中只要包含中文、Markdown、反引号或大量引号，就不要用这种内联方式。
- 对包含中文、Markdown、反引号、`\n`、反斜杠或多层引号的内容，优先使用 `base64` 负载或 `scp` 传输；避免让 PowerShell 参与正文转义。
- 禁止使用 PowerShell 双引号字符串把整段正文直接包进 `ssh` 命令，尤其是正文里含 `` ` ``、`\n`、中文或 YAML/Markdown 时；这类写法容易把反引号解释成转义并污染文件内容。
- 需要在远端改现有文本文件时，优先上传完整新文件后再 `mv` 覆盖，或在远端运行简短的 `python3`/`perl` 脚本处理；不要在本地命令串里手写复杂替换逻辑。
- 每次远端写回后必须立即验证：至少检查 `wc -l`、`sed -n '1,40p'` 或 `tail -n 20`；如果是脚本，还要补 `bash -n`/解释器语法检查，确认没有出现整文件一行、空文件或乱码。
