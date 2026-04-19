# reports

该目录只保留可复用的报告源文件和小型结果摘要，也就是仓内约定的 `md/tex/assets` 集合。

## 目录分工

- `analysis/`: 可复用分析、阶段性复盘、专题技术结论
- `closeout/`: 实验 family / infra 收口、gate 结论、最终 closeout
- `governance/`: 仓库治理、清理、审计、patch ledger、结构性收口
- `handoff/`: 会话接续、仓库快照、handoff 材料
- `external/`: 对外咨询包、外部模型回流、外部技术路线评审
- `paper/`: 论文主稿、提纲、LaTeX 报告与 paper-facing supporting materials
- `assets/`: 小型 `csv/json/png/txt` 与辅助脚本

根目录默认只保留本 README；不要再往 `reports/` 根目录直接堆新的 loose `md/tex` 文件。

## 应版本化的内容

- `*.md`、`*.tex`
- `assets/` 下的小型 `csv`、`json`、`png`、`txt`
- 为了生成图表或复核结论而保留的轻量脚本与配置

## 默认不进 Git 的内容

- LaTeX / PDF 构建产物：`*.pdf`、`*.aux`、`*.log`、`*.out`、`*.toc`、`*.fls`、`*.fdb_latexmk`、`*.xdv`
- 临时编译残留，例如 `texput.log`

## 放置规则

- 需要长期复用的完整分析，放 `analysis/`；`TRACKING.md` 只保留状态、下一步和一句话结论。
- family 或 infra 的最终 verdict，放 `closeout/`。
- 仓库结构、路径、审计、cleanup、patch ledger，放 `governance/`。
- 接续说明和冻结快照，放 `handoff/`。
- 对外咨询包、外部评审和模型回流，放 `external/`。
- 论文主稿与 `.tex` supporting 文档，放 `paper/`。

## 例外规则

- 如果 `reports/assets/` 下需要保留一个“看起来像日志”的文件作为脚本真实输入，可以保留，但必须让相邻脚本或报告正文能看出它的用途，避免和 LaTeX build 日志混淆。
- 如果某个 `reports/external/<bundle_dir>/` 明确是对外咨询包、handoff bundle 或冻结代码/配置快照，也允许整体版本化；前提是目录内有 `README` 或上传顺序说明，能解释为什么这些代码/配置快照属于报告证据链的一部分。
