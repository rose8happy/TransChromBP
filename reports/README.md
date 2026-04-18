# reports

该目录只保留可复用的报告源文件和小型结果摘要，也就是仓内约定的 `md/tex/assets` 集合。

## 应版本化的内容

- `*.md`、`*.tex`
- `assets/` 下的小型 `csv`、`json`、`png`、`txt`
- 为了生成图表或复核结论而保留的轻量脚本与配置

## 默认不进 Git 的内容

- `reports/` 根目录的 LaTeX / PDF 构建产物：`*.pdf`、`*.aux`、`*.log`、`*.out`、`*.toc`、`*.fls`、`*.fdb_latexmk`、`*.xdv`
- 临时编译残留，例如 `texput.log`

## 例外规则

- 如果 `reports/assets/` 下需要保留一个“看起来像日志”的文件作为脚本真实输入，可以保留，但必须让相邻脚本或报告正文能看出它的用途。
- 需要长期复用的完整分析，放独立报告；`TRACKING.md` 只保留状态、下一步和一句话结论。
