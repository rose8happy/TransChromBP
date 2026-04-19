# 仓库归档卫生整理（2026-04-17）

## 目标

在不改动训练语义和实验配置的前提下，先把主档案仓里最容易制造混乱的部分收紧：live 文档口径、`reports/` 目录约束、以及误留在工作树里的报告构建副产物。

## 本轮实际整理

1. 把 `README.md` 和 `DEVELOPMENT.md` 里仍残留的官方 ChromBPNet 根路径统一到当前已验证的 `/data1/zhoujiazhen/bylw_atac/chromBPNet`。
2. 新增 `reports/README.md`，明确 `reports/` 只保留可复用的 `md/tex/assets` 源文件，小型摘要允许版本化，根目录 LaTeX / PDF 构建产物默认本地自留。
3. 为根目录报告构建产物补了 `.gitignore` 规则，并清理掉当前误留在工作树里的 `reports/*.pdf`、`*.aux`、`*.log`、`*.out`、`*.toc`、`*.fls`、`*.fdb_latexmk`、`*.xdv`。
4. 逐个判定根目录 6 份参考材料后，做了明确分流：
   - `2024.12.25.630221v2.full.pdf`、`s41586-025-10014-0.pdf`、`libssl1.1_1.1.1f-1ubuntu2_amd64.deb` 归类为本地自留 lookup material，迁入 `references/local-only/`。
   - `bias_training_curve.svg`、`bias_training_curve_detailed.svg`、`bias_training_log.csv` 因为被 `reports/analysis/chrombpnet_tutorial_bias_training_analysis_20260123.md` 显式引用且体积很小，归类为正式小型报告资产，迁入 `reports/assets/chrombpnet_bias_training_20260123/`。

## 本轮刻意不做的事

- 不批量改写历史 `docs/plan/` 文档里曾经讨论过的 `chrombpnet_official` alias。那部分属于历史设计上下文；只要 live 文档已经统一、且历史文档内部明确写出“alias 尚未验证”，就先不做大面积重写。
- 不改训练代码、launcher、远端运行目录或实验 manifest。
- 不继续扩大这轮整理范围；像 `bpnet_model_annotated.py`、`visualize_bpnet.py`、`REPRODUCTION_NOTES.md`、`TRAINING_ANALYSIS.md` 这类“位置可能不理想但仍有真实内容”的顶层文件，暂时只记录为下一轮候选，不在这轮顺手搬动。

## 下一步建议

如果还要继续做仓库整理，最值得单开一轮的不是再扫一遍 `.gitignore`，而是继续处理剩余顶层单文件的结构归属：

- 像 `bpnet_model_annotated.py`、`visualize_bpnet.py` 更像脚本或代码摘录，应该判断是否迁到 `scripts/` / `vendor/` / `docs/learning/`。
- 像 `TRAINING_ANALYSIS.md`、`REPRODUCTION_NOTES.md` 更像报告或长说明，应该判断是否迁到 `reports/` / `docs/`。
