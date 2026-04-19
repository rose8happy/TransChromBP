# 工作目录整理计划（2026-04-05）

## 结论

需要整理，但当前不是磁盘空间告急，而是“主档案分叉 + 运行目录与源码目录边界不清”已经开始抬高后续出错概率。

- 本机 WSL 空间充足：`/` 使用率约 `19%`。
- 6000 `/data1` 空间充足：使用率约 `8%`；`ntv2_residual_short10_testfull_sidecar_20260405_170520` 已于 `2026-04-05 17:22 CST` 结束，当前已不存在运行中 sidecar 对归档动作的硬限制。
- 6002 `/` 也未到危险线：使用率约 `41%`。
- 真正需要优先收口的是“哪份代码/配置/状态记录才是唯一真源”。

## 已执行首轮清理（2026-04-05 17:19-17:35 CST）

- 本机根目录已移除可重建的 LaTeX 构建残留，减少顶层噪音。
- `tmp_remote_edit/` 已从约 `1.6M` 收口到约 `896K`，删除了重复下载脚本、临时 tarball、copy 文件和 `__pycache__`。
- 6000/6002 `TransChromBP` repo 内的 `.pytest_cache`、`__pycache__`、`.checked` 与历史 `.bak` 已清理。
- 6000 远端 `TRACKING` 入口已收成单一入口：
  - live 文件同步到 `/data1/zhoujiazhen/bylw_atac/chromBPNet/TRACKING.md`
  - `/data1/zhoujiazhen/bylw_atac/TRACKING.md` 与 `/data1/zhoujiazhen/bylw_atac/TransChromBP/TRACKING.md` 已改成只读指路说明
- 当前尚未触碰 `outputs/checkpoints`、`outputs/foundation_cache`、`outputs/metrics` 等正式结果目录。

## 已执行第二轮归档（2026-04-05 17:35-17:50 CST）

- 6000 已建立 foundation side quest 归档根：
  - `/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/archive/foundation_sidequest`
- 已迁入该归档根的 top-level 结果：
  - `foundation_cache`
  - `genos_cache`
  - `genos_cached_probes`
  - `genos_cached_probes_20260323`
  - `genos_feasibility`
  - `nt_v2_probe`
- 已迁入该归档根的 checkpoint family：
  - `caduceus_*`
  - `genos_*`
  - `ntv2_residual_short10_*`
- 上述 6000 路径全部保留原路径软链接，保证既有报告与脚本引用不直接失效。

- 6002 已建立历史 checkpoint 归档根：
  - `/home/zhengwei/bylw_atac/TransChromBP/outputs/archive/historical_6002/checkpoints`
- 6002 现有历史 checkpoint 已整体迁入 archive；只对文档明确引用的两条路径保留原路径软链接：
  - `gm12878_v2_6002_pilot20_s1234_20260324_105529`
  - `v2fix_20260323_attnpool_s1234`

- 当前结果：
  - 6000 foundation side quest 结果已有明确 archive 落点，规模约 `32G`
  - 6002 `outputs/checkpoints/` 顶层已不再堆满历史目录，只保留两条兼容软链接

## 已执行第三轮小清理与复核（2026-04-06 14:38-14:45 CST）

- 本机已删除 3 个明显垃圾/一次性残留：
  - `s41586-025-10014-0.pdf:Zone.Identifier`
  - `tmp_remote_edit/transchrombp/training/train_ddp.py.checked`
  - `tmp_remote_edit/ntv2_gate_smoke_after_idle_20260405.sh`
- 对 `tmp_remote_edit/` 做了再次审计：
  - `tmp_remote_edit/transchrombp/` 当前剩余文件已经全部能在 `vendor/transchrombp/transchrombp/` 找到对应副本；
  - 当前不宜直接把 `tmp_remote_edit/` 全部清成只剩 README，因为 `tmp_remote_edit/TransChromBP/scripts/` 里仍有 11 个历史 helper script 还没迁回正式命名空间。
- 6002 再次确认了 runtime mirror 角色：
  - `outputs/checkpoints/` 顶层只剩两条兼容软链接；
  - 其余历史 checkpoint 已在 `outputs/archive/historical_6002/checkpoints/`。
- 6000 `outputs/` 顶层当前是“active paper-facing 目录 + foundation sidequest 软链接归档”混合态；
  - 本轮没有继续大规模搬移 `checkpoints/`，避免破坏既有报告、日志和脚本对旧路径的依赖。

## 现状审计

| 位置 | 主要现状 | 判断 |
|---|---|---|
| 本机 `/home/zhengwei/project/python/TransChromBP` | 仓库约 `2.7G`，大头是 `.venv 2.2G`、`.venv-report 490M`、`.git 412M`；根目录当前剩余的是被 ignore 的参考 PDF/安装包；`tmp_remote_edit/` 约 `900K`，其中 `transchrombp/` 已与 `vendor/` 对齐，但 `TransChromBP/scripts/` 仍保留 11 个未迁移的历史 helper script。 | 不是容量问题；真正待收口的是 staging 边界和“哪些历史脚本要迁回正式命名空间”。 |
| 6000 `/data1/zhoujiazhen/bylw_atac` | 大目录主要是 `foundation_models 312G`、`ATACseq 156G`、`chrombpnet_datasets 76G`；根目录 `TRACKING.md` 已改成只读入口说明。 | 容量不紧；状态入口已收口，但后续必须继续避免双写。 |
| 6000 `/data1/zhoujiazhen/bylw_atac/TransChromBP` | 工作树约 `88G`，几乎都在 `outputs/`：其中 `checkpoints 72G`、`foundation_cache 9.8G`、`genos_cache 4.3G`、`preprocessing 1.4G`；目录内是 git repo，但当前状态不是“少量未提交修改”，而是“初始 scaffold 被 staged，实际实验文件大批量 untracked”。目前 foundation side quest 归档层已落地。 | 仍是最需要继续整理的点，但已经从“完全混放”推进到“已有 archive 根，可继续做 active/archive 分层”。 |
| 6002 `/home/zhengwei/bylw_atac` | `TransChromBP 30G`、`chrombpnet_datasets 60G`、`logs 166M`；临时目录都不大。 | 容量不紧。 |
| 6002 `/home/zhengwei/bylw_atac/TransChromBP` | 目录内没有 `.git`，当前更像运行镜像/手工同步副本；空间主要在 `outputs/checkpoints 22G` 和 `outputs/preprocessing 7.6G`。 | 适合当 runtime mirror，不适合继续当“主档案”。 |

## 关键风险

1. `TRACKING.md` 单入口已经收口，但后续必须持续维持这一约束。若再在 `/data1/zhoujiazhen/bylw_atac/` 或 `/data1/zhoujiazhen/bylw_atac/TransChromBP/` 直接追加状态，会再次把 live 状态写散。
2. 代码真源不唯一。本机仓库、6000 `TransChromBP` git 工作树、6002 无 git 副本三处都承载了“当前模型代码”的不同切片。
3. 远端 `TransChromBP` 把源码、配置、运行时缓存、日志、checkpoint 长期混放在一个树下，短期可用，长期不利于归档和复现。
4. 本机根目录虽然大多已被 `.gitignore` 覆盖，但仍堆着历史 PDF、LaTeX 构建文件和 `tmp_remote_edit` staging 产物，降低可读性。

## 建议动作顺序

### P0：已完成的前置收口

1. `ntv2_residual_short10_testfull_sidecar_20260405_170520` 已完成判读与文档回写，归类为 `fail-or-unsafe`。
2. 当前默认真源已收口为本机 `/home/zhengwei/project/python/TransChromBP`；远端 `TransChromBP` 按 runtime / output workdir 使用。
3. 6000 的 live `TRACKING.md` 已统一到 `/data1/zhoujiazhen/bylw_atac/chromBPNet/TRACKING.md`；根目录和 `TransChromBP` 内的 `TRACKING.md` 均已改成指路说明。

### P1：本周内完成的结构整理

1. 为远端 `TransChromBP/outputs/` 建立明确的归档层。
建议至少拆成 `active/`、`archive/` 或按实验线/日期分层，优先整理 `checkpoints/` 与 `foundation_cache/`。
2. 对 6002 明确角色。
如果它只是单卡实验镜像，就不要再把它当主代码仓；必要时保留 `outputs/`，代码改回本机统一维护后再同步。
3. 清理本机 staging 噪音。
先处理 `tmp_remote_edit/TransChromBP/scripts/` 中仍只存在于 staging 的 11 个 helper script，决定“迁回正式目录 / 单独归档 / 明确废弃”；在这一步完成前，不要把整个 `tmp_remote_edit/` 粗暴清空。随后再删除已经与 `vendor/` 完全重合的 `tmp_remote_edit/transchrombp/` 副本。

### P2：不急，但值得顺手处理

1. 把本机根目录的论文 PDF/LaTeX 构建产物收拢到更明确的位置，避免继续混在仓库根目录。
2. 统一远端日志命名和归档习惯，避免 `launch.log`、`best_test.log`、`cache_train.log` 长期散在同一层。

## 当前不建议做的事

- 不要在尚未选定真源前直接大规模清空 6000/6002 `TransChromBP/outputs/checkpoints/`。
- 不要继续默认把“远端当前目录”当作文档、代码、状态的同一入口。
- 不要在还没处理那 11 个历史 helper script 前，直接把 `tmp_remote_edit/TransChromBP/scripts/` 整树删除。
