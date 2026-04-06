# ChromBPNet 官方依赖台账（2026-04-06）

## 1. 当前总契约

- canonical official root：`/data1/zhoujiazhen/bylw_atac/chrombpnet_official`
- active workflow 的默认官方源码根已经切到 6000 的 external official root；在后续删除本地官方 payload 之前，本仓树里仍可能保留历史 `chrombpnet/` 文件，但 operator 不应把它当作运行时权威来源。
- 本文只记录 externalization 之后仍然依赖的官方运行时文件族，方便 operator 在 6000/6002 上核对路径、补丁来源和故障定位。

上面的根路径与仓库角色收口，来自 [chrombpnet_official_externalization_design_20260406.md](../docs/plan/chrombpnet_official_externalization_design_20260406.md) 的目标态定义：6000 `/data1/zhoujiazhen/bylw_atac/chrombpnet_official` 是唯一官方参考仓；本台账只固定运行时权威来源，不宣称本仓的历史官方 payload 已经物理删除。

## 2. 仍在使用的官方文件族

| 官方文件 / 文件族 | 当前用途 | canonical runtime path | 当前必须保留的 patch / delta | 为什么现在还重要 | 现有证据 |
| --- | --- | --- | --- | --- | --- |
| `chrombpnet/training/predict.py` | strict-compare official arm、`select_best_epoch.py` 的逐 checkpoint 外部评估、`run_paper_aligned_fast_1seed.sh` 的 official metrics 采集 | `/data1/zhoujiazhen/bylw_atac/chrombpnet_official/chrombpnet/training/predict.py` | 必须保留当前 `--split` 支持，否则 selector 会重新退回 test leakage；当前 official metrics JSON contract 也仍由这份入口驱动 | strict compare 的 official 侧评估仍直接调用这份入口；selector 是否在 `valid` split 选 best、以及输出 JSON 是否带上当前官方指标字段，都取决于这份脚本 | [strict_compare_code_audit_20260327.md](strict_compare_code_audit_20260327.md) 记录了 selector 泄漏问题与 `predict.py` 增加 `--split` 的修复；[select_best_epoch.py](../scripts/paper_aligned_repro/select_best_epoch.py) 与 [run_paper_aligned_fast_1seed.sh](../scripts/paper_aligned_repro/run_paper_aligned_fast_1seed.sh) 现在都显式从 `CHROMBPNET_OFFICIAL_ROOT` 解析它 |
| `chrombpnet/training/metrics.py` | official `predict.py` 运行时的 profile/count/classification metrics 生成 | `/data1/zhoujiazhen/bylw_atac/chrombpnet_official/chrombpnet/training/metrics.py` | 必须保留当前 classification 指标字段与 `metrics_only` 语义；如果回退到旧版，official compare 会再次出现缺字段或额外 PNG 副产物 | 我们当前引用的 official 指标口径仍由这份文件定义；尤其是 `metrics_only` 行为、JSD 实现细节和 classification 字段补齐，都依赖远端官方仓中的这份实现 | [strict_compare_code_audit_20260327.md](strict_compare_code_audit_20260327.md) 记录了 `metrics_only` 下仍写 PNG 的问题并定位到 `metrics.py`；[strict_compare_followup_assessment_20260327.md](strict_compare_followup_assessment_20260327.md) 记录了 official JSD 计算仍来自 `metrics.py`；[tutorial_L3_shared_region_closure_20260330.md](tutorial_L3_shared_region_closure_20260330.md) 记录了 6000 远端 `predict.py` / `metrics.py` 版本落后会导致 official classification 指标缺失 |
| `chrombpnet/helpers/make_gc_matched_negatives/*` | dataset prep 的 GC helper 来源、tutorial `step3_get_background_regions.sh` 的背景区间生成、6002 数据准备时从 6000 staged copy helper | `/data1/zhoujiazhen/bylw_atac/chrombpnet_official/chrombpnet/helpers/make_gc_matched_negatives/` | 必须把这组 helper 视为 external authority：6000 直接从官方根调用，6002 只允许 staged copy 这三个 helper，tutorial step3 也只认 `MAKE_GC_MATCHED_NEGATIVES_SCRIPT` 或 `CHROMBPNET_OFFICIAL_ROOT` | active workflow 已不再把本仓里的这套 helper 当成权威来源；6000 dataset prep 现在直接从官方根调用，6002 也仍要从 6000 官方根拉取 `get_gc_content.py` / `get_gc_matched_negatives.py` / `get_genomewide_gc_bins.py` | [chrombpnet_official_externalization_design_20260406.md](../docs/plan/chrombpnet_official_externalization_design_20260406.md) 把这组 helper 列为外置化后必须保留的 active official family；[run_remote_chrombpnet_dataset_prep.sh](../scripts/run_remote_chrombpnet_dataset_prep.sh)、[start_6000_chrombpnet_dataset_prep.sh](../scripts/start_6000_chrombpnet_dataset_prep.sh)、[start_6002_chrombpnet_dataset_prep.sh](../scripts/start_6002_chrombpnet_dataset_prep.sh) 与 [step3_get_background_regions.sh](../workflows/tutorial/step3_get_background_regions.sh) 都把它作为运行时 helper 来源 |

## 3. 操作含义

- 需要查官方实现时，不再回本仓找 `chrombpnet/...`；直接去 6000 的 `/data1/zhoujiazhen/bylw_atac/chrombpnet_official`。
- 需要跑 official strict compare、外部 best-epoch selector 或 tutorial 背景区间生成时，先确认 `CHROMBPNET_OFFICIAL_ROOT=/data1/zhoujiazhen/bylw_atac/chrombpnet_official`。
- 若 official 侧再次出现“缺字段 / 指标不一致 / GC helper 找不到”，优先检查的是这三组 canonical 文件，而不是默认回本仓里的历史 `chrombpnet/` payload 找运行时权威实现。
