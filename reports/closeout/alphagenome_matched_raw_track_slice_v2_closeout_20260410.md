# AlphaGenome Matched Raw-Track Slice v2 Closeout（2026-04-10）

## Verdict

- `pass`
- 这是一次 **technical / external-coordinate sidecar closeout** 的通过，说明 AlphaGenome 的 matched raw-track 链路已在 `16` 个位点的小面板上稳定跑通。
- 该结论**不表示** AlphaGenome 在模型质量上优于当前 teacher，也**不把**这条 sidecar 升级成 benchmark，更**不改变** A6000 formal gate 的判读语义。

## Run Identity

- run name：`alphagenome_matched_raw_track_slice_v2_20260410`
- machine / env：`6000 / alphagenome env`
- runtime repo：`/data1/zhoujiazhen/bylw_atac/TransChromBP`
- log：`/data1/zhoujiazhen/bylw_atac/TransChromBP/logs/alphagenome_matched_raw_track_slice_v2_20260410.log`
- output dir：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/alphagenome_pilot/alphagenome_matched_raw_track_slice_v2_20260410`
- regions csv：`scripts/alphagenome_pilot/regions_k562_tutorial_matched_panel_v2.csv`
- panel size：`16`
- output_type：`ATAC`
- ontology_term：`EFO:0002067`
- target_width：`1000`
- aggregate：`mean`

## Evidence

### 1. 运行已完成，且 closeout 所需的结束标记齐全

- 日志尾部以成功写盘标记结束：
  - `Wrote summary: outputs/alphagenome_pilot/alphagenome_matched_raw_track_slice_v2_20260410/summary.csv`
  - `Wrote metadata: outputs/alphagenome_pilot/alphagenome_matched_raw_track_slice_v2_20260410/region_metadata.jsonl`
  - `Wrote merged comparison: outputs/alphagenome_pilot/alphagenome_matched_raw_track_slice_v2_20260410/merged_locus_totals.csv`

### 2. 必需产物齐全

输出目录内已确认存在以下 closeout 必需文件：

- `summary.csv`
- `region_metadata.jsonl`
- `profiles/`
- `run_meta.json`
- `merged_locus_totals.csv`

其中 `merged_locus_totals.csv` 已成功生成，说明 AlphaGenome summary 与本地 locus totals 的合并表已形成。

### 3. 16 个 loci 全部完成，且每个位点都保留了 1 条可用的 ATAC track

- `summary.csv` 含 `16` 条数据行，对应 `16` 个目标 loci 全部完成。
- 每一行都满足 `num_tracks_used=1`，说明每个位点都保留了 `1` 条可比较的 ATAC track。
- `region_metadata.jsonl` 的头部复核显示，每个抽样位点都满足：
  - `tracks_after_ontology_filter=1`
  - `tracks_after_strand_filter=1`
  - `used_unstranded_only=true`

这说明本次 v2 面板没有出现“大部分位点因 ontology / strand filter 失效”的问题，满足 sidecar technical gate 的最低要求。

### 4. 本次 closeout 的含义仅限 external coordinate，不延伸到模型优劣

本次 closeout 只确认：

- AlphaGenome API / track matching / ontology 过滤 / profile 落盘 / merged comparison 这条链路，在 `16` 个 matched loci 上稳定可复用。
- 当前 v2 面板足以作为一个更完整的外部坐标 sidecar 保留。

本次 closeout 不确认：

- AlphaGenome 在这 `16` 个 loci 上“优于”当前模型。
- 这次 sidecar 应升级成大 benchmark。
- A6000 formal gate 应因为这次 sidecar `pass` 而提前判正、判负，或改变 compare target。

## Gate Decision

- 最终判定：`pass`
- 通过含义：`alphagenome_matched_raw_track_slice_v2_20260410` 已满足本次 technical / external-coordinate sidecar closeout 的要求。
- 当前落点：把 v2 结果作为已完成的外部坐标对照保留，不占 `6000` active slot。
- 明确边界：
  - 不把这次 `pass` 改写成模型质量胜利。
  - 不把这次 v2 结果自动升级成 benchmark。
  - 不把这次 v2 closeout 用来改写仍在运行中的 A6000 formal gate 语义。

## Closeout Conclusion

`alphagenome_matched_raw_track_slice_v2_20260410` 已完成并通过本次 technical / external-coordinate sidecar closeout。证据链清楚地表明：`16` 个目标 loci 全部成功落盘、每个位点都保留了 `1` 条可用的 `ATAC` track、`merged_locus_totals.csv` 已干净形成。该结果的定位到此为止：它是一个已完成的 AlphaGenome 外部坐标 sidecar，不是 benchmark，也不对 A6000 formal gate 的最终 verdict 做任何提前暗示。
