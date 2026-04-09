# AlphaGenome Matched Raw-Track Slice v1 Closeout（2026-04-10）

## Verdict

- `pass`
- 这是一次 **technical/alignment gate** 的通过，不是性能胜负结论，也不构成“AlphaGenome 优于当前模型”的证据。

## Run Identity

- run name：`alphagenome_matched_raw_track_slice_v1_20260410`
- machine：`6000 / A6000 x2`（本次任务性质是独立 API / 推理 pilot，不是双卡训练）
- runtime repo：`/data1/zhoujiazhen/bylw_atac/TransChromBP`
- log：`/data1/zhoujiazhen/bylw_atac/TransChromBP/logs/alphagenome_matched_raw_track_slice_v1_20260410.log`
- output dir：`/data1/zhoujiazhen/bylw_atac/TransChromBP/outputs/alphagenome_pilot/alphagenome_matched_raw_track_slice_v1_20260410`
- start time：`2026-04-10 00:22:40 CST`
- end time：`2026-04-10 00:22:55 CST`
  - 说明：结束时间按最终日志与产物文件的 `mtime` 近似，足以用于本次 closeout 口径。

## Evidence

### 1. 运行已完成，而不是 launch-blocked / fail

- closeout 复核时已无匹配中的运行进程。
- 日志尾部只出现成功写盘标记：
  - `Wrote summary`
  - `Wrote metadata`
  - `Wrote merged comparison`
- 同轮复核未命中 `Traceback`、`ERROR`、`failed` 等失败标记。

### 2. 必需产物齐全

输出目录内已确认存在以下文件 / 目录：

- `summary.csv`
- `region_metadata.jsonl`
- `profiles/`
- `run_meta.json`
- `merged_locus_totals.csv`

进一步复核显示：

- `profiles/*.npz` 数量为 `4`
- `summary.csv` 有 `4` 行 locus 结果
- `region_metadata.jsonl` 有 `4` 行 metadata
- `merged_locus_totals.csv` 已成功生成并与本地 locus totals 合并

### 3. 4 个固定 loci 全部成功

`summary.csv` 中的 4 个 label 全部存在，且与设计一致：

- `peak_high`
- `peak_mid`
- `nonpeak_high`
- `nonpeak_mid`

每个位点的 `num_tracks_used` 都为 `1`，说明 ontology/filter 后每个位点都保留了至少 1 条可比较的 `ATAC` track，满足本次 gate 的最低要求。

本轮 metadata 复核显示 4 个位点保留的 track 名称均为：

- `EFO:0002067 ATAC-seq`

### 4. 对齐产物已形成，但不做性能优劣宣称

`merged_locus_totals.csv` 已成功生成，说明 AlphaGenome summary 与本地已有 locus totals 已形成统一对照表。

本次 closeout 只确认：

- API / 取数 / ontology 过滤 / profile 落盘 / merge 这条链路已稳定跑通
- 4 个固定 loci 都能形成可读的 matched raw-track coordinate

本次 closeout 不确认：

- AlphaGenome 在更大样本上优于当前模型
- AlphaGenome 应升级为大 benchmark 或主线评测
- 当前模型需要因为这次 pilot 而改变训练方向

## Gate Decision

- 最终判定：`pass`
- 通过含义：`AlphaGenome matched raw-track slice v1` 已满足技术/对齐 gate，可以作为一个窄而可复用的 external coordinate 基线保留。
- 唯一允许的下一步：若继续该方向，只能扩到一个稍大的 matched panel，规模约 `12-20 loci`。
- 明确禁止：
  - 把这次 `pass` 改写成“模型质量胜利”
  - 直接扩成大 benchmark
  - 在没有新设计的前提下无边界补跑

## Closeout Conclusion

`alphagenome_matched_raw_track_slice_v1_20260410` 已完成并通过 technical/alignment gate。对 `6000` 而言，这条 pilot 已经从“候选 backlog”转成“已完成的窄 external-coordinate 基线”；如果继续推进，下一步只能是 `12-20 loci` 的小幅扩面，而不是更大的 benchmark 或训练线。
