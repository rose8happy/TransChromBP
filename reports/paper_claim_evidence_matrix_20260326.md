# TransChromBP 论文 Claim-Evidence Matrix（2026-03-26）

配套数字总表见 [paper_metric_source_table_20260326.csv](/home/zhengwei/project/python/chromBPNet/reports/assets/paper_metric_source_table_20260326.csv)。

## 1. 可写入主文的核心结论

| Claim ID | 建议写法 | 允许强度 | 主要证据 | 关键数字 | 写作注意 |
|---|---|---|---|---|---|
| C1 | Transformer 在 bias-factorized base-resolution ATAC 建模中带来真实收益 | `shows` | `V2-full` vs `V2-noTF` 三 seed matched ablation | `0.31419±0.00012 / 0.83676±0.00531` vs `0.32393±0.00021 / 0.82503±0.00489`，且 gap 始终 `~1e-4` | 这条只回答“有没有收益”，不要和 shortcut 复核混写 |
| C2 | `full/debiased gap` 是必要的 bias reliance 诊断口径 | `shows` | A/C/noTF/corrected B 四条 clean test | A=`0.00147`，C=`0.00972 / 0.00997`（双 seed），`noTF_deb2=0.00487`，`noTF_sg1_deb0=0.02682 / 0.01679`，corrected B=`0.00268` | 重点是“单看 full 不够”，不是“已证明某个强机制” |
| C3 | 在当前受控矩阵中，`debiased profile supervision` 比 `stop-gradient` 更关键 | `suggests` | A vs C（含 C 双 seed） | A=`TF + sg=false + deb2` 仍很干净；C=`TF + sg=true + deb0` 在两 seed 上都出现轻度 gap | 只能写成当前矩阵下的结论，不要上升为普适定律 |
| C4 | 最终默认模型 `corrected B` 没有明显 bias reliance，可作为 paper-facing 主模型 | `shows` | corrected B clean test + B 两个 held-out seed | corrected B gap=`0.00268`；`B_s42/B_s1234` 为 `0.31467/0.85030` 与 `0.31453/0.84885` | 这里要把“最终模型”和“机制诊断 scaffold”分开写 |
| C5 | `center pool` 是当前最稳的默认 readout 方案 | `shows` | `B_s42`、`B_s1234`、`F_s42`、`F_s1234` | `B` 两 seed 均稳定；`F_s1234` profile 退化到 `0.42691` | 不要再写成 `B/F` 并列主候选 |
| C6 | 模型在独立 GM12878 和独立 K562 数据上仍能稳定工作 | `shows` | 6002 独立真实数据 held-out | GM12878=`0.42265 / 0.80396`；K562=`0.61235 / 0.85857` | 必须明确这两条是“独立真实数据验证”，不是 tutorial 主 benchmark |
| C7 | 当前 foundation-model 接入线未显示 clean gain，应作为 appendix / future work 的受控负结果 | `suggests` | Genos 融合主表 + Caduceus matched A/B + NT v2 两轮 residual held-out gate | Genos `G1/G2/P2` 均未同时提升两项主指标；Caduceus 为 `near-null / marginal positive`；NT v2 coarse residual held-out 为 `0.3560 / 0.7729`，后续 bins16 center-aligned residual 仍只有 `0.3588 / 0.7516`，且相对 matched `short10 no-foundation control` (`0.3193 / 0.8298`) 继续落后 | 当前已否定 summary 注入、token-fusion、coarse residual 与 bins16 center-aligned residual 这些已实测 recipe；不外推到未实测路线 |
| C8 | 在 matched 的 `L3` shared-region system compare 上，TransChromBP 仍稳定优于 official ChromBPNet | `shows` | `tutorial_L3_shared_region_closure_20260330.md` + source table | official=`0.33853/0.33989/0.69958`；Trans=`0.31319/0.31547/0.84016`；分类=`0.82295/0.83148/0.75248` vs `0.87861/0.87530/0.80991` | 只能写成 `shared-region system comparison`，不能升格成 architecture-only 归因 |

## 2. 只能弱化、不能强写的结论

| Claim ID | 原始强写法 | 现在允许的写法 | 为什么要降级 | 若想升格，最值钱的补证据 |
|---|---|---|---|---|
| W1 | 我们发现并证明了强 `Profile Shortcut` | 我们在早期探索中观察到 bias reliance 风险，但 clean revalidation 未复现强 shortcut | 四条 clean test 里 C 为 `0.00972 / 0.00997`（双 seed），仍属于轻度 | 不建议为此专门补旧线 |
| W2 | `stop-gradient` 单独修复了问题 | `stop-gradient` 是辅助性的 bias-isolation 设计，而 `deb2` 更像主效应 | A/C 对照不支持“sg 单独决定”，且 C 第二个 seed 已复现 | clean-matrix 主线已够；后续更值得补的是 tutorial L2 selector / stability follow-up |
| W3 | Transformer 更容易出现 bias reliance | 当前 clean matrix 不支持该表述，甚至更接近相反结论 | matched 的 `notf_sg1_deb0` 双 seed gap=`0.02682 / 0.01679`，两者都高于 TF 的 `0.00972 / 0.00997` | 如 reviewer 仍追问，优先补 shared-region compare，而不是继续堆 shortcut seed |
| W4 | conv-only 不会出现 shortcut | 当前证据与该表述直接矛盾 | conv-only 的 `notf_sg1_deb0` 仍是全矩阵最高风险格，只是第二个 seed 未越过 `0.02` 门槛 | 同上 |
| W5 | 显式 bias factorization 是当前 paper-facing recipe 的绝对必要组件 | 当前 single-seed no-bias 补证不支持该强写；更准确的写法是“在当前 recipe 下，去掉 bias branch 未显示清晰收益差异” | no-bias `epoch 26` held-out test=`0.31496 / 0.31715 / 0.84978`，与 corrected-B / center-pool 单 seed 和双 seed 结果同档 | 若 reviewer 坚持，可再补第二个 seed 或 matched shared-region no-bias；否则维持 supplementary 边界结果 |

## 3. 不应出现在主文中的结论

| Claim ID | 不应写的内容 | 问题 |
|---|---|---|
| N1 | “这是 Transformer 特有现象” | 当前 matched control 与此相反：`notf_sg1_deb0` 的双 seed gap 都高于 TF |
| N2 | “纯卷积架构中不会出现 Profile Shortcut” | 当前 conv-only 的 `notf_sg1_deb0` 恰好是全矩阵最高风险格 |
| N3 | “stop-gradient 单独修复了该问题” | 与 A/C clean matrix 不一致 |
| N4 | “我们已经明确证明了强 shortcut 机制” | 当前 clean revalidation 只支持轻度 bias reliance 说法 |

## 4. 章节落地建议

| 章节 | 应保留的主结论 | 应删除/改写的内容 |
|---|---|---|
| 摘要 | bias-safe Transformer、full/debiased 诊断、`deb2` 更关键、corrected B 稳定 | 强 shortcut、Transformer 特有、纯卷积不会出现 |
| 引言/贡献点 | 模型、诊断框架、受控复核、final readout、foundation side quest 负结果 | “首次发现强失效模式”式表述 |
| Shortcut 节 | 改名为 “Bias reliance 诊断与受控复核” | 旧的强机制解释、scale 轨迹主叙事、旧版 `2x2` 析因主结论 |
| 讨论 | 诊断框架的意义、当前证据边界、architecture sensitivity 未定 | “普适架构风险已被证明” |
| 结论 | 安全接入 Transformer + 可复用诊断 + 最终模型稳健；foundation 仅保留为受控 side quest | 强 shortcut 发现和单一修复论 |

## 5. 当前最值得推进的收口动作

1. 维持双语主稿对 `C7/C8` 的同步口径
   `L3 shared-region` 继续作为外部 system compare 主证据；foundation 线继续只写成 appendix / future work 的受控负结果，不回退到“Genos 单线负结果”。

2. 维持 `AUROC/AUPRC/F1` 的 supplementary 定位
   这组分类指标已经闭环并进入稿件，但仍不建议挤进主表，以免稀释 profile/count 主证据。

3. 用户提供型元数据独立处理，不反向影响主结论
   作者单位、致谢、最终对外 release 说明属于投稿前元数据，不应阻塞 claim-evidence 与主文论证的收口。
