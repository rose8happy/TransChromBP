# Profile Shortcut 假设复核计划（2026-03-25）

> 状态说明（2026-04-06）：本文件保留为 2026-03-25 当时的机制复核计划与执行记录，不再代表当前权威论文口径。后续 clean-matrix 双 seed 结果已由 `reports/analysis/profile_shortcut_revalidation_summary_20260326.md` 与 `reports/paper/paper_claim_evidence_matrix_20260326.md` 收口；当前稳定结论是“Transformer 特有 shortcut / 特有 gap 机制”不成立，风险更接近 supervision-driven 的 bias reliance。下文出现的 “Transformer 更容易触发/放大” 等措辞，应按当时待检验假设理解，不能再直接当成当前结论引用。

## 1. 为什么现在必须重做

当前 `reports/paper/transchrombp_paper_cn_v1.tex` 里对 “Profile Shortcut” 的叙述过强，至少有两层证据缺口：

1. `fixcw -> teacher_v2` 不是单变量变化。
   同时改了：
   - `bias_branch.profile_pool_factor: 0 -> 32`
   - `fusion.profile_bias_stop_gradient: false -> true`
   - `loss.debiased_profile_weight: 0.0 -> 2.0`

2. 最近 `v2fix` 读出头实验里被当成“只改 count head”的配置并不干净。
   `v2fix_center_pool.yaml`、`v2fix_attn_pool.yaml`、`v2fix_profile_refine.yaml`、`v2fix_pscale01.yaml`、`v2fix_groupnorm.yaml`
   都漏掉了 `local_tower` 和 `bias_branch.profile_pool_factor`，按当前代码会分别回退到：
   - `local_tower.n_dil_layers = 0`
   - `bias_branch.profile_pool_factor = 0`

3. 2026-03-25 代码复核又发现一个实现级 blocker：
   `heads.count_pool_mode` 在模型前向里没有真正生效，count 分支一直在走
   `encoded.mean(dim=1)`。因此历史上的 `B/F` 结果不能再被视为对
   “center pool / attention pool” 本身的干净验证，必须在修复实现后重跑
   corrected `B` / corrected `F`。

因此：

- `B=center pool` 目前仍然可以作为“结果不错的经验候选”，
- 但它不适合直接拿来当“Profile Shortcut 机制已经被证明”的实验底座。

当前最适合做机制复核的底座是 clean `teacher_v2` scaffold：

- 保留 `local_tower`
- 保留 `profile_pool_factor = 32`
- 只针对 `stop-gradient / debiased supervision / Transformer on-off` 做受控改动

## 2. 这次要回答的三个问题

### H1：`debiased_profile_weight=2.0` 单独是否已经足够抑制 shortcut？

对照：

- 现有参考：`ablation_tf_20260318_full_s42`
  条件 = `TF + sg=true + deb2`
- 新跑：`TF + sg=false + deb2`

判读：

- 如果新跑仍出现明显 full/debiased gap，说明 `debiased_profile_weight` 不足以单独消除 shortcut，`stop-gradient` 仍是关键变量。
- 如果 gap 仍接近 0，则说明先前把收益主要归因给 `stop-gradient` 的写法过强。

### H2：`stop-gradient` 单独是否已经足够抑制 shortcut？

对照：

- 现有参考：`ablation_tf_20260318_full_s42`
  条件 = `TF + sg=true + deb2`
- 新跑：`TF + sg=true + deb0`

判读：

- 如果 `sg=true + deb0` 仍能保持极小 gap，说明 `stop-gradient` 本身足以压住 shortcut，`debiased_profile_weight=2.0` 更多是在帮助训练质量。
- 如果 gap 明显放大，则当前 paper 应改成“stop-gradient + explicit debiased supervision 共同修复”，而不是只强调前者。

### H3：conv-only 在当前 scaffold 下会不会也出现类似 shortcut？

对照：

- 现有参考：`ablation_tf_20260318_noTF_s42`
  条件 = `noTF + sg=true + deb2`
- 新跑：`noTF + sg=false + deb2`

判读：

- 如果 noTF 在关掉 `sg` 后仍几乎没有 gap，则“Transformer 更容易触发 shortcut”这一说法得到支持。
- 如果 noTF 也出现明显 gap，则论文里“纯卷积不会出现 shortcut”的强断言必须删除，最多只能保留为“Transformer 下更容易放大”。

## 3. 第一批三卡矩阵

统一原则：

- scaffold 固定为 `teacher_v2` 语义
- seed 先固定为 `42`
- bias checkpoint 统一复用 `ablation_tf_20260318_shared_bias/best.pt`
- 训练模板统一复用 `train_ablation_v2_main.yaml`
  这样可以直接和 `ablation_tf_20260318_full_s42 / noTF_s42` 做 matched 对照

### 3.1 新跑的三条

| run 角色 | 结构/损失 | 目的 | 计划 GPU |
|---|---|---|---|
| `tf_sg0_deb2` | `TF + sg=false + deb2` | 检验 `deb2` 是否足以替代 `sg` | 6000 GPU0 |
| `tf_sg1_deb0` | `TF + sg=true + deb0` | 检验 `sg` 是否足以替代 `deb2` | 6000 GPU1 |
| `notf_sg0_deb2` | `noTF + sg=false + deb2` | 检验 conv-only 是否也会出现 gap | 6002 GPU0 |

### 3.2 直接复用的现有参考

| 参考 run | 条件 | 用途 |
|---|---|---|
| `ablation_tf_20260318_full_s42` | `TF + sg=true + deb2` | clean safe baseline |
| `ablation_tf_20260318_noTF_s42` | `noTF + sg=true + deb2` | clean noTF safe baseline |

## 4. 关键指标与判读门槛

主看 5 个量：

1. `peak.profile_target_jsd_full_mean`
2. `peak.profile_target_jsd_debiased_mean`
3. `full - debiased gap`
4. `effective_profile_scale`
5. `profile_bias_rms_over_signal_rms`

建议判读：

- `gap <= 0.005`
  视为“没有看到实质 shortcut”
- `0.005 < gap <= 0.02`
  视为“有轻度 bias reliance，需要结合 scale/rms 再判断”
- `gap > 0.02`
  视为“shortcut 明显存在，不能再写成修复已被稳健证明”

辅助看：

- `count_pearson_full`
- `count_pearson_debiased`
- `best_epoch`

如果 `gap` 变大同时 `effective_profile_scale` 和 `profile_bias_rms_over_signal_rms` 也同步升高，才更符合“shortcut/依赖 bias profile”而不是普通训练波动。

## 5. 为什么第一批不直接用当前 B=center pool

原因不是它分数不够好，而是它现在不够“干净”：

- `v2fix_center_pool.yaml` 漏掉了 `local_tower`
- 同时漏掉了 `profile_pool_factor`

这会把“readout head 改动”和“局部塔/偏差分辨率变化”混在一起。
所以它适合继续当性能候选，不适合当机制判定底座。

如果第一批 clean matrix 明确支持某个方向，第二批再回到 corrected center-pool scaffold 做 paper-facing 复验更合理。
当前这一步已经从“建议”升级为“必须”：先修复 `count_pool_mode` 实现，再跑 corrected `B`。

## 5.1 2026-03-25 当日执行更新

- `count_pool_mode` 已在 `src/transchrombp/models/transchrombp.py` 真正接入前向：
  - `full` = 全长 mean pool
  - `center` = 仅对中心 `output_len` 窗口 mean pool
  - `attention` = learnable attention pooling
- 6002 最小 probe 已直接验证该逻辑不是“只读到配置值”：
  - 同一份“中心 1000bp = 1、两侧 = 0”的假特征上
  - `full` pooled mean = `0.4730`
  - `center` pooled mean = `1.0`
  - `attention` 路径可正常前向
- 远端同步时还发现两个实现级不同步问题：
  - `src/transchrombp/models/genos_adapter.py` 在远端缺失，新的 `transchrombp.py` import 会直接失败
  - 6002 的 `src/transchrombp/data/real_data.py` 还是旧签名，不接受 `genos_cache_dir`
- 上述文件已同步到 6000/6002。corrected `B` 第一次启动因旧 `real_data.py` 秒挂，不记作实验失败；已在同步后重启。
- corrected `B` 当前正式 run：
  - `profile_shortcut_20260325_tf_center_sg1_deb2_s42_6002_retry1`
  - 机器/GPU：6002 / RTX 3080
  - 语义：`teacher_v2 scaffold + center pool + sg=true + deb2`
  - 日志：`/home/zhengwei/bylw_atac/logs/profile_shortcut_20260325_tf_center_sg1_deb2_s42_6002_retry1.log`
  - 启动时间：2026-03-25 12:51 CST
  - 当前状态：已稳定进入 epoch 1
  - 粗略结束时间：2026-03-26 03:30 CST 左右（若 early-stop 会更早）

## 6. 第二批触发条件

只有第一批出现下列情况之一时，才建议追加：

1. `tf_sg0_deb2` 与 `tf_sg1_deb0` 都接近 safe baseline，分不出主因
2. `notf_sg0_deb2` 结果介于“明显 gap”和“几乎无 gap”之间
3. 三条 run 的 best epoch/曲线形态差异过大，怀疑单 seed 偶然性

那时再补 `seed=1234` 的 matched repeat。

## 6.1 2026-03-26 追加的最高优先级 matched control

在 clean matrix 四条 test 已经补齐之后，当前最值得追加的单条实验不再是多 seed 扩全，
而是：

- `noTF + sg=true + deb0`

原因：

1. 它是与 `C = TF + sg=true + deb0` 最 matched 的架构对照。
2. 它直接回答：当显式 debiased profile supervision 被拿掉时，轻度 bias reliance
   是否也会出现在 conv-only scaffold 中。
3. 这条结果会直接决定论文讨论里能否继续保留任何“architecture sensitivity”
   的保守措辞。

新的判读规则：

- 如果 `notf_sg1_deb0` 仍保持 `gap <= 0.005`
  则当前证据会更接近“弱监督下 Transformer 比 conv-only 更容易放大 bias reliance”，
  但正文依然最多写成 `suggests`。
- 如果 `notf_sg1_deb0` 也达到和 `C` 接近的 gap
  则论文里所有与 “Transformer 更敏感” 相关的表述都应删除，只保留
  “缺少 debiased supervision 会增加 bias reliance”。
- 如果结果落在中间区间（如 `0.005-0.02`）
  则把这条作为“当前无法确认 architecture-specific”的证据，而不是继续扩全矩阵。

计划执行：

| run 角色 | 结构/损失 | 目的 | 优先 GPU |
|---|---|---|---|
| `notf_sg1_deb0` | `noTF + sg=true + deb0` | 判断轻度 bias reliance 是否也会在 conv-only 下出现 | 6002 GPU0 |

执行更新（2026-03-26 01:15 CST）：

- 已在 6002 / RTX 3080 后台启动：
  - run name: `profile_shortcut_20260326_notf_sg1_deb0_s42_6002`
  - launcher log: `/home/zhengwei/bylw_atac/logs/profile_shortcut_20260326_notf_sg1_deb0_s42_6002.launch.log`
  - train log: `/home/zhengwei/bylw_atac/logs/profile_shortcut_20260326_notf_sg1_deb0_s42_6002.log`
  - 关键配置差异：`sequence_encoder.enabled=false`、`fusion.profile_bias_stop_gradient=true`、`loss.debiased_profile_weight=0.0`
- 启动后已确认不是空跑：
  - 3080 占用约 `1385 MiB / 97%`
  - 主日志已推进到 `epoch 1 step 1260/14761`
- 按当前吞吐粗估：
  - 每 epoch 约 `6.5` 分钟
  - 若跑满 `30 epoch`，约在 `2026-03-26 04:30 CST` 左右结束
  - 若较早触发 early-stop，会明显早于该时间

## 6.2 2026-03-26 并行启动的 sidecar repeat

在用户明确同意“不要让 A6000 闲着”的前提下，当前唯一值得并行补的 sidecar
不是扩全矩阵，而是：

- `tf_sg1_deb0` 的第二个 seed

原因：

1. 这条不会改变主判断路径；主判断仍由 `notf_sg1_deb0` 决定。
2. 但如果 `notf_sg1_deb0` 最终支持“可能存在弱的 Transformer sensitivity”，
   那么 `tf_sg1_deb0` 的第二个 seed 可以立刻告诉我们 C 线的轻度 gap
   是稳定现象还是单 seed 偶然。
3. 相比再开新的结构线，这条 repeat 的信息利用率更高，也更直接服务于论文措辞强度。

执行更新（2026-03-26 01:25 CST）：

- 已在 6000 / A6000 GPU0 后台启动：
  - run name: `profile_shortcut_20260326_tf_sg1_deb0_s1234_6000`
  - launcher log: `/data1/zhoujiazhen/bylw_atac/logs/profile_shortcut_20260326_tf_sg1_deb0_s1234_6000.launch.log`
  - train log: `/data1/zhoujiazhen/bylw_atac/logs/profile_shortcut_20260326_tf_sg1_deb0_s1234_6000.log`
  - 关键配置差异：`sequence_encoder.enabled=true`、`profile_bias_stop_gradient=true`、`debiased_profile_weight=0.0`、`seed=1234`
- 启动后已确认不是空跑：
  - A6000 GPU0 占用约 `14701 MiB / 99%`
  - `nvidia-smi pmon` 显示本进程 `PID 16297` 与另一位用户的双卡任务共享 GPU0
  - 主日志已推进到 `epoch 1 step 120/14761`
- 由于当前 GPU0 处于共享状态，吞吐明显慢于无 contention 情况：
  - 当前粗估每 epoch 约 `30-35` 分钟
  - 若跑满 `30 epoch`，约在 `2026-03-26 16:00-18:00 CST` 结束
  - 若中途 early-stop，会更早收口

## 7. 当时的临时论文改写建议（现已被后续结论覆盖）

以下内容是 2026-03-25 在新证据出来前的临时降级建议，仅保留作历史参考：

- 不要写“纯卷积架构中不会出现”
- 改写成“现有证据表明，Transformer scaffold 下更容易暴露出显著的 full/debiased gap；这一现象是否为 Transformer 特有，仍需 matched conv-only 对照进一步确认”

以及：

- 不要把 `fixcw -> v2` 的收益直接写成 `stop-gradient` 单独造成
- 改写成“v2 方案在 stop-gradient、bias profile resolution bottleneck、debiased profile supervision 共同作用下显著缩小了 gap；各因素的独立贡献仍需进一步拆分”
