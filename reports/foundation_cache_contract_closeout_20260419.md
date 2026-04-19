# Foundation Cache Contract Closeout（2026-04-19）

## 一句话结论

`autonomy/20260406-structure` 的两条 branch-only 提交已经有了足够明确的 canonical closeout：`24649aa` 负责 foundation cache contract refactor，`764b8c0` 负责 alignment regression guard。现在这条 infra line 可以从“依赖 mounted worktree 浏览”转成“通过 closeout + snapshot tag 恢复”。

## 1. 这条 infra line 固定了什么

- branch: `autonomy/20260406-structure`
- snapshot tag: `snapshot/autonomy-20260406-structure/20260406`
- closeout tag: `closeout/foundation_cache_contract/20260419`

### `24649aa` `wip: snapshot foundation cache contract`

核心内容：

- `vendor/transchrombp/transchrombp/utils/foundation_contract.py`
- `tests/test_foundation_contract.py`
- `vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py`
- `vendor/transchrombp/transchrombp/scripts/build_foundation_cache.py`
- `vendor/transchrombp/transchrombp/scripts/run_ntv2_residual_short10.sh`
- `vendor/transchrombp/transchrombp/training/train_ddp.py`
- `docs/plan/foundation_cache_contract_refactor_20260406.md`

它解决的问题是把 foundation cache 的合同层从 `train_ddp.py` / `evaluate_checkpoint.py` 的重复字符串逻辑里抽出来，统一 feature contract、cache request 与 config validation。

### `764b8c0` `wip: add foundation cache alignment regression`

核心内容：

- `tests/test_foundation_cache_alignment.py`

它把随后发现的 alignment regression 固化成显式回归守卫，避免 cache builder、训练和 held-out evaluator 的 split/seed/manifest 语义再漂移。

## 2. 与现有 master-side 报告的关系

- `reports/project_plan_code_review_20260405.md` 解释了这条问题线最初是如何暴露出来的。
- `reports/foundation_cache_contract_snapshot_20260406.md` 是上一轮的保留说明，回答“为什么还不能直接删掉这条 branch/worktree”。
- 本文负责最终 closeout：明确两条 branch-only 提交分别保存了什么、为什么现在可以卸载 mounted worktree、以及如何恢复现场。

## 3. 为什么没有把这条 branch 的代码直接回流 master

- 这条线本质是 infra snapshot，不是当前 live experiment 主线。
- branch 上的内容仍包含当时的局部实现选择与测试切片；它们适合作为恢复入口保留，但不适合在没有显式结构性任务的前提下直接混回当前 canonical trunk。
- 当前更高价值的动作是把解释层补齐，让维护者不再需要依赖 mounted worktree 才能理解这条 branch。

## 4. 当前已经满足的 gate

- branch clean
- snapshot tag 已在
- formal closeout 已在
- registry / canonical report 已能解释这条 line 的职责和恢复入口

因此现在可以安全地：

- 卸载 mounted worktree
- 保留 branch `autonomy/20260406-structure`
- 保留 snapshot / closeout tags

## 5. 恢复方式

如需重新查看这条 infra line：

1. 先读：
   - `reports/project_plan_code_review_20260405.md`
   - `reports/foundation_cache_contract_snapshot_20260406.md`
   - 本文
2. 再从以下入口恢复代码现场：
   - branch `autonomy/20260406-structure`
   - `snapshot/autonomy-20260406-structure/20260406`
   - `closeout/foundation_cache_contract/20260419`

当前不建议继续为了“本地可点开浏览”而保留 mounted worktree。
