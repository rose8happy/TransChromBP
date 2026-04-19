# Foundation Cache Contract Snapshot（2026-04-06）

## 一句话结论

`autonomy/20260406-structure` 不是 live experiment 线，而是一条仍保留独有提交的 infra snapshot branch；当前不适合直接删 branch 或卸载 worktree，因为它承载的 foundation contract / alignment regression 还没有别的 canonical 报告把这两条提交解释清楚。

## 1. 当前 branch 上固定了什么

- `24649aa` `wip: snapshot foundation cache contract`
  - 关键内容：
    - `vendor/transchrombp/transchrombp/utils/foundation_contract.py`
    - `tests/test_foundation_contract.py`
    - `build_foundation_cache.py` / `evaluate_checkpoint.py` / `train_ddp.py` 上的 contract 接口调整
- `764b8c0` `wip: add foundation cache alignment regression`
  - 关键内容：
    - `tests/test_foundation_cache_alignment.py`

这两条提交共同表达的是：foundation cache 构建、训练、held-out 评估之间的 contract 曾被单独抽象并补了回归守卫，但目前尚未回流成一份更正式的 infra closeout。

## 2. 为什么这条 branch 现在还要保留

- 它相对 `master` 仍有独有提交，不是“仅靠 snapshot tag 就能完全解释”的空壳 branch。
- registry 需要一份能回答“这条 infra branch 到底保存了什么、为什么还挂着”的 canonical 文档。
- 在没有更强 closeout 前，直接卸载 worktree 只会把问题从“目录杂乱”变成“档案解释断层”。

## 3. 当前可接受的读取方式

- 状态索引：`docs/experiments/registry.md`
- 历史审查背景：`reports/project_plan_code_review_20260405.md`
- 冻结 tag：`snapshot/autonomy-20260406-structure/20260406`

## 4. 下一步建议

如果要继续收这条 infra 线，先做二选一：

1. 把 `24649aa` / `764b8c0` 的关键实现或测试择净回流 `master`，再卸载 worktree。
2. 或者补一份更完整的 infra closeout，明确哪些文件只需靠 snapshot 保留、哪些必须继续可浏览。

在此之前，继续保留这条 mounted worktree 是合理的。
