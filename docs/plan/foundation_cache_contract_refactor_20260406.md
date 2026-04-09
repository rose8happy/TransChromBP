# foundation cache contract refactor（2026-04-06）

## 目标

把 cached foundation feature 的“合同层”从 `train_ddp.py` / `evaluate_checkpoint.py` 里的重复字符串逻辑抽成共享 helper，先统一三件事：

1. `foundation_model.feature_name` / `feature_layout` / `feature_tokens` / `hidden_size` 的解释；
2. `data.foundation_cache_features` 与 `foundation_cache_dir` 的校验；
3. `foundation_<feature_name>` batch key 和 `layer_07__bins4_mean` 这种 feature 名称的解析规则。

这一步故意不改训练数学，也不重写 launcher 或 cache builder。

## 本轮变更

- 新增 `vendor/transchrombp/transchrombp/utils/foundation_contract.py`
  - `resolve_foundation_contract()`
  - `infer_required_foundation_cache_features()`
  - `infer_foundation_cache_build_request()`
  - `validate_foundation_cache_config()`
  - `validate_foundation_cache_build_request()`
  - `parse_foundation_feature_name()`
  - `foundation_batch_key()`
- `train_ddp.py`
  - `extract_foundation_feature_kwargs()` 改为通过共享 helper 解释 feature contract
  - 训练前的 foundation cache 校验改为复用 `validate_foundation_cache_config()`
- `evaluate_checkpoint.py`
  - held-out 评估前复用同一份 foundation cache 校验
- `build_foundation_cache.py`
  - 新增可选 `--model_config`
  - 在真正开跑前校验当前 `layers / feature_types` 是否覆盖 model config 所需 feature
  - `split / max_records / region_source / seed` 解析改走共享 helper
- `run_ntv2_residual_short10.sh`
  - 从 `MODEL_CONFIG` 自动推导 `CACHE_LAYERS` / `CACHE_FEATURE_TYPES`
  - 若手工 override 的 cache request 不覆盖必需 feature，直接在 preflight 阶段报错
  - 不再把 `7,14 + global_mean bins4_mean` 写死在脚本里
- `evaluate_checkpoint.py`
  - held-out dataset 的 `split / max_records / region_source / seed` 解析也改走同一份 helper
- 新增 `tests/test_foundation_contract.py`
  - 覆盖现有 NT v2 residual / cross-attention config
  - 覆盖一个 synthetic bins16 residual config
  - 覆盖缺少 `foundation_cache_features` 时的报错
  - 覆盖 cache build request 的自动推导与校验

## 快速验证

- `PYTHONPATH=$PWD/vendor/transchrombp python3 -m unittest tests.test_foundation_contract -v`
- `PYTHONPATH=$PWD/vendor/transchrombp python3 -m py_compile vendor/transchrombp/transchrombp/utils/foundation_contract.py vendor/transchrombp/transchrombp/training/train_ddp.py vendor/transchrombp/transchrombp/evaluation/evaluate_checkpoint.py`
- `PYTHONPATH=$PWD/vendor/transchrombp python3 -m py_compile vendor/transchrombp/transchrombp/scripts/build_foundation_cache.py`
- `bash -n vendor/transchrombp/transchrombp/scripts/run_ntv2_residual_short10.sh`

## 当前刻意没做的事

- 还没有尝试泛化到 `genos_cached` 或 `caduceus`。

## 下一步

如果继续沿这个分支整理，最值得做的是第二阶段：

1. 给 launcher 增加 manifest 级 preflight，提前报错“已有 cache 与当前 `model_config` / train semantics 不一致”；
2. 再把同类合同层推广到 `genos_cached` 与其他 foundation launcher，减少同样的重复字符串逻辑。
