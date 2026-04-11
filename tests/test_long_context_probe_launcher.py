import re
from pathlib import Path

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ARCHIVE_ROOT = PROJECT_ROOT / "vendor/transchrombp"
BASE_DATA_CONFIG = ARCHIVE_ROOT / "configs/data/data_tutorial_canonical_v1.yaml"
BASE_MODEL_CONFIG = ARCHIVE_ROOT / "configs/model/transchrombp_teacher_v2_center_pool.yaml"
BASE_TRAIN_CONFIG = ARCHIVE_ROOT / "configs/train/train_tutorial_teacher_v2_readout_short10.yaml"
DATA_CONFIG = ARCHIVE_ROOT / "configs/data/data_tutorial_canonical_v1_longctx4096.yaml"
MODEL_CONFIG = ARCHIVE_ROOT / "configs/model/transchrombp_teacher_v2_center_pool_longctx4096.yaml"
TRAIN_CONFIG = ARCHIVE_ROOT / "configs/train/train_tutorial_teacher_v2_longctx4096_short10.yaml"
LAUNCHER = ARCHIVE_ROOT / "scripts/run_teacher_v2_longctx4096_probe.sh"

FORBIDDEN_DECODER_KEYS = {
    "decoder",
    "decoders",
    "profile_decoder",
    "count_decoder",
    "multi_scale_decoder",
}


def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def flatten_dict(tree: dict, prefix: str = "") -> dict[str, object]:
    flat: dict[str, object] = {}
    for key, value in tree.items():
        path = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            flat.update(flatten_dict(value, path))
        else:
            flat[path] = value
    return flat


def diff_paths(lhs: dict, rhs: dict) -> set[str]:
    lhs_flat = flatten_dict(lhs)
    rhs_flat = flatten_dict(rhs)
    changed = {key for key in lhs_flat.keys() | rhs_flat.keys() if lhs_flat.get(key) != rhs_flat.get(key)}
    return changed


def extract_train_invocation_block(launcher_text: str) -> str:
    match = re.search(
        r'exec python -m transchrombp\.training\.train_ddp \\\n(?:.+\\\n)+\s+"\$@"',
        launcher_text,
    )
    assert match, "missing train_ddp exec block"
    return match.group(0)


def test_data_config_only_changes_long_context_input_len() -> None:
    base_data_cfg = load_yaml(BASE_DATA_CONFIG)
    data_cfg = load_yaml(DATA_CONFIG)

    assert diff_paths(base_data_cfg, data_cfg) == {"window.input_len"}
    assert data_cfg["window"]["input_len"] == 4096


def test_model_config_preserves_corrected_b_center_pool_semantics() -> None:
    base_model_cfg = load_yaml(BASE_MODEL_CONFIG)
    model_cfg = load_yaml(MODEL_CONFIG)

    assert diff_paths(base_model_cfg, model_cfg) == {"model_name", "sequence_encoder.max_len"}
    assert model_cfg["model_name"] == "transchrombp_teacher_v2_center_pool_longctx4096"
    assert model_cfg["sequence_encoder"]["max_len"] == 4096
    assert FORBIDDEN_DECODER_KEYS.isdisjoint(model_cfg.keys())


def test_train_config_only_changes_long_context_probe_knobs() -> None:
    base_train_cfg = load_yaml(BASE_TRAIN_CONFIG)
    train_cfg = load_yaml(TRAIN_CONFIG)

    assert diff_paths(base_train_cfg, train_cfg) == {
        "data.config_path",
        "data.input_len",
        "data.max_seq_len",
        "logging.run_name",
    }
    assert train_cfg["data"]["config_path"] == "configs/data/data_tutorial_canonical_v1_longctx4096.yaml"
    assert train_cfg["data"]["input_len"] == 4096
    assert train_cfg["data"]["max_seq_len"] == 4096
    assert train_cfg["logging"]["run_name"] == "teacher_v2_center_pool_longctx4096_short10"


def test_launcher_uses_isolated_output_dir_and_binds_probe_configs() -> None:
    launcher_text = LAUNCHER.read_text(encoding="utf-8")
    train_invocation = extract_train_invocation_block(launcher_text)

    assert 'DATA_CONFIG="${DATA_CONFIG:-${ROOT_DIR}/configs/data/data_tutorial_canonical_v1_longctx4096.yaml}"' in launcher_text
    assert 'MODEL_CONFIG="${MODEL_CONFIG:-${ROOT_DIR}/configs/model/transchrombp_teacher_v2_center_pool_longctx4096.yaml}"' in launcher_text
    assert 'TRAIN_CONFIG="${TRAIN_CONFIG:-${ROOT_DIR}/configs/train/train_tutorial_teacher_v2_longctx4096_short10.yaml}"' in launcher_text
    assert 'OUTPUT_DIR="${OUTPUT_DIR:-${ROOT_DIR}/outputs}"' in launcher_text
    assert '--data-config "${DATA_CONFIG}"' in train_invocation
    assert '--model-config "${MODEL_CONFIG}"' in train_invocation
    assert '--train-config "${TRAIN_CONFIG}"' in train_invocation
    assert '--output-dir "${OUTPUT_DIR}"' in train_invocation
