#!/usr/bin/env bash
# NT v2 cached residual-head pilot launcher.
# Sequence:
#   1. build dataset-aligned NT v2 cache (train + valid)
#   2. optionally run unified valid probe against baseline checkpoint
#   3. run corrected-B + NT v2 residual-head short10 training
#   4. (default) build held-out test cache with matched nonpeak_ratio and run
#      evaluate_checkpoint to produce test-full gate metrics

set -euo pipefail

if [ -z "${BASH_VERSION:-}" ]; then
    echo "[error] Please run this script with bash." >&2
    exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TRAIN_GPU_IDS="${TRAIN_GPU_IDS:-${GPU_IDS:-${1:-${GPU_ID:-0,1}}}}"
CACHE_GPU_IDS="${CACHE_GPU_IDS:-${TRAIN_GPU_IDS}}"
CACHE_GPU_ID="${CACHE_GPU_ID:-${CACHE_GPU_IDS%%,*}}"
NPROC_PER_NODE="${NPROC_PER_NODE:-2}"
CACHE_NPROC_PER_NODE="${CACHE_NPROC_PER_NODE:-${NPROC_PER_NODE}}"
HELDOUT_CACHE_GPU_IDS="${HELDOUT_CACHE_GPU_IDS:-${CACHE_GPU_IDS}}"
HELDOUT_CACHE_GPU_ID="${HELDOUT_CACHE_GPU_ID:-${HELDOUT_CACHE_GPU_IDS%%,*}}"
HELDOUT_CACHE_NPROC_PER_NODE="${HELDOUT_CACHE_NPROC_PER_NODE:-${CACHE_NPROC_PER_NODE}}"
MASTER_ADDR="${MASTER_ADDR:-127.0.0.1}"
MASTER_PORT_BASE="${MASTER_PORT_BASE:-29820}"
VENV_DIR="${VENV_DIR:-/data1/zhoujiazhen/bylw_atac/.venvs/genos-1.2b}"
OUTPUT_BASE="${OUTPUT_BASE:-${ROOT_DIR}/outputs}"
FOUNDATION_CACHE_DIR="${FOUNDATION_CACHE_DIR:-${OUTPUT_BASE}/foundation_cache/ntv2_tutorial_canonical_v1}"
NT_MODEL_DIR="${NT_MODEL_DIR:-/data1/zhoujiazhen/bylw_atac/foundation_models/nucleotide_transformer/nt-v2-500m-multi-species}"
TRAIN_CONFIG="${TRAIN_CONFIG:-${ROOT_DIR}/configs/train/train_tutorial_foundation_short10.yaml}"
DATA_CONFIG="${DATA_CONFIG:-${ROOT_DIR}/configs/data/data_tutorial_canonical_v1.yaml}"
MODEL_CONFIG="${MODEL_CONFIG:-${ROOT_DIR}/configs/model/transchrombp_teacher_v2_center_pool_ntv2_residual.yaml}"
BASELINE_CHECKPOINT="${BASELINE_CHECKPOINT:-}"
BASELINE_MODEL_CONFIG="${BASELINE_MODEL_CONFIG:-${ROOT_DIR}/configs/model/transchrombp_teacher_v2_center_pool.yaml}"
PYTHON_BIN="${PYTHON_BIN:-python}"
CACHE_SPLITS="${CACHE_SPLITS:-train valid}"
CACHE_FEATURE_TYPES="${CACHE_FEATURE_TYPES:-global_mean bins4_mean}"
CACHE_RECORD_SPLITS="${CACHE_RECORD_SPLITS:-valid}"
CACHE_LAYERS="${CACHE_LAYERS:-7,14}"
CACHE_BATCH_SIZE="${CACHE_BATCH_SIZE:-8}"
CACHE_DTYPE="${CACHE_DTYPE:-float16}"
RUN_HELDOUT_GATE="${RUN_HELDOUT_GATE:-true}"
HELDOUT_SPLIT="${HELDOUT_SPLIT:-test}"
HELDOUT_NONPEAK_RATIO="${HELDOUT_NONPEAK_RATIO:-1.0}"
HELDOUT_REGION_SOURCE="${HELDOUT_REGION_SOURCE:-both}"
HELDOUT_MAX_REGIONS="${HELDOUT_MAX_REGIONS:-0}"
HELDOUT_OUTPUT="${HELDOUT_OUTPUT:-}"
EVAL_GPU_ID="${EVAL_GPU_ID:-${TRAIN_GPU_IDS%%,*}}"
TRAIN_DRY_RUN_STEPS="${TRAIN_DRY_RUN_STEPS:-0}"

if [ -z "${BATCH_SIZE_PER_GPU:-}" ]; then
    if [ "${NPROC_PER_NODE}" -gt 1 ]; then
        BATCH_SIZE_PER_GPU=10
    else
        BATCH_SIZE_PER_GPU=20
    fi
fi
if [ -z "${NUM_WORKERS:-}" ]; then
    if [ "${NPROC_PER_NODE}" -gt 1 ]; then
        NUM_WORKERS=4
    else
        NUM_WORKERS=2
    fi
fi
if [ -z "${DDP_FIND_UNUSED_PARAMETERS:-}" ]; then
    DDP_FIND_UNUSED_PARAMETERS=true
fi
RUN_HELDOUT_GATE="$(printf '%s' "${RUN_HELDOUT_GATE}" | tr '[:upper:]' '[:lower:]')"
if [ "${RUN_HELDOUT_GATE}" != "true" ] && [ "${RUN_HELDOUT_GATE}" != "false" ]; then
    echo "[error] RUN_HELDOUT_GATE must be true or false, got: ${RUN_HELDOUT_GATE}" >&2
    exit 1
fi

read -r -a CACHE_SPLIT_ARGS <<< "${CACHE_SPLITS}"
read -r -a CACHE_FEATURE_TYPE_ARGS <<< "${CACHE_FEATURE_TYPES}"
read -r -a CACHE_RECORD_SPLIT_ARGS <<< "${CACHE_RECORD_SPLITS}"
if [ "${#CACHE_SPLIT_ARGS[@]}" -eq 0 ]; then
    echo "[error] CACHE_SPLITS is empty" >&2
    exit 1
fi
if [ "${#CACHE_FEATURE_TYPE_ARGS[@]}" -eq 0 ]; then
    echo "[error] CACHE_FEATURE_TYPES is empty" >&2
    exit 1
fi

export PATH="${VENV_DIR}/bin:${PATH:-}"
if [ -d "${ROOT_DIR}/src/transchrombp" ]; then
    PACKAGE_IMPORT_ROOT="${ROOT_DIR}/src"
elif [ -d "${ROOT_DIR}/../transchrombp" ]; then
    PACKAGE_IMPORT_ROOT="$(cd "${ROOT_DIR}/.." && pwd)"
else
    echo "[error] Could not locate transchrombp package import root from ${ROOT_DIR}" >&2
    exit 1
fi
export PYTHONPATH="${PACKAGE_IMPORT_ROOT}:${PYTHONPATH:-}"

require_file() {
    local path="$1"
    local label="$2"
    if [ ! -f "${path}" ]; then
        echo "[error] Missing ${label}: ${path}" >&2
        exit 1
    fi
}

require_dir() {
    local path="$1"
    local label="$2"
    if [ ! -d "${path}" ]; then
        echo "[error] Missing ${label}: ${path}" >&2
        exit 1
    fi
}

count_csv_items() {
    local csv="$1"
    if [ -z "${csv}" ]; then
        echo 0
        return 0
    fi
    awk -F',' '{print NF}' <<<"${csv}"
}

launch_train() {
    local master_port="$1"
    shift
    export CUDA_VISIBLE_DEVICES="${TRAIN_GPU_IDS}"
    if [ "${NPROC_PER_NODE}" -gt 1 ]; then
        torchrun \
            --nnodes=1 \
            --nproc_per_node="${NPROC_PER_NODE}" \
            --master_addr="${MASTER_ADDR}" \
            --master_port="${master_port}" \
            -m transchrombp.training.train_ddp \
            "$@"
    else
        "${PYTHON_BIN}" -m transchrombp.training.train_ddp "$@"
    fi
}

require_file "${TRAIN_CONFIG}" "train config"
require_file "${DATA_CONFIG}" "data config"
require_file "${MODEL_CONFIG}" "model config"
require_dir "${NT_MODEL_DIR}" "NT v2 model directory"

TRAIN_GPU_COUNT="$(count_csv_items "${TRAIN_GPU_IDS}")"
if [ "${NPROC_PER_NODE}" -gt "${TRAIN_GPU_COUNT}" ]; then
    echo "[error] NPROC_PER_NODE=${NPROC_PER_NODE} exceeds visible TRAIN_GPU_IDS=${TRAIN_GPU_IDS}" >&2
    exit 1
fi
CACHE_GPU_COUNT="$(count_csv_items "${CACHE_GPU_IDS}")"
if [ "${CACHE_NPROC_PER_NODE}" -gt "${CACHE_GPU_COUNT}" ]; then
    echo "[error] CACHE_NPROC_PER_NODE=${CACHE_NPROC_PER_NODE} exceeds visible CACHE_GPU_IDS=${CACHE_GPU_IDS}" >&2
    exit 1
fi
HELDOUT_CACHE_GPU_COUNT="$(count_csv_items "${HELDOUT_CACHE_GPU_IDS}")"
if [ "${HELDOUT_CACHE_NPROC_PER_NODE}" -gt "${HELDOUT_CACHE_GPU_COUNT}" ]; then
    echo "[error] HELDOUT_CACHE_NPROC_PER_NODE=${HELDOUT_CACHE_NPROC_PER_NODE} exceeds visible HELDOUT_CACHE_GPU_IDS=${HELDOUT_CACHE_GPU_IDS}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_BASE}" "${FOUNDATION_CACHE_DIR}"

RUNTIME_DIR="${OUTPUT_BASE}/runtime/ntv2_residual_short10"
mkdir -p "${RUNTIME_DIR}"
EFFECTIVE_TRAIN_CONFIG="${RUNTIME_DIR}/train_tutorial_foundation_short10_runtime.yaml"
"${PYTHON_BIN}" - "${TRAIN_CONFIG}" "${EFFECTIVE_TRAIN_CONFIG}" "${BATCH_SIZE_PER_GPU}" "${NUM_WORKERS}" "${DDP_FIND_UNUSED_PARAMETERS}" <<'PY'
import sys
from pathlib import Path
import yaml

src, dst, batch_size, num_workers, find_unused = sys.argv[1:6]
with open(src, "r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
cfg.setdefault("data", {})["batch_size_per_gpu"] = int(batch_size)
cfg["data"]["num_workers"] = int(num_workers)
cfg.setdefault("trainer", {})["find_unused_parameters"] = str(find_unused).lower() == "true"
Path(dst).write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
PY

echo "=============================================="
echo "NT v2 Residual Short10"
echo "ROOT_DIR: ${ROOT_DIR}"
echo "CACHE_GPU_IDS: ${CACHE_GPU_IDS}"
echo "CACHE_GPU_ID: ${CACHE_GPU_ID} (used only when CACHE_NPROC_PER_NODE=1)"
echo "CACHE_NPROC_PER_NODE: ${CACHE_NPROC_PER_NODE}"
echo "TRAIN_GPU_IDS: ${TRAIN_GPU_IDS} (preferred dual-GPU training)"
echo "NPROC_PER_NODE: ${NPROC_PER_NODE}"
echo "MASTER_PORT_BASE: ${MASTER_PORT_BASE}"
echo "Cache dir: ${FOUNDATION_CACHE_DIR}"
echo "NT model dir: ${NT_MODEL_DIR}"
echo "Train config: ${TRAIN_CONFIG}"
echo "Effective train config: ${EFFECTIVE_TRAIN_CONFIG}"
echo "Model config: ${MODEL_CONFIG}"
echo "Data config: ${DATA_CONFIG}"
echo "BATCH_SIZE_PER_GPU: ${BATCH_SIZE_PER_GPU}"
echo "NUM_WORKERS: ${NUM_WORKERS}"
echo "GLOBAL_BATCH: $((BATCH_SIZE_PER_GPU * NPROC_PER_NODE))"
echo "DDP_FIND_UNUSED_PARAMETERS: ${DDP_FIND_UNUSED_PARAMETERS}"
echo "CACHE_SPLITS: ${CACHE_SPLITS}"
echo "CACHE_FEATURE_TYPES: ${CACHE_FEATURE_TYPES}"
echo "RUN_HELDOUT_GATE: ${RUN_HELDOUT_GATE}"
echo "HELDOUT_SPLIT: ${HELDOUT_SPLIT}"
echo "HELDOUT_NONPEAK_RATIO: ${HELDOUT_NONPEAK_RATIO}"
echo "HELDOUT_MAX_REGIONS: ${HELDOUT_MAX_REGIONS}"
echo "EVAL_GPU_ID: ${EVAL_GPU_ID}"
echo "TRAIN_DRY_RUN_STEPS: ${TRAIN_DRY_RUN_STEPS}"
echo "=============================================="

cache_manifest_missing=0
for split in "${CACHE_SPLIT_ARGS[@]}"; do
    if [ ! -f "${FOUNDATION_CACHE_DIR}/manifest_${split}.json" ]; then
        cache_manifest_missing=1
        break
    fi
done

if [ "${cache_manifest_missing}" -eq 1 ]; then
    echo ""
    cache_cmd=(
        "${ROOT_DIR}/scripts/build_foundation_cache.py"
        --data_config "${DATA_CONFIG}"
        --train_config "${EFFECTIVE_TRAIN_CONFIG}"
        --output_dir "${FOUNDATION_CACHE_DIR}"
        --backend nt_v2
        --model_dir "${NT_MODEL_DIR}"
        --splits "${CACHE_SPLIT_ARGS[@]}"
        --layers "${CACHE_LAYERS}"
        --feature_types "${CACHE_FEATURE_TYPE_ARGS[@]}"
        --batch_size "${CACHE_BATCH_SIZE}"
        --dtype "${CACHE_DTYPE}"
    )
    if [ "${#CACHE_RECORD_SPLIT_ARGS[@]}" -gt 0 ]; then
        cache_cmd+=(--record_splits "${CACHE_RECORD_SPLIT_ARGS[@]}")
    fi
    if [ "${CACHE_NPROC_PER_NODE}" -gt 1 ]; then
        echo "[Step 0] Building dataset-aligned NT v2 cache on CACHE_GPU_IDS=${CACHE_GPU_IDS} with ${CACHE_NPROC_PER_NODE} workers (splits: ${CACHE_SPLITS})..."
        CUDA_VISIBLE_DEVICES="${CACHE_GPU_IDS}" \
        torchrun \
            --nnodes=1 \
            --nproc_per_node="${CACHE_NPROC_PER_NODE}" \
            --master_addr="${MASTER_ADDR}" \
            --master_port="$((MASTER_PORT_BASE - 1))" \
            "${cache_cmd[@]}"
    else
        echo "[Step 0] Building dataset-aligned NT v2 cache on single GPU ${CACHE_GPU_ID} (splits: ${CACHE_SPLITS})..."
        CUDA_VISIBLE_DEVICES="${CACHE_GPU_ID}" "${PYTHON_BIN}" "${cache_cmd[@]}"
    fi
else
    echo "[Step 0] Cache already exists, skipping."
fi

if [ -n "${BASELINE_CHECKPOINT}" ] && [ -f "${BASELINE_CHECKPOINT}" ] && [ -f "${FOUNDATION_CACHE_DIR}/records_valid.jsonl" ]; then
    echo ""
    echo "[Step 1] Running unified valid probe..."
    "${PYTHON_BIN}" "${ROOT_DIR}/../../scripts/foundation_model_probe.py" \
        --records-jsonl "${FOUNDATION_CACHE_DIR}/records_valid.jsonl" \
        --cache-dir "${FOUNDATION_CACHE_DIR}" \
        --split valid \
        --baseline-checkpoint "${BASELINE_CHECKPOINT}" \
        --baseline-model-config "${BASELINE_MODEL_CONFIG}" \
        --transchrombp-root "${ROOT_DIR}" \
        --output-dir "${OUTPUT_BASE}/foundation_probes/ntv2_tutorial_valid"
else
    echo "[Step 1] Probe skipped. Set BASELINE_CHECKPOINT and ensure records_valid.jsonl exists."
fi

echo ""
echo "[Step 2] Launching residual-head short10 training on TRAIN_GPU_IDS=${TRAIN_GPU_IDS}..."
RUN_NAME="${RUN_NAME:-ntv2_residual_short10_s42}"
train_args=(
    --train-config "${EFFECTIVE_TRAIN_CONFIG}"
    --model-config "${MODEL_CONFIG}"
    --data-config "${DATA_CONFIG}"
    --run-name "${RUN_NAME}"
    --output-dir "${OUTPUT_BASE}"
    --foundation-cache-dir "${FOUNDATION_CACHE_DIR}"
    --batch-size-per-gpu "${BATCH_SIZE_PER_GPU}"
)
if [ "${TRAIN_DRY_RUN_STEPS}" -gt 0 ]; then
    train_args+=(--dry-run-steps "${TRAIN_DRY_RUN_STEPS}")
fi
launch_train "${MASTER_PORT_BASE}" "${train_args[@]}"

if [ "${RUN_HELDOUT_GATE}" = "true" ]; then
    echo ""
    EFFECTIVE_HELDOUT_CONFIG="${RUNTIME_DIR}/train_tutorial_foundation_short10_testfull_runtime.yaml"
    "${PYTHON_BIN}" - "${EFFECTIVE_TRAIN_CONFIG}" "${EFFECTIVE_HELDOUT_CONFIG}" "${HELDOUT_NONPEAK_RATIO}" "${HELDOUT_REGION_SOURCE}" "${HELDOUT_SPLIT}" "${HELDOUT_MAX_REGIONS}" <<'PY'
import sys
from pathlib import Path
import yaml

src, dst, nonpeak_ratio, region_source, split_name, max_regions = sys.argv[1:7]
with open(src, "r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
data_cfg = cfg.setdefault("data", {})
data_cfg["nonpeak_ratio"] = float(nonpeak_ratio)
data_cfg[f"{split_name}_region_source"] = region_source
data_cfg[f"max_{split_name}_regions"] = int(max_regions)
Path(dst).write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
PY

    heldout_manifest="${FOUNDATION_CACHE_DIR}/manifest_${HELDOUT_SPLIT}.json"
    heldout_manifest_matches="false"
    if [ -f "${heldout_manifest}" ]; then
        heldout_manifest_matches="$("${PYTHON_BIN}" - "${heldout_manifest}" "${EFFECTIVE_HELDOUT_CONFIG}" "${DATA_CONFIG}" "${HELDOUT_NONPEAK_RATIO}" "${HELDOUT_REGION_SOURCE}" "${HELDOUT_MAX_REGIONS}" "${HELDOUT_SPLIT}" <<'PY'
import json
import sys
from pathlib import Path

manifest_path, heldout_cfg, data_cfg, nonpeak_ratio, region_source, max_regions, split_name = sys.argv[1:8]
try:
    payload = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
except Exception:
    print("false")
    raise SystemExit(0)
want_train = str(Path(heldout_cfg).resolve())
want_data = str(Path(data_cfg).resolve())
got_train = str(Path(payload.get("train_config_path", "")).resolve())
got_data = str(Path(payload.get("data_config_path", "")).resolve())
want_nonpeak_ratio = float(nonpeak_ratio)
want_region_source = str(region_source)
want_max_records = int(max_regions)
got_split = str(payload.get("dataset_split", payload.get("split", ""))).strip().lower()
got_nonpeak_ratio = float(payload.get("nonpeak_ratio", -1.0))
got_region_source = str(payload.get("region_source", ""))
got_max_records = int(payload.get("max_records", -1))
matches = (
    got_train == want_train
    and got_data == want_data
    and got_split == str(split_name).strip().lower()
    and got_nonpeak_ratio == want_nonpeak_ratio
    and got_region_source == want_region_source
    and got_max_records == want_max_records
)
print("true" if matches else "false")
PY
)"
    fi

    if [ "${heldout_manifest_matches}" != "true" ]; then
        echo "[Step 3] Building held-out cache for split=${HELDOUT_SPLIT} (nonpeak_ratio=${HELDOUT_NONPEAK_RATIO})..."
        heldout_cache_cmd=(
            "${ROOT_DIR}/scripts/build_foundation_cache.py"
            --data_config "${DATA_CONFIG}"
            --train_config "${EFFECTIVE_HELDOUT_CONFIG}"
            --output_dir "${FOUNDATION_CACHE_DIR}"
            --backend nt_v2
            --model_dir "${NT_MODEL_DIR}"
            --splits "${HELDOUT_SPLIT}"
            --layers "${CACHE_LAYERS}"
            --feature_types "${CACHE_FEATURE_TYPE_ARGS[@]}"
            --batch_size "${CACHE_BATCH_SIZE}"
            --dtype "${CACHE_DTYPE}"
        )
        if [ "${HELDOUT_CACHE_NPROC_PER_NODE}" -gt 1 ]; then
            CUDA_VISIBLE_DEVICES="${HELDOUT_CACHE_GPU_IDS}" \
            torchrun \
                --nnodes=1 \
                --nproc_per_node="${HELDOUT_CACHE_NPROC_PER_NODE}" \
                --master_addr="${MASTER_ADDR}" \
                --master_port="$((MASTER_PORT_BASE + 1))" \
                "${heldout_cache_cmd[@]}"
        else
            CUDA_VISIBLE_DEVICES="${HELDOUT_CACHE_GPU_ID}" "${PYTHON_BIN}" "${heldout_cache_cmd[@]}"
        fi
    else
        echo "[Step 3] Held-out cache manifest already matches current held-out runtime config, skipping rebuild."
    fi

    BEST_CKPT="${OUTPUT_BASE}/checkpoints/${RUN_NAME}/best.pt"
    if [ ! -f "${BEST_CKPT}" ]; then
        echo "[error] Missing best checkpoint: ${BEST_CKPT}" >&2
        exit 1
    fi

    mkdir -p "${OUTPUT_BASE}/metrics"
    HELDOUT_TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
    HELDOUT_OUTPUT_PATH="${HELDOUT_OUTPUT:-${OUTPUT_BASE}/metrics/${RUN_NAME}_best_${HELDOUT_SPLIT}_full_${HELDOUT_TIMESTAMP}.json}"

    echo "[Step 4] Running held-out gate evaluation on GPU ${EVAL_GPU_ID}..."
    eval_cmd=(
        "${PYTHON_BIN}" -m transchrombp.evaluation.evaluate_checkpoint
        --checkpoint "${BEST_CKPT}"
        --split "${HELDOUT_SPLIT}"
        --data-config "${DATA_CONFIG}"
        --nonpeak-ratio "${HELDOUT_NONPEAK_RATIO}"
        --region-source "${HELDOUT_REGION_SOURCE}"
        --batch-size "${BATCH_SIZE_PER_GPU}"
        --num-workers "${NUM_WORKERS}"
        --output "${HELDOUT_OUTPUT_PATH}"
        --log-every 100
    )
    if [ "${HELDOUT_MAX_REGIONS}" -gt 0 ]; then
        eval_cmd+=(--max-regions "${HELDOUT_MAX_REGIONS}")
    fi
    CUDA_VISIBLE_DEVICES="${EVAL_GPU_ID}" "${eval_cmd[@]}"

    echo "[Step 4] Held-out metrics written: ${HELDOUT_OUTPUT_PATH}"
fi

echo ""
echo "=============================================="
echo "NT v2 residual short10 completed."
echo "Results in: ${OUTPUT_BASE}/checkpoints, ${OUTPUT_BASE}/logs, ${OUTPUT_BASE}/metrics"
echo "=============================================="
