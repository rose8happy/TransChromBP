#!/usr/bin/env bash
set -euo pipefail

ROOT_BASE="${ROOT_BASE:-/data1/zhoujiazhen/bylw_atac}"
VENV_DIR="${VENV_DIR:-${ROOT_BASE}/.venvs/caduceus-ps}"
MODEL_DIR="${MODEL_DIR:-${ROOT_BASE}/foundation_models/caduceus/caduceus-ps_seqlen-131k_d_model-256_n_layer-16}"
REPO_DIR="${REPO_DIR:-${ROOT_BASE}/foundation_models/caduceus/caduceus_repo}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
TORCH_INDEX_URL="${TORCH_INDEX_URL:-https://download.pytorch.org/whl/cu121}"
CADUCEUS_SAFETENSORS_URL="${CADUCEUS_SAFETENSORS_URL:-}"

DO_VENV=1
DO_DOWNLOAD=1
DO_SMOKE=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-venv) DO_VENV=0; shift ;;
        --skip-download) DO_DOWNLOAD=0; shift ;;
        --skip-smoke) DO_SMOKE=0; shift ;;
        *)
            echo "[error] Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$1"; }

mkdir -p "$(dirname "${VENV_DIR}")" "$(dirname "${MODEL_DIR}")"

write_model_metadata() {
    mkdir -p "${MODEL_DIR}"

    cat > "${MODEL_DIR}/config.json" <<'EOF'
{
  "architectures": [
    "CaduceusForMaskedLM"
  ],
  "auto_map": {
    "AutoConfig": "configuration_caduceus.CaduceusConfig",
    "AutoModel": "modeling_caduceus.Caduceus",
    "AutoModelForMaskedLM": "modeling_caduceus.CaduceusForMaskedLM",
    "AutoModelForSequenceClassification": "modeling_caduceus.CaduceusForSequenceClassification"
  },
  "bidirectional": true,
  "bidirectional_strategy": "add",
  "bidirectional_weight_tie": true,
  "complement_map": {
    "0": 0,
    "1": 1,
    "2": 2,
    "3": 3,
    "4": 4,
    "5": 5,
    "6": 6,
    "7": 10,
    "8": 9,
    "9": 8,
    "10": 7,
    "11": 11,
    "12": 12,
    "13": 13,
    "14": 14,
    "15": 15
  },
  "d_model": 256,
  "fused_add_norm": true,
  "initializer_cfg": {
    "initializer_range": 0.02,
    "n_residuals_per_layer": 1,
    "rescale_prenorm_residual": true
  },
  "model_type": "caduceus",
  "n_layer": 16,
  "norm_epsilon": 1e-05,
  "pad_vocab_size_multiple": 8,
  "rcps": true,
  "residual_in_fp32": false,
  "rms_norm": true,
  "ssm_cfg": {
    "bias": false,
    "conv_bias": true,
    "d_conv": 4,
    "d_state": 16,
    "dt_init": "random",
    "dt_init_floor": 0.0001,
    "dt_max": 0.1,
    "dt_min": 0.001,
    "dt_rank": "auto",
    "dt_scale": 1.0,
    "expand": 2,
    "use_fast_path": true
  },
  "torch_dtype": "float32",
  "transformers_version": "4.38.1",
  "vocab_size": 16
}
EOF

    cat > "${MODEL_DIR}/tokenizer_config.json" <<'EOF'
{
  "add_prefix_space": false,
  "added_tokens_decoder": {
    "0": {
      "content": "[CLS]",
      "lstrip": false,
      "normalized": false,
      "rstrip": false,
      "single_word": false,
      "special": true
    },
    "1": {
      "content": "[SEP]",
      "lstrip": false,
      "normalized": false,
      "rstrip": false,
      "single_word": false,
      "special": true
    },
    "2": {
      "content": "[BOS]",
      "lstrip": false,
      "normalized": false,
      "rstrip": false,
      "single_word": false,
      "special": true
    },
    "3": {
      "content": "[MASK]",
      "lstrip": false,
      "normalized": false,
      "rstrip": false,
      "single_word": false,
      "special": true
    },
    "4": {
      "content": "[PAD]",
      "lstrip": false,
      "normalized": false,
      "rstrip": false,
      "single_word": false,
      "special": true
    },
    "6": {
      "content": "[UNK]",
      "lstrip": false,
      "normalized": false,
      "rstrip": false,
      "single_word": false,
      "special": true
    }
  },
  "auto_map": {
    "AutoTokenizer": [
      "tokenization_caduceus.CaduceusTokenizer",
      null
    ]
  },
  "bos_token": "[BOS]",
  "clean_up_tokenization_spaces": true,
  "cls_token": "[CLS]",
  "eos_token": "[SEP]",
  "mask_token": "[MASK]",
  "model_max_length": 131072,
  "pad_token": "[PAD]",
  "padding_side": "left",
  "sep_token": "[SEP]",
  "tokenizer_class": "CaduceusTokenizer",
  "unk_token": "[UNK]"
}
EOF

    cat > "${MODEL_DIR}/special_tokens_map.json" <<'EOF'
{
  "bos_token": "[BOS]",
  "cls_token": "[CLS]",
  "eos_token": "[SEP]",
  "mask_token": "[MASK]",
  "pad_token": "[PAD]",
  "sep_token": "[SEP]",
  "unk_token": "[UNK]"
}
EOF
}

copy_repo_python_files() {
    for file in configuration_caduceus.py modeling_caduceus.py modeling_rcps.py tokenization_caduceus.py; do
        cp "${REPO_DIR}/caduceus/${file}" "${MODEL_DIR}/${file}"
    done
}

if [[ "${DO_VENV}" == "1" ]]; then
    if [[ ! -d "${VENV_DIR}" ]]; then
        log "[env] create venv ${VENV_DIR}"
        "${PYTHON_BIN}" -m venv "${VENV_DIR}"
    fi

    # shellcheck source=/dev/null
    source "${VENV_DIR}/bin/activate"

    log "[env] upgrade pip/setuptools/wheel"
    python -m pip install --upgrade pip setuptools wheel packaging ninja

    log "[env] pin numpy<2 for torch/mamba compatibility"
    python -m pip install --upgrade "numpy<2"

    log "[env] install torch 2.2.0 cu121"
    python -m pip install \
        --index-url "${TORCH_INDEX_URL}" \
        torch==2.2.0 \
        torchvision==0.17.0 \
        torchaudio==2.2.0

    log "[env] install transformer/hf runtime deps"
    python -m pip install \
        huggingface-hub==0.24.7 \
        transformers==4.38.1 \
        safetensors==0.4.2 \
        sentencepiece==0.2.0 \
        einops==0.7.0 \
        biopython==1.81

    log "[env] install TransChromBP data/runtime deps"
    python -m pip install \
        PyYAML==6.0.2 \
        pyBigWig==0.3.24 \
        pyfaidx==0.8.1.3

    log "[env] install causal-conv1d + mamba-ssm"
    python -m pip install --no-build-isolation causal-conv1d==1.2.0.post2
    python -m pip install --no-build-isolation mamba-ssm==1.2.0.post1
else
    if [[ ! -x "${VENV_DIR}/bin/python" ]]; then
        echo "[error] Missing existing venv: ${VENV_DIR}" >&2
        exit 1
    fi
    # shellcheck source=/dev/null
    source "${VENV_DIR}/bin/activate"
fi

if [[ "${DO_DOWNLOAD}" == "1" ]]; then
    if [[ ! -d "${REPO_DIR}/.git" ]]; then
        log "[repo] clone official Caduceus repo -> ${REPO_DIR}"
        git clone https://github.com/kuleshov-group/caduceus.git "${REPO_DIR}"
    else
        log "[repo] update official Caduceus repo"
        git -C "${REPO_DIR}" pull --ff-only
    fi

    if [[ -n "${CADUCEUS_SAFETENSORS_URL}" ]]; then
        log "[direct] assemble local model directory from repo + embedded metadata"
        write_model_metadata
        copy_repo_python_files

        log "[direct] download model.safetensors from resolved URL -> ${MODEL_DIR}"
        curl --fail --location --retry 3 --retry-delay 2 \
            --output "${MODEL_DIR}/model.safetensors" \
            "${CADUCEUS_SAFETENSORS_URL}"
    else
        log "[hf] download Caduceus-PS snapshot -> ${MODEL_DIR}"
        python - <<'PY'
from huggingface_hub import snapshot_download

snapshot_download(
    repo_id="kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16",
    local_dir="/data1/zhoujiazhen/bylw_atac/foundation_models/caduceus/caduceus-ps_seqlen-131k_d_model-256_n_layer-16",
    local_dir_use_symlinks=False,
)
PY
    fi
fi

if [[ "${DO_SMOKE}" == "1" ]]; then
    log "[smoke] load tokenizer/model and extract hidden states"
    python - <<'PY'
import torch
from transformers import AutoModelForMaskedLM, AutoTokenizer

model_dir = "/data1/zhoujiazhen/bylw_atac/foundation_models/caduceus/caduceus-ps_seqlen-131k_d_model-256_n_layer-16"
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
tokenizer = AutoTokenizer.from_pretrained(model_dir, local_files_only=True, trust_remote_code=True)
model = AutoModelForMaskedLM.from_pretrained(
    model_dir,
    local_files_only=True,
    trust_remote_code=True,
    torch_dtype=torch.bfloat16,
)

seq = "ACGT" * 32
inputs = tokenizer([seq], return_tensors="pt", padding=True, add_special_tokens=False)
model = model.eval()

print("device", device.type)
print("tokenizer", tokenizer.__class__.__name__)
print("model", model.__class__.__name__)
print("input_ids", tuple(inputs["input_ids"].shape))

if device.type == "cuda":
    model = model.to(device)
    inputs = {k: v.to(device) for k, v in inputs.items()}
    with torch.no_grad():
        outputs = model(**inputs, output_hidden_states=True)
    print("hidden_last", tuple(outputs.hidden_states[-1].shape))
else:
    print("hidden_last", "skipped_cpu_forward")

PY
fi

log "[done] Caduceus-PS assets ready"
