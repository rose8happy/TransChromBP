#!/usr/bin/env bash
set -euo pipefail

REMOTE_HOST="zhengwei@127.0.0.1"
REMOTE_PORT="6002"
REMOTE_KEY="/home/zhengwei/.ssh/codex_6002_ed25519"
REMOTE_ROOT="${REMOTE_ROOT:-/home/zhengwei/bylw_atac}"
MODE="${1:-core}"
RUN_TAG="${RUN_TAG:-transchrombp_official_downloads_$(date +%Y%m%d_%H%M%S)}"
REMOTE_JOB_DIR="$REMOTE_ROOT/.codex_jobs"
REMOTE_SCRIPT="$REMOTE_JOB_DIR/${RUN_TAG}.sh"
REMOTE_LOG="$REMOTE_ROOT/logs/${RUN_TAG}.log"

case "$MODE" in
  core|raw|all) ;;
  *)
    echo "usage: $0 [core|raw|all]" >&2
    exit 1
    ;;
esac

ssh_base=(
  ssh
  -i "$REMOTE_KEY"
  -p "$REMOTE_PORT"
  "$REMOTE_HOST"
)

"${ssh_base[@]}" "mkdir -p '$REMOTE_JOB_DIR' '$REMOTE_ROOT/logs' '$REMOTE_ROOT/chrombpnet_tutorial/data'"

"${ssh_base[@]}" "cat > '$REMOTE_SCRIPT.tmp'" <<'REMOTE_SCRIPT'
#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:?missing root}"
MODE="${2:-core}"
DATA_DIR="$ROOT/chrombpnet_tutorial/data"
FOLDS_DIR="$DATA_DIR/folds"
BIAS_DIR="$DATA_DIR/bias_models"

mkdir -p "$DATA_DIR" "$FOLDS_DIR" "$BIAS_DIR"

timestamp() {
  date '+%Y-%m-%d %H:%M:%S'
}

log() {
  printf '[%s] %s\n' "$(timestamp)" "$*"
}

download_to() {
  local url="$1"
  local out="$2"
  local part="${out}.part"

  mkdir -p "$(dirname "$out")"

  if [ -s "$out" ] && [ ! -e "$part" ]; then
    log "skip existing $out"
    return 0
  fi

  log "download $url -> $out"
  wget -c --tries=5 --waitretry=5 -O "$part" "$url"
  mv "$part" "$out"
}

extract_zip() {
  local zip_path="$1"
  local out_dir="$2"

  mkdir -p "$out_dir"
  if command -v unzip >/dev/null 2>&1; then
    unzip -oq "$zip_path" -d "$out_dir"
  else
    python3 - "$zip_path" "$out_dir" <<'PY'
import pathlib
import sys
import zipfile

zip_path = pathlib.Path(sys.argv[1])
out_dir = pathlib.Path(sys.argv[2])
out_dir.mkdir(parents=True, exist_ok=True)
with zipfile.ZipFile(zip_path) as zf:
    zf.extractall(out_dir)
PY
  fi
}

log "start official downloads mode=$MODE root=$ROOT"

download_to \
  "https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz" \
  "$DATA_DIR/hg38.fa.gz"

hg38_target="$ROOT/hg38.fa"
hg38_min_bytes=3000000000
hg38_size=0
if [ -e "$hg38_target" ]; then
  hg38_size="$(stat -c %s "$hg38_target" 2>/dev/null || echo 0)"
fi

if [ "$hg38_size" -lt "$hg38_min_bytes" ]; then
  log "build $DATA_DIR/hg38.fa.gz -> $hg38_target (current_size=$hg38_size)"
  gzip -dc "$DATA_DIR/hg38.fa.gz" > "$ROOT/hg38.fa.tmp"
  mv "$ROOT/hg38.fa.tmp" "$ROOT/hg38.fa"
else
  log "skip existing $hg38_target size=$hg38_size"
fi

download_to \
  "https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv" \
  "$DATA_DIR/hg38.chrom.sizes"

download_to \
  "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz" \
  "$DATA_DIR/blacklist.bed.gz"

download_to \
  "https://www.encodeproject.org/files/ENCFF333TAT/@@download/ENCFF333TAT.bed.gz" \
  "$DATA_DIR/overlap.bed.gz"

download_to \
  "https://storage.googleapis.com/chrombpnet_data/input_files/ENCSR868FGK_merged.bam" \
  "$DATA_DIR/ENCSR868FGK_merged.bam"

download_to \
  "https://storage.googleapis.com/chrombpnet_data/input_files/ENCSR868FGK_relaxed_peaks_no_blacklist.bed" \
  "$DATA_DIR/ENCSR868FGK_relaxed_peaks_no_blacklist.bed"

download_to \
  "https://storage.googleapis.com/chrombpnet_data/input_files/ENCSR868FGK_nonpeaks_no_blacklist.bed" \
  "$DATA_DIR/ENCSR868FGK_nonpeaks_no_blacklist.bed"

download_to \
  "https://zenodo.org/records/7443683/files/folds.zip?download=1" \
  "$DATA_DIR/folds.zip"
log "extract $DATA_DIR/folds.zip -> $FOLDS_DIR"
extract_zip "$DATA_DIR/folds.zip" "$FOLDS_DIR"

download_to \
  "https://zenodo.org/records/7443683/files/bias_models.zip?download=1" \
  "$DATA_DIR/bias_models.zip"
log "extract $DATA_DIR/bias_models.zip -> $BIAS_DIR"
extract_zip "$DATA_DIR/bias_models.zip" "$BIAS_DIR"

if [ "$MODE" = "raw" ] || [ "$MODE" = "all" ]; then
  download_to \
    "https://www.encodeproject.org/files/ENCFF077FBI/@@download/ENCFF077FBI.bam" \
    "$DATA_DIR/rep1.bam"
  download_to \
    "https://www.encodeproject.org/files/ENCFF128WZG/@@download/ENCFF128WZG.bam" \
    "$DATA_DIR/rep2.bam"
  download_to \
    "https://www.encodeproject.org/files/ENCFF534DCE/@@download/ENCFF534DCE.bam" \
    "$DATA_DIR/rep3.bam"
fi

log "download summary"
du -sh "$ROOT/hg38.fa" "$DATA_DIR" 2>/dev/null || true
find "$DATA_DIR" -maxdepth 2 -type f | sort
log "finished official downloads"
REMOTE_SCRIPT

"${ssh_base[@]}" "mv '$REMOTE_SCRIPT.tmp' '$REMOTE_SCRIPT' && chmod +x '$REMOTE_SCRIPT'"

pid="$("${ssh_base[@]}" "nohup bash '$REMOTE_SCRIPT' '$REMOTE_ROOT' '$MODE' > '$REMOTE_LOG' 2>&1 < /dev/null & echo \$!")"

echo "started pid=$pid"
echo "remote_script=$REMOTE_SCRIPT"
echo "remote_log=$REMOTE_LOG"
echo "check: ssh -i $REMOTE_KEY -p $REMOTE_PORT $REMOTE_HOST 'tail -n 40 $REMOTE_LOG'"
