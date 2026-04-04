#!/usr/bin/env bash
set -euo pipefail

REMOTE_HOST="zhengwei@127.0.0.1"
REMOTE_PORT="6002"
REMOTE_KEY="/home/zhengwei/.ssh/codex_6002_ed25519"
REMOTE_ROOT="${REMOTE_ROOT:-/home/zhengwei/bylw_atac}"
MODE="${1:-cell_lines}"
RUN_TAG="${RUN_TAG:-dataset_cache_downloads_$(date +%Y%m%d_%H%M%S)}"
REMOTE_JOB_DIR="$REMOTE_ROOT/.codex_jobs"
REMOTE_SCRIPT="$REMOTE_JOB_DIR/${RUN_TAG}.sh"
REMOTE_LOG="$REMOTE_ROOT/logs/${RUN_TAG}.log"

case "$MODE" in
  cell_lines|peaks_only|all) ;;
  *)
    echo "usage: $0 [cell_lines|peaks_only|all]" >&2
    exit 1
    ;;
esac

ssh_base=(
  ssh
  -i "$REMOTE_KEY"
  -p "$REMOTE_PORT"
  "$REMOTE_HOST"
)

"${ssh_base[@]}" "mkdir -p '$REMOTE_JOB_DIR' '$REMOTE_ROOT/logs' '$REMOTE_ROOT/chrombpnet_datasets' '$REMOTE_ROOT/ATACseq/peaks'"

"${ssh_base[@]}" "cat > '$REMOTE_SCRIPT.tmp'" <<'REMOTE_SCRIPT'
#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:?missing root}"
MODE="${2:-cell_lines}"

CELL_DIR="$ROOT/chrombpnet_datasets"
PEAK_DIR="$ROOT/ATACseq/peaks"

mkdir -p "$CELL_DIR/GM12878" "$CELL_DIR/K562" "$PEAK_DIR"

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

download_cell_lines() {
  download_to "https://www.encodeproject.org/files/ENCFF962FMH/@@download/ENCFF962FMH.bam" "$CELL_DIR/GM12878/rep1.bam"
  download_to "https://www.encodeproject.org/files/ENCFF981FXV/@@download/ENCFF981FXV.bam" "$CELL_DIR/GM12878/rep2.bam"
  download_to "https://www.encodeproject.org/files/ENCFF440GRZ/@@download/ENCFF440GRZ.bam" "$CELL_DIR/GM12878/rep3.bam"
  download_to "https://www.encodeproject.org/files/ENCFF748UZH/@@download/ENCFF748UZH.bed.gz" "$CELL_DIR/GM12878/overlap.bed.gz"

  download_to "https://www.encodeproject.org/files/ENCFF852RRY/@@download/ENCFF852RRY.bam" "$CELL_DIR/K562/rep1.bam"
  download_to "https://www.encodeproject.org/files/ENCFF557RSA/@@download/ENCFF557RSA.bam" "$CELL_DIR/K562/rep2.bam"
  download_to "https://www.encodeproject.org/files/ENCFF444DLQ/@@download/ENCFF444DLQ.bed.gz" "$CELL_DIR/K562/overlap.bed.gz"
}

download_peak_pool() {
  download_to "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5214nnn/GSM5214550/suppl/GSM5214550_ENCFF203QDA_IDR_ranked_peaks_GRCh38.bed.gz" "$PEAK_DIR/GSM5214550_ENCFF203QDA_IDR_ranked_peaks_GRCh38.bed.gz"
  download_to "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5214nnn/GSM5214756/suppl/GSM5214756_ENCFF124CDP_IDR_ranked_peaks_GRCh38.bed.gz" "$PEAK_DIR/GSM5214756_ENCFF124CDP_IDR_ranked_peaks_GRCh38.bed.gz"
  download_to "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE267154&format=file&file=GSE267154_overlap.optimal_peak.narrowPeak.gz" "$PEAK_DIR/GSE267154_overlap.optimal_peak.narrowPeak.gz"
  download_to "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5214nnn/GSM5214186/suppl/GSM5214186_ENCFF204FWC_IDR_thresholded_peaks_GRCh38.bigBed" "$PEAK_DIR/GSM5214186_ENCFF204FWC_IDR_thresholded_peaks_GRCh38.bigBed"
  download_to "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5214nnn/GSM5214097/suppl/GSM5214097_ENCFF149QZK_IDR_ranked_peaks_GRCh38.bed.gz" "$PEAK_DIR/GSM5214097_ENCFF149QZK_IDR_ranked_peaks_GRCh38.bed.gz"
  download_to "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5214nnn/GSM5214058/suppl/GSM5214058_ENCFF428ROW_IDR_ranked_peaks_GRCh38.bed.gz" "$PEAK_DIR/GSM5214058_ENCFF428ROW_IDR_ranked_peaks_GRCh38.bed.gz"
}

log "start dataset cache downloads mode=$MODE root=$ROOT"

if [ "$MODE" = "cell_lines" ] || [ "$MODE" = "all" ]; then
  download_cell_lines
fi

if [ "$MODE" = "peaks_only" ] || [ "$MODE" = "all" ]; then
  download_peak_pool
fi

log "download summary"
du -sh "$CELL_DIR/GM12878" "$CELL_DIR/K562" "$PEAK_DIR" 2>/dev/null || true
log "finished dataset cache downloads"
REMOTE_SCRIPT

"${ssh_base[@]}" "mv '$REMOTE_SCRIPT.tmp' '$REMOTE_SCRIPT' && chmod +x '$REMOTE_SCRIPT'"

pid="$("${ssh_base[@]}" "nohup bash '$REMOTE_SCRIPT' '$REMOTE_ROOT' '$MODE' > '$REMOTE_LOG' 2>&1 < /dev/null & echo \$!")"

echo "started pid=$pid"
echo "remote_script=$REMOTE_SCRIPT"
echo "remote_log=$REMOTE_LOG"
echo "check: ssh -i $REMOTE_KEY -p $REMOTE_PORT $REMOTE_HOST 'tail -n 40 $REMOTE_LOG'"
