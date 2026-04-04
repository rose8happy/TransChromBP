#!/usr/bin/env bash
set -euo pipefail

REMOTE_HOST="zhengwei@127.0.0.1"
REMOTE_PORT="6002"
REMOTE_KEY="/home/zhengwei/.ssh/codex_6002_ed25519"
REMOTE_ROOT="${REMOTE_ROOT:-/home/zhengwei/bylw_atac}"
RUN_TAG="${RUN_TAG:-tutorial_raw_prep_$(date +%Y%m%d_%H%M%S)}"
REMOTE_JOB_DIR="$REMOTE_ROOT/.codex_jobs"
REMOTE_SCRIPT="$REMOTE_JOB_DIR/${RUN_TAG}.sh"
REMOTE_LOG="$REMOTE_ROOT/logs/${RUN_TAG}.log"

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
DATA_DIR="$ROOT/chrombpnet_tutorial/data"
ENV_DIR="$ROOT/.mamba/envs/transchrombp"

mkdir -p "$DATA_DIR"
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"

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

  if [ -s "$out" ] && samtools quickcheck "$out" >/dev/null 2>&1; then
    log "skip existing valid $out"
    return 0
  fi

  rm -f "$part"
  log "download $url -> $out"
  wget -c --tries=5 --waitretry=5 -O "$part" "$url"
  mv "$part" "$out"
  samtools quickcheck "$out"
}

log "start tutorial raw prep root=$ROOT"

download_to "https://www.encodeproject.org/files/ENCFF077FBI/@@download/ENCFF077FBI.bam" "$DATA_DIR/rep1.bam"
download_to "https://www.encodeproject.org/files/ENCFF128WZG/@@download/ENCFF128WZG.bam" "$DATA_DIR/rep2.bam"
download_to "https://www.encodeproject.org/files/ENCFF534DCE/@@download/ENCFF534DCE.bam" "$DATA_DIR/rep3.bam"

log "merge and sort tutorial BAMs"
rm -f "$DATA_DIR/merged_unsorted.bam.tmp" "$DATA_DIR/merged.bam.tmp" "$DATA_DIR/merged.bam.tmp.bai"
samtools merge -f "$DATA_DIR/merged_unsorted.bam.tmp" "$DATA_DIR/rep1.bam" "$DATA_DIR/rep2.bam" "$DATA_DIR/rep3.bam"
samtools sort -@4 "$DATA_DIR/merged_unsorted.bam.tmp" -o "$DATA_DIR/merged.bam.tmp"
samtools index "$DATA_DIR/merged.bam.tmp"
mv "$DATA_DIR/merged.bam.tmp" "$DATA_DIR/merged.bam"
mv "$DATA_DIR/merged.bam.tmp.bai" "$DATA_DIR/merged.bam.bai"
rm -f "$DATA_DIR/merged_unsorted.bam.tmp"
samtools quickcheck "$DATA_DIR/merged.bam"

log "tutorial raw prep summary"
ls -lh "$DATA_DIR/rep1.bam" "$DATA_DIR/rep2.bam" "$DATA_DIR/rep3.bam" "$DATA_DIR/merged.bam" "$DATA_DIR/merged.bam.bai"
log "finished tutorial raw prep"
REMOTE_SCRIPT

"${ssh_base[@]}" "mv '$REMOTE_SCRIPT.tmp' '$REMOTE_SCRIPT' && chmod +x '$REMOTE_SCRIPT'"

pid="$("${ssh_base[@]}" "nohup bash '$REMOTE_SCRIPT' '$REMOTE_ROOT' > '$REMOTE_LOG' 2>&1 < /dev/null & echo \$!")"

echo "started pid=$pid"
echo "remote_script=$REMOTE_SCRIPT"
echo "remote_log=$REMOTE_LOG"
echo "check: ssh -i $REMOTE_KEY -p $REMOTE_PORT $REMOTE_HOST 'tail -n 40 $REMOTE_LOG'"
