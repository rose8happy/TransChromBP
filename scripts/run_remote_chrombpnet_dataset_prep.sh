#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash run_remote_chrombpnet_dataset_prep.sh \
    --root /path/to/bylw_atac \
    --env-dir /path/to/env \
    [--python-bin /path/to/python] \
    [--datasets GM12878,K562] \
    [--threads 4] \
    [--nice-level 10] \
    [--inputlen 2114] \
    [--stride 1000]
EOF
}

ROOT=""
ENV_DIR=""
PYTHON_BIN=""
DATASETS="GM12878,K562"
THREADS="${THREADS:-4}"
NICE_LEVEL="${NICE_LEVEL:-10}"
INPUTLEN="${INPUTLEN:-2114}"
STRIDE="${STRIDE:-1000}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --root)
      ROOT="$2"
      shift 2
      ;;
    --env-dir)
      ENV_DIR="$2"
      shift 2
      ;;
    --python-bin)
      PYTHON_BIN="$2"
      shift 2
      ;;
    --datasets)
      DATASETS="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --nice-level)
      NICE_LEVEL="$2"
      shift 2
      ;;
    --inputlen)
      INPUTLEN="$2"
      shift 2
      ;;
    --stride)
      STRIDE="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

ROOT="${ROOT:?missing --root}"
ENV_DIR="${ENV_DIR:?missing --env-dir}"
PYTHON_BIN="${PYTHON_BIN:-$ENV_DIR/bin/python}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GC_BINS_PY="$SCRIPT_DIR/get_genomewide_gc_bins.py"
GC_CONTENT_PY="$SCRIPT_DIR/get_gc_content.py"
GC_MATCH_PY="$SCRIPT_DIR/get_gc_matched_negatives.py"

DATA_ROOT="$ROOT/chrombpnet_datasets"
TUTORIAL_DATA_DIR="${TUTORIAL_DATA_DIR:-$ROOT/chrombpnet_tutorial/data}"
REF_DIR="${REF_DIR:-$ROOT/chrombpnet_refs}"
GENOME="${GENOME:-$ROOT/hg38.fa}"
CHROM_SIZES="${CHROM_SIZES:-$TUTORIAL_DATA_DIR/hg38.chrom.sizes}"
BLACKLIST_GZ="${BLACKLIST_GZ:-$TUTORIAL_DATA_DIR/blacklist.bed.gz}"
FOLDS_JSON="${FOLDS_JSON:-$TUTORIAL_DATA_DIR/folds.json}"
GC_PREFIX="$REF_DIR/genomewide_gc_hg38_stride_${STRIDE}_inputlen_${INPUTLEN}"
GC_BED="${GC_PREFIX}.bed"
BLACKLIST_BED="$REF_DIR/blacklist.bed"
TMP_ROOT="${TMP_ROOT:-$ROOT/tmp_chrombpnet_dataset_prep}"

export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export MPLBACKEND=Agg

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

fail() {
  log "ERROR: $*"
  exit 1
}

require_cmd() {
  local cmd="$1"
  command -v "$cmd" >/dev/null 2>&1 || fail "command not found: $cmd"
}

require_file() {
  local path="$1"
  [[ -f "$path" ]] || fail "missing file: $path"
}

run_low() {
  log "run: $*"
  nice -n "$NICE_LEVEL" "$@"
}

ensure_link() {
  local target="$1"
  local link_path="$2"

  if [[ -e "$link_path" ]]; then
    return 0
  fi

  ln -s "$target" "$link_path"
}

decompress_if_needed() {
  local src="$1"
  local dst="$2"
  local dst_tmp="${dst}.tmp"

  if [[ -s "$dst" ]]; then
    return 0
  fi

  mkdir -p "$(dirname "$dst")"
  rm -f "$dst_tmp"

  if [[ "$src" == *.gz ]]; then
    log "decompress: $src -> $dst"
    gzip -dc "$src" > "$dst_tmp"
  else
    log "copy: $src -> $dst"
    cp "$src" "$dst_tmp"
  fi

  mv "$dst_tmp" "$dst"
}

prepare_references() {
  mkdir -p "$REF_DIR" "$TMP_ROOT"

  require_file "$GC_BINS_PY"
  require_file "$GC_CONTENT_PY"
  require_file "$GC_MATCH_PY"
  require_file "$GENOME"
  require_file "$CHROM_SIZES"
  require_file "$BLACKLIST_GZ"
  require_file "$FOLDS_JSON"

  require_cmd samtools
  require_cmd bedtools
  require_cmd bedGraphToBigWig
  require_cmd gzip
  require_cmd awk
  require_cmd sort

  "$PYTHON_BIN" - <<'PY'
import importlib
for module in ["pyBigWig", "pyfaidx", "pandas", "numpy", "matplotlib", "tqdm"]:
    importlib.import_module(module)
PY

  if [[ ! -s "${GENOME}.fai" ]]; then
    run_low samtools faidx "$GENOME"
  fi

  if [[ ! -s "$BLACKLIST_BED" ]]; then
    decompress_if_needed "$BLACKLIST_GZ" "$BLACKLIST_BED"
  fi

  if [[ ! -s "$GC_BED" ]]; then
    run_low "$PYTHON_BIN" "$GC_BINS_PY" \
      -g "$GENOME" \
      -o "$GC_PREFIX" \
      -f "$INPUTLEN" \
      -s "$STRIDE"
  else
    log "skip existing genomewide gc bins: $GC_BED"
  fi
}

prepare_dataset_links() {
  local dataset_dir="$1"
  ensure_link "$GENOME" "$dataset_dir/hg38.fa"
  ensure_link "${GENOME}.fai" "$dataset_dir/hg38.fa.fai"
  ensure_link "$CHROM_SIZES" "$dataset_dir/hg38.chrom.sizes"
  ensure_link "$BLACKLIST_GZ" "$dataset_dir/blacklist.bed.gz"
  ensure_link "$FOLDS_JSON" "$dataset_dir/folds.json"
}

prepare_peaks_bed() {
  local dataset_dir="$1"
  local prep_dir="$2"
  local peaks_src_gz="$dataset_dir/overlap.bed.gz"
  local peaks_src_plain="$dataset_dir/overlap.bed"
  local peaks_dst="$prep_dir/inputs/overlap.bed"

  if [[ -s "$peaks_src_gz" ]]; then
    decompress_if_needed "$peaks_src_gz" "$peaks_dst"
  elif [[ -s "$peaks_src_plain" ]]; then
    decompress_if_needed "$peaks_src_plain" "$peaks_dst"
  else
    fail "missing overlap peaks for $dataset_dir"
  fi
}

prepare_merged_bam() {
  local dataset_dir="$1"
  local merged_bam="$dataset_dir/merged.bam"
  local merged_tmp="${merged_bam}.tmp"
  local merged_bai="${merged_bam}.bai"
  local bam_count=0
  local bam_inputs=()

  mapfile -t bam_inputs < <(find "$dataset_dir" -maxdepth 1 -type f -name 'rep*.bam' | sort -V)
  bam_count="${#bam_inputs[@]}"
  [[ "$bam_count" -gt 0 ]] || fail "no rep*.bam found under $dataset_dir"

  if [[ ! -s "$merged_bam" ]]; then
    rm -f "$merged_tmp"
    if [[ "$bam_count" -eq 1 ]]; then
      log "single input bam, copying to merged.bam"
      cp "${bam_inputs[0]}" "$merged_tmp"
    else
      run_low samtools merge -@ "$THREADS" -f "$merged_tmp" "${bam_inputs[@]}"
    fi
    mv "$merged_tmp" "$merged_bam"
  else
    log "skip existing merged bam: $merged_bam"
  fi

  run_low samtools quickcheck "$merged_bam"

  if [[ ! -s "$merged_bai" ]]; then
    run_low samtools index -@ "$THREADS" "$merged_bam"
  else
    log "skip existing bam index: $merged_bai"
  fi
}

prepare_bigwig() {
  local dataset_dir="$1"
  local prep_dir="$2"
  local merged_bam="$dataset_dir/merged.bam"
  local bigwig_dir="$prep_dir/bigwig"
  local tmp_dir="$TMP_ROOT/$(basename "$dataset_dir")/bigwig_sort"
  local out_prefix="$bigwig_dir/merged"
  local bedgraph="${out_prefix}.unstranded.bedGraph"
  local bedgraph_tmp="${bedgraph}.tmp"
  local bigwig="${out_prefix}_unstranded.bw"

  mkdir -p "$bigwig_dir" "$tmp_dir"

  if [[ -s "$bigwig" ]]; then
    log "skip existing bigwig: $bigwig"
    return 0
  fi

  rm -f "$bedgraph_tmp"
  log "build bigwig from $merged_bam"
  nice -n "$NICE_LEVEL" bash -c '
    bedtools bamtobed -i "$1" \
      | awk -v OFS="\t" '"'"'{if ($6=="+"){print $1,$2,$3,$4,$5,$6} else if ($6=="-"){print $1,$2,$3+1,$4,$5,$6}}'"'"' \
      | LC_ALL=C sort -T "$2" -k1,1 \
      | bedtools genomecov -bg -5 -i stdin -g "$3" \
      | LC_ALL=C sort -T "$2" -k1,1 -k2,2n > "$4"
  ' _ "$merged_bam" "$tmp_dir" "$CHROM_SIZES" "$bedgraph_tmp"

  mv "$bedgraph_tmp" "$bedgraph"
  run_low bedGraphToBigWig "$bedgraph" "$CHROM_SIZES" "$bigwig"
}

prepare_nonpeaks() {
  local prep_dir="$1"
  local peaks_bed="$prep_dir/inputs/overlap.bed"
  local bg_dir="$prep_dir/background/tutorial_folds"
  local final_nonpeaks="$bg_dir/nonpeaks_tutorial_folds.bed"
  local flank_size=$((INPUTLEN / 2))

  mkdir -p "$bg_dir"

  if [[ -s "$final_nonpeaks" ]]; then
    log "skip existing nonpeaks: $final_nonpeaks"
    return 0
  fi

  if [[ ! -s "$bg_dir/exclude.bed" ]]; then
    nice -n "$NICE_LEVEL" bash -c '
      bedtools slop -i "$1" -g "$2" -b "$3" | cut -f1-3 > "$4"
    ' _ "$BLACKLIST_BED" "$CHROM_SIZES" "$flank_size" "$bg_dir/blacklist_slop.bed"

    nice -n "$NICE_LEVEL" bash -c '
      bedtools slop -i "$1" -g "$2" -b "$3" | cut -f1-3 > "$4"
    ' _ "$peaks_bed" "$CHROM_SIZES" "$flank_size" "$bg_dir/peaks_slop.bed"

    cat "$bg_dir/blacklist_slop.bed" "$bg_dir/peaks_slop.bed" \
      | bedtools sort \
      | bedtools merge -i stdin > "$bg_dir/exclude.bed.tmp"

    mv "$bg_dir/exclude.bed.tmp" "$bg_dir/exclude.bed"
    rm -f "$bg_dir/blacklist_slop.bed" "$bg_dir/peaks_slop.bed"
  else
    log "skip existing exclude bed: $bg_dir/exclude.bed"
  fi

  if [[ ! -s "$bg_dir/foreground.gc.bed" ]]; then
    run_low "$PYTHON_BIN" "$GC_CONTENT_PY" \
      -i "$peaks_bed" \
      -c "$CHROM_SIZES" \
      -g "$GENOME" \
      -op "$bg_dir/foreground.gc" \
      -il "$INPUTLEN"
  else
    log "skip existing foreground gc bed: $bg_dir/foreground.gc.bed"
  fi

  if [[ ! -s "$bg_dir/candidate.negatives.bed" ]]; then
    nice -n "$NICE_LEVEL" bash -c '
      bedtools intersect -v -a "$1" -b "$2" > "$3"
    ' _ "$GC_BED" "$bg_dir/exclude.bed" "$bg_dir/candidate.negatives.bed.tmp"
    mv "$bg_dir/candidate.negatives.bed.tmp" "$bg_dir/candidate.negatives.bed"
  else
    log "skip existing candidate negatives: $bg_dir/candidate.negatives.bed"
  fi

  if [[ ! -s "$bg_dir/negatives.bed" ]]; then
    run_low "$PYTHON_BIN" "$GC_MATCH_PY" \
      -c "$bg_dir/candidate.negatives.bed" \
      -f "$bg_dir/foreground.gc.bed" \
      -o "$bg_dir/negatives" \
      -fl "$FOLDS_JSON" \
      -npr 2
  else
    log "skip existing negatives: $bg_dir/negatives.bed"
  fi

  awk -v OFS="\t" -v summit="$flank_size" \
    '{print $1, $2, $3, ".", ".", ".", ".", ".", ".", summit}' \
    "$bg_dir/negatives.bed" > "${final_nonpeaks}.tmp"
  mv "${final_nonpeaks}.tmp" "$final_nonpeaks"
}

prepare_dataset() {
  local dataset="$1"
  local dataset_dir="$DATA_ROOT/$dataset"
  local prep_dir="$dataset_dir/prep_v1"

  [[ -d "$dataset_dir" ]] || fail "missing dataset dir: $dataset_dir"

  log "========== dataset: $dataset =========="
  mkdir -p "$prep_dir/inputs"

  prepare_dataset_links "$dataset_dir"
  prepare_peaks_bed "$dataset_dir" "$prep_dir"
  prepare_merged_bam "$dataset_dir"
  prepare_bigwig "$dataset_dir" "$prep_dir"
  prepare_nonpeaks "$prep_dir"

  ls -lh \
    "$dataset_dir/merged.bam" \
    "$dataset_dir/merged.bam.bai" \
    "$prep_dir/bigwig/merged_unstranded.bw" \
    "$prep_dir/background/tutorial_folds/nonpeaks_tutorial_folds.bed"
}

main() {
  local dataset_spec
  local dataset
  local dataset_list=()

  prepare_references

  dataset_spec="${DATASETS// /,}"
  IFS=',' read -r -a dataset_list <<< "$dataset_spec"

  for dataset in "${dataset_list[@]}"; do
    [[ -n "$dataset" ]] || continue
    prepare_dataset "$dataset"
  done

  log "all requested datasets prepared"
}

main "$@"
