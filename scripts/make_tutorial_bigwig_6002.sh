#!/usr/bin/env bash
set -euo pipefail

ENV_DIR="${ENV_DIR:-/home/zhengwei/bylw_atac/.mamba/envs/transchrombp}"
ROOT_DIR="${ROOT_DIR:-/home/zhengwei/bylw_atac/TransChromBP}"
DATA_DIR="${DATA_DIR:-/home/zhengwei/bylw_atac/chrombpnet_tutorial/data}"
OUT_DIR="${OUT_DIR:-$ROOT_DIR/outputs/preprocessing/tutorial_canonical_v1/step2_bigwig}"
TMP_DIR="${TMP_DIR:-/home/zhengwei/bylw_atac/tmp_bigwig_sort}"

BAM="${BAM:-$DATA_DIR/merged.bam}"
CHROM="${CHROM:-$DATA_DIR/hg38.chrom.sizes}"
OUT_PREFIX="${OUT_PREFIX:-$OUT_DIR/merged}"
BEDGRAPH="${BEDGRAPH:-$OUT_PREFIX.unstranded.bedGraph}"

mkdir -p "$OUT_DIR" "$TMP_DIR"

export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"

samtools quickcheck "$BAM"

rm -f "$BEDGRAPH.tmp"

bedtools bamtobed -i "$BAM" \
  | awk -v OFS="\t" '{if ($6=="+"){print $1,$2,$3,$4,$5,$6} else if ($6=="-"){print $1,$2,$3+1,$4,$5,$6}}' \
  | LC_COLLATE=C sort -T "$TMP_DIR" -k1,1 \
  | bedtools genomecov -bg -5 -i stdin -g "$CHROM" \
  | LC_COLLATE=C sort -T "$TMP_DIR" -k1,1 -k2,2n > "$BEDGRAPH.tmp"

mv "$BEDGRAPH.tmp" "$BEDGRAPH"
bedGraphToBigWig "$BEDGRAPH" "$CHROM" "${OUT_PREFIX}_unstranded.bw"

ls -lh "$BEDGRAPH" "${OUT_PREFIX}_unstranded.bw"
