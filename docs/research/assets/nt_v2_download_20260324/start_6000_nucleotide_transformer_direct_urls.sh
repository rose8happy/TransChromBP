#!/usr/bin/env bash
set -euo pipefail

ROOT=/data1/zhoujiazhen/bylw_atac
MODEL_DIR="$ROOT/foundation_models/nucleotide_transformer/nt-v2-500m-multi-species"
MANIFEST="$ROOT/foundation_models/logs/nt_v2_500m_resolved_urls.tsv"

mkdir -p "$MODEL_DIR"

if [[ ! -f "$MANIFEST" ]]; then
	echo "missing manifest: $MANIFEST" >&2
	exit 1
fi

download_one() {
	local rel="$1"
	local url="$2"
	local out="$MODEL_DIR/$rel"
	local part="${out}.part"

	mkdir -p "$(dirname "$out")"

	if [[ -s "$out" ]]; then
		echo "skip existing $rel"
		return 0
	fi

	echo "[$(date '+%F %T')] downloading $rel"
	wget --continue --tries=5 --output-document "$part" "$url"
	mv "$part" "$out"
	du -sh "$out"
}

while IFS=$'\t' read -r rel url; do
	case "$rel" in
		pytorch_model.bin|jax_model/pytree_ckpt.joblib)
			download_one "$rel" "$url"
			;;
		*)
			;;
	esac
done < "$MANIFEST"

echo "[$(date '+%F %T')] direct-url download finished"
