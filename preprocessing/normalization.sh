#!/usr/bin/env bash
# Scaled coverage tracks (bigWig) from deduplicated BAMs with deepTools bamCoverage.
# Usage: bash normalization.sh <project_dir> <scale_table> [bin_size]   (default bin_size: 10)
# scale_table: tab-separated, column 1 = BAM path or file name, column 2 = scale factor
#   passed to --scaleFactor as is (raw counts are multiplied by it, e.g. 1e6 / (library size x TMM)).
#   Lines starting with '#' and a header whose 2nd column is not numeric are ignored.
# Out: normalization/bs<bin_size>/<sample>.TMM.bigwig
source "$(dirname "$0")/common.sh"
init_work_dir "${1:-}"
table="${2:?second argument must be the scale-factor table}"
bin="${3:-10}"
require_files "$table"
require_tools bamCoverage
out_dir="$WORK_DIR/normalization/bs$bin"
mkdir -p "$out_dir"

shopt -s nullglob
bams=("$WORK_DIR"/align/chip/dedup/*.dedup.bam)
((${#bams[@]})) || die "no BAM files in $WORK_DIR/align/chip/dedup"

for bam in "${bams[@]}"; do
  file=$(basename "$bam")
  sample=${file%.dedup.bam}
  # Exact match on the file name, so RING1_BF_rep1 never matches RING1_BF_rep10 or a control.
  scale=$(awk -F'\t' -v f="$file" '!/^#/ {n=split($1,p,"/"); if (p[n]==f) {print $2; exit}}' "$table")
  [[ "$scale" =~ ^[0-9.eE+-]+$ ]] || die "no numeric scale factor for $file in $table"
  out="$out_dir/$sample.TMM.bigwig"
  echo "$sample: bamCoverage, bin $bin, scale $scale"
  is_done "$out" && continue
  bamCoverage -p "$THREADS" -b "$bam" -o "$out" --binSize "$bin" \
    --normalizeUsing None --scaleFactor "$scale" > "$LOGS/$sample.bamCoverage.log" 2>&1
done
