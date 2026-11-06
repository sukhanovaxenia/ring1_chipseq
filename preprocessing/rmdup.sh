#!/usr/bin/env bash
# Remove PCR duplicates with Picard MarkDuplicates.
# Usage: bash rmdup.sh <project_dir>
# In:  <project_dir>/align/chip/*.bam (already deduplicated *.dedup.bam are ignored)
# Out: align/chip/dedup/<sample>.dedup.bam(.bai), logs/<sample>.dup_metrics.txt, logs/<sample>.rmdup.log
source "$(dirname "$0")/common.sh"
init_work_dir "${1:-}"
require_tools picard samtools
out_dir="$WORK_DIR/align/chip/dedup"
mkdir -p "$out_dir"

shopt -s nullglob
bams=()
for f in "$WORK_DIR"/align/chip/*.bam; do [[ "$f" == *.dedup.bam ]] || bams+=("$f"); done
((${#bams[@]})) || die "no BAM files in $WORK_DIR/align/chip"

for bam in "${bams[@]}"; do
  sample=$(basename "$bam" .bam)
  out="$out_dir/$sample.dedup.bam"
  echo "$sample: removing duplicates"
  is_done "$out" && continue
  picard MarkDuplicates --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER coordinate \
    -I "$bam" -O "$out" -M "$LOGS/$sample.dup_metrics.txt" \
    > "$LOGS/$sample.rmdup.log" 2>&1
  # samtools index gives the <name>.bam.bai name that every downstream tool expects.
  samtools index -@ "$THREADS" "$out"
done
