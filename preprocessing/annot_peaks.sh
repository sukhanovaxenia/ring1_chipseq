#!/usr/bin/env bash
# Annotate bdgdiff regions with HOMER annotatePeaks.pl.
# Usage: bash annot_peaks.sh <project_dir> pooled|reps [name_pattern]
#   pooled: peaks/diff/*.bed      -> annotation/<name>.annotated.bed,      peaks/stats/<name>.stats.txt
#   reps:   peaks/reps/diff/*.bed -> annotation/reps/<name>.annotated.bed, peaks/reps/stats/<name>.stats.txt
#   name_pattern (optional): annotate only files whose name contains it, e.g. BF_unique
source "$(dirname "$0")/common.sh"
init_work_dir "${1:-}"
mode="${2:?second argument must be pooled or reps}"
pattern="${3:-}"
require_vars HOMER_GENOME
require_tools annotatePeaks.pl
case "$mode" in
  pooled) in_dir="$WORK_DIR/peaks/diff";      out_dir="$WORK_DIR/annotation";      stats="$WORK_DIR/peaks/stats" ;;
  reps)   in_dir="$WORK_DIR/peaks/reps/diff"; out_dir="$WORK_DIR/annotation/reps"; stats="$WORK_DIR/peaks/reps/stats" ;;
  *) die "mode must be pooled or reps, got: $mode" ;;
esac
mkdir -p "$out_dir" "$stats"

shopt -s nullglob
n=0
for bed in "$in_dir"/*.bed; do
  name=$(basename "$bed" .bed)
  [[ -z "$pattern" || "$name" == *"$pattern"* ]] || continue
  n=$((n + 1))
  out="$out_dir/$name.annotated.bed"
  echo "$name: HOMER annotation"
  is_done "$out" && continue
  annotatePeaks.pl "$bed" "$HOMER_GENOME" -cpu "$THREADS" -annStats "$stats/$name.stats.txt" \
    > "$out" 2> "$LOGS/$name.homer.log"
done
((n)) || die "no BED files matched in $in_dir"
