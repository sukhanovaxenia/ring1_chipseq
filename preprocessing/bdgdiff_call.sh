#!/usr/bin/env bash
# Differential binding between the two conditions with MACS2 bdgdiff.
# Usage: bash bdgdiff_call.sh <project_dir> pooled|reps
#   pooled: peaks/<P>_<c1|c2>_*        -> peaks/diff/<P>_{c1_unique,c2_unique,c1_c2_common}.bed
#   reps:   peaks/reps/<P>_<c>_<rep>_*  -> peaks/reps/diff/<P>_<rep>_{...}.bed (rep N vs rep N)
# Depths (--d1/--d2) are read from each *_peaks.xls: min(treatment, control) tags after filtering.
# BED3 copies for deepTools are written to <diff dir>/bed3/.
source "$(dirname "$0")/common.sh"
init_work_dir "${1:-}"
mode="${2:?second argument must be pooled or reps}"
require_tools macs2
((${#CONDITIONS[@]} == 2)) || die "bdgdiff needs exactly two CONDITIONS"
c1=${CONDITIONS[0]}; c2=${CONDITIONS[1]}

run_bdgdiff() {   # <prefix1> <prefix2> <output prefix> <out dir>
  local p1="$1" p2="$2" tag="$3" out="$4" d1 d2 f
  require_files "${p1}_treat_pileup.bdg" "${p1}_control_lambda.bdg" "${p1}_peaks.xls" \
                "${p2}_treat_pileup.bdg" "${p2}_control_lambda.bdg" "${p2}_peaks.xls"
  d1=$(macs_depth "${p1}_peaks.xls")
  d2=$(macs_depth "${p2}_peaks.xls")
  echo "$tag: $c1 vs $c2, depth $d1 / $d2"
  mkdir -p "$out/bed3"
  if ! is_done "$out/${tag}_${c1}_${c2}_common.bed"; then
    macs2 bdgdiff \
      --t1 "${p1}_treat_pileup.bdg" --c1 "${p1}_control_lambda.bdg" \
      --t2 "${p2}_treat_pileup.bdg" --c2 "${p2}_control_lambda.bdg" \
      --d1 "$d1" --d2 "$d2" -g 100 -l 200 --cutoff 3 \
      -o "$out/${tag}_${c1}_unique.bed" "$out/${tag}_${c2}_unique.bed" "$out/${tag}_${c1}_${c2}_common.bed" \
      > "$LOGS/$tag.bdgdiff.log" 2>&1
  fi
  for f in "$out/${tag}_${c1}_unique.bed" "$out/${tag}_${c2}_unique.bed" "$out/${tag}_${c1}_${c2}_common.bed"; do
    to_bed3 "$f" "$out/bed3/$(basename "$f")"
    printf '  %-40s %6d regions\n' "$(basename "$f")" "$(wc -l < "$out/bed3/$(basename "$f")")"
  done
}

case "$mode" in
  pooled)
    run_bdgdiff "$WORK_DIR/peaks/${SAMPLE_PREFIX}_$c1" "$WORK_DIR/peaks/${SAMPLE_PREFIX}_$c2" \
      "$SAMPLE_PREFIX" "$WORK_DIR/peaks/diff" ;;
  reps)
    for rep in "${REPLICATES[@]}"; do
      run_bdgdiff "$WORK_DIR/peaks/reps/$(chip_name "$c1" "$rep")" "$WORK_DIR/peaks/reps/$(chip_name "$c2" "$rep")" \
        "${SAMPLE_PREFIX}_$rep" "$WORK_DIR/peaks/reps/diff"
    done ;;
  *) die "mode must be pooled or reps, got: $mode" ;;
esac
