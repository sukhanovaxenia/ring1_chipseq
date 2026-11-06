#!/usr/bin/env bash
# Per-replicate differential binding (BF vs HS) with MACS2 bdgdiff.
# Usage: bash scripts/bdgdiff_reps.sh <project_root>
# Expects MACS2 callpeak -B outputs in <project_root>/peaks/reps/:
#   RING1_{BF,HS}_rep{1,2}_{treat_pileup.bdg,control_lambda.bdg,peaks.xls}
set -euo pipefail

WORK_DIR="${1:?Usage: bdgdiff_reps.sh <project_root>}"
REPS="$WORK_DIR/peaks/reps"
OUT="$REPS/diff"
mkdir -p "$OUT" "$WORK_DIR/logs"

# Sequencing depth for bdgdiff: the smaller of treatment and control
# tag counts after filtering, as reported by MACS2 in *_peaks.xls.
depth() {
  local xls="$1" t c
  t=$(grep -m1 '# tags after filtering in treatment' "$xls" | awk -F': ' '{print $2}')
  c=$(grep -m1 '# tags after filtering in control' "$xls" | awk -F': ' '{print $2}')
  echo $(( t < c ? t : c ))
}

for rep in rep1 rep2; do
  bf="$REPS/RING1_BF_${rep}"
  hs="$REPS/RING1_HS_${rep}"
  for f in "${bf}_treat_pileup.bdg" "${bf}_control_lambda.bdg" "${bf}_peaks.xls" \
           "${hs}_treat_pileup.bdg" "${hs}_control_lambda.bdg" "${hs}_peaks.xls"; do
    [[ -s "$f" ]] || { echo "Missing or empty: $f" >&2; exit 1; }
  done
  d1=$(depth "${bf}_peaks.xls")
  d2=$(depth "${hs}_peaks.xls")
  echo "${rep}: depth BF=${d1}, HS=${d2}"
  macs2 bdgdiff \
    --t1 "${bf}_treat_pileup.bdg" --c1 "${bf}_control_lambda.bdg" \
    --t2 "${hs}_treat_pileup.bdg" --c2 "${hs}_control_lambda.bdg" \
    --d1 "$d1" --d2 "$d2" -g 100 -l 200 --cutoff 3 \
    -o "$OUT/RING1_${rep}_BF_unique.bed" "$OUT/RING1_${rep}_HS_unique.bed" "$OUT/RING1_${rep}_BF_HS_common.bed" \
    > "$WORK_DIR/logs/RING1_${rep}.bdgdiff.log" 2>&1
done
