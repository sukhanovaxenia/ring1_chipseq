#!/usr/bin/env bash
# MACS2 peak calling on pooled replicates of each condition (with -B bedGraph output).
# Usage: bash peak_calling.sh <project_dir> [condition ...]   (default: all CONDITIONS)
# In:  align/chip/dedup/<SAMPLE_PREFIX>_<cond>_<rep>.dedup.bam and the matching controls
# Out: peaks/<SAMPLE_PREFIX>_<cond>_{peaks.narrowPeak,peaks.xls,treat_pileup.bdg,control_lambda.bdg,...}
source "$(dirname "$0")/common.sh"
init_work_dir "${1:-}"; shift
pick_conditions "$@"
require_tools macs2
mkdir -p "$WORK_DIR/peaks"
bam_dir="$WORK_DIR/align/chip/dedup"

for cond in "${SELECTED[@]}"; do
  name="${SAMPLE_PREFIX}_${cond}"
  chip=(); ctrl=()
  for rep in "${REPLICATES[@]}"; do
    chip+=("$bam_dir/$(chip_name "$cond" "$rep").dedup.bam")
    ctrl+=("$bam_dir/$(control_name "$cond" "$rep").dedup.bam")
  done
  require_files "${chip[@]}" "${ctrl[@]}"
  echo "$name: pooled peak calling"
  echo "  ChIP:    ${chip[*]##*/}"
  echo "  control: ${ctrl[*]##*/}"
  is_done "$WORK_DIR/peaks/${name}_peaks.narrowPeak" && continue
  macs2 callpeak -t "${chip[@]}" -c "${ctrl[@]}" -f BAM -g hs -n "$name" \
    --outdir "$WORK_DIR/peaks" -B --cutoff-analysis \
    > "$LOGS/$name.macs2.log" 2>&1
done
