#!/usr/bin/env bash
# MACS2 peak calling for each replicate against its own control (with -B bedGraph output).
# Usage: bash peak_calling_reps.sh <project_dir> [condition ...]   (default: all CONDITIONS)
# In:  align/chip/dedup/<SAMPLE_PREFIX>_<cond>_<rep>.dedup.bam and the matching controls
# Out: peaks/reps/<SAMPLE_PREFIX>_<cond>_<rep>_{peaks.narrowPeak,peaks.xls,treat_pileup.bdg,control_lambda.bdg,...}
source "$(dirname "$0")/common.sh"
init_work_dir "${1:-}"; shift
pick_conditions "$@"
require_tools macs2
out_dir="$WORK_DIR/peaks/reps"
mkdir -p "$out_dir"
bam_dir="$WORK_DIR/align/chip/dedup"

for cond in "${SELECTED[@]}"; do
  for rep in "${REPLICATES[@]}"; do
    name=$(chip_name "$cond" "$rep")
    chip="$bam_dir/$name.dedup.bam"
    ctrl="$bam_dir/$(control_name "$cond" "$rep").dedup.bam"
    require_files "$chip" "$ctrl"
    echo "$name: peak calling against $(basename "$ctrl")"
    is_done "$out_dir/${name}_peaks.narrowPeak" && continue
    macs2 callpeak -t "$chip" -c "$ctrl" -f BAM -g hs -n "$name" \
      --outdir "$out_dir" -B > "$LOGS/$name.macs2.log" 2>&1
  done
done
