#!/usr/bin/env bash
# Adapter and quality trimming of single-end ChIP-seq reads with Trimmomatic.
# Usage: bash trimming.sh <project_dir>
# In:  <project_dir>/raw/chip/*.fastq.gz (or *.fq.gz)
# Out: <project_dir>/trim/<sample>.trim.fastq.gz, logs/<sample>.trim.log
source "$(dirname "$0")/common.sh"
init_work_dir "${1:-}"
require_vars ADAPTERS
require_files "$ADAPTERS"
require_tools trimmomatic
mkdir -p "$WORK_DIR/trim"

shopt -s nullglob
raw=("$WORK_DIR"/raw/chip/*.fastq.gz "$WORK_DIR"/raw/chip/*.fq.gz)
((${#raw[@]})) || die "no FASTQ files in $WORK_DIR/raw/chip"

for file in "${raw[@]}"; do
  name=$(basename "$file")
  name=${name%.fastq.gz}; name=${name%.fq.gz}
  sample=${name%"${RAW_CHIP_SUFFIX:-}"}
  out="$WORK_DIR/trim/$sample.trim.fastq.gz"
  echo "$sample: trimming"
  is_done "$out" && continue
  trimmomatic SE -threads "$THREADS" "$file" "$out" \
    "ILLUMINACLIP:$ADAPTERS:2:30:10" AVGQUAL:28 TRAILING:26 MINLEN:25 \
    > "$LOGS/$sample.trim.log" 2>&1
done
