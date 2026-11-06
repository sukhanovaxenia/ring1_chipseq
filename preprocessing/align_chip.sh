#!/usr/bin/env bash
# Align trimmed single-end ChIP-seq reads with bowtie2 and sort/index the BAM.
# Usage: bash align_chip.sh <project_dir>
# In:  <project_dir>/trim/*.trim.fastq.gz
# Out: align/chip/<sample>.bam(.bai);
#      logs/<sample>.bowtie2.log (alignment summary and timing), logs/<sample>.flagstat.txt
source "$(dirname "$0")/common.sh"
init_work_dir "${1:-}"
require_vars BOWTIE2_INDEX
require_tools bowtie2 samtools
[[ -s "$BOWTIE2_INDEX.1.bt2" || -s "$BOWTIE2_INDEX.1.bt2l" ]] || die "bowtie2 index not found: $BOWTIE2_INDEX"
mkdir -p "$WORK_DIR/align/chip"

shopt -s nullglob
reads=("$WORK_DIR"/trim/*.trim.fastq.gz)
((${#reads[@]})) || die "no trimmed reads in $WORK_DIR/trim"

for file in "${reads[@]}"; do
  sample=$(basename "$file" .trim.fastq.gz)
  bam="$WORK_DIR/align/chip/$sample.bam"
  echo "$sample: aligning"
  is_done "$bam" && continue
  # bowtie2 writes its summary to stderr -> own log; --no-unal drops unaligned reads.
  bowtie2 -p "$THREADS" --very-sensitive --no-unal -t -x "$BOWTIE2_INDEX" -U "$file" \
      2> "$LOGS/$sample.bowtie2.log" \
    | samtools view -@ "$THREADS" -b - \
    | samtools sort -@ "$THREADS" -o "$bam" - 2> "$LOGS/$sample.sort.log"
  samtools index -@ "$THREADS" "$bam"
  samtools flagstat -@ "$THREADS" "$bam" > "$LOGS/$sample.flagstat.txt"
  grep -m1 'overall alignment rate' "$LOGS/$sample.bowtie2.log" || true
done
