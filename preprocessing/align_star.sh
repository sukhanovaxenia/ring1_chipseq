#!/usr/bin/env bash
# Align RNA-seq reads with STAR and count reads per gene (for integration with ChIP-seq).
# Usage: bash align_star.sh <project_dir>
# In:  <project_dir>/raw/rna/*.fastq.gz (or *.fq.gz)
# Out: align/star/<sample>.Aligned.sortedByCoord.out.bam(.bai), <sample>.ReadsPerGene.out.tab,
#      STAR's own Log.final.out next to them
source "$(dirname "$0")/common.sh"
init_work_dir "${1:-}"
require_vars STAR_INDEX
[[ -d "$STAR_INDEX" ]] || die "STAR index not found: $STAR_INDEX"
require_tools STAR samtools
mkdir -p "$WORK_DIR/align/star"

shopt -s nullglob
raw=("$WORK_DIR"/raw/rna/*.fastq.gz "$WORK_DIR"/raw/rna/*.fq.gz)
((${#raw[@]})) || die "no FASTQ files in $WORK_DIR/raw/rna"

for file in "${raw[@]}"; do
  name=$(basename "$file")
  name=${name%.fastq.gz}; name=${name%.fq.gz}
  sample=${name%"${RAW_RNA_SUFFIX:-}"}
  prefix="$WORK_DIR/align/star/$sample."
  bam="${prefix}Aligned.sortedByCoord.out.bam"
  echo "$sample: STAR alignment"
  is_done "$bam" && continue
  STAR --runThreadN "$THREADS" --runMode alignReads --genomeDir "$STAR_INDEX" \
    --readFilesIn "$file" --readFilesCommand zcat --quantMode GeneCounts \
    --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN "$THREADS" \
    --outFileNamePrefix "$prefix" > "$LOGS/$sample.star.log" 2>&1
  samtools index "$bam"
  grep -m1 'Uniquely mapped reads %' "${prefix}Log.final.out" || true
done
