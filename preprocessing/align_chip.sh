#!/bin/bash

WORK_DIR="$1"

shift #remove WORK_DIR from the argument list

for file in $(ls $WORK_DIR/trim/*trim.fastq.gz); do
	# prefix=$(echo $file | sed 's/test\///g' | sed 's/\.trim\.fastq\.gz//g')
	filename=$(basename $file)
	prefix="${filename%.trim.fastq.gz}"
	
	echo "$prefix is being aligned"
	bowtie2 -p 8 --no-unal -t --quiet --very-sensitive -x /home/intens2026/data/ref/hg38 -U $WORK_DIR/trim/$prefix.trim.fastq.gz | samtools view -@ 8 -F 4 -Sb | samtools sort -@ 8 -o $WORK_DIR/align/chip/$prefix.bam > $WORK_DIR/logs/$prefix.align.log 2>&1
	samtools index $WORK_DIR/align/chip/$prefix.bam
done
