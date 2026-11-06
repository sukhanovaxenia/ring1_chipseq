#!/bin/bash

for file in $(ls raw/rna/*gz); do
	prefix=$(echo $file | sed 's/raw\/rna\///g' | sed 's/_RNASeq_biological_replicate_1_and_2_Homo_sapiens_RNA-Seq\.fastq\.gz//g')
	echo "$prefix is being aligned with STAR"
	STAR --runThreadN 12 --runMode alignReads --genomeDir ../ref/hg38/STAR_genome/ --quantMode GeneCounts \
        --outFileNamePrefix align/star/$prefix. --readFilesIn $file --readFilesCommand zcat \
	--outBAMsortingThreadN 10 --outSAMtype BAM SortedByCoordinate
done
