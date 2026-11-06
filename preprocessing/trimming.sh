#!/bin/bash

for file in $(ls raw/chip/*gz); do
	prefix=$(echo $file | sed 's/raw\/chip\///g' | sed 's/_ChipSeq_biological_replicate_1_and_2_Homo_sapiens_ChIP-Seq\.fastq\.gz//g')
	echo "$prefix is being trimmed"
	trimmomatic SE -threads 6 $file trimmed/$prefix.trimmed.fq.gz \
		ILLUMINACLIP:~/ref/Sequencing_adaptors.fasta:2:30:10 \
		AVGQUAL:28 TRAILING:26 MINLEN:25
done
