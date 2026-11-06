#!/bin/bash

WORK_DIR="$1"

shift #remove WORK_DIR from the argument list


for file in $(ls $WORK_DIR/align/chip/*.bam); do
	# prefix=$(echo $file | sed 's/\.\///g' | sed 's/\.bam//g')
	filename=$(basename ${file})
	prefix="${filename%.bam}"
	
	echo "$prefix is being deduped"
	echo "Output dir is $WORKDIR/align/chip"
	picard MarkDuplicates --REMOVE_DUPLICATES true \
	--CREATE_INDEX true --ASSUME_SORT_ORDER coordinate \
	-I $file \
	-M $prefix.metric.txt \
	-O $WORK_DIR/align/chip/${prefix}.dedup.bam \
	> $WORK_DIR/logs/${prefix}.rmdup.log 2>&1
done
