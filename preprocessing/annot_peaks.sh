#!/bin/bash

WORK_DIR="$1"

shift;

for file in $(ls $WORK_DIR/peaks/reps/diff/*.bed | grep BF); do
	filename=$(basename $file)
	prefix=$(echo $filename | sed 's/\.bed//g')
	echo "$prefix is being annotated with homer"
	annotatePeaks.pl $file hg38 -cpu 6 -annStats $WORK_DIR/peaks/reps/stats/$prefix.stats.txt > $WORK_DIR/annotation/reps/$prefix.annotated.bed
done
