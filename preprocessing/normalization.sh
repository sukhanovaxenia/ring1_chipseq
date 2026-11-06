#!/bin/bash

WORK_DIR="$1"
NORM="$2"

shift #remove WORK_DIR from the argument list

for i in $WORK_DIR/align/chip/dedup/*dedup.bam
do
	## extract norm.factor from the file we created in R for the current sample:
	NormFactor=$(grep "${i}" $NORM | cut -f2)

	# prefix=$(echo $i | sed 's/align\/chip\///g' | sed 's/\.dedup\.bam//g') would leave WORK_DIR prefix unparsed
    filename=$(basename $i)

	prefix="${filename%.dedup.bam}"
	
	echo "$prefix is being normalized" 

	## use bamCoverage:
	bamCoverage -p 8 -b $i -o $WORK_DIR/normalization/$prefix.TMM.bigwig -bs 1 --scaleFactor ${NormFactor} > $WORK_DIR/logs/$prefix.norm.log 2>&1

done
