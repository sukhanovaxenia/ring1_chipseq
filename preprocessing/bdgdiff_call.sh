#!/bin/bash

WORK_DIR="$1"

shift #remove WORK_DIR from the argument list

for i in $(ls $WORK_DIR/peaks/reps/*.bdg | grep 'treat' | grep 'HS' )
do
	filename=$(basename $i)
	prefix=$(echo $filename | cut -d_ -f1)
	suffix=$(echo $filename | cut -d. -f1 | cut -d_ -f3)

	echo "Compare states for ${prefix}_${suffix}"
	macs2 bdgdiff --t1 $WORK_DIR/peaks/reps/${prefix}_BF_${suffix}.dedup_treat_pileup.bdg --t2 $WORK_DIR/peaks/reps/${prefix}_HS_${suffix}_treat_pileup.bdg \
		--c1 $WORK_DIR/peaks/reps/${prefix}_BF_${suffix}.dedup_control_lambda.bdg --c2 $WORK_DIR/peaks/reps/${prefix}_HS_${suffix}_control_lambda.bdg \
		-o $WORK_DIR/peaks/reps/diff/${prefix}_${suffix}_BF_unique.bed $WORK_DIR/peaks/reps/diff/${prefix}_${suffix}_HS_unique.bed $WORK_DIR/peaks/reps/diff/${prefix}_${suffix}_BF_HS_common.bed \
		--d1 41338352 --d2 56570853 > $WORK_DIR/logs/${prefix}_${suffix}.bdgdiff.log 2>&1
done