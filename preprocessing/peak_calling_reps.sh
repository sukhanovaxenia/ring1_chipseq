#!/bin/bash

WORK_DIR="$1"
CTRL_SUFFIX="$2"
shift #remove WORK_DIR from the argument list

#calling peaks for BF
for i in $(ls $WORK_DIR/align/chip/dedup/*.dedup.bam | grep 'HS' | grep -E -v 'control' | grep 'rep1')
do
	filename=$(basename $i)
	prefix="${filename%.dedup.bam}"
	# another approach to filer input|control is
	# [[ "$prefix" =~ input|control ]] && continue
	echo $filename
	control=$(echo "$prefix" | sed -E 's/_(BF|HS)_rep([12])/_'${CTRL_SUFFIX}'_\1_rep\2/')

	echo "Calling ${prefix} against ${control}"

	macs2 callpeak --bdg --outdir $WORK_DIR/peaks/reps/ -n $prefix -f BAM -g hs \
		-t $WORK_DIR/align/chip/dedup/$prefix.dedup.bam \
		-c $WORK_DIR/align/chip/dedup/$control.dedup.bam > $WORK_DIR/logs/$prefix.macs2.log 2>&1
done
	# macs2 callpeak --bdg --outdir peaks/reps/ -n RING1_BF_rep2 -f BAM -g hs \
	# 	-t align/chip/RING1_BF.dedup.bam \
	# 	-c align/chip/RING1_input_BF.dedup.bam 2> logs/BF2_macs2.log

	#calling peaks for HS
	# macs2 callpeak --bdg --outdir peaks/reps/ -n RING1_HS_rep1 -f BAM -g hs \
	# 	-t align/chip/RING1_HS.dedup.bam \
	# 	-c align/chip/RING_input_HS.dedup.bam 2> logs/HS1_macs2.log

	# macs2 callpeak --bdg --outdir peaks/reps -n RING1_HS_rep2 -f BAM -g hs \
	# 		-t align/chip/RING1_HS.dedup.bam \
	# 		-c align/chip/RING_input_HS.dedup.bam 2> logs/HS2_macs2.log
