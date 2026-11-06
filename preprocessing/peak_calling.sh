#!/bin/bash

WORK_DIR="$1"
CTRL_SUFFIX="$2"
shift #remove WORK_DIR from the argument list

#calling peaks for BF
for i in $(ls $WORK_DIR/align/chip/*.dedup.bam | grep 'rep1' | grep -E -v 'input|control' | grep 'HS')
do
	filename=$(basename $i)
	prefix="${filename%_rep1.dedup.bam}"
	# control=$(echo "$prefix" | sed -E 's/1_(BF|HS)_rep([12])/_input_\1_rep\2/')
	if [[ "$prefix" == *"BF"* ]]; then
		control=$(echo $prefix | sed -E 's/1_(BF|HS)/1_'$CTRL_SUFFIX'\_\1/g')
	else
		control=$(echo $prefix | sed -E 's/1_(BF|HS)/_'$CTRL_SUFFIX'\_\1/g')
	fi

	echo "Calling united peaks for $prefix against $control"

	macs2 callpeak --bdg --outdir $WORK_DIR/peaks/ -n $prefix -f BAM -g hs \
		-t $WORK_DIR/align/chip/${prefix}_rep1.dedup.bam $WORK_DIR/align/chip/${prefix}_rep2.dedup.bam \
		-c $WORK_DIR/align/chip/${control}_rep1.dedup.bam $WORK_DIR/align/chip/${control}_rep2.dedup.bam \
		--cutoff-analysis > $WORK_DIR/logs/${prefix}_rep1.macs2_norep.log 2>&1
done
#calling peaks for HS
# macs2 callpeak --bdg --outdir peaks/ -n RING1_HS -f BAM -g hs \
# 	-t align/chip/SRR1055707_GSM1295091_RING1_HS.dedup.bam align/chip/SRR1055706_GSM1295091_RING1_HS.dedup.bam \
# 	-c align/chip/SRR1055711_GSM1295093_RING_input_HS.dedup.bam align/chip/SRR1055710_GSM1295093_RING_input_HS.dedup.bam \
# 	--cutoff-analysis 2> HS1_macs2_norep.log
