WORK_DIR="$1"

shift #remove WORK_DIR from the argument list

for i in $(ls $WORK_DIR/peaks/*.bdg | grep 'treat' | grep 'BF')
do
	filename=$(basename $i)
	prefix=$(echo $filename | cut -d_ -f1)

	echo "Compare states for ${prefix}"
	macs2 bdgdiff --t1 $WORK_DIR/peaks/${prefix}_BF_treat_pileup.bdg --t2 $WORK_DIR/peaks/${prefix}_HS_treat_pileup.bdg \
		--c1 $WORK_DIR/peaks/${prefix}_BF_control_lambda.bdg --c2 $WORK_DIR/peaks/${prefix}_HS_control_lambda.bdg \
		-o $WORK_DIR/peaks/diff/${prefix}_BF_unique.bed $WORK_DIR/peaks/diff/${prefix}_HS_unique.bed $WORK_DIR/peaks/diff/${prefix}_BF_HS_common.bed \
		--d1 41338352 --d2 56570853 > $WORK_DIR/logs/${prefix}.bdgdiff.log 2>&1
done