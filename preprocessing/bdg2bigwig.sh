WORK_DIR="$1"
shift #remove WORK_DIR from the argument list

#calling peaks for BF
for i in $(ls $WORK_DIR/peaks/reps/*.bdg)
do
	filename=$(basename $i)
	prefix="${filename%.bdg}"
	# another approach to filer input|control is
	# [[ "$prefix" =~ input|control ]] && continue
	echo $filename

    echo "Coverting ${prefix} BDG file..."

    sort -k1,1 -k2,2n $WORK_DIR/peaks/reps/$prefix.bdg > $WORK_DIR/peaks/reps/$prefix.sorted.bdg 2> $WORK_DIR/logs/$prefix.bdg_sort.log
    bedGraphToBigWig $WORK_DIR/peaks/reps/$prefix.sorted.bdg /home/intens2026/data/ref/hg38.chr.sizes $WORK_DIR/peaks/reps/$prefix.bw 2>&1 $WORK_DIR/logs/$prefix.bdg2bw.log
done