
WORK_DIR="$1"

shift;

# for file in $WORK_DIR/normalization/bin/*.binned.bigwig; do
# 	filename=$(basefile $file)
# 	prefix=$(echo $filename | sed s'/\.binned\.bigwig//g')
# 	annotatePeaks.pl $WORK_DIR/peaks/$prefix.wig hg38 > $WORK_DIR/stats/norm/bin/$prefix.stats.txt
# done

# for file in $WORK_DIR/normalization/*.binned.bigwig; do
# 	filename=$(basefile $file)
# 	prefix=$(echo $filename | sed s'/\.bigwig//g')
# 	annotatePeaks.pl $WORK_DIR/peaks/$prefix.wig hg38 > $WORK_DIR/stats/norm/$prefix.stats.txt
# done

for file in $WORK_DIR/peaks/diff/*.bed; do
	filename=$(basefile $file)
	prefix=$(echo $filename | sed s'/\.bed//g')
	annotatePeaks.pl $file hg38 -cpu 6 -annStats $WORK_DIR/peaks/stats/$prefix.stats.txt > $WORK_DIR/annotation/$prefix.annotated.bed
done