#!/bin/bash

WORK_DIR="$1"

shift;

computeMatrix reference-point --referencePoint TSS \
-b 1000 -a 1000 \
-R $WORK_DIR/peaks/diff/bed3/RING1_BF_unique.bed \
-S $WORK_DIR/normalization/bin/RING1_BF_rep*.binned.TMM.bigwig \
--skipZeros \
-o $WORK_DIR/visualization/matrix/BF_TSS.gz \
-p 6 \
--outFileSortedRegions $WORK_DIR/visualization/regionsRING1BF_TSS.bed \
> $WORK_DIR/logs/RING1_BF.comp_matrix.log 2>&1

computeMatrix reference-point --referencePoint TSS \
-b 1000 -a 1000 \
-R $WORK_DIR/peaks/diff/bed3/RING1_HS_unique.bed \
-S $WORK_DIR/normalization/bin/RING1_HS_rep*.binned.TMM.bigwig \
--skipZeros \
-o $WORK_DIR/visualization/matrix/HS_TSS.gz \
-p 6 \
--outFileSortedRegions $WORK_DIR/visualization/regionsRING1HS_TSS.bed \
> $WORK_DIR/logs/RING1_HS.comp_matrix.log 2>&1

# now visualize

plotProfile -m $WORK_DIR/visualization/matrix/BF_TSS.gz \
-out $WORK_DIR/visualization/figures/TSS_BF_profile.png \
--perGroup \
--colors green purple \
--plotTitle "" --samplesLabel "Rep1" "Rep2" \
--refPointLabel "TSS" \
-T "RING1 BF read density" \
-z "" \
> $WORK_DIR/logs/TSS_BF_profile.log 2>&1

plotProfile -m $WORK_DIR/visualization/matrix/HS_TSS.gz \
-out $WORK_DIR/visualization/figures/TSS_HS_profile.png \
--perGroup \
--colors green purple \
--plotTitle "" --samplesLabel "Rep1" "Rep2" \
--refPointLabel "TSS" \
-T "RING1 HS read density" \
-z "" \
> $WORK_DIR/logs/TSS_HF_profile.log 2>&1

# heatmap with the same info

plotHeatmap -m $WORK_DIR/visualization/matrix/BF_TSS.gz \
-out $WORK_DIR/visualization/figures/TSS_BF_profile-heatmap.png \
--colorMap RdBu \
--zMin -2 --zMax 2 \
> $WORK_DIR/logs/TSS_BF_profile-heatmap.log 2>&1

plotHeatmap -m $WORK_DIR/visualization/matrix/HS_TSS.gz \
-out $WORK_DIR/visualization/figures/TSS_HS_profile-heatmap.png \
--colorMap RdBu \
--zMin -2 --zMax 2 \
> $WORK_DIR/logs/TSS_HS_profile-heatmap.log 2>&1