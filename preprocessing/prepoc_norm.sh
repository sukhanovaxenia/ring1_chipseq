#!/bin/bash

# preprocessing raw chip counts for edgeR scaling

bedtools

## BED (0-based) to SAF (1-based):
awk 'OFS="\t" {print $1"_"$2+1"_"$3, $1, $2+1, $3, "+"}' peaks/reference_peaks.bed > reference_peaks.saf

## quantify bam files over peaks:
featureCounts -a reference_peaks.saf -F SAF -o reference_peaks.countmatrix align/chip/*.bam
