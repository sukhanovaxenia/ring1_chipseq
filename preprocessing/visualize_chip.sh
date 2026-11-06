#!/usr/bin/env bash
# deepTools profiles and heatmaps of ChIP signal around the TSS of the pooled bdgdiff regions.
# Usage: bash visualize_chip.sh <project_dir> [bigwig_dir] [flank_bp]
#   bigwig_dir default: normalization/bs10 (output of normalization.sh); flank default: 1000
# In:  peaks/diff/bed3/<P>_<cond>_unique.bed (from bdgdiff_call.sh pooled),
#      <bigwig_dir>/<P>_<cond>_rep*.TMM.bigwig
# Out: visualization/matrix/<cond>_TSS.gz, visualization/regions<P><cond>_TSS.bed,
#      visualization/figures/TSS_<cond>_profile.png, TSS_<cond>_profile-heatmap.png
source "$(dirname "$0")/common.sh"
init_work_dir "${1:-}"
bw_dir="${2:-$WORK_DIR/normalization/bs10}"
flank="${3:-1000}"
require_tools computeMatrix plotProfile plotHeatmap
vis="$WORK_DIR/visualization"
mkdir -p "$vis/matrix" "$vis/figures"

for cond in "${CONDITIONS[@]}"; do
  regions="$WORK_DIR/peaks/diff/bed3/${SAMPLE_PREFIX}_${cond}_unique.bed"
  require_files "$regions"
  bws=(); labels=()
  for rep in "${REPLICATES[@]}"; do
    bws+=("$bw_dir/$(chip_name "$cond" "$rep").TMM.bigwig"); labels+=("$rep")
  done
  require_files "${bws[@]}"
  matrix="$vis/matrix/${cond}_TSS.gz"
  echo "$cond: matrix and plots (TSS +/- $flank bp)"

  computeMatrix reference-point --referencePoint TSS -b "$flank" -a "$flank" \
    -R "$regions" -S "${bws[@]}" --samplesLabel "${labels[@]}" --skipZeros -p "$THREADS" \
    -o "$matrix" --outFileSortedRegions "$vis/regions${SAMPLE_PREFIX}${cond}_TSS.bed" \
    > "$LOGS/${SAMPLE_PREFIX}_${cond}.computeMatrix.log" 2>&1

  plotProfile -m "$matrix" -out "$vis/figures/TSS_${cond}_profile.png" \
    --perGroup --colors green purple --refPointLabel TSS \
    --regionsLabel "${SAMPLE_PREFIX}_${cond}_unique" \
    --plotTitle "${SAMPLE_PREFIX} ${cond} read density" \
    > "$LOGS/TSS_${cond}_profile.log" 2>&1

  # Signal is non-negative: a sequential colour map starting at zero uses the whole scale.
  plotHeatmap -m "$matrix" -out "$vis/figures/TSS_${cond}_profile-heatmap.png" \
    --colorMap Reds --zMin 0 --refPointLabel TSS \
    --regionsLabel "${SAMPLE_PREFIX}_${cond}_unique" \
    > "$LOGS/TSS_${cond}_profile-heatmap.log" 2>&1
done
