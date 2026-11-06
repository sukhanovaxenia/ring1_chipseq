#!/usr/bin/env bash
# Convert MACS2 bedGraph files (pileup, lambda) to bigWig for genome browsers.
# Usage: bash bdg2bigwig.sh <project_dir> [bdg_dir]   (default bdg_dir: peaks/reps)
# Out: <bdg_dir>/<name>.bw, logs/<name>.bdg2bw.log
source "$(dirname "$0")/common.sh"
init_work_dir "${1:-}"
bdg_dir="${2:-$WORK_DIR/peaks/reps}"
require_vars CHROM_SIZES
require_files "$CHROM_SIZES"
require_tools bedClip bedGraphToBigWig

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

shopt -s nullglob
files=("$bdg_dir"/*.bdg)
((${#files[@]})) || die "no .bdg files in $bdg_dir"

for bdg in "${files[@]}"; do
  name=$(basename "$bdg" .bdg)
  [[ "$name" == *.sorted ]] && continue          # leftovers of the old version
  out="$bdg_dir/$name.bw"
  echo "$name: bedGraph -> bigWig"
  is_done "$out" && continue
  # MACS2 intervals can run past chromosome ends: clip before conversion.
  LC_ALL=C sort -k1,1 -k2,2n "$bdg" > "$tmp/$name.sorted.bdg"
  bedClip "$tmp/$name.sorted.bdg" "$CHROM_SIZES" "$tmp/$name.clip.bdg"
  bedGraphToBigWig "$tmp/$name.clip.bdg" "$CHROM_SIZES" "$out" > "$LOGS/$name.bdg2bw.log" 2>&1
  rm -f "$tmp/$name".*.bdg
done
