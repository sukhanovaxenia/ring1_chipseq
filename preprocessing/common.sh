# shellcheck shell=bash
# Shared settings and helpers for the preprocessing scripts. Source it, do not run it.
# Requires bash >= 4 (associative arrays).
set -euo pipefail

PRE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
if [[ ! -f "$PRE_DIR/paths.sh" ]]; then
  echo "Create $PRE_DIR/paths.sh from paths.example.sh first" >&2
  exit 1
fi
# shellcheck source=paths.example.sh
source "$PRE_DIR/paths.sh"
THREADS=${THREADS:-8}
FORCE=${FORCE:-0}            # FORCE=1 recomputes outputs that already exist

die() { echo "ERROR: $*" >&2; exit 1; }

# Abort unless every named variable is set and non-empty.
require_vars() {
  local v
  for v in "$@"; do [[ -n "${!v:-}" ]] || die "variable $v is not set (see paths.sh)"; done
}

# Abort unless every file exists and is non-empty.
require_files() {
  local f
  for f in "$@"; do [[ -s "$f" ]] || die "missing or empty file: $f"; done
}

# Abort unless every command is on PATH.
require_tools() {
  local t
  for t in "$@"; do command -v "$t" >/dev/null || die "$t not found in PATH"; done
}

# True if the output exists and FORCE is not set; prints a note.
is_done() {
  if [[ -s "$1" && "$FORCE" != 1 ]]; then
    echo "  skip: $1 exists (FORCE=1 to redo)"
    return 0
  fi
  return 1
}

# First positional argument: project directory; creates the logs folder.
init_work_dir() {
  WORK_DIR="${1:?first argument must be the project directory}"
  [[ -d "$WORK_DIR" ]] || die "not a directory: $WORK_DIR"
  WORK_DIR=$(cd "$WORK_DIR" && pwd)
  LOGS="$WORK_DIR/logs"
  mkdir -p "$LOGS"
}

# Conditions to process: those given on the command line, otherwise all from paths.sh.
pick_conditions() {
  if (($#)); then SELECTED=("$@"); else SELECTED=("${CONDITIONS[@]}"); fi
  local c
  for c in "${SELECTED[@]}"; do
    [[ -n "${CONTROL_PREFIX[$c]:-}" ]] || die "no CONTROL_PREFIX for condition $c"
  done
}

chip_name()    { echo "${SAMPLE_PREFIX}_$1_$2"; }        # condition, replicate
control_name() { echo "${CONTROL_PREFIX[$1]}_$2"; }      # condition, replicate

# Sequencing depth for bdgdiff: the smaller of treatment and control tag counts
# after filtering, as reported by MACS2 in *_peaks.xls.
macs_depth() {
  local xls="$1" t c
  t=$(grep -m1 '# tags after filtering in treatment' "$xls" | awk -F': ' '{print $2}')
  c=$(grep -m1 '# tags after filtering in control' "$xls" | awk -F': ' '{print $2}')
  [[ "$t" =~ ^[0-9]+$ && "$c" =~ ^[0-9]+$ ]] || die "cannot read tag counts from $xls"
  echo $(( t < c ? t : c ))
}

# Strip bdgdiff/MACS2 BED files to BED3 (drops track lines), for deepTools.
to_bed3() {
  awk 'BEGIN{OFS="\t"} !/^(track|browser|#)/ && NF>=3 {print $1,$2,$3}' "$1" > "$2"
}
