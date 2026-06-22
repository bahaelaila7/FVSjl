#!/usr/bin/env bash
# mask_natural.sh — produce NATURAL-ONLY copies of every scenario by stripping all
# management/disturbance keywords (cuts, fire/FFE block, pests). The remaining run is
# pure natural dynamics (growth/mortality/density/regen) for the congruence sweep.
# Usage: mask_natural.sh <scenario-dir> <out-dir>
set -euo pipefail
SRC="${1:-$(dirname "$0")/scenarios}"
OUT="${2:-/tmp/nat_keys}"
rm -rf "$OUT"; mkdir -p "$OUT"

# cut/management single-line keywords to drop
CUTKW='THINBTA|THINATA|THINBBA|THINABA|THINDBH|THINPRSC|THINSDI|THINHT|THINCC|THINAUTO|THINMIST|THINRDEN|SETPTHIN|THINPT|SALVAGE|YARDLOSS|SPECPREF|CUTEFF|TCONDMLT|SPGROUP|SPLEAVE|LEAVESP|PRUNE'
# extension blocks opened by these and closed by END (fire/FFE, etc.)
EXTKW='FMIN|ESTB|FFE'

for k in "$SRC"/*.key; do
  n="$(basename "$k" .key)"
  awk -v cut="^($CUTKW)" -v ext="^($EXTKW)" '
    BEGIN{skip=0}
    skip==1 && $0 ~ /^END/ {skip=0; next}
    skip==1 {next}
    $0 ~ ext {skip=1; next}
    $0 ~ cut {next}
    {print}
  ' "$k" > "$OUT/$n.key"
  [ -f "$SRC/$n.tre" ] && cp "$SRC/$n.tre" "$OUT/$n.tre"
done
echo "masked $(ls "$OUT"/*.key | wc -l) keys → $OUT"
