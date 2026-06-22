#!/usr/bin/env bash
# gen_natural_coverage.sh — natural-process COVERAGE scenarios that the dense snt01-
# stocking stands never exercise: understocked stands that reach the low-density /
# no-mortality (MORTS t<=t55d10) and sparse-crown (CROWN relsdi<=1) branches. Built by
# subsampling snt01.tre. (Coverage confirmed via the MCOV branch probe; both are
# bit-congruent to Oracle A — see NATURAL_SWEEP.md.)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; OUT="$HERE/scenarios"
SRC="${FVS_SN_KEY:-/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key}"
TRE="${SRC%.key}.tre"
base="$OUT/_natbase.key"; awk '/^PROCESS/{print "PROCESS"; print "STOP"; exit} /NOAUTOES/{print} {if(!/^PROCESS/)print}' "$SRC" > "$base" 2>/dev/null || awk '/^PROCESS/{print "PROCESS"; print "STOP"; exit}{print}' "$SRC" > "$base"
cp "$base" "$OUT/sparse_lowdens.key"; awk 'NR%6==1' "$TRE" > "$OUT/sparse_lowdens.tre"   # low density → background/no-mortality
cp "$base" "$OUT/sparse_min.key";     awk 'NR==1||NR==15' "$TRE" > "$OUT/sparse_min.tre" # very sparse → CROWN relsdi<=1
rm -f "$base"
echo "natural-coverage scenarios: sparse_lowdens, sparse_min"
