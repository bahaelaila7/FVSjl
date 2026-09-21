#!/usr/bin/env bash
# regime_sweep.sh — run the full regime + extension matrix for one variant, then triage + sign-tally each.
#
# Usage: regime_sweep.sh <VARIANT> <none_list> <regime_list> <outdir> [maxjobs]
#   none      runs on <none_list> (the big stratified sample, e.g. extract_sample.jl <V> 400)
#   every other regime runs on <regime_list> (a smaller stratified sample, e.g. extract_sample.jl <V> 120)
# Writes <outdir>/<regime>.csv (ledger), <outdir>/<regime>.log, <outdir>/<regime>.tally (signed population
# tally) and <outdir>/<regime>.triage. ⚠ Warm the FVSjl precompile cache BEFORE calling this (e.g. run the
# gate first): concurrent jobs that all recompile race the shared cache and corrupt it.
set -uo pipefail
V=$1; NONE_LIST=$2; REG_LIST=$3; OUT=$4; MAXJ=${5:-5}
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../../.." && pwd)
mkdir -p "$OUT"
REGIMES=(none thinbba salvage plant simfire mistletoe climate rootdis econ cover)

run_one() {
  local r=$1 list=$2
  ( cd "$ROOT" && LEDGER="$OUT/$r.csv" julia --project=. test/harness/fia/ledger_fia.jl "$list" "$V" "$r" \
      > "$OUT/$r.log" 2>&1
    bash test/harness/fia/triage_ledger.sh "$OUT/$r.csv" > "$OUT/$r.triage" 2>&1
    VARIANT=$V REGIME=$r julia --project=. test/harness/fia/sign_tally.jl "$list" > "$OUT/$r.tally" 2>&1
    echo "$r DONE $(date +%H:%M:%S)" >> "$OUT/_progress" )
}

for r in "${REGIMES[@]}"; do
  [ -f "$OUT/$r.tally" ] && { echo "skip $r (done)"; continue; }
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJ" ]; do sleep 5; done
  list=$REG_LIST; [ "$r" = none ] && list=$NONE_LIST
  echo "$r START $(date +%H:%M:%S)" >> "$OUT/_progress"
  run_one "$r" "$list" &
done
wait
echo "ALL DONE $(date +%H:%M:%S)" >> "$OUT/_progress"
