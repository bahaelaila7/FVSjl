#!/usr/bin/env bash
# slow_tier.sh — SLOW TIER of the tiered suite (not in Pkg.test; run before a master merge and nightly).
#   1. GOLDEN DRIFT: regenerate every committed tiered fixture from the CURRENT (freshly relinked) live oracles into a
#      temp dir and diff against the committed goldens. Any difference = the oracle changed (relink, FVS source edit,
#      build-flag change) ⇒ investigate before trusting the fast tier. PROVENANCE.toml sha256 changes are reported.
#   2. POPULATION: run the regime matrix sweep (regime_sweep.sh: ledger + triage + signed tally per regime) for each
#      variant on its 400/120 stratified samples; print the tallies.
# Usage: test/harness/tiered/slow_tier.sh [VARIANT ...]   (default: every variant with a committed fixture)
# Env: SLOW_OUT (default /workspace/.tiered_slow/<date>), SLOW_SKIP_POP=1 (drift check only), MAXJ (sweep jobs, 5)
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../../.." && pwd)
FX=$ROOT/test/fixtures/tiered
OUT=${SLOW_OUT:-/workspace/.tiered_slow/$(date +%Y%m%d_%H%M%S)}; mkdir -p "$OUT"
VS=("$@"); [ ${#VS[@]} -eq 0 ] && VS=($(ls "$FX" | grep -v '\.' | tr a-z A-Z))
echo "slow tier → $OUT   variants: ${VS[*]}"
status=0
for V in "${VS[@]}"; do
  v=$(echo "$V" | tr A-Z a-z)
  K=$(grep -E '^k = ' "$FX/$v/PROVENANCE.toml" | awk '{print $3}')
  echo "== [$V] regenerate goldens (K=$K) from the current oracle"
  (cd "$ROOT" && julia --project=. test/harness/tiered/make_fixtures.jl "$V" "$K" "$OUT/regen" > "$OUT/regen_$v.log" 2>&1)
  # compare everything except the provenance timestamps/commits
  d=$(diff -rq -x PROVENANCE.toml "$FX/$v" "$OUT/regen/$v" | grep -v "SNAPSHOT.tsv" || true)
  if [ -n "$d" ]; then
    echo "   GOLDEN DRIFT ($(echo "$d" | wc -l) files differ) — the oracle output changed:"; echo "$d" | head -20 | sed 's/^/     /'
    status=1
  else
    echo "   goldens identical"
  fi
  diff <(grep -E 'oracle_sha256|oracle =' "$FX/$v/PROVENANCE.toml") <(grep -E 'oracle_sha256|oracle =' "$OUT/regen/$v/PROVENANCE.toml") >/dev/null \
    || echo "   NOTE: oracle binary changed (sha256) since the goldens were generated"
  if [ "${SLOW_SKIP_POP:-0}" != "1" ]; then
    S400=$ROOT/test/harness/fia/${v}_sample.txt; S120=$ROOT/test/harness/fia/${v}_sample120.txt
    [ -f "$S400" ] || (cd "$ROOT" && julia --project=. test/harness/fia/extract_sample.jl "$V" 400 "$S400" >/dev/null 2>&1)
    [ -f "$S120" ] || (cd "$ROOT" && julia --project=. test/harness/fia/extract_sample.jl "$V" 120 "$S120" >/dev/null 2>&1)
    echo "== [$V] population regime sweep → $OUT/pop_$v"
    bash "$ROOT/test/harness/fia/regime_sweep.sh" "$V" "$S400" "$S120" "$OUT/pop_$v" "${MAXJ:-5}"
    for t in "$OUT/pop_$v"/*.tally; do
      printf "   %-10s " "$(basename "$t" .tally)"; sed -n '/SIGN-TALLY/,/^$/p' "$t" | grep -E '^(TPA|BA|TCuFt)' | \
        awk '{printf "%s %s/%s | ", $1, substr($3,index($3,"=")+1), substr($4,index($4,"=")+1)}'; echo
    done
  fi
done
exit $status
