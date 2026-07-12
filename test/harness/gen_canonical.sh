#!/usr/bin/env bash
# gen_canonical.sh — copy the FVS reference multi-stand keys (net01/cst01/lst01) into
# test/fixtures/canonical/ and capture their live golden .sum, for test_canonical_multistand.jl.
#
# These reference keys bundle realistic management scenarios (control / thinning / shelterwood
# +ECON / FFE fire / bare-ground plant) across 5 stands each — the strongest per-scenario
# drop-in gate. Golden sums come from the freshly-relinked live binaries (/workspace/FVSjl/tmp/oracles/FVSxx_new).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
FIX="$ROOT/fixtures/canonical"; mkdir -p "$FIX"
TESTS=/workspace/ForestVegetationSimulator/tests

gen() { # <variant-dir> <stem> <bin>
  local vd="$1" stem="$2" bin="$3"
  cp "$TESTS/$vd/$stem.key" "$FIX/$stem.key"
  cp "$TESTS/$vd/$stem.tre" "$FIX/$stem.tre"
  ( cd "$FIX" && "$bin" --keywordfile="$stem.key" >/dev/null 2>&1 || true; cp "$stem.sum" "$stem.live.sum" )
  echo "  $stem: $(awk '/^-999/{n++}END{print n}' "$FIX/$stem.sum") stands, $(grep -cE '^(19|20)[0-9][0-9]' "$FIX/$stem.sum") rows"
}

gen FVSsn snt01 /workspace/FVSjl/tmp/oracles/FVSsn_new
gen FVSne net01 /workspace/FVSjl/tmp/oracles/FVSne_new
gen FVScs cst01 /workspace/FVSjl/tmp/oracles/FVScs_new
gen FVSls lst01 /workspace/FVSjl/tmp/oracles/FVSls_new
find "$FIX" -maxdepth 1 -type f ! -name '*.key' ! -name '*.tre' ! -name '*.live.sum' -delete
echo "canonical multi-stand fixtures + golden sums written to $FIX"
