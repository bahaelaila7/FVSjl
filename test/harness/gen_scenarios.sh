#!/usr/bin/env bash
# gen_scenarios.sh — derive a matrix of Southern-variant test scenarios from the
# known-valid snt01.key, varying the dimensions we want broader coverage on:
# forest type / ecological unit, site species + index, fire, and cycle count.
#
# Each scenario is a single-stand .key (first stand of snt01, up to the first
# PROCESS) with targeted keyword edits. Starting from a valid key keeps every
# generated input legal FVS, so the Fortran oracle accepts them. The companion
# inline TREEDATA is preserved (so species mix is shared) except for the
# species-mix scenarios, which swap the leading species code.
#
# Output: test/harness/scenarios/*.key  (run through the three-way harness next).

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="${FVS_SN_KEY:-/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key}"
OUT="$HERE/scenarios"
mkdir -p "$OUT"
rm -f "$OUT"/*.key

# First stand only: everything up to and including the first TREEDATA block's
# data + ECHOSUM, then a terminating PROCESS/STOP. Take lines up to first PROCESS.
base="$OUT/_base.key"
awk '/^PROCESS/{print "PROCESS"; print "STOP"; exit} {print}' "$SRC" > "$base"

emit() {  # emit <name> <sed-script>
    local name="$1" script="$2"
    sed "$script" "$base" > "$OUT/$name.key"
    echo "  generated $name.key"
}

# 1. baseline (snt01 stand 1, unmodified)
cp "$base" "$OUT/s00_baseline.key"; echo "  generated s00_baseline.key"

# 2. site-index variations (drives growth/mortality across the site range)
emit s01_site_low   's/^SITECODE.*/SITECODE          63      40./'
emit s02_site_high  's/^SITECODE.*/SITECODE          63      90./'

# 3. site-species variations (different master-group → different SITSET fan-out)
emit s03_sitesp_lp  's/^SITECODE.*/SITECODE          13      70./'   # loblolly pine
emit s04_sitesp_yp  's/^SITECODE.*/SITECODE          45      75./'   # yellow-poplar

# 4. ecological-unit / forest variations (different physiographic + FORTYP path)
emit s05_ecounit_m221 's/231Dd /M221   /'
emit s06_ecounit_232  's/231Dd /232    /'
emit s07_forest_808   's/^STDINFO        80106/STDINFO        80806/'

# 5. cycle-count variations
emit s08_cyc3  's/^NUMCYCLE.*/NUMCYCLE         3.0/'
emit s09_cyc20 's/^NUMCYCLE.*/NUMCYCLE        20.0/'

# 6. fire scenario (FFE simulated wildfire mid-run)
sed 's/^NUMCYCLE.*/NUMCYCLE         5.0/' "$base" | \
  awk '/^PROCESS/{print "FMIN"; print "SIMFIRE"; print "END"; print; next} {print}' \
  > "$OUT/s10_fire.key"; echo "  generated s10_fire.key"

# 7. management: a thinning from below
awk '/^PROCESS/{print "THINBTA            2010.       0.        5.       999.       0.       1."; print; next} {print}' \
  "$base" > "$OUT/s11_thinbta.key"; echo "  generated s11_thinbta.key"

rm -f "$base"

# Every scenario shares snt01's tree list; copy the companion .tre next to each key
# (FVSjl/FVSjulia/Fortran all read <stem>.tre from the run dir).
TRE="${SRC%.key}.tre"
if [ -f "$TRE" ]; then
    for k in "$OUT"/*.key; do cp "$TRE" "${k%.key}.tre"; done
fi
echo "Scenarios in $OUT:"; ls "$OUT"/*.key | wc -l
