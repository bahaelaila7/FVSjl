#!/usr/bin/env bash
# gen_species_scenarios.sh — broad species + fire coverage for the harness.
#
# SPECIES: overwrite snt01.tre's species field (cols 34-36) with each of a wide
# representative set of SN species (pines, oaks, hickories, maples, gums, ashes,
# etc.), yielding homogeneous stands the FIA forest-type key classifies across many
# types. This exercises per-species growth/translation AND a broad spread of FORTYP.
#
# FIRE: several SIMFIRE variants (timing, repeat) on the baseline composition.
#
# Realism doesn't matter — only that FVSjulia and the live Fortran agree (oracle
# validation), and that FVSjl's FORTYP/init match across the matrix.

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/scenarios"
SRCKEY="${FVS_SN_KEY:-/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key}"
SRCTRE="${SRCKEY%.key}.tre"
mkdir -p "$OUT"
rm -f "$OUT"/sp_*.key "$OUT"/sp_*.tre "$OUT"/fire_*.key "$OUT"/fire_*.tre

base="$OUT/_spbase.key"
awk '/^PROCESS/{print "PROCESS"; print "STOP"; exit} {print}' "$SRCKEY" > "$base"

# Representative SN species (2-char alpha) across the major groups.
SPECIES="LP:lobpine SP:shleafpine SA:slashpine LL:longleafpine VP:virginiapine \
WP:whitepine PP:pitchpine SR:sprucepine PU:sandpine PD:pondpine \
WO:whiteoak SO:soredoak SK:scarletoak CO:chestnutoak RO:redoak PO:postoak \
BO:blackoak LO:liveoak WK:wateroak SN:swampchestnutoak CK:chinkapinoak \
HI:hickory RM:redmaple SM:sugarmaple YP:yelpoplar SU:sweetgum AB:beech \
BG:blackgum WA:whiteash SY:sycamore BC:blackcherry HB:hackberry BY:baldcypress \
JU:redcedar EL:elm BD:basswood"

for spec in $SPECIES; do
    code="${spec%%:*}"; label="${spec##*:}"
    name="sp_${label}"
    cp "$base" "$OUT/$name.key"
    awk -v sp="$code " '{print substr($0,1,33) sp substr($0,37)}' "$SRCTRE" > "$OUT/$name.tre"
done

# Fire scenarios (SIMFIRE) — scheduled wildfire at a year, varying timing/severity.
# Field layout (10-col, from FVSak/akt01.key): year, flame/intensity, fuel model, %area.
mkfire() {  # mkfire <name> <year> <intensity> <pctarea>
    local name="$1" yr="$2" inten="$3" area="$4"
    local sim
    sim="$(printf 'SIMFIRE   %10d%10.2f%10d%10.1f' "$yr" "$inten" 1 "$area")"
    sed 's/^NUMCYCLE.*/NUMCYCLE         5.0/' "$base" | \
      awk -v s="$sim" '/^PROCESS/{print "FMIN"; print s; print "END"; print; next} {print}' \
      > "$OUT/$name.key"
    cp "$SRCTRE" "$OUT/$name.tre"
}
mkfire fire_early 2000 10.00 50.0
mkfire fire_mid   2010 20.00 75.0
mkfire fire_late  2020  6.00 30.0

rm -f "$base"
echo "species scenarios: $(ls "$OUT"/sp_*.key | wc -l), fire scenarios: $(ls "$OUT"/fire_*.key | wc -l)"
