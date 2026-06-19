#!/usr/bin/env bash
# gen_species_scenarios.sh — species-varied scenarios to exercise FORTYP across
# different forest types. Each overwrites the species field (cols 34-36) of
# snt01.tre with a single species, yielding a homogeneous stand the FIA
# forest-type key classifies differently. These give the FORTYP port real
# validation targets (the snt01-derived scenarios all share one tree list → all
# FORTYP=520). Realism doesn't matter — only that FVSjulia and Fortran agree.

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/scenarios"
SRCKEY="${FVS_SN_KEY:-/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key}"
SRCTRE="${SRCKEY%.key}.tre"
mkdir -p "$OUT"

# single-stand key template (first stand of snt01)
base="$OUT/_spbase.key"
awk '/^PROCESS/{print "PROCESS"; print "STOP"; exit} {print}' "$SRCKEY" > "$base"

# species → short label (2-char SN alpha codes)
for spec in "LP:lobpine" "SP:shleafpine" "WO:whiteoak" "YP:yelpoplar" "RM:redmaple" "SA:slashpine"; do
    code="${spec%%:*}"; label="${spec##*:}"
    name="sp_${label}"
    cp "$base" "$OUT/$name.key"
    awk -v sp="$code " '{print substr($0,1,33) sp substr($0,37)}' "$SRCTRE" > "$OUT/$name.tre"
    echo "  generated $name (species=$code)"
done
rm -f "$base"
echo "species scenarios: $(ls "$OUT"/sp_*.key | wc -l)"
