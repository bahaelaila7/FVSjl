#!/usr/bin/env bash
# gen_species_scenarios.sh — homogeneous stands for ALL 90 SN species + fire variants.
#
# SPECIES (all_<CODE>): overwrite snt01.tre's species field (cols 34-36) with each
# of the 90 SN alpha codes, so every species' growth/crown/mortality coefficient row
# and the FORTYP classification are exercised. The 90 species span 45 distinct FIA
# forest types. (7 codes Fortran itself can't run as a pure stand → those SKIP.)
#
# FIRE (fire_*): SIMFIRE variants — timing, intensity, fuel model, salvage, repeat —
# exercising the FFE codepaths. Field layout (FVSak): year, intensity, fuel, %area.
#
# Realism doesn't matter — only that FVSjulia/Fortran agree and FVSjl matches.

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/scenarios"; mkdir -p "$OUT"
SRCKEY="${FVS_SN_KEY:-/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key}"
TRE="${SRCKEY%.key}.tre"
JULIA="${JULIA:-julia}"
FVSJL="${FVSJL:-/workspace/FVSjl}"
rm -f "$OUT"/all_*.key "$OUT"/all_*.tre "$OUT"/fire_*.key "$OUT"/fire_*.tre

base="$OUT/_spbase.key"
awk '/^PROCESS/{print "PROCESS"; print "STOP"; exit} {print}' "$SRCKEY" > "$base"

# all 90 SN alpha codes, from the loaded coefficients
CODES="$("$JULIA" --project="$FVSJL" -e 'using FVSjl; print(join([strip(x) for x in FVSjl.coefficients(Southern()).code_alpha]," "))')"
for code in $CODES; do
    cp "$base" "$OUT/all_${code}.key"
    awk -v sp="$code " '{print substr($0,1,33) sp substr($0,37)}' "$TRE" > "$OUT/all_${code}.tre"
done

# fire variants (SIMFIRE field layout: year, intensity, fuel-model, %area)
simfire() { printf 'SIMFIRE   %10d%10.2f%10d%10.1f' "$1" "$2" "$3" "$4"; }
firekey() {  # firekey <name> <ncycle> <body-awk-block>
    local name="$1" nc="$2" blk="$3"
    sed "s/^NUMCYCLE.*/NUMCYCLE         ${nc}.0/" "$base" | \
      awk -v b="$blk" '/^PROCESS/{print "FMIN"; printf "%s", b; print "END"; print; next} {print}' \
      > "$OUT/$name.key"
    cp "$TRE" "$OUT/$name.tre"
}
firekey fire_early   5 "$(simfire 2000 10 1 50)\n"
firekey fire_mid     5 "$(simfire 2010 20 1 75)\n"
firekey fire_late    5 "$(simfire 2020  6 1 30)\n"
firekey fire_fuel2   5 "$(simfire 2005 12 2 50)\n"
firekey fire_fuel9   5 "$(simfire 2005 12 9 50)\n"
firekey fire_fuel11  5 "$(simfire 2005 12 11 50)\n"
firekey fire_salvage 5 "$(simfire 2010 15 4 60)\nSALVAGE        2010.\n"
firekey fire_repeat  8 "$(simfire 2000 10 1 40)\n$(simfire 2020 20 1 60)\n"

rm -f "$base"
echo "species: $(ls "$OUT"/all_*.key | wc -l), fire: $(ls "$OUT"/fire_*.key | wc -l)"
