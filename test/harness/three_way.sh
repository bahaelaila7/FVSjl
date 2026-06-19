#!/usr/bin/env bash
# three_way.sh — the full Fortran ↔ FVSjulia ↔ FVSjl harness over the scenario set.
#
# Leg 1+2 (validate_oracle.sh): Fortran FVSsn vs FVSjulia (Oracle A), full .sum —
#   proves the oracle is trustworthy for each scenario.
# Leg 3 (here): FVSjl cycle-0 stand summary vs the Fortran .sum's first data row
#   (TPA BA SDI CCF TopHt QMD). Until FVSjl has the .sum writer (C5), cycle-0 is the
#   slice we can three-way; it exercises init / NOTRE / stand statistics on every
#   scenario and surfaces divergences early.
#
# Usage: three_way.sh [scenario-dir]
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
SCEN="${1:-$HERE/scenarios}"
FVSJULIA="${FVSJULIA:-/workspace/FVSjulia}"
FTBASE="$FVSJULIA/tests/fortran_baseline.sh"
JULIA="${JULIA:-julia}"

# Leg 1+2: oracle validation
echo "### Leg 1+2: Fortran vs FVSjulia (oracle validation) ###"
bash "$HERE/validate_oracle.sh" "$SCEN"
echo
echo "### Leg 3: FVSjl cycle-0 stand summary vs Fortran .sum row 1 ###"
echo "                    TPA   BA  SDI  CCF  TopHt QMD"

pass=0; fail=0; skip=0
for key in "$SCEN"/*.key; do
    name="$(basename "$key" .key)"
    ftdir="$(mktemp -d)"
    if ! bash "$FTBASE" "$key" "$ftdir" >/dev/null 2>&1 || [ ! -f "$ftdir/$name.sum" ]; then
        printf "  %-18s SKIP (no Fortran .sum)\n" "$name"; skip=$((skip+1)); rm -rf "$ftdir"; continue
    fi
    # first data row (skip -999), fields: yr age TPA BA SDI CCF TopHt QMD
    row="$(grep -v -- '-999' "$ftdir/$name.sum" | head -1)"
    ref="$(echo "$row" | awk '{printf "%d %d %d %d %d %.1f", $3,$4,$5,$6,$7,$8}')"
    rm -rf "$ftdir"

    got="$("$JULIA" --project="$PWD" "$HERE/fvsjl_cycle0.jl" "$key" 2>/dev/null | tail -1)"
    # compare TPA/BA/SDI/CCF/TopHt within ±1 (rounding), QMD within 0.1
    ok=1
    read -r a1 a2 a3 a4 a5 a6 <<<"$ref"; read -r b1 b2 b3 b4 b5 b6 <<<"$got"
    for pair in "$a1 $b1" "$a2 $b2" "$a3 $b3" "$a4 $b4" "$a5 $b5"; do
        set -- $pair; d=$(( $1 - $2 )); [ "${d#-}" -le 1 ] || ok=0
    done
    awk "BEGIN{exit !(($a6-$b6)<0.11 && ($b6-$a6)<0.11)}" || ok=0
    if [ "$ok" -eq 1 ]; then
        printf "  %-18s OK    ft=[%s] fvsjl=[%s]\n" "$name" "$ref" "$got"; pass=$((pass+1))
    else
        printf "  %-18s DIFF  ft=[%s] fvsjl=[%s]\n" "$name" "$ref" "$got"; fail=$((fail+1))
    fi
done
echo "------------------------------------------------------------"
echo "FVSjl cycle-0 vs Fortran: OK=$pass DIFF=$fail SKIP=$skip"
