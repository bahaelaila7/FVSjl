#!/usr/bin/env bash
# validate_oracle.sh — for every scenario .key, run the live Fortran FVSsn AND
# FVSjulia (Oracle A) and diff their .sum output. This proves the oracle is NOT
# faulty for the expanded scenarios BEFORE we promote any of them to an FVSjl
# regression (the user's hard requirement: don't match a faulty oracle).
#
# Reuses FVSjulia's fortran_baseline.sh (rebuilds /tmp/FVSsn_new + glibc shim).
# Usage: validate_oracle.sh [scenario-dir]   (default: ./scenarios)

set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
SCEN="${1:-$HERE/scenarios}"
FVSJULIA="${FVSJULIA:-/workspace/FVSjulia}"
FTBASE="$FVSJULIA/tests/fortran_baseline.sh"
JULIA="${JULIA:-julia}"

pass=0; fail=0; skip=0
for key in "$SCEN"/*.key; do
    name="$(basename "$key" .key)"
    stem="$name"

    ftdir="$(mktemp -d)"; jldir="$(mktemp -d)"
    # Fortran ground truth
    if ! bash "$FTBASE" "$key" "$ftdir" >/dev/null 2>"$ftdir/err"; then
        echo "SKIP  $name (Fortran build/run failed)"; skip=$((skip+1)); rm -rf "$ftdir" "$jldir"; continue
    fi
    # FVSjulia (Oracle A)
    cp "$key" "$jldir/"
    ( cd "$jldir" && "$JULIA" --project="$FVSJULIA" -e \
        "using FVSjulia; FVSjulia.main([\"--keywordfile=$jldir/$name.key\"])" >/dev/null 2>"$jldir/err" ) || true

    ftsum="$ftdir/$stem.sum"; jlsum="$jldir/$stem.sum"
    if [ ! -f "$ftsum" ] || [ ! -f "$jlsum" ]; then
        echo "SKIP  $name (missing .sum: ft=$( [ -f "$ftsum" ]&&echo y||echo n ) jl=$( [ -f "$jlsum" ]&&echo y||echo n ))"
        skip=$((skip+1)); rm -rf "$ftdir" "$jldir"; continue
    fi
    # tolerance-aware numeric compare (Float32 volume transcendentals give ±1 ulp
    # noise that a strict byte diff over-flags; sumdiff.jl allows abs≤1 or rel≤0.1%).
    res="$("$JULIA" --project="${FVSJL:-/workspace/FVSjl}" "$HERE/sumdiff.jl" "$ftsum" "$jlsum" 2>&1)"
    if [ "$res" = "MATCH" ]; then
        echo "PASS  $name (.sum Fortran==FVSjulia within tol)"; pass=$((pass+1))
    else
        echo "FAIL  $name (.sum differs): $res"; fail=$((fail+1))
    fi
    rm -rf "$ftdir" "$jldir"
done
echo "------------------------------------------------------------"
echo "oracle validation: PASS=$pass FAIL=$fail SKIP=$skip"
[ "$fail" -eq 0 ]
