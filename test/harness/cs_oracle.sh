#!/usr/bin/env bash
# cs_oracle.sh — produce live FVS Northeast ground-truth output for a .key file.
#
# The NE port's oracle. The shipped bin/FVScs fails on this box's GLIBC (needs
# 2.38/2.43), so — exactly like the SN harness (fortran_baseline.sh) — we relink
# a working binary (/tmp/FVScs_new) from the resolved .o files in the NE build dir
# plus the glibc shim, then run it on a keyfile in an isolated dir.
#
# Usage:  cs_oracle.sh <keyfile> <outdir>      → prints <outdir>, writes <stem>.sum etc.
#
# NOTE: tests/FVScs/net01.sum.save is STALE (built from source dated 20250930; the
# current tree is 20260401). Cycle-0 still matches bit-exact, but cycle-1+ diverge
# from the .save due to model updates — validate cycle-1+ against THIS live binary,
# not the committed .save.
set -euo pipefail

BUILDDIR="${FVSNE_BUILDDIR:-/workspace/ForestVegetationSimulator/bin/FVScs_buildDir}"
BIN="${FVSNE_BIN:-/tmp/FVScs_new}"
SHIM=/tmp/glibc_shim.o

_ensure_binary() {
    [ -x "$BIN" ] && return 0
    if [ ! -f "$SHIM" ]; then
        cat > /tmp/glibc_shim.c <<'EOF'
#define _GNU_SOURCE
#include <stdio.h>
#include <stdarg.h>
int __isoc23_sscanf(const char *s, const char *f, ...){va_list a;va_start(a,f);int r=vsscanf(s,f,a);va_end(a);return r;}
int __isoc99_sscanf(const char *s, const char *f, ...){va_list a;va_start(a,f);int r=vsscanf(s,f,a);va_end(a);return r;}
EOF
        gcc -c -O2 /tmp/glibc_shim.c -o "$SHIM"
    fi
    local n; n=$(ls "$BUILDDIR"/*.o 2>/dev/null | wc -l)
    [ "$n" -lt 100 ] && { echo "ERROR: $BUILDDIR has only $n .o files" >&2; return 1; }
    ( cd "$BUILDDIR" && gfortran -o "$BIN" *.o "$SHIM" -lpthread -ldl )
}

main() {
    local key="$1" out="$2"
    _ensure_binary
    mkdir -p "$out"
    local kbase kdir; kbase=$(basename "$key"); kdir=$(cd "$(dirname "$key")" && pwd)
    cp "$key" "$out/$kbase"
    for ext in tre trl chp sng; do
        [ -f "$kdir/${kbase%.key}.$ext" ] && cp "$kdir/${kbase%.key}.$ext" "$out/" || true
    done
    ( cd "$out" && "$BIN" --keywordfile="$out/$kbase" >/dev/null 2>&1 ) || true
    echo "$out"
}
main "$@"
