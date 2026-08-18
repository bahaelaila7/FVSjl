#!/bin/bash
# Relink FVSie_dfb from the ALREADY-built ieobj/ (no recompile), so an instrumented
# .o copied over its pristine name survives. $1 optional output suffix.
set -e
D=/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb
BD=/workspace/ForestVegetationSimulator/bin/FVSie_buildDir
G16=/workspace/.iework/g16obj
SHIM=/workspace/.crwork/isoc23_shim.o
OBJ=$D/ieobj
OUT=$D/FVSie_dfb${1:-}
BASEOBJ=$(ls $G16/*.o | grep -v '/exdfb.o')
COBJ=""; for c in "$BD"/*.c; do COBJ="$COBJ ${c%.c}.o"; done
gfortran-16 -o "$OUT" $BASEOBJ $OBJ/*.o $COBJ "$SHIM" -lpthread -ldl 2>$D/ielink.err && echo "LINK_OK $OUT" || { echo LINK_FAIL; tail -5 $D/ielink.err; }
