#!/bin/bash
# Relink FVSie_dfb<SUFFIX> from clean ieobj/ swapping in one or more instrumented
# objects. Usage: relink_ovr.sh <suffix> <instr_obj1> [instr_obj2 ...]
# Each instr_obj basename (e.g. dfbgo.o) replaces the pristine ieobj/<name>.
set -e
D=/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb
BD=/workspace/ForestVegetationSimulator/bin/FVSie_buildDir
G16=/workspace/.iework/g16obj
SHIM=/workspace/.crwork/isoc23_shim.o
OBJ=$D/ieobj
SUF=$1; shift
OUT=$D/FVSie_dfb${SUF}
EXCL=""
for o in "$@"; do EXCL="$EXCL -e /$(basename "$o")\$"; done
BASEOBJ=$(ls $G16/*.o | grep -v '/exdfb.o')
if [ -n "$EXCL" ]; then DFBSET=$(ls $OBJ/*.o | grep -v $EXCL); else DFBSET=$(ls $OBJ/*.o); fi
COBJ=""; for c in "$BD"/*.c; do COBJ="$COBJ ${c%.c}.o"; done
gfortran-16 -o "$OUT" $BASEOBJ $DFBSET "$@" $COBJ "$SHIM" -lpthread -ldl 2>$D/ovr.err \
  && echo "LINK_OK $OUT" || { echo LINK_FAIL; tail -4 $D/ovr.err; }
