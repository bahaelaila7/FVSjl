#!/bin/bash
# Relink FVSbc WITH the real Douglas-fir Beetle model (dfb/*.o) in place of the exdfb.o stub.
# Optionally swap ONE instrumented dfb object: pass its .o path as $1.
set -e
BD=/workspace/ForestVegetationSimulator/bin/FVSbc_buildDir
SHIM=/workspace/.crwork/isoc23_shim.o
STUB=/workspace/.bcwork/bc_stubs.o
DFBOBJ=/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb/obj
OUT=/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb/FVSbc_dfb
# all buildDir objects except the exdfb stub
BASEOBJ=$(ls $BD/*.o | grep -v "/exdfb.o")
# dfb objects; if an instrumented override is given, drop the matching base name from the dfb set
if [ -n "$1" ]; then
  ovr="$1"; ovbase=$(basename "$ovr")
  DFBSET=$(ls $DFBOBJ/*.o | grep -v "/$ovbase")
  gfortran-16 -o "$OUT" $BASEOBJ $DFBSET "$ovr" "$STUB" "$SHIM"
else
  DFBSET=$(ls $DFBOBJ/*.o)
  gfortran-16 -o "$OUT" $BASEOBJ $DFBSET "$STUB" "$SHIM"
fi
echo "linked $OUT"
