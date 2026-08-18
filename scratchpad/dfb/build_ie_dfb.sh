#!/bin/bash
# Build the 20 real DFB objects for IE and relink FVSie WITH the DFB model
# (real dfb/*.o in place of the exdfb.o no-op stub). Reads pristine sources only;
# writes only into scratchpad. Optional $1 = an instrumented dfb .o to swap in.
set -e
D=/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb
BD=/workspace/ForestVegetationSimulator/bin/FVSie_buildDir
SRC=/workspace/ForestVegetationSimulator/dfb
G16=/workspace/.iework/g16obj
SHIM=/workspace/.crwork/isoc23_shim.o
OBJ=$D/ieobj
OUT=$D/FVSie_dfb
mkdir -p "$OBJ"; cd "$OBJ"; rm -f build.log
# the 20 IE DFB objects: all routines + dfblkdie (IE block data)
FILES="dfbdam dfbdbh dfbdrv dfber dfbgo dfbhed dfbin dfbind dfbint dfbinv dfbmod dfbmrt dfbout dfbprb dfbran dfbsch dfbtab dfbwin dfblkdie"
for f in $FILES; do
  gfortran-16 -std=legacy -w -fno-automatic -O0 -c -I"$BD" -I"$SRC/common" "$SRC/$f.f" -o "$f.o" 2>>build.log || echo "FAIL $f" | tee -a build.log
done
# txnote (referenced by dfbin.f via TXNOTE)
gfortran-16 -std=legacy -w -fno-automatic -O0 -c -I"$BD" /workspace/ForestVegetationSimulator/wsbwe/txnote.f -o txnote.o 2>>build.log || echo "FAIL txnote"
echo "=== dfb objects built: $(ls *.o | wc -l) ==="
# base objects = full IE g16 build minus exdfb.o
BASEOBJ=$(ls $G16/*.o | grep -v '/exdfb.o')
COBJ=""; for c in "$BD"/*.c; do COBJ="$COBJ ${c%.c}.o"; done
if [ -n "${1:-}" ]; then
  ovr="$1"; ovbase=$(basename "$ovr")
  DFBSET=$(ls $OBJ/*.o | grep -v "/$ovbase")
  gfortran-16 -o "$OUT" $BASEOBJ $DFBSET "$ovr" $COBJ "$SHIM" -lpthread -ldl 2>$D/ielink.err && echo "LINK_OK $OUT" || { echo LINK_FAIL; grep -i 'undefined reference' $D/ielink.err | sed 's/.*to //' | sort -u | head; }
else
  gfortran-16 -o "$OUT" $BASEOBJ $OBJ/*.o $COBJ "$SHIM" -lpthread -ldl 2>$D/ielink.err && echo "LINK_OK $OUT" || { echo LINK_FAIL; grep -i 'undefined reference' $D/ielink.err | sed 's/.*to //' | sort -u | head; }
fi
