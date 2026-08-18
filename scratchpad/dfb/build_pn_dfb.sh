#!/bin/bash
# Build the DFB objects for PN (dfblkdpn block data, IDFSPC=16) and relink FVSpn_dfb.
set -e
D=/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb
BD=/workspace/ForestVegetationSimulator/bin/FVSpn_buildDir
SRC=/workspace/ForestVegetationSimulator/dfb
G16=/workspace/.pnwork/g16obj
SHIM=/workspace/.crwork/isoc23_shim.o
OBJ=$D/pnobj
OUT=$D/FVSpn_dfb
mkdir -p "$OBJ"; cd "$OBJ"; rm -f build.log
FILES="dfbdam dfbdbh dfbdrv dfber dfbgo dfbhed dfbin dfbind dfbint dfbinv dfbmod dfbmrt dfbout dfbprb dfbran dfbsch dfbtab dfbwin dfblkdpn"
for f in $FILES; do
  gfortran-16 -std=legacy -w -fno-automatic -O0 -c -I"$BD" -I"$SRC/common" "$SRC/$f.f" -o "$f.o" 2>>build.log || echo "FAIL $f" | tee -a build.log
done
gfortran-16 -std=legacy -w -fno-automatic -O0 -c -I"$BD" /workspace/ForestVegetationSimulator/wsbwe/txnote.f -o txnote.o 2>>build.log || echo "FAIL txnote"
echo "=== dfb objects built: $(ls *.o | wc -l) ==="
BASEOBJ=$(ls $G16/*.o | grep -v '/exdfb.o')
COBJ=""; for c in "$BD"/*.c; do COBJ="$COBJ ${c%.c}.o"; done
gfortran-16 -o "$OUT" $BASEOBJ $OBJ/*.o $COBJ "$SHIM" -lpthread -ldl 2>$D/pnlink.err \
  && echo "LINK_OK $OUT" || { echo LINK_FAIL; grep -i 'undefined reference' $D/pnlink.err | sed 's/.*to //' | sort -u | head; }
