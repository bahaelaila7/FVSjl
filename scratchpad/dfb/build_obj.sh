#!/bin/bash
# Compile DFB Fortran sources into objects for a BC (Douglas-fir Beetle) oracle relink.
# Reads pristine sources; writes only into scratchpad obj/. Does NOT touch pristine tree.
set -e
OBJ=/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb/obj
mkdir -p "$OBJ"
cd "$OBJ"
rm -f build.log
BD=/workspace/ForestVegetationSimulator/bin/FVSbc_buildDir
SRC=/workspace/ForestVegetationSimulator/dfb
FILES="dfbdrv dfbmod dfbmrt dfbdbh dfbind dfbinv dfbint dfbran dfbin dfbgo dfber dfbprb dfbsch dfbwin dfbtab dfbdam dfbout dfbhed dfbincr dfblkdbc"
for f in $FILES; do
  if gfortran-16 -std=legacy -w -fno-automatic -c -I"$BD" -I"$SRC/common" "$SRC/$f.f" -o "$f.o" 2>> build.log; then
    :
  else
    echo "FAIL $f" | tee -a build.log
  fi
done
echo "=== objects built ==="
ls *.o | wc -l
echo "=== build.log ==="
cat build.log
