#!/bin/bash
# FVSon_g16 full rebuild — self-consistent g16 ABI, mirrors .emwork/build_g16.sh.
# Adds the 4 DB-writer sources missing from the ON buildDir/sourceList.
set -u
SRC=/workspace/ForestVegetationSimulator/bin/FVSon_buildDir
OBJ=/workspace/.onwork/g16obj
SHIM=/workspace/.onwork/isoc23_shim.o
OUT=/workspace/.onwork/FVSon_g16
FC=gfortran-16
FLAGS="-c -std=legacy -w -fno-automatic -finit-local-zero -O0"
mkdir -p "$OBJ"; cd "$SRC" || exit 2
nok=0; nfb=0
for f in *.f *.for; do
  b=$(basename "$f" | sed 's/\.[^.]*$//')
  if $FC $FLAGS "$f" -o "$OBJ/$b.o" 2>/dev/null; then nok=$((nok+1))
  elif [ -f "$b.o" ]; then cp "$b.o" "$OBJ/$b.o"; nfb=$((nfb+1)); fi
done
echo "fortran: compiled_ok=$nok fell_back=$nfb"
# 4 DB-writer sources missing from ON buildDir (compile from canonical trees)
$FC $FLAGS -I"$SRC" /workspace/ForestVegetationSimulator/dbsqlite/dbs_fiavbc_trls.f   -o "$OBJ/dbs_fiavbc_trls.o"
$FC $FLAGS -I"$SRC" /workspace/ForestVegetationSimulator/dbsqlite/dbs_fiavbc_cutlst.f -o "$OBJ/dbs_fiavbc_cutlst.o"
$FC $FLAGS -I"$SRC" /workspace/ForestVegetationSimulator/dbsqlite/dbs_fiavbc_atrtls.f -o "$OBJ/dbs_fiavbc_atrtls.o"
$FC $FLAGS -I"$SRC" /workspace/ForestVegetationSimulator/vdbsqlite/dbsreference.f      -o "$OBJ/dbsreference.o"
# C + C++ objects: compile FROM SOURCE with gcc-16/g++-16 for one consistent ABI
CDEF="-DANSI -DCMPgcc -w -O0 -I$SRC"
ncok=0
for c in "$SRC"/*.c; do
  b=$(basename "$c" .c)
  gcc-16 -c $CDEF "$c" -o "$OBJ/$b.o" 2>/dev/null && ncok=$((ncok+1))
done
for c in "$SRC"/*.cpp; do
  b=$(basename "$c" .cpp)
  g++-16 -c $CDEF "$c" -o "$OBJ/$b.o" 2>/dev/null && ncok=$((ncok+1))
done
echo "c/c++: compiled_ok=$ncok"
cd "$OBJ" || exit 3
if $FC -o "$OUT" *.o "$SHIM" -lpthread -ldl -lstdc++ 2>/tmp/claude-1000/g16link_on.err; then
  echo "LINK_OK $OUT"; ls -la "$OUT"
else
  echo "LINK_FAILED:"; grep -i 'undefined reference' /tmp/claude-1000/g16link_on.err | sed 's/.*to //' | sort -u | head
fi
