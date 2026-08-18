#!/bin/bash
# Relink the ON oracle into /workspace/.onwork.
#   relink_on.sh clean                     -> pristine /workspace/.onwork/FVSon_clean
#   relink_on.sh <name> /path/to/patched.o -> FVSon_<name> with one instrumented object swapped
# The 4 EXTRA objects (dbs_fiavbc_*, dbsreference) are DB-writer routines missing from
# the ON buildDir/sourceList; compiled from canonical /dbsqlite + /vdbsqlite by build_missing.sh.
set -e
BD=/workspace/ForestVegetationSimulator/bin/FVSon_buildDir
OW=/workspace/.onwork
SHIM=$OW/isoc23_shim.o
EXTRA="$OW/dbs_fiavbc_trls.o $OW/dbs_fiavbc_cutlst.o $OW/dbs_fiavbc_atrtls.o $OW/dbsreference.o"
OUT=$OW/FVSon_$1
if [ "$1" = "clean" ]; then
  gfortran-16 -o "$OUT" $(ls $BD/*.o) $EXTRA "$SHIM"
else
  OBJ="$2"; base=$(basename "$OBJ")
  gfortran-16 -o "$OUT" "$OBJ" $(ls $BD/*.o | grep -v "/$base") $EXTRA "$SHIM"
fi
echo "linked $OUT"
