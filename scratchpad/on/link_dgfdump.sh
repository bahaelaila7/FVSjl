#!/bin/bash
set -e
OBJ=/workspace/.onwork/g16obj
SHIM=/workspace/.onwork/isoc23_shim.o
OUT=/workspace/.onwork/FVSon_dgfdump
DUMP=/workspace/.onwork/dgf_dump.o
cd "$OBJ"
gfortran-16 -o "$OUT" $DUMP $(ls *.o | grep -v '^dgf.o$') "$SHIM" -lpthread -ldl -lstdc++ 2>/tmp/claude-1000/dgfdlink.err
echo "linked $OUT"
