#!/bin/bash
# Scope which variants link real wpbr (br*.o) vs the exbrus.f no-op stub.
FVS=/workspace/ForestVegetationSimulator
for d in "$FVS"/bin/FVS*_buildDir; do
  v=$(basename "$d" | sed 's/FVS//; s/_buildDir//')
  n=$(ls "$d"/br*.o 2>/dev/null | wc -l)
  e=$(ls "$d"/exbrus.o 2>/dev/null | wc -l)
  echo "$v: br_objs=$n exbrus_stub=$e"
done
