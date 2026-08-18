#!/bin/bash
FVS=/workspace/ForestVegetationSimulator
echo "=== which br*.o are present (ie) ==="
ls "$FVS"/bin/FVSie_buildDir/br*.o 2>/dev/null | xargs -n1 basename
echo "=== count of REAL wpbr objs (brcank/brdam/brin/brinit) per variant ==="
for d in "$FVS"/bin/FVS*_buildDir; do
  v=$(basename "$d" | sed 's/FVS//; s/_buildDir//')
  c=0
  for o in brcank brdam brin brinit brgi brcgro brial brblkd; do
    [ -f "$d/$o.o" ] && c=$((c+1))
  done
  echo "$v: real_model_objs=$c/8"
done
