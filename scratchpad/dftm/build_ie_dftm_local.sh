#!/bin/bash
# FVSie_dftm relink — IE object set with base/exdftm.f no-op stub replaced by real dftm/*.o.
# Optional first arg = directory of *.f overrides (instrumented dftm objects) to compile LAST.
set -u
SRC=/workspace/ForestVegetationSimulator/bin/FVSie_buildDir
DFTM=/workspace/ForestVegetationSimulator/dftm
WORK=/workspace/FVSjl/scratchpad/dftm
OBJ=$WORK/ieobj
SHIM=/workspace/.crwork/isoc23_shim.o
OVR="${1:-}"
OUT="${2:-$WORK/FVSie_dftm}"
FC=gfortran-16
FLAGS="-c -std=legacy -w -fno-automatic -O0 -I$DFTM -I$SRC"
mkdir -p "$OBJ"; cd "$SRC" || exit 2
nok=0; nfb=0
# Only compile buildDir once (reuse cached .o on re-runs unless missing)
if [ ! -f "$OBJ/.builddir_done" ]; then
  for f in *.f *.for; do
    b=$(basename "$f" | sed 's/\.[^.]*$//')
    [ "$b" = "exdftm" ] && continue
    if $FC $FLAGS "$f" -o "$OBJ/$b.o" 2>/dev/null; then nok=$((nok+1))
    elif [ -f "$b.o" ]; then cp "$b.o" "$OBJ/$b.o"; nfb=$((nfb+1)); fi
  done
  touch "$OBJ/.builddir_done"
fi
ndf=0
cd "$DFTM" || exit 2
for f in *.f; do
  b=$(basename "$f" .f)
  case "$b" in tminitec|tminitem|tminitso|tminittt) continue;; esac
  if $FC $FLAGS "$f" -o "$OBJ/$b.o" 2>/tmp/claude-1000/dftm_${b}.err; then ndf=$((ndf+1))
  else echo "DFTM COMPILE FAIL: $b"; head -3 /tmp/claude-1000/dftm_${b}.err; fi
done
$FC $FLAGS /workspace/ForestVegetationSimulator/base/exppe.f -o "$OBJ/exppe.o" 2>/dev/null && echo "exppe stub added"
# Instrumented overrides compiled last (win over stock dftm .o)
if [ -n "$OVR" ] && [ -d "$OVR" ]; then
  for f in "$OVR"/*.f; do
    [ -e "$f" ] || continue
    b=$(basename "$f" .f)
    if $FC $FLAGS "$f" -o "$OBJ/$b.o" 2>/tmp/claude-1000/dftm_ovr_${b}.err; then echo "OVERRIDE $b"
    else echo "OVERRIDE FAIL: $b"; head -5 /tmp/claude-1000/dftm_ovr_${b}.err; fi
  done
fi
echo "buildDir: compiled_ok=$nok fell_back=$nfb ; dftm: compiled=$ndf"
COBJ=""; for c in "$SRC"/*.c; do COBJ="$COBJ ${c%.c}.o"; done
cd "$OBJ" || exit 3
if $FC -o "$OUT" *.o $COBJ "$SHIM" -lpthread -ldl 2>/tmp/claude-1000/ie_dftm_link.err; then
  echo "LINK_OK $OUT"; ls -la "$OUT"
else
  echo "LINK_FAILED:"; grep -i 'undefined\|multiple def' /tmp/claude-1000/ie_dftm_link.err | sort -u | head -20
fi
