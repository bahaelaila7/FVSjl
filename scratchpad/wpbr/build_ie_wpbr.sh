#!/bin/bash
# FVSie_wpbr relink — IE object set with base/exbrus.f no-op stub replaced by real
# wpbr/*.o (IE block data brblkdie.f only; brblkd/cr/so excluded — all share
# BLOCK DATA BRBLKD). Optional first arg = dir of instrumented *.f overrides.
set -u
SRC=/workspace/ForestVegetationSimulator/bin/FVSie_buildDir
WPBR=/workspace/ForestVegetationSimulator/wpbr
WORK=/workspace/FVSjl/.claude/worktrees/agent-a7dcc56dae713b70f/scratchpad/wpbr
OBJ=$WORK/ieobj
SHIM=/workspace/.crwork/isoc23_shim.o
OVR="${1:-}"
OUT="${2:-$WORK/FVSie_wpbr}"
FC=gfortran-16
FLAGS="-c -std=legacy -w -fno-automatic -O0 -I$WPBR -I$SRC"
mkdir -p "$OBJ"; cd "$SRC" || exit 2
nok=0; nfb=0
if [ ! -f "$OBJ/.builddir_done" ]; then
  for f in *.f *.for; do
    b=$(basename "$f" | sed 's/\.[^.]*$//')
    [ "$b" = "exbrus" ] && continue
    if $FC $FLAGS "$f" -o "$OBJ/$b.o" 2>/dev/null; then nok=$((nok+1))
    elif [ -f "$b.o" ]; then cp "$b.o" "$OBJ/$b.o"; nfb=$((nfb+1)); fi
  done
  touch "$OBJ/.builddir_done"
fi
nbr=0
cd "$WPBR" || exit 2
for f in *.f; do
  b=$(basename "$f" .f)
  # IE variant: keep brblkdie only; drop the base/CR/SO block data (dup BLOCK DATA BRBLKD)
  case "$b" in brblkd|brblkdcr|brblkdso|exbrus) continue;; esac
  if $FC $FLAGS "$f" -o "$OBJ/$b.o" 2>/tmp/claude-1000/wpbr_${b}.err; then nbr=$((nbr+1))
  else echo "WPBR COMPILE FAIL: $b"; head -3 /tmp/claude-1000/wpbr_${b}.err; fi
done
# missing-dep stubs seen with the insect ports
$FC $FLAGS /workspace/ForestVegetationSimulator/base/exppe.f -o "$OBJ/exppe.o" 2>/dev/null && echo "exppe stub added"
if [ -n "$OVR" ] && [ -d "$OVR" ]; then
  for f in "$OVR"/*.f; do
    [ -e "$f" ] || continue
    b=$(basename "$f" .f)
    if $FC $FLAGS "$f" -o "$OBJ/$b.o" 2>/tmp/claude-1000/wpbr_ovr_${b}.err; then echo "OVERRIDE $b"
    else echo "OVERRIDE FAIL: $b"; head -5 /tmp/claude-1000/wpbr_ovr_${b}.err; fi
  done
fi
echo "buildDir: compiled_ok=$nok fell_back=$nfb ; wpbr: compiled=$nbr"
COBJ=""; for c in "$SRC"/*.c; do COBJ="$COBJ ${c%.c}.o"; done
cd "$OBJ" || exit 3
if $FC -o "$OUT" *.o $COBJ "$SHIM" -lpthread -ldl 2>/tmp/claude-1000/ie_wpbr_link.err; then
  echo "LINK_OK $OUT"; ls -la "$OUT"
else
  echo "LINK_FAILED:"; grep -i 'undefined\|multiple def' /tmp/claude-1000/ie_wpbr_link.err | sort -u | head -30
fi
