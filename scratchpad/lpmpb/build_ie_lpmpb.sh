#!/bin/bash
# FVSie_lpmpb relink — IE object set with base/exmpb.f no-op stub replaced by real lpmpb/*.o.
# Optional first arg = directory of *.f overrides (instrumented lpmpb objects) compiled LAST.
# Second arg = output binary path.
set -u
SRC=/workspace/ForestVegetationSimulator/bin/FVSie_buildDir
LPMPB=/workspace/ForestVegetationSimulator/lpmpb
WORK=/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/lpmpb
OBJ=$WORK/ieobj
SHIM=/workspace/.crwork/isoc23_shim.o
OVR="${1:-}"
OUT="${2:-$WORK/FVSie_lpmpb}"
FC=gfortran-16
FLAGS="-c -std=legacy -w -fno-automatic -O0 -I$LPMPB -I$SRC"
mkdir -p "$OBJ"; cd "$SRC" || exit 2
nok=0; nfb=0
if [ ! -f "$OBJ/.builddir_done" ]; then
  for f in *.f *.for; do
    b=$(basename "$f" | sed 's/\.[^.]*$//')
    [ "$b" = "exmpb" ] && continue
    if $FC $FLAGS "$f" -o "$OBJ/$b.o" 2>/dev/null; then nok=$((nok+1))
    elif [ -f "$b.o" ]; then cp "$b.o" "$OBJ/$b.o"; nfb=$((nfb+1)); fi
  done
  touch "$OBJ/.builddir_done"
fi
nlp=0
cd "$LPMPB" || exit 2
for f in *.f; do
  b=$(basename "$f" .f)
  # keep ONLY the IE block data (mpblkdie); drop generic + other-variant block data
  case "$b" in mpblkd|mpblkdbm|mpblkdci|mpblkdcr|mpblkdec|mpblkdem|mpblkdso|mpblkdtt|mpblkdut) continue;; esac
  if $FC $FLAGS "$f" -o "$OBJ/$b.o" 2>/tmp/claude-1000/lpmpb_${b}.err; then nlp=$((nlp+1))
  else echo "LPMPB COMPILE FAIL: $b"; head -4 /tmp/claude-1000/lpmpb_${b}.err; fi
done
# PPE stub (in case any landscape symbol is referenced)
$FC $FLAGS /workspace/ForestVegetationSimulator/base/exppe.f -o "$OBJ/exppe.o" 2>/dev/null && echo "exppe stub added"
# TXNOTE (wsbwe) — ??X.exe warning flag setter referenced by MPBIN
$FC $FLAGS /workspace/ForestVegetationSimulator/wsbwe/txnote.f -o "$OBJ/txnote.o" 2>/dev/null && echo "txnote added"
if [ -n "$OVR" ] && [ -d "$OVR" ]; then
  for f in "$OVR"/*.f; do
    [ -e "$f" ] || continue
    b=$(basename "$f" .f)
    if $FC $FLAGS "$f" -o "$OBJ/$b.o" 2>/tmp/claude-1000/lpmpb_ovr_${b}.err; then echo "OVERRIDE $b"
    else echo "OVERRIDE FAIL: $b"; head -6 /tmp/claude-1000/lpmpb_ovr_${b}.err; fi
  done
fi
echo "buildDir: compiled_ok=$nok fell_back=$nfb ; lpmpb: compiled=$nlp"
COBJ=""; for c in "$SRC"/*.c; do COBJ="$COBJ ${c%.c}.o"; done
cd "$OBJ" || exit 3
if $FC -o "$OUT" *.o $COBJ "$SHIM" -lpthread -ldl 2>/tmp/claude-1000/ie_lpmpb_link.err; then
  echo "LINK_OK $OUT"; ls -la "$OUT"
else
  echo "LINK_FAILED:"; grep -i 'undefined\|multiple def' /tmp/claude-1000/ie_lpmpb_link.err | sort -u | head -30
fi
