#!/bin/bash
# FVSbc_wsbwe relink — CI object set + real wsbwe/*.o with CI host coeffs
# (bwebkci/bwebmsci included; EM/other-variant + generic + PPE + txnote excluded).
set -u
SRC=/workspace/ForestVegetationSimulator/bin/FVSbc_buildDir
BWE=/workspace/ForestVegetationSimulator/wsbwe
WORK=/workspace/.bcwork/wsbwe
OBJ=$WORK/bcobj
SHIM=/workspace/.crwork/isoc23_shim.o
OVR="${1:-}"
OUT="${2:-$WORK/FVSbc_wsbwe}"
FC=gfortran-16
FLAGS="-c -std=legacy -w -fno-automatic -O0 -I$BWE -I$SRC"
mkdir -p "$OBJ"; cd "$SRC" || exit 2
nok=0; nfb=0
if [ ! -f "$OBJ/.builddir_done" ]; then
  for f in *.f *.for; do
    [ -e "$f" ] || continue
    b=$(basename "$f" | sed 's/\.[^.]*$//')
    [ "$b" = "exbudl" ] && continue
    if $FC $FLAGS "$f" -o "$OBJ/$b.o" 2>/dev/null; then nok=$((nok+1))
    elif [ -f "$b.o" ]; then cp "$b.o" "$OBJ/$b.o"; nfb=$((nfb+1)); fi
  done
  touch "$OBJ/.builddir_done"
fi
nbw=0
cd "$BWE" || exit 2
for f in *.f; do
  b=$(basename "$f" .f)
  case "$b" in
    bwebk|bwebkci|bwebkbm|bwebkem|bwebkec|bwebkso|bwebktt) continue;;
    bwebms|bwebmsci|bwebmsbm|bwebmsem|bwebmsec|bwebmsso|bwebmstt) continue;;
    bweppatv|bweppgt|bwepppt) continue;;
    txnote) continue;;
  esac
  if $FC $FLAGS "$f" -o "$OBJ/$b.o" 2>/tmp/claude-1000/wsbwe_bc_${b}.err; then nbw=$((nbw+1))
  else echo "WSBWE COMPILE FAIL: $b"; head -3 /tmp/claude-1000/wsbwe_bc_${b}.err; fi
done
if [ -n "$OVR" ] && [ -d "$OVR" ]; then
  for f in "$OVR"/*.f; do [ -e "$f" ] || continue; b=$(basename "$f" .f)
    $FC $FLAGS "$f" -o "$OBJ/$b.o" 2>/tmp/claude-1000/wsbwe_bc_ovr_${b}.err && echo "OVERRIDE $b" || { echo "OVR FAIL $b"; head -5 /tmp/claude-1000/wsbwe_bc_ovr_${b}.err; }
  done
fi
echo "buildDir: ok=$nok fb=$nfb ; wsbwe: compiled=$nbw"
COBJ=""; for c in "$SRC"/*.c; do [ -e "$c" ] || continue; COBJ="$COBJ ${c%.c}.o"; done
cd "$OBJ" || exit 3
if $FC -o "$OUT" *.o $COBJ "$SHIM" -lpthread -ldl 2>/tmp/claude-1000/bc_wsbwe_link.err; then
  echo "LINK_OK $OUT"; ls -la "$OUT" | awk '{print $5,$9}'
else
  echo "LINK_FAILED:"; grep -i 'undefined\|multiple def' /tmp/claude-1000/bc_wsbwe_link.err | sort -u | head -20
fi
