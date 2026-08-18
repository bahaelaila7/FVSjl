#!/bin/bash
# FVSem_wsbwe relink — EM object set with base/exbudl.f no-op stub replaced by
# the real wsbwe/*.o (EM host coefficients bwebkem/bwebmsem). Mirrors
# scratchpad/dftm/build_ie_dftm.sh. WSBWE is stand-level (driver from base
# GRINCR/GRADD), so this yields a runnable stand oracle (unlike PPE-gated WWPB).
#
# Authoritative file set = canada/bin/FVSbcc_notReady.txt wsbwe block, with the
# BC host coeffs swapped to EM (bwebkbc->bwebkem, bwebmsbc->bwebmsem). Excludes:
# generic bwebk/bwebms, other-variant bwebk*/bwebms*, the 3 PPE handoff files
# (bweppatv/bweppgt/bwepppt — nothing in base calls them), and txnote.
#
# Optional arg1 = dir of *.f instrumented overrides (compiled LAST, win over stock).
set -u
SRC=/workspace/ForestVegetationSimulator/bin/FVSem_buildDir
BWE=/workspace/ForestVegetationSimulator/wsbwe
WORK=/workspace/FVSjl/.claude/worktrees/agent-aef698467e4369552/scratchpad/wsbwe
OBJ=$WORK/emobj
SHIM=/workspace/.crwork/isoc23_shim.o
OVR="${1:-}"
OUT="${2:-$WORK/FVSem_wsbwe}"
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
# wsbwe file set (EM). Skip generic/other-variant/PPE/txnote.
nbw=0
cd "$BWE" || exit 2
for f in *.f; do
  b=$(basename "$f" .f)
  case "$b" in
    bwebk|bwebkbc|bwebkbm|bwebkci|bwebkec|bwebkso|bwebktt) continue;;
    bwebms|bwebmsbc|bwebmsbm|bwebmsci|bwebmsec|bwebmsso|bwebmstt) continue;;
    bweppatv|bweppgt|bwepppt) continue;;
    txnote) continue;;
  esac
  if $FC $FLAGS "$f" -o "$OBJ/$b.o" 2>/tmp/claude-1000/wsbwe_${b}.err; then nbw=$((nbw+1))
  else echo "WSBWE COMPILE FAIL: $b"; head -3 /tmp/claude-1000/wsbwe_${b}.err; fi
done
# instrumented overrides last
if [ -n "$OVR" ] && [ -d "$OVR" ]; then
  for f in "$OVR"/*.f; do
    [ -e "$f" ] || continue
    b=$(basename "$f" .f)
    if $FC $FLAGS "$f" -o "$OBJ/$b.o" 2>/tmp/claude-1000/wsbwe_ovr_${b}.err; then echo "OVERRIDE $b"
    else echo "OVERRIDE FAIL: $b"; head -5 /tmp/claude-1000/wsbwe_ovr_${b}.err; fi
  done
fi
echo "buildDir: compiled_ok=$nok fell_back=$nfb ; wsbwe: compiled=$nbw"
COBJ=""; for c in "$SRC"/*.c; do [ -e "$c" ] || continue; COBJ="$COBJ ${c%.c}.o"; done
cd "$OBJ" || exit 3
if $FC -o "$OUT" *.o $COBJ "$SHIM" -lpthread -ldl 2>/tmp/claude-1000/em_wsbwe_link.err; then
  echo "LINK_OK $OUT"; ls -la "$OUT"
else
  echo "LINK_FAILED:"; grep -i 'undefined\|multiple def' /tmp/claude-1000/em_wsbwe_link.err | sort -u | head -30
fi
