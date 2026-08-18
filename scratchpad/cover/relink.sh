#!/bin/bash
set -u
cd /workspace/FVSjl/.claude/worktrees/agent-ab2fafff0c8d72cb5/scratchpad/cover || exit 2
OBJ=/workspace/.emwork/g16obj
SRC=/workspace/ForestVegetationSimulator/bin/FVSem_buildDir
SHIM=/workspace/.crwork/isoc23_shim.o
OBJS=$(ls "$OBJ"/*.o | grep -v '/cvcw.o$')
COBJ=""
for c in "$SRC"/*.c; do COBJ="$COBJ ${c%.c}.o"; done
gfortran-16 -o FVSem_g16cov $OBJS cvcw_instr.o $COBJ "$SHIM" -lpthread -ldl 2>link.err
echo "link exit=$?"
grep -i undefined link.err | head
ls -la FVSem_g16cov 2>/dev/null
