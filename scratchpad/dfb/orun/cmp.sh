#!/bin/bash
cd /tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb/orun
ORACLE=/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb/FVSie_dfb
rm -f *.sum
printf 'ctrl.key\ntrees.tre\n' | "$ORACLE" >/dev/null 2>&1
cp ctrl.sum ctrl_clean.sum
rm -f *.sum
printf 'dfb5.key\ntrees.tre\n' | "$ORACLE" >/dev/null 2>&1
cp dfb5.sum dfb_on.sum
echo "=== CONTROL (no DFB) ==="
grep -v '^-999' ctrl_clean.sum
echo "=== DFB ON (MANSTART+MANSCHED 2) ==="
grep -v '^-999' dfb_on.sum
