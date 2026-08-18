#!/bin/bash
# $1 = binary suffix (e.g. _instr), $2 = keyfile, $3 = tree file
cd /tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb/orun
ORACLE=/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb/FVSie_dfb$1
rm -f *.sum
printf '%s\n%s\n' "$2" "$3" | "$ORACLE"
echo "EXIT=$?"
