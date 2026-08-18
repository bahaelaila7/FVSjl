#!/bin/bash
cd /tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb/orun
ORACLE=/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb/FVSie_dfb
rm -f *.sum *.out fort.* DFBOUT
printf '%s\n' "$1" | "$ORACLE"
echo "EXIT=$?"
ls -la *.sum *.out 2>/dev/null
