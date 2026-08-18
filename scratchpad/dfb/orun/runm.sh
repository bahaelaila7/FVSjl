#!/bin/bash
# $1=suffix $2=keyfile $3=treefile $4=mainout
cd /tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb/orun
ORACLE=/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb/FVSie_dfb$1
rm -f *.sum "$4"
printf '%s\n%s\n%s\n' "$2" "$3" "$4" | "$ORACLE" 2>dfbgo_err.txt
echo "EXIT=$?"
