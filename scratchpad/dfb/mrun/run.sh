#!/bin/bash
# Usage: run.sh <binary-suffix> <keybase>   (keyfile=<keybase>.key, tre=<keybase>.tre)
# Produces <keybase>.sum and captures stderr (the g16 dump) to <keybase>.<suffix>.dump
cd /tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb/mrun
SUF=$1; KB=$2
ORACLE=/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb/FVSie_dfb${SUF}
cp -f iet01.tre "${KB}.tre"
rm -f "${KB}.sum" "${KB}.out" DFBOUT fort.*
printf '%s\n' "${KB}.key" | "$ORACLE" 2>"${KB}${SUF}.dump" >/dev/null
echo "EXIT=$? -> ${KB}.sum $(wc -l < ${KB}.sum 2>/dev/null) lines, dump $(wc -l < ${KB}${SUF}.dump) lines"
