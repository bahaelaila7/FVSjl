#!/bin/bash
cd /tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb/pnrun
ORACLE=/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb/FVSpn_dfb
rm -f pnctrl.sum pndfb.sum
printf 'pnctrl.key\n' | "$ORACLE" >/dev/null 2>pnc.err; echo "ctrl exit=$?"
printf 'pndfb.key\n'  | "$ORACLE" >/dev/null 2>pnd.err; echo "dfb  exit=$?"
echo "year  ctrlMOR  dfbMOR"
grep -v '^-999' pnctrl.sum | awk '{print $1, $25}' > /tmp/pnc.txt
grep -v '^-999' pndfb.sum  | awk '{print $25}'      > /tmp/pnd.txt
paste /tmp/pnc.txt /tmp/pnd.txt
