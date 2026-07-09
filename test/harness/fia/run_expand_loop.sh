#!/usr/bin/env bash
# run_expand_loop.sh — continuous driver for the full-population FIA coverage sweep.
# Repeatedly runs one batch-cycle (SN first, then NE→CS→LS) until EITHER the dig-queue reaches DIGCAP (~200
# dig-worthy discrepancies → PAUSE for the caller to root-cause/fix) OR all variants are exhausted.
# Cursor-based & checkpointed every cycle ⇒ safely resumable. Runs in the background across many cycles.
# Env: BATCH (default 2000), DIGCAP (default 200).
set -u
cd /workspace/FVSjl
export BATCH=${BATCH:-2000}
DIGCAP=${DIGCAP:-200}
DIGQ=docs/fia_dig_queue.csv
while true; do
  out=$(bash test/harness/fia/run_expand_cycle.sh 2>&1)
  echo "$out"
  echo "$out" | grep -q "ALL_VARIANTS_EXHAUSTED" && { echo "LOOP_DONE all exhausted"; break; }
  echo "$out" | grep -q "RUN FAILED" && { echo "LOOP_HALT run failed"; break; }
  dign=0; [ -f $DIGQ ] && dign=$(($(wc -l < $DIGQ)-1))
  if [ "$dign" -ge "$DIGCAP" ]; then echo "DIG_PAUSE digq=$dign (>= $DIGCAP) — pausing sweep for dig/fix"; break; fi
done
