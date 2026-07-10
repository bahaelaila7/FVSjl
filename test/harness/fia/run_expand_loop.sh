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
SC=${SWEEP_WORK:-/workspace/FVSjl/.sweep_work}   # persistent /workspace volume (see run_expand_cycle.sh)
mkdir -p $SC/expand
LOGF=$SC/expand/cycle_last.log
while true; do
  # Write the cycle output to a FILE (not `$(...)`) — the ledger/run logs can contain null bytes on empty
  # strata, which command substitution mangles (the earlier spurious LOOP_HALT). grep on the file is null-safe.
  bash test/harness/fia/run_expand_cycle.sh > $LOGF 2>&1
  cat $LOGF
  grep -qa "ALL_VARIANTS_EXHAUSTED" $LOGF && { echo "LOOP_DONE all exhausted"; break; }
  grep -qa "RUN FAILED" $LOGF && { echo "LOOP_HALT run failed"; break; }
  dign=0; [ -f $DIGQ ] && dign=$(($(wc -l < $DIGQ)-1))
  if [ "$dign" -ge "$DIGCAP" ]; then echo "DIG_PAUSE digq=$dign (>= $DIGCAP) — pausing sweep for dig/fix"; break; fi
done
