#!/usr/bin/env bash
# resume_sweep.sh — one-command recovery of the FIA/FVS coverage sweep after a container restart.
#
# A restart wipes the EPHEMERAL overlay (/ and /tmp) but preserves the PERSISTENT btrfs volume (/workspace).
# So the verification STATE is safe — the per-stand results (data/fia_sweep.db), the progress cursor
# (test/harness/fia/expand/<v>.cursor, also mirrored in the DB progress table), and the working dir
# (.sweep_work) all live on /workspace. The only casualties are REGENERABLE: the live-FVS oracle binary
# (/workspace/FVSjl/tmp/oracles/FVSsn_new) and transient batch files. This script rebuilds the oracle if missing, reconciles the cursor
# (max of the file and the DB mirror), reports state, and resumes the loop. Nothing is re-swept from scratch.
#
#   bash test/harness/fia/resume_sweep.sh [SN|NE|CS|LS]   # default SN; env BATCH/CYCLE_TO/DIGCAP as usual
set -u
cd /workspace/FVSjl
V=${1:-SN}; vl=$(echo "$V" | tr 'A-Z' 'a-z')
DB=/workspace/FVSjl/data/fia_sweep.db
JL="julia --project=."
declare -A ORACLE=( [SN]=/workspace/FVSjl/tmp/oracles/FVSsn_new [NE]=/workspace/FVSjl/tmp/oracles/FVSne_new [CS]=/workspace/FVSjl/tmp/oracles/FVScs_new [LS]=/workspace/FVSjl/tmp/oracles/FVSls_new )
BIN=${ORACLE[$V]}

echo "== resume_sweep $V =="

# 0. Julia depot — ~/.julia is on the ephemeral overlay, so a restart can EVICT installed packages (SQLite/etc.),
#    which silently breaks expand_batch/ledger_fia (they import SQLite) and stalls the sweep with "DONE". Re-download
#    at the manifest-pinned versions (no version change; floor-safe) if a core dep is missing.
if ! julia --project=. -e 'import SQLite' >/dev/null 2>&1; then
  echo "julia depot missing packages — Pkg.instantiate() ..."
  julia --project=. -e 'using Pkg; Pkg.instantiate()' 2>&1 | tail -3
fi

# 1. Oracle binary — regenerable; relink if the ephemeral /tmp copy was lost on restart.
if [ ! -x "$BIN" ]; then
  echo "oracle $BIN absent — relinking via test/harness/${vl}_oracle.sh ..."
  # sn_oracle.sh (and siblings) relink $BIN as a side effect of running any keyfile; use a shipped example.
  bash test/harness/${vl}_oracle.sh examples/fia/163384065010854.key >/dev/null 2>&1 || true
  [ -x "$BIN" ] && echo "  relinked ok ($(du -h "$BIN" | cut -f1))" || echo "  WARN: relink failed — run test/harness/${vl}_oracle.sh manually"
else
  echo "oracle $BIN present."
fi

# 2. Reconcile the cursor: max(working file, DB mirror). Either could be stale after an abrupt stop; the max is
#    always safe (re-sweeping is idempotent — the DB upserts by CN).
CURF=test/harness/fia/expand/$vl.cursor
fc=0; [ -f "$CURF" ] && fc=$(cat "$CURF")
dc=$($JL test/harness/fia/sweep_db.jl getcursor "$DB" "$V" 2>/dev/null); dc=${dc:-0}
cur=$fc; [ "$dc" -gt "$cur" ] && cur=$dc
mkdir -p "$(dirname "$CURF")"; echo "$cur" > "$CURF"
$JL test/harness/fia/sweep_db.jl setcursor "$DB" "$V" "$cur" >/dev/null 2>&1
echo "cursor reconciled: file=$fc db=$dc -> $cur"

# 3. Report durable state.
echo "-- durable coverage (data/fia_sweep.db) --"
$JL test/harness/fia/sweep_db.jl stats "$DB" "$V" 2>/dev/null | sed 's/^/  /'
nd=$($JL test/harness/fia/sweep_db.jl digs "$DB" "$V" 2>/dev/null | grep -c '^'"$V"',' || true)
echo "  needs_dig worklist: $nd  (list: $JL test/harness/fia/sweep_db.jl digs $DB $V)"

# 4. Resume the sweep (single foreground loop; Ctrl-C or DIGCAP to stop).
echo "-- resuming loop (BATCH=${BATCH:-100} CYCLE_TO=${CYCLE_TO:-480} DIGCAP=${DIGCAP:-200}) --"
exec bash test/harness/fia/run_expand_loop.sh
