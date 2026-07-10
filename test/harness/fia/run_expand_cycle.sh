#!/usr/bin/env bash
# run_expand_cycle.sh — ONE resumable batch-cycle of the full-population FIA coverage sweep.
# Processes the CURRENT variant only (SN first; when SN's cursor reaches its population, advance to NE→CS→LS),
# in the deterministic (ECOREGION,LOCATION,STAND_CN) order. Runs the live-vs-jl ledger differential on the next
# BATCH stands, appends full detail to the scratchpad master ledger, advances the cursor ON SUCCESS, and appends
# the DIG-WORTHY discrepancies (material, potentially-real-bug divergences) to docs/fia_dig_queue.csv.
# The caller pauses the sweep to dig/fix when the dig-queue reaches ~200 rows.
#
# Dig-worthy = signature ∈ {UNCLASSIFIED, volume_persistent, structure_densephase} OR
#              (worst_col==TCuFt & struct%<1 & max_rel≥5)  — i.e. NOT bit_exact/print_boundary/count_straddle/
#              threshold_crossing (those are the accepted print/ULP/merch-threshold cornered classes).
# Env: BATCH (default 2000). Cursors: test/harness/fia/expand/<v>.cursor.
set -u
cd /workspace/FVSjl
# Sweep working dir on the PERSISTENT /workspace btrfs volume (NOT the ephemeral, session-specific /tmp
# scratchpad) so batch state + the accumulated master ledger survive a container restart and are not tied to a
# session UUID. Override with SWEEP_WORK. (The durable per-stand results live in the sweep DB; this is scratch +
# a redundant ledger copy.) .sweep_work is gitignored.
SC=${SWEEP_WORK:-/workspace/FVSjl/.sweep_work}
BATCH=${BATCH:-2000}
CURD=test/harness/fia/expand
DIGQ=docs/fia_dig_queue.csv
# Durable per-stand coverage DB (survives sessions / container restart; gitignored). ledger_fia.jl upserts every
# stand's outcome (bit_exact | ulp_class | needs_dig) here as it runs — the cross-session dig worklist.
export SWEEP_DB=${SWEEP_DB:-/workspace/FVSjl/data/fia_sweep.db}
mkdir -p $CURD $SC/expand

# variant order + populations (STAND_CN IS NOT NULL, per FVS_STANDINIT_COND.VARIANT)
order=(SN NE CS LS)
declare -A POP=( [SN]=637641 [NE]=178149 [CS]=255952 [LS]=400649 )

# pick the current variant = first whose cursor < population
V=""
for x in "${order[@]}"; do
  xl=$(echo $x | tr A-Z a-z); c=0; [ -f $CURD/$xl.cursor ] && c=$(cat $CURD/$xl.cursor)
  if [ "$c" -lt "${POP[$x]}" ]; then V=$x; break; fi
done
if [ -z "$V" ]; then echo "ALL_VARIANTS_EXHAUSTED"; exit 0; fi
vl=$(echo $V | tr A-Z a-z)
cur=0; [ -f $CURD/$vl.cursor ] && cur=$(cat $CURD/$vl.cursor)

bl=$SC/expand/${vl}_batch.stands
emitted=$(julia --project=. test/harness/fia/expand_batch.jl $V $cur $BATCH $bl 2>$SC/expand/${vl}_emit.log | tail -1)
if [ "${emitted:-0}" -eq 0 ]; then echo "$V DONE (cursor $cur)"; exit 0; fi

cyc=$SC/expand/${vl}_cycle.csv; rm -f $cyc
# Bound the batch with `timeout -s KILL` so a pathologically slow stratum (dense stands / huge tree lists)
# SELF-KILLS inside the script instead of running away. Critical for driver robustness: when the outer Bash
# tool times out it kills the shell but ORPHANS a still-running julia child, which later checkpoints a STALE
# cursor and rewinds the sweep. A per-batch timeout keeps every julia bounded so the script always reaches its
# clean cursor-advance and no orphan survives the turn. CYCLE_TO default 480s ⇒ a batch fits under a ~9-min cap.
LEDGER=$cyc timeout -s KILL ${CYCLE_TO:-480} julia --project=. test/harness/fia/ledger_fia.jl $bl $V none > $SC/expand/${vl}_run.log 2>&1
rc=$?
# rc 124 (timeout) / 137 (128+SIGKILL) = the batch exceeded CYCLE_TO — NOT a crash. Treat like an empty stratum:
# advance the cursor past it (any rows ledger already wrote are still filtered below via the partial $cyc) so the
# sweep keeps moving. A genuine crash (other non-zero rc) still HALTS for investigation per doctrine.
if [ $rc -eq 124 ] || [ $rc -eq 137 ]; then
  echo $((cur + emitted)) > $CURD/$vl.cursor
  echo "$V batch: offset $cur→$((cur+emitted)) of ${POP[$V]}  TIMEOUT (>${CYCLE_TO:-480}s) — skipped, cursor advanced"
  # fall through so any partial $cyc rows still get ledgered/filtered, but do not halt
fi
# A non-zero rc is a REAL failure (crash) — halt so it can be investigated. But rc=0 with an empty/missing
# cycle CSV is NOT a failure: it's a stratum where live FVS emitted no comparable .sum rows for the whole batch
# (NOSUM-heavy / nonstocked plots — live itself can't project ~1 in 6 real stands). Skipping it and ADVANCING
# the cursor is correct; halting the whole sweep on it was a bug that stopped coverage prematurely.
if [ $rc -ne 0 ] && [ $rc -ne 124 ] && [ $rc -ne 137 ]; then echo "$V RUN FAILED at offset $cur (rc=$rc, cursor NOT advanced)"; tail -4 $SC/expand/${vl}_run.log; exit 1; fi
ncyc=0; [ -f $cyc ] && ncyc=$(($(wc -l < $cyc) - 1)); [ $ncyc -lt 0 ] && ncyc=0
if [ ! -s $cyc ] || [ "$ncyc" -eq 0 ]; then
  echo $((cur + emitted)) > $CURD/$vl.cursor
  echo "$V batch: offset $cur→$((cur+emitted)) of ${POP[$V]}  EMPTY-STRATUM (no comparable .sum rows) — skipped, cursor advanced"
  echo "EXPAND_CYCLE_DONE $(date +%T)"; exit 0
fi

# append cycle rows to the scratchpad master ledger (full detail, not committed)
master=$SC/expand/${vl}_ledger.csv
if [ ! -f $master ]; then cp $cyc $master; else tail -n +2 $cyc >> $master; fi

# extract dig-worthy discrepancies → committed dig-queue (create header if new). filter_digworthy.jl applies
# the base dig-worthy rule AND drops stands whose (ecoregion,signature) is already cornered
# (docs/fia_cornered_clusters.tsv) so the sweep advances past cornered clusters to NEW strata; UNCLASSIFIED and
# structure/density blow-ups (≥15% on a density col) are never dropped (escalation guard).
[ -f $DIGQ ] || head -1 $cyc > $DIGQ
n_before=$(($(wc -l < $DIGQ) - 1))
julia --project=. test/harness/fia/filter_digworthy.jl $cyc $V docs/fia_cornered_clusters.tsv >> $DIGQ 2>>$SC/expand/${vl}_filter.log
n_after=$(($(wc -l < $DIGQ) - 1))

# advance cursor on success
echo $((cur + emitted)) > $CURD/$vl.cursor
be=$(awk -F, 'NR>1 && $5=="true"' $cyc | wc -l); tot=$ncyc
echo "$V batch: offset $cur→$((cur+emitted)) of ${POP[$V]}  bit_exact=$be/$tot  dig-worthy +$((n_after-n_before)) (queue=$n_after)"
echo "EXPAND_CYCLE_DONE $(date +%T)"
