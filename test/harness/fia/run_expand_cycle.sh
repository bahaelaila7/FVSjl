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
SC=/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad
BATCH=${BATCH:-2000}
CURD=test/harness/fia/expand
DIGQ=docs/fia_dig_queue.csv
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
LEDGER=$cyc julia --project=. test/harness/fia/ledger_fia.jl $bl $V none > $SC/expand/${vl}_run.log 2>&1
rc=$?
if [ $rc -ne 0 ] || [ ! -f $cyc ]; then echo "$V RUN FAILED at offset $cur (rc=$rc, cursor NOT advanced)"; tail -4 $SC/expand/${vl}_run.log; exit 1; fi

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
be=$(awk -F, 'NR>1 && $5=="true"' $cyc | wc -l); tot=$(($(wc -l < $cyc)-1))
echo "$V batch: offset $cur→$((cur+emitted)) of ${POP[$V]}  bit_exact=$be/$tot  dig-worthy +$((n_after-n_before)) (queue=$n_after)"
echo "EXPAND_CYCLE_DONE $(date +%T)"
