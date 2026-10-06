#!/usr/bin/env bash
# lane.sh <lane> <V:lo:hi> [<V:lo:hi> ...]  — resumable east full-population lane. Cursors in lanes/<lane>.<V>.cursor
# Env: BATCH(1000) DIGCAP(100) CYCLE_TO(2400)
set -u
W=/workspace/.wt-fiaeast; E=/workspace/.sweep_east
source /workspace/env.sh; cd $W
lane=$1; shift
BATCH=${BATCH:-1000}; DIGCAP=${DIGCAP:-100}; CYCLE_TO=${CYCLE_TO:-2400}
export SWEEP_DB=$E/fia_sweep_east_2026-10.$lane.db SKIP_DONE=1
DIGQ=$E/dig_queue.csv
for spec in "$@"; do
  IFS=: read V lo hi <<< "$spec"; vl=$(echo $V|tr A-Z a-z)
  cf=$E/lanes/$lane.$vl.cursor; [ -f $cf ] || echo $lo > $cf
  while true; do
    cur=$(cat $cf); [ "$cur" -ge "$hi" ] && { echo "$lane $V DONE $cur"; break; }
    dign=0; [ -f $DIGQ ] && dign=$(($(wc -l < $DIGQ)-1))
    [ "$dign" -ge "$DIGCAP" ] && { echo "DIG_PAUSE $lane digq=$dign"; exit 3; }
    n=$BATCH; [ $((cur+n)) -gt $hi ] && n=$((hi-cur))
    bl=$E/lanes/$lane.batch.stands; cyc=$E/lanes/$lane.cycle.csv; rm -f $cyc
    sed -n "$((cur+1)),$((cur+n))p" $E/pop_$vl.txt | awk -v v=$V '{print $1"\t"v}' > $bl
    LEDGER=$cyc timeout -s KILL $CYCLE_TO nice -n 10 julia --project=. test/harness/fia/ledger_fia.jl $bl $V none > $E/lanes/$lane.run.log 2>&1
    rc=$?
    if [ $rc -ne 0 ] && [ $rc -ne 137 ]; then echo "$lane RUN FAILED rc=$rc at $V $cur"; tail -4 $E/lanes/$lane.run.log; exit 1; fi
    [ $rc -eq 137 ] && echo "$V,$cur,$((cur+n)),timeout,$CYCLE_TO" >> $E/skipped_ranges.csv
    if [ -s $cyc ]; then
      m=$E/lanes/$lane.$vl.ledger.csv; if [ ! -f $m ]; then cp $cyc $m; else tail -n +2 $cyc >> $m; fi
      ( flock 9; [ -s $DIGQ ] || head -1 $cyc > $DIGQ
        julia --project=. test/harness/fia/filter_digworthy.jl $cyc $V docs/fia_cornered_clusters.tsv >> $DIGQ 2>>$E/lanes/$lane.filter.log ) 9>$E/dig.lock
    fi
    echo $((cur+n)) > $cf
    echo "$(date +%T) $lane $V $cur→$((cur+n))/$hi $(tail -1 $E/lanes/$lane.run.log | cut -c1-110)" >> $E/lanes/$lane.progress
  done
done
echo "LANE_DONE $lane"
