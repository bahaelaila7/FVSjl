#!/usr/bin/env bash
# coverage_breadth.sh — rank western variants by BREADTH of FVSjl source covered (USER 2026-09-05:
# order the enriched regime-matrix sweep widest-breadth-first, population=tiebreak). For each variant:
# build a small INDEXED sub-DB, run coverage_breadth.jl across the regime set under --code-coverage=user,
# count DISTINCT covered src lines + files, record. Output: /tmp/cov_breadth_results.tsv (variant lines files).
set -u
cd /workspace/FVSjl
export JULIA_DEPOT_PATH=/workspace/.julia_depot
NST=${COV_NSTANDS:-3}
REGIMES=${REGIMES:-none,plant,thinbba,salvage,simfire}
VARIANTS=${VARIANTS:-"CR EM UT AK SO IE WS CI EC CA WC PN BM TT NC"}
OUT=/tmp/cov_breadth_results.tsv
: > $OUT
echo -e "variant\tpop_rank\tcov_lines\tcov_files\truns_ok" | tee -a $OUT
for V in $VARIANTS; do
  vl=$(echo $V | tr 'A-Z' 'a-z')
  db=/tmp/cov_${vl}.db
  # build a small indexed sub-DB (stratified) from the master
  julia --project=. test/harness/fia/build_subdb.jl $V $NST $db > /tmp/cov_build_${vl}.log 2>&1
  if [ ! -s "$db" ]; then echo -e "$V\t-\tBUILD_FAIL\t-\t-" | tee -a $OUT; continue; fi
  find src -name '*.jl.*.cov' -delete 2>/dev/null
  COV_NSTANDS=$NST julia --project=. --code-coverage=user test/harness/fia/coverage_breadth.jl $V $db $REGIMES > /tmp/cov_run_${vl}.log 2>&1
  ok=$(grep -oE 'ok=[0-9]+' /tmp/cov_run_${vl}.log | tail -1 | cut -d= -f2)
  lines=$(cat $(find src -name '*.jl.*.cov' 2>/dev/null) 2>/dev/null | grep -cE '^ *[1-9]')
  files=$(find src -name '*.jl.*.cov' 2>/dev/null | wc -l)
  echo -e "$V\t-\t${lines:-0}\t${files:-0}\t${ok:-0}" | tee -a $OUT
  find src -name '*.jl.*.cov' -delete 2>/dev/null
  rm -f $db
done
echo "=== RANKED BY BREADTH (cov_lines desc) ===" | tee -a $OUT
sort -t$'\t' -k3 -rn $OUT | grep -vE '^variant|BUILD_FAIL|RANKED' | tee -a $OUT
