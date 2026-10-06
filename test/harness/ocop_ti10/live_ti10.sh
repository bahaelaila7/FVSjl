#!/bin/bash
# live_ti10.sh V STAND REGIME  (V = OC|OP) — run the LIVE private tiered oracle (/workspace/.{oc,op}work/tiered/FVS{oc,op}_tiered, the
# crash-fixed relink of the stock objects) on the committed tiered fixture key at TIMEINT 10 (FINT=10 != the 5-yr YR), 4 cycles,
# with TREELIDB. Output: $OUT/V_STAND_REGIME/live/{out_live.db,s.sum,s.key,...}. Then build_ti10_fixture.jl turns it into goldens.
set -e
V=$1; ST=$2; RG=$3; OUT=${OUT:-/workspace/.ocop4/ti10}
R=$(cd "$(dirname "$0")/../.." && pwd)                      # .../test
FX=$R/fixtures/tiered/$(echo $V | tr A-Z a-z)
D=$OUT/${V}_${ST}_${RG}; rm -rf $D; mkdir -p $D/live
cp $FX/stands.db $D/live/
sed 's/^out.db/out_live.db/; s/^NUMCYCLE         5.0/NUMCYCLE         4.0\nTIMEINT           0.0      10.0\nTREELIST           0/' $FX/${ST}_${RG}.key | awk '!(/^TREELIST/ && seen++ )' > $D/live/s.key
grep -q TREELIDB $D/live/s.key || sed -i 's/^SUMMARY$/SUMMARY\nTREELIDB/' $D/live/s.key
BIN=/workspace/.$(echo $V | tr A-Z a-z)work/tiered/FVS$(echo $V | tr A-Z a-z)_tiered
(cd $D/live && $BIN --keywordfile=s.key >/dev/null 2>&1)
ls $D/live
