#!/usr/bin/env bash
# triage_ledger.sh — isolate REAL-BUG candidates from a fia_ledger CSV produced by ledger_fia.jl.
#
# The ledger's `signature` column is a deterministic bucket over MEASURED facts. Most non-bit-exact
# stands fall into KNOWN-CORNERED buckets (count_straddle, threshold_crossing, print_boundary,
# volume_persistent, structure_densephase) — documented named primitives (compounded-ULP volume drift,
# self-thinning count-straddle, dense-phase growth-ranking), NOT bugs. Those buckets were both-sides-traced
# in the SN 5000 pass and cornered; the ledger tracks their COUNTS so a later fix can be diffed for flips.
#
# This script separates two things:
#   (A) SIGNATURE HISTOGRAM  — population of every bucket (the flip-detection baseline).
#   (B) REAL-BUG CANDIDATES  — the TIGHT filter that isolates a NEW logic bug from the cornered noise. The
#       discriminator (learned from the SN FORKOD zero-volume bug, slice 41) is a BASE-cubic-volume
#       divergence with CLEAN structure and LARGE magnitude — i.e. not a merch/board threshold-col straddle
#       and not explained by a structure move:
#         1. UNCLASSIFIED                                  — matched no deterministic rule at all.
#         2. worst_col=TCuFt & struct%<1 & max_rel>=5%     — total-cubic diverges while structure is EXACT
#                                                            (the zero-vol / vol-eq-misassignment signature).
#       TCuFt (base cubic) is the discriminator: every tree with DBH>0 contributes, so a material TCuFt gap
#       with clean structure means the volume EQUATION is wrong (not a merch/board threshold-col straddle,
#       not a structure move propagating into volume). The merch cols (MCuFt/SCuFt/BdFt) legitimately zero
#       out below merch size — a 100% there is a cornered threshold_crossing, converges, ratio exact — so a
#       generic "big vol%" filter is NOISE. The struct%<1 gate excludes structure_densephase (structure
#       moved ⇒ volume moving with it is expected). Validated: 0 candidates on the post-fix SN 5000 ledger
#       (would have flagged the FORKOD zero-vol bug pre-fix). Widen only if a flip is suspected.
#
# Usage: triage_ledger.sh <ledger.csv> [...more csvs]
# Columns: 1 variant,2 regime,3 cn,4 n_cycles,5 bit_exact,6 div_cols,7 worst_col,8 worst_cycle,
#          9 max_rel_pct,10 max_abs_diff,11 struct_max_rel_pct,12 vol_max_rel_pct,13 density_bitexact,
#          14 converges,15 signature
set -euo pipefail
[ $# -ge 1 ] || { echo "usage: triage_ledger.sh <ledger.csv> [...]" >&2; exit 2; }

CAND='NR>1 && ( $15=="UNCLASSIFIED" || \
                ($7=="TCuFt" && $11+0<1.0 && $9+0>=5.0) )'

for f in "$@"; do
  [ -f "$f" ] || { echo "MISSING: $f" >&2; continue; }
  echo "================================================================"
  echo "LEDGER: $f"
  total=$(($(wc -l < "$f") - 1))
  echo "rows (stands): $total"
  echo "--- signature histogram (cornered-class populations; flip baseline) ---"
  awk -F, 'NR>1{c[$15]++} END{for(s in c) printf "  %-22s %d\n", s, c[s]}' "$f" | sort -k2 -rn
  echo "--- REAL-BUG CANDIDATES (tight filter; need both-sides trace) ---"
  awk -F, "$CAND {
             print \"  [\"\$15\"] cn=\"\$3\" worst=\"\$7\"@\"\$8\" max_rel=\"\$9\"% max_abs=\"\$10 \
                   \" struct%=\"\$11\" vol%=\"\$12\" dens_be=\"\$13\" conv=\"\$14 }" "$f"
  n=$(awk -F, "$CAND" "$f" | wc -l)
  echo "  => $n real-bug candidate(s)"
done
