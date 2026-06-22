#!/usr/bin/env bash
# gen_estab_scenarios.sh — regeneration/establishment COVERAGE scenarios. The only
# establishment exercise is snt01's 5th stand ("BARE GROUND PLANT": NOTREES + ESTAB
# block). Extract it standalone (bare_plant) and a NATURAL-keyword variant
# (bare_natural, the natural-process form). Both regenerate to 800 TPA at cyc1 in
# Oracle A. Targets for the ESTAB port. (No .tre — NOTREES.)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; OUT="$HERE/scenarios"
SRC="${FVS_SN_KEY:-/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key}"
# the BARE GROUND stand = the keyword block from its STDIDENT to STOP
awk '/BARE GROUND/{f=1} f{print} f&&/^STOP/{exit}' "$SRC" | \
  awk 'BEGIN{print "STDIDENT"} {print}' > "$OUT/bare_plant.key"
# NATURAL variant: overwrite cols 1-7 (keeps field columns), retitle
awk '/^PLANT /{print "NATURAL" substr($0,8); next} /BARE GROUND PLANT/{print "         BARE GROUND NATURAL"; next} {print}' \
  "$OUT/bare_plant.key" > "$OUT/bare_natural.key"
echo "establishment scenarios: bare_plant, bare_natural"
