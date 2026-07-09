#!/usr/bin/env bash
# =============================================================================
# export_and_run.sh — export FVS-ready FIA plots (by CN) to standalone files, then run.
#
#   cd FVSjl && bash examples/fia/export_and_run.sh  [path/to/fia.db]  [CN1,CN2,…]
#
# Given an FVS-ready FIA SQLite database (the FVS_STANDINIT_COND / FVS_TREEINIT_COND tables)
# and a list of stand CNs, this shows the whole round trip:
#   1. EXPORT   CNs -> standalone <CN>.key + <CN>.tre     (legacy, no DB needed to run)
#   2. EXPORT   CNs -> standalone <CN>.yaml + <CN>.csv    (modern, self-describing variant)
#   3. RUN      each exported stand (prints its .sum)
#   4. CONVERT  the exported files between forms with fvsjl-translate
#   5. VALIDATE the export vs the DATABASE reader (fidelity check)
#
# This folder already ships one exported stand (CN 163384065010854) so you can inspect the
# format WITHOUT a database. Pass a real FIA .db (e.g. SQLITE_FIADB_ENTIRE.db) to regenerate.
# Full reference: docs/TOOLS.md.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/../.."          # → the FVSjl project root

JL="julia --project=."
DB="${1:-}"                          # optional path to an FVS-ready FIA SQLite DB
CNS="${2:-163384065010854}"          # comma-separated CN list (defaults to the shipped one)
OUT=$(mktemp -d)

echo "### 0. The stand shipped in examples/fia/ (no DB needed) ###############"
echo "# a materialized standalone .key (STDINFO/DESIGN/… from FVS_STANDINIT_COND):"
cat examples/fia/163384065010854.key
echo "# run it directly (reads 163384065010854.tre beside it):"
$JL bin/fvsjl-run.jl examples/fia/163384065010854.key | sed -n '1,4p'

if [ -z "$DB" ] || [ ! -f "$DB" ]; then
  echo
  echo "### (pass a real FVS-ready FIA .db as arg 1 to run steps 1–5) ##########"
  echo "#   bash examples/fia/export_and_run.sh /path/to/SQLITE_FIADB_ENTIRE.db 163384065010854,..."
  exit 0
fi

echo "### 1. EXPORT CNs -> legacy .key + .tre ###############################"
$JL bin/fvsjl-fia-export.jl "$DB" "$CNS" "$OUT/key" --format key

echo "### 2. EXPORT CNs -> modern .yaml + .csv (self-describing variant) #####"
$JL bin/fvsjl-fia-export.jl "$DB" "$CNS" "$OUT/yaml" --format yaml

FIRST="${CNS%%,*}"                   # the first CN in the list
echo "### 3. RUN an exported stand (no DB needed) ###########################"
$JL bin/fvsjl-run.jl "$OUT/key/$FIRST.key"  | sed -n '1,4p'
# the .yaml carries its own variant:, so it runs with no --variant flag:
$JL bin/fvsjl-run.jl "$OUT/yaml/$FIRST.yaml" | sed -n '1,4p'

echo "### 4. CONVERT an exported file between forms #########################"
$JL bin/fvsjl-translate.jl "$OUT/key/$FIRST.key" "$OUT/$FIRST.roundtrip.yaml"
echo "# wrote $OUT/$FIRST.roundtrip.yaml"

echo "### 5. VALIDATE the export against the DATABASE reader ################"
$JL bin/fvsjl-fia-export.jl "$DB" "$FIRST" "$OUT/val" --validate | grep -E "validate"

echo "### done. Exported/converted files under: $OUT"
