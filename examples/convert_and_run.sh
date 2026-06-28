#!/usr/bin/env bash
# =============================================================================
# convert_and_run.sh — run the FVSjl examples and convert between every form.
#
#   cd FVSjl && bash examples/convert_and_run.sh
#
# Shows the four things you can do with a stand:
#   1. RUN it           (.key / .yaml — the engine reads either; .tre / .csv trees)
#   2. CONVERT keywords  .key  <-> .yaml
#   3. CONVERT trees     .tre  <-> .csv
#   4. UNRAVEL the semantic YAML (format: fvs-stand) -> .key for stock FVS
#
# Full reference: docs/FORMATS.md (formats) and docs/KEYWORDS.md (every keyword).
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."          # → the FVSjl project root

JL="julia --project=."
TR="$JL bin/fvsjl-translate.jl"
OUT=$(mktemp -d)
# The SN tree FORMAT used by these stands (needed only for .tre <-> .csv on a stand
# whose TREEFMT is non-default — pass it as the 3rd arg to the translator).
FMT="(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,T52,I2,T66,5I1,T54,7I1,T75,F3.0)"

echo "### 1. RUN a stand (prints the .sum) ###################################"
# Legacy .key (reads thinba.tre beside it):
$JL -e 'using FVSjl; print(run_keyfile("examples/thinba/thinba.key"))'
# Semantic YAML + CSV trees (fully modern; reads thinsdi.csv beside it):
$JL -e 'using FVSjl; print(run_keyfile("examples/semantic/thinsdi.yaml"))'
# Multi-stand semantic file → three .sum blocks:
$JL -e 'using FVSjl; print(run_keyfile("examples/semantic/multistand.yaml"))'

echo "### 2. CONVERT keywords  .key <-> .yaml ###############################"
# .key → order-aware keyword-stream YAML (add --flat for the legacy flat list):
$TR examples/thinba/thinba.key "$OUT/thinba.yaml"
# .yaml → .key (feed stock FVS); works for BOTH yaml flavors (stream + semantic):
$TR "$OUT/thinba.yaml"              "$OUT/thinba.key"
$TR examples/semantic/thinsdi.yaml "$OUT/thinsdi_from_semantic.key"

echo "### 3. CONVERT trees  .tre <-> .csv ###################################"
# .tre → named-column CSV (pass the stand's FORMAT for a non-default layout):
$TR examples/thinba/thinba.tre "$OUT/thinba.csv" "$FMT"
# .csv → .tre (back to fixed-column for stock FVS):
$TR "$OUT/thinba.csv"          "$OUT/thinba.tre" "$FMT"

echo "### 4. UNRAVEL semantic YAML → .key (canonical FVS order) #############"
$TR examples/semantic/multistand.yaml "$OUT/multistand.key"
echo "----- examples/semantic/multistand.yaml unravels to: -----"
sed -n '1,16p' "$OUT/multistand.key"

echo "### done — artifacts in $OUT ##########################################"
