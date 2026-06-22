#!/usr/bin/env bash
# gen_cut_scenarios.sh — coverage scenarios for the CUTS methods/modifiers that the
# decision flow flags as UNTESTED (SPECPREF / THINPRSC / SALVAGE / YARDLOSS and the
# other thin methods). The headline snt01.key/sn.key exercise these only in their
# thinned stands (2-4), which the suite never validated; these are single-stand keys
# (one keyword each) so the .sum diff vs Oracle A isolates the one semantic. Syntax
# is copied verbatim from sn.key/snt01.key where available (authoritative).
#
# Usage: gen_cut_scenarios.sh   → writes scenarios/cut_*.{key,tre}
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/scenarios"; mkdir -p "$OUT"
SRCKEY="${FVS_SN_KEY:-/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key}"
TRE="${SRCKEY%.key}.tre"
rm -f "$OUT"/cut_*.key "$OUT"/cut_*.tre

# single-stand base = stand-1 keywords up to (not incl.) the first PROCESS.
base="$OUT/_cutbase.txt"
awk '/^PROCESS/{exit} {print}' "$SRCKEY" > "$base"

mkcut() {  # mkcut name  <keyword line>...
  local name="$1"; shift
  { cat "$base"; printf '%s\n' "$@"; echo "PROCESS"; echo "STOP"; } > "$OUT/$name.key"
  cp "$TRE" "$OUT/$name.tre"
}

# --- cut MODIFIER + method combos (syntax verbatim from sn.key) ---
# SPECPREF reorders the THINBTA cut by species preference (sp1 low pref, sp26 high).
mkcut cut_specpref \
  "SPECPREF      2000.0       1.0     999.0" \
  "SPECPREF      2000.0      26.0    9999.0" \
  "THINBTA       2000.0     157.0"
# THINPRSC prescription thin.
mkcut cut_thinprsc "THINPRSC      2000.0     0.999"
# YARDLOSS volume-accounting modifier on a THINDBH cut.
mkcut cut_yardloss \
  "YARDLOSS                    .5        .7        .5        1." \
  "THINDBH       2000.0                  3."
# THINSDI: thin to a residual SDI at a year.
mkcut cut_thinsdi "THINSDI       2000.0     250.0"

# NOTE — these need a fixture/syntax before they exercise their semantic on snt01:
#   SALVAGE  → no-op here (snt01 has no dead/damaged trees; needs a dead-tree .tre).
#   THINHT / THINAUTO → my OPNEW field guess produced no cut in the oracle; get the
#   exact field layout from initre.f label_4050 / label_3300 before adding.
# Until then they are NOT generated (a no-op scenario falsely looks "passing").

rm -f "$base"
echo "cut scenarios: $(ls "$OUT"/cut_*.key | wc -l)"
ls "$OUT"/cut_*.key | xargs -n1 basename
