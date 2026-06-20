#!/usr/bin/env bash
# gen_codepath_scenarios.sh — scenarios chosen to hit distinct engine codepaths
# that homogeneous single-species stands can't reach.
#
# mix_*    : 2-species compositions reaching the oak-pine (okpn) and conifer
#            (sfhp) forest-type groups, so ALL 8 _dgf_forest_group branches are hit.
# sitset_* : vary SITECODE across the 9 site-index master-group representatives →
#            each SITSET A/B/C/D site-index codepath.
# spctrn_* : foreign (non-SN) species codes → the 562-row species-translation
#            crosswalk codepath.

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/scenarios"; mkdir -p "$OUT"
SRCKEY="${FVS_SN_KEY:-/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key}"
TRE="${SRCKEY%.key}.tre"
rm -f "$OUT"/mix_*.key "$OUT"/mix_*.tre "$OUT"/sitset_*.key "$OUT"/sitset_*.tre "$OUT"/spctrn_*.key "$OUT"/spctrn_*.tre

base="$OUT/_cpbase.key"
awk '/^PROCESS/{print "PROCESS"; print "STOP"; exit} {print}' "$SRCKEY" > "$base"

# --- mixes: A/B by record parity (50/50) or A/A/B (66/33) ---
mkmix() {  # mkmix name A B pattern(2|3)
    local name="$1" A="$2" B="$3" pat="$4"
    cp "$base" "$OUT/$name.key"
    if [ "$pat" = 2 ]; then
        awk -v a="$A " -v b="$B " '{print substr($0,1,33) (NR%2==0?a:b) substr($0,37)}' "$TRE" > "$OUT/$name.tre"
    else
        awk -v a="$A " -v b="$B " '{print substr($0,1,33) (NR%3==0?b:a) substr($0,37)}' "$TRE" > "$OUT/$name.tre"
    fi
}
mkmix mix_lp_wo LP WO 2; mkmix mix_sp_so SP SO 2; mkmix mix_ll_so LL SO 3
mkmix mix_lp_hi LP HI 2; mkmix mix_vp_wo VP WO 3; mkmix mix_lp_rm LP RM 2
mkmix mix_sa_yp SA YP 2; mkmix mix_sp_hi SP HI 3; mkmix mix_ll_sa LL SA 2
mkmix mix_sa_ll3 SA LL 3; mkmix mix_ll_lp LL LP 2; mkmix mix_sr_sa SR SA 2
mkmix mix_wp_hm WP HM 2

# --- SITSET: SITECODE = each master-group representative species index ---
for pair in 5:g1sp 64:g2so 45:g3yp 12:g4wp 14:g5vp 65:g6sk 74:g7co 10:g8pp; do
    sidx="${pair%%:*}"; lbl="${pair##*:}"
    awk -v s="$sidx" '/^SITECODE/{printf "SITECODE  %10d      60.\n", s; next} {print}' "$base" > "$OUT/sitset_$lbl.key"
    cp "$TRE" "$OUT/sitset_$lbl.tre"
done

# --- SPCTRN: foreign (non-SN) species codes ---
for spec in BF:balsamfir JP:jackpine NS:norwayspruce QA:quakingaspen; do
    code="${spec%%:*}"; lbl="${spec##*:}"
    cp "$base" "$OUT/spctrn_$lbl.key"
    awk -v sp="$code " '{print substr($0,1,33) sp substr($0,37)}' "$TRE" > "$OUT/spctrn_$lbl.tre"
done

rm -f "$base"
echo "mix: $(ls "$OUT"/mix_*.key|wc -l), sitset: $(ls "$OUT"/sitset_*.key|wc -l), spctrn: $(ls "$OUT"/spctrn_*.key|wc -l)"
