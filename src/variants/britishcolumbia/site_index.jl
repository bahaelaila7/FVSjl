# =============================================================================
# site_index.jl (britishcolumbia) — BC site/habitat (bc/habtyp.f + bc/becset.f + bc/sitset.f). Chunk 2.
# BC's full site model parses a BEC (Biogeoclimatic) STRING → {Region,Zone,SubZone} → PrettyName →
# KODTYP → ITYPE + zone-keyed BAMAX (becset.f, ~1000 lines). PORTED HERE: the DG-relevant core —
# the KODTYP→ITYPE table (becset.f:1116, NI-habitat codes = IE MTYPE) + SDIDEF=BAMAX/(0.5454154·PMSDIU/100)
# (Stage, like IE). The BEC-STRING PARSER (becset PrettyName→KODTYP) is DEFERRED: synthetic BC stands
# supply a numeric KODTYP (STDINFO habitat) + BAMAX, which is all the DG/mortality need.
# =============================================================================

# becset.f:1116 SELECT CASE(KODTYP)→ITYPE. KODTYP are NI-habitat codes (= subset of IE MTYPE).
const BC_KODTYP_ITYPE = Dict{Int,Int}(130=>1, 250=>4, 310=>7, 320=>8, 420=>10, 470=>11,
                                      620=>19, 640=>20, 660=>21, 670=>22, 730=>27)

"bc/becset.f: KODTYP (NI-habitat code) → ITYPE (1..30). DEFAULT 21."
bc_kodtyp_itype(kodtyp::Integer) = get(BC_KODTYP_ITYPE, Int(kodtyp), 21)

function bc_site_index_setup!(s::StandState)
    p = s.plot
    kodtyp = Int(p.habitat_code)
    itype = kodtyp > 0 ? bc_kodtyp_itype(kodtyp) : max(Int(p.habitat_input), 1)
    (itype < 1 || itype > 30) && (itype = 21)
    p.habitat_input = Int32(itype)                            # DG reads MAPHAB(ITYPE) (IE structure)
    p.forest_idx <= 0 && (p.forest_idx = Int32(4))    # bc/grinit.f:201 IFOR=4; bc/forkod.f never changes it
    # SDImax (Stage): SDIDEF = BAMAX/(0.5454154·PMSDIU/100). BAMAX from user (control.ba_max) or the
    # BC site-series default (sitset.f SELECT CASE Zone/SubZone/Series, in m²/ha → convert to ft²/ac).
    # ⚠ full sitset table is chunk-2 TODO; ICHmw2/01 (all_BC, the hardcoded zone) = 89 m²/ha (oracle .out
    # "MAXIMUM BASAL AREA FOR ICHmw2/01 IS SET BY DEFAULT TO 89.0 SQ M/HA"). 89/0.2295643 = 387.69 ft²/ac.
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 85.0f0
    z, sz, ser, _ = bc_becset(s)                              # BEC zone (default ICHmw2/01); sitset BAMAX table
    bamax = s.control.ba_max > 0f0 ? s.control.ba_max : bc_sitset_bamax(z, sz, bc_iseries(ser))
    @inbounds for sp in 1:15
        p.sp_sdi_def[sp] <= 0f0 && (p.sp_sdi_def[sp] = bamax / (0.5454154f0 * (pmsdiu / 100f0)))
    end
    return s
end

site_setup!(s::StandState, ::BritishColumbia) = bc_site_index_setup!(s)

# canada/bc cratet.f:70-110 — at CRATET (after SITSET and NOTRE) the PY (10) site index goes to a 50-YEAR age base
# (Alexander, Tackle & Dahms 1967, RM-29): TEMCCF = Σ CCFCAL(mode 1)·PROB over records 1..IREC1, floored at 125;
# SITEAR(10) = 9.89311 − 0.19177·50 + 0.00124·50² − 0.00082·(TEMCCF−125)·SI + 0.01387·50·SI − 0.0000455·50²·SI.
# The same block as tt|ci|ut|ie/cratet.f (BC keeps only CASE(10)); run from setup_growth! on the NOTRE-expanded PROB.
function bc_cratet_site_adjust!(s::StandState)
    p, t = s.plot, s.trees
    temccf = 0f0
    @inbounds for i in 1:t.n
        temccf += bc_tree_ccf(Int(t.species[i]), t.dbh[i], t.tpa[i])
    end
    temccf < 125f0 && (temccf = 125f0)
    si = p.sp_site_index[10]
    p.sp_site_index[10] = 9.89311f0 - 0.19177f0 * 50f0 + 0.00124f0 * 50f0^2 -
        0.00082f0 * (temccf - 125f0) * si + 0.01387f0 * 50f0 * si - 0.0000455f0 * 50f0^2 * si
    return s
end
