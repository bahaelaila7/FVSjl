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
    p.forest_idx <= 0 && (p.forest_idx = Int32(1))
    # SDImax (Stage): SDIDEF = BAMAX/(0.5454154·PMSDIU/100). BAMAX from user (control.ba_max) or the
    # BC site-series default (sitset.f SELECT CASE Zone/SubZone/Series, in m²/ha → convert to ft²/ac).
    # ⚠ full sitset table is chunk-2 TODO; ICHmw2/01 (all_BC, the hardcoded zone) = 89 m²/ha (oracle .out
    # "MAXIMUM BASAL AREA FOR ICHmw2/01 IS SET BY DEFAULT TO 89.0 SQ M/HA"). 89/0.2295643 = 387.69 ft²/ac.
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 85.0f0
    bamax = s.control.ba_max > 0f0 ? s.control.ba_max : 89.0f0 / BC_FT2pACRtoM2pHA
    @inbounds for sp in 1:15
        p.sp_sdi_def[sp] <= 0f0 && (p.sp_sdi_def[sp] = bamax / (0.5454154f0 * (pmsdiu / 100f0)))
    end
    return s
end

site_setup!(s::StandState, ::BritishColumbia) = bc_site_index_setup!(s)
