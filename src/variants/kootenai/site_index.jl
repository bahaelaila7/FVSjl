# =============================================================================
# site_index.jl (kootenai) — KT habitat-type mapping (habtyp.f) + SDImax defaults (sitset.f).
#
# Uses the validated tables in habitat_tables.jl. Two outputs feed later chunks:
#   KKTYPE (1..175, "Kootenai habitat type") -> chunk-3 DG (MAPHAB(KKTYPE,ISPC)->DGHAB). Stored in
#     plot.habitat_code (no shared engine code reads it; only KT DG will).
#   ITYPE (1..30, "Inland Empire habitat type") -> sitset BAMAX/SDIDEF. Stored in plot.habitat_input.
# VALIDATED end-to-end vs live FVSkt (stand 753200841290487, input PV_CODE=531):
#   KKTYPE(KOTHAB=531) -> ITYPE=14 -> IE rep 530 -> BAMAXA[14]=440 -> SDIDEF=949.
# =============================================================================

# kt/habtyp.f: input habitat code KODTYP -> (KKTYPE, ITYPE). Returns (kktype, itype).
# KKTYPE: first K in 1..175 with KODTYP < KOTHAB[K], minus 1 (0 -> 97 sentinel). Then the represented
# Kootenai code KOTHAB[KKTYPE] is mapped to the IE ITYPE via the JTYPE(95)/KTYPE(95) search.
function kt_habtyp(kodtyp_in::Integer)
    kk = findfirst(k -> kodtyp_in < KT_KOTHAB[k], 1:175)
    kktype = kk === nothing ? 97 : kk - 1
    kktype == 0 && (kktype = 97)
    # KODTYP now the represented Kootenai code (sentinel 97 has no KOTHAB entry -> keep input)
    kod = (1 <= kktype <= 175) ? Int(KT_KOTHAB[kktype]) : Int(kodtyp_in)
    jj = findfirst(k -> kod < KT_JTYPE[k], 1:95)
    jk = jj === nothing ? 96 : jj
    itype = Int(KT_KTYPE[jk - 1])
    return kktype, itype
end

function kt_site_index_setup!(s::StandState)
    p = s.plot
    kodtyp_in = Int(p.habitat_code)
    if kodtyp_in > 0
        kktype, itype = kt_habtyp(kodtyp_in)
        p.habitat_code  = Int32(kktype)   # KKTYPE for chunk-3 DG (MAPHAB)
        p.habitat_input = Int32(itype)    # ITYPE for sitset BAMAX/SDIDEF
    else
        itype = Int(p.habitat_input)      # already an ITYPE, or 0
    end
    # sitset.f: BAMAX defaults to BAMAXA[ITYPE]; SDIDEF = BAMAX/(0.5454154*(PMSDIU/100)).
    bamax = s.control.ba_max
    if bamax <= 0f0 && 1 <= itype <= 30
        bamax = KT_BAMAXA[itype]
    end
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0
    @inbounds for i in 1:nspecies(s.variant)
        if p.sp_sdi_def[i] <= 0f0 && bamax > 0f0
            p.sp_sdi_def[i] = bamax / (0.5454154f0 * pmsdiu)
        end
    end
    return s
end

site_setup!(s::StandState, ::Kootenai) = kt_site_index_setup!(s)
