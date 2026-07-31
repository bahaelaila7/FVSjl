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

# kt/forkod.f: translate the user forest LOCATION code (KODFOR) into KT's location subscript KOTFOR
# (1..10, default 8 "to avoid blowups in REGENT & DGF" per forkod.f:125). KOTFOR feeds the DG/REGENT
# MAPLOC(KOTFOR,ISPC) location-class lookup; stored in p.forest_idx (kt_dgcons! reads it there).
function kt_forkod!(p)
    kodfor = Int(p.user_forest_code)
    # reservation pseudo-code crosswalk (forkod.f:56-71)
    kodfor == 8109 && (kodfor = 11300000)   # Kootenai off-res trust -> Kaniksu 113
    kodfor == 8133 && (kodfor = 11000000)   # Flathead reservation   -> Flathead 110
    kodfor == 8137 && (kodfor = 11800000)   # Coeur d'Alene res.     -> St. Joe  118
    # pad shorter forest location codes with trailing zeros (forkod.f:78-93)
    if 100 <= kodfor <= 9999999
        if     kodfor <= 999;    kodfor *= 100000
        elseif kodfor <= 9999;   kodfor *= 10000
        elseif kodfor <= 99999;  kodfor *= 1000
        elseif kodfor <= 999999; kodfor *= 100
        else;                    kodfor *= 10
        end
    end
    # KOTFOR translation (forkod.f:125-166) — sequential IFs, later matches override earlier
    kotfor = 8
    inter = kodfor - 10000000
    inter < 0 && (inter = 0)
    ifore = inter ÷ 100000
    idist = inter ÷ 1000
    icomp = inter - idist * 1000
    (ifore == 10 || ifore == 0) && (kotfor = 8)
    ifore == 14 && (kotfor = 7)
    idist == 1402 && (kotfor = 2)
    idist == 1403 && (kotfor = 3)
    idist == 406  && (kotfor = 7)
    (idist == 1404 || idist == 407) && (kotfor = 4)
    idist == 1401 && (kotfor = (16 <= icomp <= 27) ? 1 : 2)
    idist == 1405 && (kotfor = ((1 <= icomp <= 4) || (8 <= icomp <= 19) || icomp == 27) ? 5 : 10)
    idist == 1406 && (kotfor = (1 <= icomp <= 4) ? 6 : 9)
    p.forest_idx = Int32(kotfor)
    return p
end

function kt_site_index_setup!(s::StandState)
    p = s.plot
    kt_forkod!(p)                          # KOTFOR (forest_idx) for the DG/REGENT MAPLOC lookup
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
