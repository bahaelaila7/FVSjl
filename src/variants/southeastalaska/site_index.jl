# =============================================================================
# site_index.jl (southeastalaska) — AK site index + SDImax + forest code (ak/sitset.f + ak/forkod.f).
#
#   ak_forkod!(p)       — KODFOR → IFOR subscript (ak/forkod.f JFOR crosswalk; default Tongass 1005/IFOR 2).
#   ak_sitset!(s)       — SITEAR site-index translation (ak/sitset.f DO 10) + SDIDEF max-SDI (SDICON) load.
#   site_setup!(::SoutheastAlaska)
#
# SITEAR (ak/sitset.f): TEM=70 default (or SITEAR(ISISP) if the site species has one); ISISP defaults to 11
# (WH). For each species SITEAR(I)=SLO(I)+(TEM−SLO(ISISP))/(SHI(ISISP)−SLO(ISISP))·(SHI(I)−SLO(I)), with
# TEM floored at SLO(ISISP). VERIFIED bit-exact vs live FVSak (akt01 XSITE dump: YC=50, SS=82.5, LP=35,
# RC=57.5, WH=70, MH=42.5 — exactly the ISISP=11/TEM=70 midpoint interpolation).
#
# SDIDEF (ak/sitset.f DO 40): SDIDEF(I)=SDICON(I) (Shaw & Long) when the user set no SDIMAX/BAMAX; if a
# stand BAMAX is set, SDIDEF(I)=BAMAX/(0.5454154·PMSDIU/100). SDICAL (crown.jl) then BA-weights these into
# the per-point XMAXPT for the crown-ratio point relative density.
# =============================================================================

# ak/forkod.f DATA JFOR — accepted AK location codes; index = IFOR subscript.
const AK_JFOR = Int[1004, 1005, 703, 713, 720, 7400, 7401, 7402, 7403, 7404, 7405,
                    7406, 7407, 7408, 8134, 8135, 8112]

"ak/forkod.f — KODFOR → IFOR (1..17); default Tongass (IFOR 2 / 1005) when unrecognized or 0."
function ak_forkod!(p)
    kod = Int(p.user_forest_code)
    ifor = 2                                   # grinit default IFOR=2 (Tongass 1005)
    if kod == 1002 || kod == 1003              # old Tongass area codes → Tongass
        ifor = 2
    elseif kod == 701                          # BC/Makah combined → British Columbia
        ifor = 3
    else
        idx = findfirst(==(kod), AK_JFOR)
        idx !== nothing && (ifor = idx)        # not-found ⇒ keep default 2 (forkod ERRGRO path)
    end
    p.forest_idx = Int32(ifor)
    p.user_forest_code = Int32(AK_JFOR[ifor])  # KODFOR = JFOR(IFOR)
    return ifor
end

# ak/sitset.f DATA — SDICON (Shaw & Long max SDI), SLO/SHI (min/max site index per species).
const AK_SDICON = Float32[790., 602., 592., 387., 412., 412., 500., 654., 679., 762.,
                          682., 687., 412., 441., 441., 466., 466., 384., 562., 452.,
                          447., 447., 452.]
const AK_SLO = Float32[35., 20., 25., 25., 30., 30., 20., 40., 15., 30.,
                       35., 20., 30., 20., 75., 35., 45., 45., 45., 45.,
                       20., 30., 20.]
const AK_SHI = Float32[105., 65., 75., 55., 85., 85., 55., 125., 55., 85.,
                       105., 65., 85., 50., 135., 75., 90., 105., 85., 105.,
                       50., 70., 50.]

"ak/sitset.f — load SITEAR (site index) and SDIDEF (max SDI) for species not set by keyword."
function ak_sitset!(s::StandState)
    p = s.plot
    isisp = Int(p.site_species)
    tem = 70f0                                             # TEM = 70 default
    if isisp > 0 && p.sp_site_index[isisp] > 0f0
        tem = p.sp_site_index[isisp]
    end
    isisp == 0 && (isisp = 11)                             # sitset.f: ISISP default = 11 (WH)
    p.site_species = Int32(isisp)
    @inbounds for i in 1:23
        tem < AK_SLO[isisp] && (tem = AK_SLO[isisp])       # NB: floors TEM inside the loop (matches Fortran DO 10)
        p.sp_site_index[i] <= 0f0 &&
            (p.sp_site_index[i] = AK_SLO[i] +
                (tem - AK_SLO[isisp]) / (AK_SHI[isisp] - AK_SLO[isisp]) * (AK_SHI[i] - AK_SLO[i]))
    end
    # SDIDEF: SDICON unless a stand BAMAX is set (ak/sitset.f DO 40; PMSDIU/100=0.85 grinit default).
    bamax = s.control.ba_max
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0   # fraction (grinit PMSDIU=85 ⇒ /100)
    @inbounds for i in 1:23
        p.sp_sdi_def[i] > 0f0 && continue
        p.sp_sdi_def[i] = bamax > 0f0 ? bamax / (0.5454154f0 * pmsdiu) : AK_SDICON[i]
    end
    return s
end

function site_setup!(s::StandState, ::SoutheastAlaska)
    ak_forkod!(s.plot)
    ak_sitset!(s)
    return s
end
