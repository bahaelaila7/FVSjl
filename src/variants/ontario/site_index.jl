# =============================================================================
# ontario/site_index.jl — Ontario site-index setup (canada/on/sitset.f) + forkod.
#
# Produces the per-species SITEAR (p.sp_site_index) the Penner large-tree DGF reads as
# SIM = SITEAR(sp)·FTtoM, plus SDIDEF (p.sp_sdi_def, the SDImax defaults mortality/density
# read) and ISISP (p.site_species). Ontario's sitset is the Lake-States family fan-out:
#   • default site species = 3 (Red pine) when unset; default SI = 60 ft when unset
#     (canada/on/sitset.f:386-388).
#   • SITEAR(I) = SICOEF1(ISISP,I) + SICOEF2(ISISP,I)·SITEAR(ISISP) for every unset species
#     (DO 5), then an aspen (species 41) second pass (DO 10), then a final fallback to
#     SITEAR(ISISP) (DO 15). SICOEF1/SICOEF2 = ON_SICOEF1/2 (72×72, site_coefficients.jl).
#   • SDIDEF(I) = (BAMAX>0 ? BAMAX : BAMAX1(I)) / (0.5454154·PMSDIU/100), PMSDIU default 85
#     (grinit.f:294), BAMAX1 per-species (canada/on/sitset.f:53-64).
#   • LONT top-height branch (canada/on/sitset.f:390-477): only KODFOR∈{915,916} with a
#     NEGATIVE SITEAR(ISISP) (user entered a metric TOP HEIGHT, not a site index) and
#     OSP(ISISP)>0 — the Ontario closed-form / interpolated SI equations. INERT for a normal
#     positive site index (the ont01 default path).
#   • forkod (canada/on/forkod.f + grinit.f:240): KODFOR match against JFOR; the ON default
#     when unmatched is IFOR=9 / KODFOR=915 (grinit IFOR=9, KODFOR=0 ⇒ forkod remaps to 915);
#     Manistee 924 ⇒ IFOR=3.
#
# VALIDATED bit-exact (Float32-hex) vs instrumented FVSon_g16 on ont01: SITEAR[1..72] +
# SDIDEF[1..72] all 72/72, ISISP=3, PMSDIU=85, BAMAX=0 (scratchpad/on/validate_sitset.jl vs
# /workspace/.onwork/FVSon_sitdump). The LONT branch is a faithful transcription but has no
# top-height fixture yet ⇒ NOT dump-validated (inert for ont01 and all positive-SI stands).
# =============================================================================

# canada/on/forkod.f: accepted FVS location codes for ON (LS 902..913 + ON 915/916 + 924).
const ON_JFOR = Int32[902, 903, 904, 906, 907, 909, 910, 913, 915, 916, 924]

"""
    on_forkod!(p)

canada/on/forkod.f: map the user forest code to `IFOR` (forest_idx) and remap KODFOR
(user_forest_code) to the canonical JFOR value. An unrecognized/absent code keeps the ON
grinit default IFOR=9 ⇒ KODFOR=915 (grinit.f:240-241). Manistee (924, JFOR idx 11) folds to
Huron-Manistee (904, IFOR=3). IGL (geo_location) is set to 1 only when a valid code was found.
"""
function on_forkod!(p)
    kodfor = Int(p.user_forest_code)
    ifor = 9                                   # grinit.f:240 default (ON forest 915)
    idx = findfirst(==(Int32(kodfor)), ON_JFOR)
    useigl = idx !== nothing
    useigl && (ifor = idx)
    ifor == 11 && (ifor = 3)                   # Manistee 924 → Huron-Manistee 904 (forkod.f:64-68)
    p.forest_idx = Int32(ifor)
    p.user_forest_code = ON_JFOR[ifor]         # KODFOR = JFOR(IFOR)
    useigl && (p.geo_location = Int32(1))       # KFOR is all 1
    return ifor
end

# canada/on/sitset.f OSP(72): FVS species → Ontario SI equation index (1..17; 0 = none).
const ON_OSP = Int[
    14,15,15,15,16,13,13, 1, 2, 3,
     2, 9, 2, 2, 4, 4,17, 9, 9, 5,
     6, 6, 6, 7, 8, 9, 9, 9,10,11,
    11,11,11,11,11,11, 9, 9, 9,17,
    17,17,12, 9, 9, 9, 9, 9, 9, 9,
     9, 9, 9, 9, 9, 9, 9, 5, 9, 9,
     5, 5, 5,17,17,17, 9, 4,14,16,
    13, 2]

# canada/on/sitset.f A1..A5 / THT1..4 (17 SI-equation species) — LONT top-height branch only.
const ON_SI_A1 = Float32[0.0061, 9.0023, 0.6464, 0.2388, 0.1073, 0.1898, 0.1817, 0.1921, 0.1984, 0.1728, 0.1692, 0.5119, 40.6506, 23.7086, 22.5992, 24.4573, 28.0109]
const ON_SI_A2 = Float32[1.3539, 1.4753, 1.000, 1.1583, 1.3455, 1.2186, 1.243, 1.201, 1.2089, 1.256, 1.2648, 1.0229, 5.6605, 20.5596, 18.7782, 18.1841, 24.5365]
const ON_SI_A3 = Float32[-0.00019, 0.7996, -0.0225, -0.0102, -0.007, -0.011, -0.011, -0.01, -0.011, -0.011, -0.011, -0.0167, 1.2544, 17.3764, 15.1961, 12.4653, 12.4653]
const ON_SI_A4 = Float32[-1.0286, 0.3976, -1.1129, -1.8455, -3.3034, -2.6865, -3.0184, -2.3009, -2.4917, -3.3605, -3.4334, -1.0284, -0.1567, 14.2208, 11.5513, 6.7618, 16.9224]
const ON_SI_A5 = Float32[-0.0723, 19.26275, 0.0, -0.1883, -0.3899, -0.2717, -0.318, -0.2331, -0.2542, -0.3452, -0.3557, -0.0049, 0.0, 0.0, 0.0, 0.0, 0.0]
const ON_SI_THT1 = Float32[0,0,0,0,0,0,0,0,0,0,0,0,0, 22.56, 23.0, 24.42, 28.53]
const ON_SI_THT2 = Float32[0,0,0,0,0,0,0,0,0,0,0,0,0, 19.46, 19.78, 19.24, 25.02]
const ON_SI_THT3 = Float32[0,0,0,0,0,0,0,0,0,0,0,0,0, 16.36, 16.56, 14.06, 21.12]
const ON_SI_THT4 = Float32[0,0,0,0,0,0,0,0,0,0,0,0,0, 13.32, 13.04, 8.3, 17.3]

@inline _on_logf(x::Float32) = ccall(:logf, Float32, (Float32,), x)
@inline _on_expf(x::Float32) = ccall(:expf, Float32, (Float32,), x)
@inline _on_powf(x::Float32, y::Float32) = ccall(:powf, Float32, (Float32, Float32), x, y)

# canada/on/sitset.f:390-477 — Ontario top-height→SI equations (LONT branch). `ksp`=OSP(ISISP),
# `htneg`=SITEAR(ISISP) (negative; the entered top height in ft). Returns SIM (metres) or 0.
function _on_lont_sim(ksp::Int, htneg::Float32)::Float32
    agebh = (ksp >= 3 && ksp <= 12) ? 4.0f0 : 6.0f0
    ieqag = (ksp > 1 && ksp <= 11 && ksp != 7) ? 50 : 56
    agedif = Float32(ieqag) - agebh
    htm = -htneg * 0.3048f0                                  # FTtoM
    if ksp == 2                                              # black spruce
        (htm <= 1.3f0 || agedif <= 0f0) && return 0f0
        return ON_SI_A1[ksp] + 0.4396f0*(htm-1.3f0) + ON_SI_A2[ksp]*_on_logf(htm-1.3f0) -
               ON_SI_A3[ksp]*_on_logf(agedif) - ON_SI_A4[ksp]*_on_powf(_on_logf(agedif),2.0f0) +
               ON_SI_A5[ksp]*(htm-1.3f0)/agedif
    elseif ksp == 3                                          # tamarack
        return (4.5f0 + ON_SI_A1[ksp]*(htm*3.28f0-4.5f0) *
                _on_powf(1.0f0-_on_expf(ON_SI_A3[ksp]*agedif), ON_SI_A4[ksp]))/3.28f0
    elseif ksp == 1                                          # balsam fir
        return (4.5f0 + ON_SI_A1[ksp]*_on_powf(htm*3.28f0, ON_SI_A2[ksp]) *
                _on_powf(1.0f0-_on_expf(ON_SI_A3[ksp]*agedif),
                         ON_SI_A4[ksp]*_on_powf(htm*3.28f0, ON_SI_A5[ksp])))/3.28f0
    elseif ksp == 13                                         # white spruce (double precision KD)
        a1 = Float64(ON_SI_A1[ksp]); a2 = Float64(ON_SI_A2[ksp])
        a3 = Float64(ON_SI_A3[ksp]); a4 = Float64(ON_SI_A4[ksp]); hd = Float64(htm)
        kd = 1.0 - (1.0/(a1*(hd-1.3)^(a2+1)))^(1.0/(a3*(hd-1.3)^a4))
        return Float32(1.3 + (1.0/(a1*(hd-1.3)^a2*(1.0-kd^(Float64(agedif)/50.0))^(a3*(hd-1.3)^a4))))
    elseif ksp < 14                                          # most other defined species
        return (4.5f0 + ON_SI_A1[ksp]*_on_powf(htm*3.28f0-4.5f0, ON_SI_A2[ksp]) *
                _on_powf(1.0f0-_on_expf(ON_SI_A3[ksp]*agedif),
                         ON_SI_A4[ksp]*_on_powf(htm*3.28f0-4.5f0, ON_SI_A5[ksp])))/3.28f0
    else                                                    # ksp>=14 interpolate/extrapolate the 4 site types
        local slp, diff, simin
        if (htm >= ON_SI_THT4[ksp] && htm < ON_SI_THT3[ksp]) || htm < ON_SI_THT4[ksp]
            slp = (ON_SI_A3[ksp]-ON_SI_A4[ksp])/(ON_SI_THT3[ksp]-ON_SI_THT4[ksp]); diff = htm-ON_SI_THT4[ksp]; simin = ON_SI_A4[ksp]
        elseif htm >= ON_SI_THT3[ksp] && htm < ON_SI_THT2[ksp]
            slp = (ON_SI_A2[ksp]-ON_SI_A3[ksp])/(ON_SI_THT2[ksp]-ON_SI_THT3[ksp]); diff = htm-ON_SI_THT3[ksp]; simin = ON_SI_A3[ksp]
        else
            slp = (ON_SI_A1[ksp]-ON_SI_A2[ksp])/(ON_SI_THT1[ksp]-ON_SI_THT2[ksp]); diff = htm-ON_SI_THT2[ksp]; simin = ON_SI_A2[ksp]
        end
        return diff*slp + simin
    end
end

"""
    on_sitset!(s::StandState)

canada/on/sitset.f main body: fill SITEAR (per-species site index) + SDIDEF (per-species
SDImax). Bit-exact vs FVSon_g16 on ont01 (SITEAR/SDIDEF 72/72 Float32-hex).
"""
function on_sitset!(s::StandState)
    p = s.plot; v = s.variant
    sea = p.sp_site_index
    isisp = Int(p.site_species)
    kodfor = Int(p.user_forest_code)
    lont = (kodfor == 915 || kodfor == 916) && isisp >= 1 && isisp <= length(ON_OSP) && ON_OSP[isisp] > 0
    # default site species = red pine (3); default SI = 60 (sitset.f:386-388)
    isisp <= 0 && (isisp = 3; p.site_species = Int32(3))
    if sea[isisp] == 0f0 || (sea[isisp] < 0f0 && !lont)
        sea[isisp] = 60f0
    end
    # LONT top-height branch (sitset.f:390-477): user entered a metric top height (SITEAR<0)
    if lont && sea[isisp] < 0f0
        ksp = ON_OSP[isisp]
        sim = _on_lont_sim(ksp, sea[isisp])
        sim > 0f0 && (sea[isisp] = sim * 3.28084f0)          # MtoFT
    end
    sisp = sea[isisp]
    @inbounds for i in 1:nspecies(v)                          # DO 5
        sea[i] <= 0.0001f0 && (sea[i] = ON_SICOEF1[isisp, i] + ON_SICOEF2[isisp, i] * sisp)
    end
    if sea[41] > 0.0001f0                                     # DO 10 — aspen second pass
        s41 = sea[41]
        @inbounds for i in 1:nspecies(v)
            sea[i] <= 0.0001f0 && (sea[i] = ON_SICOEF1[41, i] + ON_SICOEF2[41, i] * s41)
        end
    end
    # SDIDEF (sitset.f:505-514) — per-species BAMAX1 (or the keyword BAMAX) / (0.5454154·PMSDIU/100)
    bamax = s.control.ba_max
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 85f0
    @inbounds for i in 1:nspecies(v)                          # DO 15
        sea[i] < 0.0001f0 && (sea[i] = sisp)
        if p.sp_sdi_def[i] <= 0f0
            num = bamax > 0f0 ? bamax : ON_BAMAX1[i]
            p.sp_sdi_def[i] = num / (0.5454154f0 * (pmsdiu / 100f0))
        end
    end
    return s
end

function on_site_index_setup!(s::StandState)
    on_forkod!(s.plot)
    on_sitset!(s)
    return s
end

site_setup!(s::StandState, ::Ontario) = on_site_index_setup!(s)
