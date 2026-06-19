# =============================================================================
# site_index.jl — fan the site-species site index out to all species (SITSET)
#
# Ported from: sn/sitset.f.
#
# SITECODE gives a site index for ONE site species; SITSET converts it (via a
# 9-group master-species system) into a site index for every species, which the
# diameter- and height-growth models need (SITEAR). Also fills default max-SDI by
# species (SDIDEF). Tables copied verbatim from sitset.f.
# =============================================================================

const SN_SDICON = Float32[655, 354, 412, 499, 490, 385, 490, 332, 398, 398,
        310, 529, 480, 499, 692, 623, 518, 371, 344, 421,
        590, 371, 371, 400, 350, 375, 276, 492, 420, 422,
        257, 147, 364, 414, 408, 423, 414, 338, 492, 430,
        155, 283, 283, 430, 478, 492, 415, 492, 492, 492,
        422, 277, 726, 430, 704, 304, 164, 492, 499, 648,
        520, 384, 361, 315, 342, 405, 326, 387, 384, 326,
        417, 336, 365, 417, 414, 342, 311, 370, 410, 343,
        447, 492, 526, 282, 263, 282, 227, 354, 492, 421]
const SN_ISNSIS = Int32[5,  6, 11, 17, 64, 35, 47, 75, 78, 15,
        16, 44, 59, 76, 45,  8, 12, 13,  2,  1,
         3,  7,  4, 14, 34, 61, 65, 87, 74, 10,
        20, 22, 24, 25, 33, 60, 62, 63, 66, 69,
        71, 73, 83]
const SN_ISNGRP = Int32[1, 1, 1, 1, 2, 2, 2, 2, 2, 3,
        3, 3, 3, 3, 3, 4, 4, 4, 5, 5,
        5, 5, 5, 5, 6, 6, 6, 6, 7, 8,
        9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
        9, 9, 9]
const SN_MAPSI = Int32[5, 5, 5, 5, 1, 1, 5, 4, 9, 8,
        1, 4, 4, 5, 3, 3, 1, 9, 9, 9,
        9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
        9, 9, 9, 6, 2, 9, 9, 9, 9, 9,
        9, 9, 9, 3, 3, 9, 2, 9, 9, 9,
        9, 9, 9, 9, 9, 9, 9, 9, 3, 9,
        6, 9, 9, 2, 6, 9, 9, 9, 9, 9,
        9, 9, 9, 7, 2, 3, 9, 2, 9, 9,
        9, 9, 9, 9, 9, 9, 6, 5, 9, 9]
const SN_MGSISP = Int32[5, 64, 45, 12, 14, 65, 74, 10, 63]
const SN_SIMIN = Float32[15, 15, 15, 35, 35, 35, 45, 45, 35, 25,
        35, 40, 40, 35, 30, 30, 35, 35, 35, 35,
        30, 35, 25, 35, 35, 15, 25, 30, 15, 15,
        15, 15, 35, 35, 35, 35, 35, 25, 15, 15,
        35, 35, 35, 30, 30, 35, 25, 35, 15, 35,
        15, 15, 30, 35, 35, 15, 15, 15, 30, 40,
        30, 35, 25, 25, 25, 30, 25, 25, 35, 25,
        35, 35, 30, 25, 25, 15, 25, 25, 30, 25,
        15, 15, 35, 35, 35, 35, 35, 15, 15, 15]
const SN_SIMAX = Float32[100,  70,  80, 100, 105, 105,  90, 125,  70,  95,
        105, 135, 125,  95, 120, 120,  90,  70,  70,  85,
        105, 100,  90,  85,  70,  40,  85,  90,  90,  40,
         45,  70,  85, 105,  95,  85, 105, 120,  50,  65,
         70,  85,  85, 125, 135, 125, 115, 125,  75, 125,
         40,  55, 105, 105,  95,  40,  70,  60, 120, 125,
         90, 105, 115, 115, 115, 125,  65,  65,  95,  65,
         95,  75, 115, 115, 115, 125,  85, 115,  65,  95,
        110,  80,  90,  90,  90,  90,  90,  55,  55,  55]

# A,B coefficients by master group (sitset.f:120).
function _sitset_ab(igrp::Int, pmom::Bool, isisp::Int)
    igrp == 1 && return pmom ? (-7.1837f0, 0.1633f0) : (-10f0, 0.2f0)
    igrp == 2 && return pmom ? (-8.6809f0, 0.1702f0) :
                 (isisp == 78 ? (-16f0, 0.2667f0) : (-12f0, 0.2f0))
    igrp == 3 && return (-4f0, 0.1f0)
    igrp == 4 && return (-9.4118f0, 0.1569f0)
    igrp == 5 && return (-9.3913f0, 0.1739f0)
    igrp == 6 && return (-10f0, 0.2f0)
    igrp == 7 && return (-8.6809f0, 0.1702f0)
    igrp == 8 && return (-7.1837f0, 0.1633f0)
    igrp == 9 && return pmom ? (-8.7442f0, 0.186f0) : (-10f0, 0.2f0)
    return (0f0, 0f0)
end

# C,D coefficients by master group index (sitset.f:156).
function _sitset_cd(i::Int, pmom::Bool, isisp::Int)
    i == 1 && return pmom ? (44f0, 6.13f0) : (50f0, 5f0)
    i == 2 && return pmom ? (51f0, 5.88f0) :
              (isisp == 78 ? (60f0, 3.75f0) : (60f0, 5f0))
    i == 3 && return (40f0, 10f0)
    i == 4 && return (60f0, 6.38f0)
    i == 5 && return (54f0, 5.75f0)
    i == 6 && return (50f0, 5f0)
    i == 7 && return (51f0, 5.88f0)
    i == 8 && return (44f0, 6.13f0)
    i == 9 && return pmom ? (47f0, 5.38f0) : (50f0, 5f0)
    return (0f0, 0f0)
end

"""
    site_index_setup!(state)

SITSET: fill `plot.sp_site_index[sp]` (SITEAR) for every species from the site
species' index, and default `plot.sp_sdi_def` (SDIDEF). Call after SITECODE/STDINFO.
"""
function site_index_setup!(s::StandState)
    p = s.plot
    sea = p.sp_site_index
    isisp = Int(p.site_species)
    imapsp = findfirst(==(Int32(isisp)), SN_ISNSIS)
    if imapsp === nothing
        isisp = 63; p.site_species = Int32(63); imapsp = 38
    end
    sea[isisp] <= 0f0 && (sea[isisp] = 70f0)
    sea[isisp] < SN_SIMIN[isisp]  && (sea[isisp] = SN_SIMIN[isisp])
    sea[isisp] >= SN_SIMAX[isisp] && (sea[isisp] = SN_SIMAX[isisp])

    rsisp = (sea[isisp] - SN_SIMIN[isisp]) / (SN_SIMAX[isisp] - SN_SIMIN[isisp])
    igrp  = Int(SN_ISNGRP[imapsp])
    pmom  = !isempty(p.eco_unit) && p.eco_unit[1] == 'M'
    a, b  = _sitset_ab(igrp, pmom, Int(p.site_species))
    imgsp  = Int(SN_MGSISP[igrp])
    mgsion = rsisp * (SN_SIMAX[imgsp] - SN_SIMIN[imgsp]) + SN_SIMIN[imgsp]
    mgspix = a + b * mgsion

    mgrsi = ntuple(9) do i
        c, d = _sitset_cd(i, pmom, Int(p.site_species))
        msi = Int(SN_MGSISP[i])
        (c + d * mgspix - SN_SIMIN[msi]) / (SN_SIMAX[msi] - SN_SIMIN[msi])
    end

    @inbounds for i in 1:MAXSP
        if sea[i] == 0f0
            mi = Int(SN_MAPSI[i])
            sea[i] = mgrsi[mi] * (SN_SIMAX[i] - SN_SIMIN[i]) + SN_SIMIN[i]
            sea[i] <= SN_SIMIN[i] && (sea[i] = SN_SIMIN[i])
            sea[i] >= SN_SIMAX[i] && (sea[i] = SN_SIMAX[i])
        end
    end

    bamax = s.control.ba_max
    pmsdiu = p.pct_sdimax_mort_hi
    @inbounds for i in 1:MAXSP
        if p.sp_sdi_def[i] <= 0f0
            p.sp_sdi_def[i] = bamax > 0f0 ?
                bamax / (0.5454154f0 * (pmsdiu / 100f0)) : SN_SDICON[i]
        end
    end
    return s
end
