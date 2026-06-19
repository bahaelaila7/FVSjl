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

# SDI/site-index defaults (sdi_max_default, site_index_min/max, site_group) and the
# site-species lookup tables (ISNSIS/ISNGRP/MGSISP) are loaded from CSV — see
# data/southern/{species_coefficients,site_species,site_master_group}.csv.

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
    isnsis = s.coef.site_species_index; isngrp = s.coef.site_master_group
    mgsisp = s.coef.master_group_rep
    simin = s.coef.species[:site_index_min]; simax = s.coef.species[:site_index_max]
    mapsi = s.coef.species[:site_group];     sdicon = s.coef.species[:sdi_max_default]
    sea = p.sp_site_index
    isisp = Int(p.site_species)
    imapsp = findfirst(==(Int32(isisp)), isnsis)
    if imapsp === nothing
        isisp = 63; p.site_species = Int32(63); imapsp = 38
    end
    sea[isisp] <= 0f0 && (sea[isisp] = 70f0)
    sea[isisp] < simin[isisp]  && (sea[isisp] = simin[isisp])
    sea[isisp] >= simax[isisp] && (sea[isisp] = simax[isisp])

    rsisp = (sea[isisp] - simin[isisp]) / (simax[isisp] - simin[isisp])
    igrp  = Int(isngrp[imapsp])
    pmom  = !isempty(p.eco_unit) && p.eco_unit[1] == 'M'
    a, b  = _sitset_ab(igrp, pmom, Int(p.site_species))
    imgsp  = Int(mgsisp[igrp])
    mgsion = rsisp * (simax[imgsp] - simin[imgsp]) + simin[imgsp]
    mgspix = a + b * mgsion

    mgrsi = ntuple(9) do i
        c, d = _sitset_cd(i, pmom, Int(p.site_species))
        msi = Int(mgsisp[i])
        (c + d * mgspix - simin[msi]) / (simax[msi] - simin[msi])
    end

    @inbounds for i in 1:MAXSP
        if sea[i] == 0f0
            mi = Int(mapsi[i])
            sea[i] = mgrsi[mi] * (simax[i] - simin[i]) + simin[i]
            sea[i] <= simin[i] && (sea[i] = simin[i])
            sea[i] >= simax[i] && (sea[i] = simax[i])
        end
    end

    bamax = s.control.ba_max
    pmsdiu = p.pct_sdimax_mort_hi
    @inbounds for i in 1:MAXSP
        if p.sp_sdi_def[i] <= 0f0
            p.sp_sdi_def[i] = bamax > 0f0 ?
                bamax / (0.5454154f0 * (pmsdiu / 100f0)) : sdicon[i]
        end
    end
    return s
end
