# =============================================================================
# diameter_growth.jl — Southern diameter-growth model (DGF + DGCONS)
#
# Ported from: sn/dgf.f.
#
# `dgcons!` builds the per-species site constant DGCON from site index, slope,
# aspect and the ecological-unit (PCOM) physiographic terms. `dgf!` then evaluates
# ln(DDS) (change in squared inside-bark diameter) for every tree into the scratch
# column `wk2`. The stochastic driver (DGDRIV: calibration + serial correlation)
# turns ln(DDS) into per-tree DG and is ported separately.
#
# Coefficient tables live in dgf_coefficients.jl. The species loop is flattened to
# a single per-tree pass (the Fortran species-chain ordering doesn't affect the
# per-tree value).
# =============================================================================

# Forest-type group of IFORTP → which categorical coefficient applies (dgf.f:453).
function _dgf_forest_group(ifortp::Integer)
    ifortp == 701 || ifortp == 801 || ifortp == 805 ? :nohd :
    ifortp in (104,105,121,124) ? :sfhp :
    ifortp in (141,142,161,162,163,164,166) ? :ylpn :
    ifortp in (165,403,404,405,406,407,409) ? :okpn :
    ifortp in (501,502,503,504,505,510,514,515) ? :upok :
    ifortp in (168,508,601,602,605,606,607,608,702,703,704,705,706,708) ? :lohd :
    ifortp in (103,167,181,401,402,506,507,511,512,513,519,520,802,803,807,809) ? :uphd :
    :none
end

# Parse the SN ecological-unit code (PCOM) → physiographic flag symbol (dgf.f:567).
function _dgf_phys_group(pcom::AbstractString)
    s   = rpad(pcom, 5)
    ch1 = s[1]; c234 = s[1:3]; c45 = s[4:5]; c4 = s[4]; c5 = s[5]
    if ch1 == 'M'
        s[2:4] == "221" && return :pm221
        s[2:4] in ("222","223") && return :pm222
        s[2:4] == "231" && return :pm231
        return :none
    end
    if c234 == "221"
        c45 in ("DD","DE") && return :s231t
        c45 in ("EJ","EG","EN") && return :p222
        return :p221
    elseif c234 in ("222","223")
        return :p222
    elseif c234 == "231"
        c4 in ('A','C','D','I') && return :s231t
        c4 in ('B','E','F','G') && return :s231l
        if c4 == 'H'
            c5 in ('A','B','C') && return :s231l
            c5 in ('D','E','F') && return :p222
            c5 in ('H','I') && return :p232
        end
        return :none
    elseif c234 == "232"
        return c45 == "FF" ? :p234 : :p232
    elseif c234 == "234"
        return :p234
    elseif c234 == "251"
        return :p255
    elseif c234 == "255"
        return c45 == "CH" ? :s231l : :p255
    elseif c234 == "411"
        return :p411
    end
    return :none
end

"""
    dgcons!(state)

Load the site-dependent DG constant `calib.dg_const[sp]` and `calib.atten[sp]`
(DGCONS, dgf.f:564). Call once per stand after site variables are known.
"""
function dgcons!(s::StandState)
    p, c = s.plot, s.calib
    phys = _dgf_phys_group(s.plot.eco_unit)
    cosa = cos(p.aspect); sina = sin(p.aspect)
    @inbounds for sp in 1:MAXSP
        c.atten[sp] = DGF_OBSERV[sp]
        base = DGF_ISIO[sp] * p.sp_site_index[sp] +
               DGF_TANS[sp] * p.slope +
               DGF_FCOS[sp] * p.slope * cosa +
               DGF_FSIN[sp] * p.slope * sina
        base += _phys_coef(phys, sp)
        c.dg_const[sp] = base
    end
    return s
end

# physiographic coefficient for the active group (one of PM*/P*/S231*).
@inline function _phys_coef(phys::Symbol, sp::Integer)
    phys === :pm221 ? DGF_PM221[sp] : phys === :pm222 ? DGF_PM222[sp] :
    phys === :pm231 ? DGF_PM231[sp] : phys === :p221  ? DGF_P221[sp]  :
    phys === :p222  ? DGF_P222[sp]  : phys === :s231t ? DGF_S231T[sp] :
    phys === :s231l ? DGF_S231L[sp] : phys === :p232  ? DGF_P232[sp]  :
    phys === :p234  ? DGF_P234[sp]  : phys === :p255  ? DGF_P255[sp]  :
    phys === :p411  ? DGF_P411[sp]  : 0f0
end

# forest-type categorical coefficient for the active group.
@inline function _ft_coef(grp::Symbol, sp::Integer)
    grp === :lohd ? DGF_FTLOHD[sp] : grp === :nohd ? DGF_FTNOHD[sp] :
    grp === :okpn ? DGF_FTOKPN[sp] : grp === :sfhp ? DGF_FTSFHP[sp] :
    grp === :uphd ? DGF_FTUPHD[sp] : grp === :upok ? DGF_FTUPOK[sp] :
    grp === :ylpn ? DGF_FTYLPN[sp] : 0f0
end

"""
    dgf!(state)

Evaluate ln(DDS) for every tree into `scratch.wk[2, i]` (DGF, dgf.f:436), using
the current diameters. `cor` is the per-species calibration correction (0 until
calibration runs). Requires `dgcons!` and AVH/BA already set.
"""
function dgf!(s::StandState)
    p, t, c = s.plot, s.trees, s.calib
    wk2 = view(s.scratch.wk, 2, :)
    ftgrp = _dgf_forest_group(p.forest_type)
    kplant = p.managed > 0 ? 1f0 : 0f0
    avh = p.avg_height
    @inbounds for i in 1:t.n
        d = t.dbh[i]
        d <= 0f0 && continue
        sp = t.species[i]
        conspp = c.dg_const[sp] + c.dg_cor[sp]

        relht = 0f0
        if avh > 0f0
            relht = t.height[i] / avh
            relht > 1.5f0 && (relht = 1.5f0)
        end
        ba_v  = p.basal_area <= 0f0 ? 25f0 : p.basal_area
        icr_i = t.crown_pct[i] <= 0 ? 25 : t.crown_pct[i]
        bal   = (1f0 - t.crown_ratio[i] / 100f0) * ba_v
        pba   = ba_v                                   # point BA (PTBAA) — stand BA until ported
        pbal  = pba * (1f0 - t.crown_ratio[i] / 100f0)
        pbal <= 0f0 && (pbal = bal)

        dds = conspp + DGF_INTERC[sp] +
              DGF_LDBH[sp]  * log(d) +
              DGF_DBH2[sp]  * d * d +
              DGF_LCRWN[sp] * log(Float32(icr_i)) +
              DGF_HREL[sp]  * relht +
              DGF_PLTB[sp]  * ba_v +
              DGF_PNTBL[sp] * pbal +
              _ft_coef(ftgrp, sp) +
              DGF_PLANT[sp] * kplant

        dds < -9.21f0 && (dds = -9.21f0)
        wk2[i] = dds
    end
    return s
end
