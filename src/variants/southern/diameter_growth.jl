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
# physiographic group symbol → its CSV coefficient-column name.
const _DG_PHYS_COL = Dict(
    :pm221 => :dg_phys_pm221, :pm222 => :dg_phys_pm222, :pm231 => :dg_phys_pm231,
    :p221  => :dg_phys_p221,  :p222  => :dg_phys_p222,  :s231t => :dg_phys_s231t,
    :s231l => :dg_phys_s231l, :p232  => :dg_phys_p232,  :p234  => :dg_phys_p234,
    :p255  => :dg_phys_p255,  :p411  => :dg_phys_p411)

# forest-type group symbol → its CSV coefficient-column name.
const _DG_FORTYPE_COL = Dict(
    :lohd => :dg_fortype_lowland_hw, :nohd => :dg_fortype_nonoak_hw,
    :okpn => :dg_fortype_oak_pine,   :sfhp => :dg_fortype_softwd_hp,
    :uphd => :dg_fortype_upland_hw,  :upok => :dg_fortype_upland_oak,
    :ylpn => :dg_fortype_yellow_pine)

function dgcons!(s::StandState)
    p, c, sd = s.plot, s.calib, s.coef.species
    prior_obs_count = sd[:dg_prior_obs_count]
    site_coef = sd[:dg_site_index]
    slope_tan = sd[:dg_slope_tan]; slope_cos = sd[:dg_slope_cos_aspect]; slope_sin = sd[:dg_slope_sin_aspect]
    phys = _dgf_phys_group(s.plot.eco_unit)
    phys_coef = haskey(_DG_PHYS_COL, phys) ? sd[_DG_PHYS_COL[phys]] : nothing
    cosa = cos(p.aspect); sina = sin(p.aspect)
    @inbounds for sp in 1:MAXSP
        c.atten[sp] = prior_obs_count[sp]
        base = site_coef[sp] * p.sp_site_index[sp] +
               slope_tan[sp] * p.slope +
               slope_cos[sp] * p.slope * cosa +
               slope_sin[sp] * p.slope * sina
        phys_coef === nothing || (base += phys_coef[sp])
        c.dg_const[sp] = base
    end
    return s
end

"""
    dgf!(state)

Evaluate ln(DDS) for every tree into `scratch.wk[2, i]` (DGF, dgf.f:436), using
the current diameters. `cor` is the per-species calibration correction (0 until
calibration runs). Requires `dgcons!` and AVH/BA already set.
"""
function dgf!(s::StandState)
    p, t, c, sd = s.plot, s.trees, s.calib, s.coef.species
    wk2 = view(s.scratch.wk, 2, :)
    intercept = sd[:dg_intercept]; ln_dbh = sd[:dg_ln_dbh]; dbh_sq = sd[:dg_dbh_squared]
    ln_crown = sd[:dg_ln_crown_pct]; rel_ht = sd[:dg_relative_height]
    stand_ba_c = sd[:dg_stand_basal_area]; point_bal = sd[:dg_point_bal]; planted = sd[:dg_planted]
    ftgrp = _dgf_forest_group(p.forest_type)
    ft_coef = haskey(_DG_FORTYPE_COL, ftgrp) ? sd[_DG_FORTYPE_COL[ftgrp]] : nothing
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
        # DGF competition term: pbal = (point BA)·(1 − PCT/100), where PCT is the
        # CROWN-modeled crown ratio (crown_ratio, set by crown_ratio!). pba=PTBAA.
        bal   = (1f0 - t.crown_ratio[i] / 100f0) * ba_v
        pba   = s.density.point_ba[t.plot_id[i]]
        pba <= 0f0 && (pba = ba_v)
        pbal  = pba * (1f0 - t.crown_ratio[i] / 100f0)
        pbal <= 0f0 && (pbal = bal)

        dds = conspp + intercept[sp] +
              ln_dbh[sp]     * log(d) +
              dbh_sq[sp]     * d * d +
              ln_crown[sp]   * log(Float32(icr_i)) +
              rel_ht[sp]     * relht +
              stand_ba_c[sp] * ba_v +
              point_bal[sp]  * pbal +
              (ft_coef === nothing ? 0f0 : ft_coef[sp]) +
              planted[sp]    * kplant

        dds < -9.21f0 && (dds = -9.21f0)
        wk2[i] = dds
    end
    return s
end

"Empirical-Bayes prior variance for DG calibration (dgdriv.f PSIGSQ)."
const DG_PSIGSQ = 0.089827273f0

"""
    calibrate_diameter_growth!(state; scale=1f0)

DGDRIV calibration pass (LSTART, dgdriv.f:150). For each species with measured
diameter growth in the input, regress the residuals of the DGF prediction against
the measured DG to get the large-tree calibration `calib.dg_cor` (COR), then shrink
it toward 0 with an empirical-Bayes weight. Deterministic (DGSD=0 ⇒ no RNG).
`scale = YR/FINT` converts the measurement period to the cycle length (1 for snt01).
Must run before `diameter_growth!` (it reads the measured DG from `trees.diam_growth`).
"""
function calibrate_diameter_growth!(s::StandState; scale::Float32 = 1f0, fnmin::Float32 = 5f0)
    t, c = s.trees, s.calib
    bark_a = s.coef.species[:bark_intercept]; bark_b = s.coef.species[:bark_slope]
    dgf!(s)                                   # WK2 = DGF prediction at WK3=DBH, COR=0
    wk2 = view(s.scratch.wk, 2, :)

    # per-species accumulators
    dn = fill(999f0, MAXSP); dx = zeros(Float32, MAXSP)
    @inbounds for i in 1:t.n
        (t.dbh[i] < 3f0 || t.diam_growth[i] <= 0f0) && continue
        sp = t.species[i]
        t.dbh[i] < dn[sp] && (dn[sp] = t.dbh[i])
        t.dbh[i] > dx[sp] && (dx[sp] = t.dbh[i])
    end
    spopn = zeros(Float32, MAXSP); spopx = zeros(Float32, MAXSP)
    dev = zeros(Float32, MAXSP); devsq = zeros(Float32, MAXSP); fn = zeros(Float32, MAXSP)
    snp = zeros(Float32, MAXSP); snx = zeros(Float32, MAXSP); sny = zeros(Float32, MAXSP)
    snxx = zeros(Float32, MAXSP); snxy = zeros(Float32, MAXSP)
    @inbounds for i in 1:t.n
        sp = t.species[i]; wk3 = t.dbh[i]; dg = t.diam_growth[i]; p = t.tpa[i]
        (wk3 < dn[sp] || wk3 > dx[sp]) && continue
        edds = exp(wk2[i]); spopn[sp] += p; spopx[sp] += edds * p
        dg <= 0f0 && continue
        bark = bark_ratio(bark_a, bark_b, sp, wk3)
        term = dg * (2f0 * bark * wk3 + dg) * scale
        term <= 0f0 && continue
        reslog = log(term) - wk2[i]
        fn[sp] += 1f0; dev[sp] += reslog; devsq[sp] += reslog^2
        snp[sp] += p; snx[sp] += p * edds; sny[sp] += p * reslog
        snxx[sp] += p * edds^2; snxy[sp] += p * reslog * edds
    end

    @inbounds for sp in 1:MAXSP
        (fn[sp] < fnmin || snp[sp] <= 0f0) && continue     # FNMIN observations to calibrate
        bpopx = spopx[sp] / spopn[sp]; bnx = snx[sp] / snp[sp]; bny = sny[sp] / snp[sp]
        csnxy = snxy[sp] - bnx * bny * snp[sp]
        csnxx = snxx[sp] - bnx * bnx * snp[sp]
        csnxx < 0f0 && continue
        slop = csnxy / csnxx
        sdpred = sqrt(csnxx / (snp[sp] * (1f0 - 1f0 / fn[sp])))
        dist = abs(bpopx - bnx) / sdpred
        regcor = bny + (bpopx - bnx) * slop
        cornew = dist > 3f0 ? bny :
                 dist <= 1f0 ? regcor :
                 bny * (dist / 2f0) + regcor * (1f0 - dist / 2f0)
        # empirical-Bayes shrinkage toward 0
        svar = devsq[sp] - dev[sp]^2 / fn[sp]
        svar_v = (svar / (fn[sp] - 1f0)) / fn[sp]
        temp = min(cornew * cornew / DG_PSIGSQ, 72f0)
        wc = 1f0 / (1f0 + exp(-0.5f0 * temp) * sqrt(svar_v / DG_PSIGSQ))
        c.dg_cor[sp] = wc * cornew
    end
    return s
end

"""
    diameter_growth!(state, ::Southern)

Variant hook: compute each tree's periodic diameter growth into `trees.diam_growth`
(DGDRIV, sn/dgdriv.f). For the deterministic case (no DGSTDEV → DGSD=0 so the
serial-correlation multiplier is 1, no growth multipliers → xdgrow=0, no
calibration → COR=0) this reduces to `DG = sqrt(d_ib² + exp(ln DDS)) − d_ib`,
bounded. (Calibration + serial correlation are added when a test stand needs them.)
"""
function diameter_growth!(s::StandState, ::Southern)
    dgf!(s)
    t = s.trees
    sd = s.coef.species
    bark_a = sd[:bark_intercept]; bark_b = sd[:bark_slope]
    dlo_v = sd[:dg_bound_dbh_lo]; dhi_v = sd[:dg_bound_dbh_hi]
    wk2 = view(s.scratch.wk, 2, :)
    @inbounds for i in 1:t.n
        sp = t.species[i]
        bark = bark_ratio(bark_a, bark_b, sp, t.dbh[i])
        d_ib = t.dbh[i] * bark
        dds  = exp(wk2[i])                                  # xdgrow = log(XDMULT)=0
        dg   = sqrt(d_ib * d_ib + dds) - d_ib
        t.diam_growth[i] = dg_bound(dlo_v, dhi_v, sp, t.dbh[i], dg, s.control.sp_size_cap)
    end
    return s
end
