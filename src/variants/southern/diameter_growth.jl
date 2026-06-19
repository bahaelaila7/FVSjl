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

DGDRIV calibration pass (LSTART, dgdriv.f:150). For each species with enough
measured diameter growth, regress the DGF residuals to get the large-tree
calibration `calib.dg_cor` (COR, empirical-Bayes shrunk) and the attenuation goal
`calib.dg_cor_goal` (WCI). Also seeds the per-tree serial-correlation residual
`trees.old_random` (OLDRN) — measured trees get their residual, calibrated species
fill the rest by regression, and uncalibrated species draw from BACHLO — and the
per-species `calib.vardg` (VARDG). Consumes the RNG, so trees are walked in FVS's
species-sorted order. `scale = YR/FINT` (1 for snt01). Run before `diameter_growth!`.
"""
function calibrate_diameter_growth!(s::StandState; scale::Float32 = 1f0, fnmin::Float32 = 5f0)
    t, c = s.trees, s.calib
    sd = s.coef.species
    bark_a = sd[:bark_intercept]; bark_b = sd[:bark_slope]; sigmar = sd[:dg_resid_sd]
    isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
    species_sort!(s)

    # Backdate diameters to the start of the measured-growth period (DENSE/LBKDEN,
    # dense.f:70-86): WK3 = sqrt(d²·r). For a tree with measured DG, r=(d−DG/bark)²/d²
    # (so WK3 = past inside-bark-adjusted dbh); unmeasured trees use the stand-average
    # ratio bagr. The calibration DGF must predict from this PAST stand state (past
    # dbh + past BA/AVH/PCT), which is what makes COR/OLDRN bit-exact.
    saved_dbh = Float32[t.dbh[i] for i in 1:t.n]
    bagr = 0f0; nb = 0f0
    @inbounds for i in 1:t.n
        g = t.diam_growth[i]; g <= 0f0 && continue
        d = t.dbh[i]; gadj = g / bark_ratio(bark_a, bark_b, t.species[i], d)
        gadj > d && continue
        bagr += 1f0 - (2f0 * d * gadj - gadj * gadj) / (d * d); nb += 1f0
    end
    nb > 0f0 && (bagr /= nb)
    @inbounds for i in 1:t.n
        d = t.dbh[i]; g = t.diam_growth[i]
        r = bagr
        if g > 0f0
            gadj = min(g / bark_ratio(bark_a, bark_b, t.species[i], d), d)
            rr = 1f0 - (2f0 * d * gadj - gadj * gadj) / (d * d)
            (gadj >= 0f0 && rr > 0f0) && (r = rr)
        end
        t.dbh[i] = sqrt(d * d * r)
    end
    compute_density!(s)                       # past-stand BA/AVH/point_ba/PCT from WK3
    dgf!(s)                                   # WK2 = DGF prediction at the PAST stand
    wk2 = view(s.scratch.wk, 2, :)

    # calibration VMLT (autcor LSTART: new=old=floor(YR))
    yr = Int(floor(s.control.year)); yr < 1 && (yr = 1)
    _, vmlt = autcor(yr, yr); c.vmlt = vmlt

    # per-species DBH range + endpoint predictions over measured trees
    dn = fill(999f0, MAXSP); dx = zeros(Float32, MAXSP)
    pn = zeros(Float32, MAXSP); px = zeros(Float32, MAXSP)
    @inbounds for i in 1:t.n
        (t.dbh[i] < 3f0 || t.diam_growth[i] <= 0f0) && continue
        sp = t.species[i]
        if t.dbh[i] < dn[sp]; dn[sp] = t.dbh[i]; pn[sp] = exp(wk2[i]); end
        if t.dbh[i] > dx[sp]; dx[sp] = t.dbh[i]; px[sp] = exp(wk2[i]); end
    end

    # per-species sums; remember each measured tree's residual
    reslog_t = zeros(Float32, t.n); measured = falses(t.n)
    spopn = zeros(Float32, MAXSP); spopx = zeros(Float32, MAXSP)
    dev = zeros(Float32, MAXSP); devsq = zeros(Float32, MAXSP); fn = zeros(Float32, MAXSP)
    snp = zeros(Float32, MAXSP); snx = zeros(Float32, MAXSP); sny = zeros(Float32, MAXSP)
    snxx = zeros(Float32, MAXSP); snxy = zeros(Float32, MAXSP)
    @inbounds for i in 1:t.n
        sp = t.species[i]; wk3 = t.dbh[i]; dg = t.diam_growth[i]; p = t.tpa[i]
        (wk3 < dn[sp] || wk3 > dx[sp]) && continue
        edds = exp(wk2[i]); spopn[sp] += p; spopx[sp] += edds * p
        dg <= 0f0 && continue
        bark = bark_ratio(bark_a, bark_b, sp, saved_dbh[i])   # bark at CURRENT dbh (dgdriv.f:225)
        term = dg * (2f0 * bark * wk3 + dg) * scale
        term <= 0f0 && continue
        reslog = log(term) - wk2[i]
        reslog_t[i] = reslog; measured[i] = true
        fn[sp] += 1f0; dev[sp] += reslog; devsq[sp] += reslog^2
        snp[sp] += p; snx[sp] += p * edds; sny[sp] += p * reslog
        snxx[sp] += p * edds^2; snxy[sp] += p * reslog * edds
    end

    # per-species COR / SIGMA / regression line / VARDG
    slop = zeros(Float32, MAXSP); bnx = zeros(Float32, MAXSP); bny = zeros(Float32, MAXSP)
    calibrated = falses(MAXSP)
    @inbounds for sp in 1:MAXSP
        c.sigma[sp] = sigmar[sp]                      # SIGMA=SIGMAR unless calibrated (dgdriv.f:196)
        if isct[sp, 1] != 0 && fn[sp] >= fnmin && snp[sp] > 0f0
            bnxv = snx[sp] / snp[sp]; bnyv = sny[sp] / snp[sp]
            csnxy = snxy[sp] - bnxv * bnyv * snp[sp]
            csnxx = snxx[sp] - bnxv * bnxv * snp[sp]
            if csnxx >= 0f0
                slp = csnxy / csnxx
                bpopx = spopx[sp] / spopn[sp]
                sdpred = sqrt(csnxx / (snp[sp] * (1f0 - 1f0 / fn[sp])))
                dist = abs(bpopx - bnxv) / sdpred
                regcor = bnyv + (bpopx - bnxv) * slp
                cornew = dist > 3f0 ? bnyv : dist <= 1f0 ? regcor :
                         bnyv * (dist / 2f0) + regcor * (1f0 - dist / 2f0)
                svar = devsq[sp] - dev[sp]^2 / fn[sp]
                svar_v = (svar / (fn[sp] - 1f0)) / fn[sp]
                temp = min(cornew * cornew / DG_PSIGSQ, 72f0)
                wc = 1f0 / (1f0 + exp(-0.5f0 * temp) * sqrt(svar_v / DG_PSIGSQ))
                corv = wc * cornew
                # out-of-range trap (cortem = exp(COR))
                if exp(corv) < 0.0821f0 || exp(corv) > 12.1825f0
                    corv = 0f0
                end
                c.dg_cor[sp] = corv
                # modified residual SD only for calibrated species (dgdriv.f:323-325)
                c.sigma[sp] = sqrt((svar + c.atten[sp] * sigmar[sp]^2) / (fn[sp] + c.atten[sp]))
                slop[sp] = slp; bnx[sp] = bnxv; bny[sp] = bnyv; calibrated[sp] = true
            end
        end
        vtemp = exp(c.sigma[sp]^2)
        c.vardg[sp] = (vtemp - 1f0) * vtemp / vmlt
    end

    # seed OLDRN in FVS species-sorted RNG order
    oldrn = t.old_random
    @inbounds for sp in 1:MAXSP
        i1 = isct[sp, 1]; i1 == 0 && continue
        i2 = isct[sp, 2]
        if calibrated[sp]
            rx = bny[sp] + (px[sp] - bnx[sp]) * slop[sp]
            rn = bny[sp] + (pn[sp] - bnx[sp]) * slop[sp]
            for k in i1:i2
                i = ind1[k]
                if measured[i]
                    oldrn[i] = reslog_t[i]
                else
                    oldrn[i] = bny[sp] + (exp(wk2[i]) - bnx[sp]) * slop[sp]
                    t.dbh[i] < dn[sp] && (oldrn[i] = rn)
                    t.dbh[i] > dx[sp] && (oldrn[i] = rx)
                end
            end
        else
            bound = DG_DGSD * c.sigma[sp]
            for k in i1:i2
                i = ind1[k]
                z = 0f0
                while true
                    z = bachlo(s.rng, 0f0, c.sigma[sp])
                    z <= bound && break
                end
                oldrn[i] = z
            end
        end
    end

    # COR attenuation goal (WCI) + clamp OLDRN to ±DGSD·SIGMA
    @inbounds for sp in 1:MAXSP
        c.dg_cor_goal[sp] = 0.5f0 * c.dg_cor[sp]
    end
    @inbounds for i in 1:t.n
        lim = DG_DGSD * c.sigma[t.species[i]]
        oldrn[i] > lim && (oldrn[i] = lim)
        oldrn[i] < -lim && (oldrn[i] = -lim)
    end

    # restore current diameters + current-stand density (the backdating was local)
    @inbounds for i in 1:t.n; t.dbh[i] = saved_dbh[i]; end
    compute_density!(s)
    return s
end

"""
    diameter_growth!(state, ::Southern; sfint=5f0)

Variant hook: compute each tree's periodic diameter growth into `trees.diam_growth`
(DGDRIV growth mode, sn/dgdriv.f). `DG = sqrt(d_ib² + exp(ln DDS)·frm) − d_ib`,
where `frm = exp(raw)`. While tripling is active (snt01 early cycles, `tripling=true`)
the serial-correlation factor is DETERMINISTIC — `raw = FM·ssigma·rhocp + corr·OLDRN`
(dgdriv.f:90,117, FM=−0.14228); once tripling stops it is the stochastic `dgscor!`
(BACHLO + AR(1)). COR is attenuated each cycle toward its goal (WCI). Trees are
walked in species-sorted order to keep any RNG draws bit-exact. `sfint = IY[icyc+1]−IY[1]`.
"""
const DG_FM = -0.14228f0      # tripling mid-record variance factor (dgdriv.f FM)

function diameter_growth!(s::StandState, ::Southern; sfint::Float32 = 5f0,
                          tripling::Bool = true)
    t, c = s.trees, s.calib
    sd = s.coef.species
    bark_a = sd[:bark_intercept]; bark_b = sd[:bark_slope]
    dlo_v = sd[:dg_bound_dbh_lo]; dhi_v = sd[:dg_bound_dbh_hi]
    isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
    oldrn = t.old_random

    # attenuate COR toward the calibration goal before predicting (dgdriv.f:76-79)
    cormlt = exp(-0.02773f0 * sfint)
    @inbounds for sp in 1:MAXSP
        c.dg_cor[sp] = c.dg_cor_goal[sp] + cormlt * c.dg_cor_goal[sp]
    end

    species_sort!(s)
    dgf!(s)
    wk2 = view(s.scratch.wk, 2, :)

    # per-cycle ARMA multipliers (cyc1: new=old=floor(YR); multi-cycle TODO)
    yr = Int(floor(s.control.year)); yr < 1 && (yr = 1)
    covmlt, vmlt = autcor(yr, yr)
    pvmlt = c.vmlt > 0f0 ? c.vmlt : vmlt
    corr = covmlt / sqrt(vmlt * pvmlt)

    @inbounds for sp in 1:MAXSP
        i1 = isct[sp, 1]; i1 == 0 && continue
        i2 = isct[sp, 2]
        vardg = c.vardg[sp]
        evarp1 = (sqrt(1f0 + 4f0 * vardg * pvmlt) + 1f0) / 2f0
        sig1   = sqrt(log(max(evarp1, 1f0 + eps(Float32))))
        evarp2 = (sqrt(1f0 + 4f0 * vardg * vmlt) + 1f0) / 2f0
        ssigma = sqrt(log(max(evarp2, 1f0 + eps(Float32))))
        rho = (sig1 > 0f0 && ssigma > 0f0) ?
              log(1f0 + corr * sqrt((evarp1 - 1f0) * (evarp2 - 1f0))) / (sig1 * ssigma) : 0f0
        rhocp = sqrt(max(1f0 - rho * rho, 0f0))
        frmbase = DG_FM * ssigma * rhocp
        for k in i1:i2
            i = ind1[k]
            bark = bark_ratio(bark_a, bark_b, sp, t.dbh[i])
            d_ib = t.dbh[i] * bark
            if tripling
                frmt = frmbase + corr * oldrn[i]           # deterministic (dgdriv.f:117)
                oldrn[i] = frmt
                frm = exp(frmt)
            else
                frm = dgscor!(s.rng, oldrn, i, ssigma, rho, rhocp, wk2[i])
            end
            dds  = exp(wk2[i])                              # xdgrow = log(XDMULT)=0
            dg   = sqrt(d_ib * d_ib + dds * frm) - d_ib
            t.diam_growth[i] = dg_bound(dlo_v, dhi_v, sp, t.dbh[i], dg, s.control.sp_size_cap)
        end
    end
    return s
end
