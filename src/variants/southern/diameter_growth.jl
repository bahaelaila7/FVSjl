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
    bark_a0 = sd[:bark_intercept]; bark_b0 = sd[:bark_slope]
    @inbounds for sp in 1:MAXSP
        c.atten[sp] = prior_obs_count[sp]
        base = site_coef[sp] * p.sp_site_index[sp] +
               slope_tan[sp] * p.slope +
               slope_cos[sp] * p.slope * cosa +
               slope_sin[sp] * p.slope * sina
        phys_coef === nothing || (base += phys_coef[sp])
        c.dg_const[sp] = base
        c.bark_a[sp] = bark_a0[sp]; c.bark_b[sp] = bark_b0[sp]   # per-stand bark copy
    end
    # Fort Bragg (IFOR==20, dgf.f:636 / bratio.f:106): override the DG attenuation
    # count for longleaf/loblolly and the special bark equations (sp 5,6,8,11,13).
    if p.forest_idx == 20
        c.atten[8] = 2056f0; c.atten[13] = 689f0
        c.bark_a[5]  = 0.1713f0;  c.bark_b[5]  = 0.87459f0
        c.bark_a[6]  = -0.26207f0; c.bark_b[6] = 0.87347f0
        c.bark_a[8]  = -0.43439f0; c.bark_b[8] = 0.91382f0
        c.bark_a[11] = -0.62033f0; c.bark_b[11]= 0.91645f0
        c.bark_a[13] = -0.4671f0;  c.bark_b[13]= 0.90198f0
    end
    # READCORD/READCORH (dgf.f:1168 / htgf.f:332): add ln(COR2)/ln(HCOR2) to DGCON/HTCON before
    # the LSTART calibration. Default off (terms = 1 ⇒ ln 1 = 0, no-op); guarded > 0 like Fortran.
    ctl = s.control
    if ctl.dg_cor2_on
        @inbounds for sp in 1:MAXSP
            ctl.dg_cor2[sp] > 0f0 && (c.dg_const[sp] += log(ctl.dg_cor2[sp]))
        end
    end
    if ctl.htg_cor2_on
        @inbounds for sp in 1:MAXSP
            ctl.htg_cor2[sp] > 0f0 && (c.htg_cor[sp] += log(ctl.htg_cor2[sp]))
        end
    end
    return s
end

"""
    dgf!(state)

Evaluate ln(DDS) for every tree into `scratch.wk[2, i]` (DGF, dgf.f:436), using
the current diameters. `cor` is the per-species calibration correction (0 until
calibration runs). Requires `dgcons!` and AVH/BA already set.
"""
function dgf!(s::StandState, ::Southern)
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

        # Fort Bragg (IFOR==20, dgf.f:515-537): special longleaf(8)/loblolly(13)
        # diameter-growth equations replace the standard DDS. dg5 is the inside-bark
        # DG; dds = ln(dg5·(2·d_ib + dg5)) re-encodes it for the sqrt growth formula.
        if p.forest_idx == 20 && (sp == 8 || sp == 13)
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d)
            dib  = d * bark
            site = p.sp_site_index[sp]
            pctf = Float32(t.crown_ratio[i])              # PCT (BA percentile)
            cr   = Float32(icr_i) / 100f0
            dg5 = sp == 8 ?
                dib * (-0.4553f0 * (0.09737f0 - exp(-0.2428f0 * d)) + 0.05574f0 * cr -
                       0.0002965f0 * ba_v - 0.00002481f0 * pba -
                       0.001192f0 * (pctf / 100f0)^(-0.9663f0) +
                       0.0010110f0 * site - 0.007711f0 * relht) :
                dib * (-0.3428f0 * (-0.1741f0 - exp(-0.1328f0 * d)) + 0.1145f0 * cr -
                       0.0001682f0 * ba_v - 0.00003978f0 * pba -
                       0.159400f0 * (pctf / 100f0)^(-0.1299f0) +
                       0.0006204f0 * site + 0.02474f0 * relht)
            dg5 < 0.01f0 && (dg5 = 0.01f0)
            dds = log(dg5 * (2f0 * dib + dg5))
        end

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
# GROWTH IDG/IHTG=1/3 (sn/cratet.f:169-180): the input DG/HTG field is a PAST (code 1) or CURRENT
# (code 3) DBH/HT measurement, not the growth increment. Convert it to the increment
# `DG = Q·(DBH − field)` with Q=+1 for code 1 (field = past, current in DBH) and Q=−1 for code 3
# (field = current, past in DBH). A field ≤0 is a MISSING measurement → −1 sentinel (every measured-
# DG consumer filters `≤0`, so the tree falls back to the stand-average backdating); the sentinel is
# essential — without it field=0 would yield the bogus increment DBH−0 = DBH. The diameter
# increment here is OUTSIDE-bark (matching cratet feeding DENSE); the outside→inside BRATIO bark
# correction is applied later, in `calibrate_diameter_growth!` (dgdriv.f:330-333). Height has no
# bark, so the HTG increment is final. Gated: a no-op for the default IDG/IHTG=0.
function apply_growth_input_types!(s::StandState)
    t = s.trees
    idg = s.control.growth_idg
    if idg == 1 || idg == 3
        q = idg == 1 ? 1f0 : -1f0
        @inbounds for i in 1:t.n
            f = t.diam_growth[i]
            t.diam_growth[i] = f <= 0f0 ? -1f0 : q * (t.dbh[i] - f)
        end
    end
    ihtg = s.control.growth_ihtg
    if ihtg == 1 || ihtg == 3
        q = ihtg == 1 ? 1f0 : -1f0
        @inbounds for i in 1:t.n
            f = t.ht_growth[i]
            t.ht_growth[i] = f <= 0f0 ? -1f0 : q * (t.height[i] - f)
        end
    end
    return
end

# Backdate live diameters to the start of the measured-growth period, IN PLACE
# (DENSE/LBKDEN, dense.f:70-128). WK3 = sqrt(d²·r): for a measured-DG tree,
# r = 1−(2·d·gadj − gadj²)/d² (the past inside-bark dbh); unmeasured trees fall back to
# the stand-average BA-growth ratio BAGR. IDG-faithful (dense.f:100-127): exclude only
# MISSING growth (IDG=1/3 missing is already −1 → a genuine 0 is KEPT; IDG=0/2 keep
# input-0 as missing); bark-divide the increment for IDG∈{0,2,3}, SKIP for IDG==1
# (dense.f:101,122). Shared by calibrate_diameter_growth! and init_crown_ratios! so the
# crown's backdated CCF and the DG calibration read ONE backdating — as FVS runs one DENSE.
function _backdate_dbh!(s::StandState)
    t = s.trees; n = t.n
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    idg = s.control.growth_idg
    ismiss = (idg == 1 || idg == 3) ? (g -> g < 0f0) : (g -> g <= 0f0)
    bagr = 0f0; nb = 0f0
    @inbounds for i in 1:n
        g = t.diam_growth[i]; ismiss(g) && continue
        d = t.dbh[i]
        gadj = idg == 1 ? g : g / bark_ratio(bark_a, bark_b, t.species[i], d)
        gadj > d && continue
        bagr += 1f0 - (2f0 * d * gadj - gadj * gadj) / (d * d); nb += 1f0
    end
    nb > 0f0 && (bagr /= nb)
    @inbounds for i in 1:n
        d = t.dbh[i]; g = t.diam_growth[i]; r = bagr
        if !ismiss(g)
            gadj = idg == 1 ? g : min(g / bark_ratio(bark_a, bark_b, t.species[i], d), d)
            rr = 1f0 - (2f0 * d * gadj - gadj * gadj) / (d * d)
            rr > 0f0 && (r = rr)
        end
        r > 0f0 && (t.dbh[i] = sqrt(d * d * r))
    end
    return s
end

function calibrate_diameter_growth!(s::StandState; scale::Float32 = 1f0, fnmin::Float32 = 5f0)
    t, c = s.trees, s.calib
    sd = s.coef.species
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b; sigmar = sd[:dg_resid_sd]
    isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
    species_sort!(s)

    # Backdate diameters to the start of the measured-growth period (DENSE/LBKDEN,
    # dense.f:70-86): WK3 = sqrt(d²·r). For a tree with measured DG, r=(d−DG/bark)²/d²
    # (so WK3 = past inside-bark-adjusted dbh); unmeasured trees use the stand-average
    # ratio bagr. The calibration DGF must predict from this PAST stand state (past
    # dbh + past BA/AVH/PCT), which is what makes COR/OLDRN bit-exact.
    saved_dbh = Float32[t.dbh[i] for i in 1:t.n]
    # NOTRE inflates DEAD-record PROB by FINT/FINTM (cycle-growth period / mortality-observation period) so the
    # recent dead are added back at the right rate to recover the BACKDATED density (notre.f:122-124). FVS keeps
    # that inflation live through calibration and "undoes" it for treelist/FFE uses (FMSSEE/PRTRLS); jl carries
    # the TRUE dead TPA everywhere, so it applies the inflation ONLY here, scoped to the two calibration density
    # passes below, and restores it before returning. Inert when FINTM=FINT (the default).
    _fintr = s.control.growth_fintm > 0f0 ? s.control.growth_fint / s.control.growth_fintm : 1f0
    _livec = t.n
    _saved_dead_tpa = (_fintr != 1f0 && t.ndead > 0) ?
        Float32[t.tpa[j] for j in (_livec + 1):(_livec + t.ndead)] : Float32[]
    if _fintr != 1f0
        @inbounds for j in (_livec + 1):(_livec + t.ndead); t.tpa[j] *= _fintr; end
    end
    # PTBAA (point basal area) in the calibration prediction uses the CURRENT-DBH
    # point BA (dgf.f:493 reads the live PTBAA the last DENSE pass filled at current
    # diameters). The percentile population is the full ITRN, so recently-dead trees
    # (history≠8) count toward PTBAA at their current dbh too (long-dead history-8 are
    # zeroed). Compute it here with the dead partition exposed at current dbh, then
    # restore. (Stand BA/AVH/PCT stay backdated; only this point-BA total is current.)
    # PTBAA uses CURRENT dbh for ALL dead trees (history 6 and 8) — a prism tree
    # contributes a fixed point-BA regardless of dbh, so even a long-dead snag counts.
    # (Only the BACKDATED percentile zeroes history-8; PTBAA does not.)
    let nlive0 = t.n
        t.n = nlive0 + t.ndead
        point_basal_area!(s)
        t.n = nlive0
    end
    cur_point_ba = copy(s.density.point_ba)
    _backdate_dbh!(s)                         # dense.f:70-128 backdating (IDG-faithful); shared w/ init_crown_ratios!
    # The backdated stand BA/AVH still include the dead trees (kept at current dbh):
    # expose the dead partition for this density pass, then restore. (PTBAA itself is
    # overridden below with the live-only current point_ba; only BA/AVH/PCT use this.)
    nlive = t.n
    saved_dead = Float32[t.dbh[j] for j in (nlive + 1):(nlive + t.ndead)]
    @inbounds for j in (nlive + 1):(nlive + t.ndead)
        t.history[j] == 8 && (t.dbh[j] = 0f0)
    end
    t.n = nlive + t.ndead
    compute_density!(s)                       # past-stand BA/AVH/point_ba/PCT
    t.n = nlive
    @inbounds for (k, j) in enumerate((nlive + 1):(nlive + t.ndead))
        t.dbh[j] = saved_dead[k]
    end
    s.density.point_ba .= cur_point_ba        # PTBAA from current DBH, live-only (above)
    # PCT (BA percentile, dense.f pass-1 PCTILE) accumulates the BACKDATED point-BA
    # weights in the CURRENT-dbh rank order (IND is fixed at setup). The percentile
    # POPULATION is the full ITRN including recently-dead trees (history≠8) at their
    # CURRENT dbh — long-dead (history 8) are in the order but zeroed. Only the live
    # trees' crown_ratio is read downstream by dgf!. compute_density! above re-sorted
    # by backdated dbh; recompute crown_ratio here.
    let nlive2 = t.n, ntot = t.n + t.ndead
        rankd = Float32[i <= nlive2 ? saved_dbh[i] : t.dbh[i] for i in 1:ntot]   # current dbh
        wk5   = Float32[i <= nlive2 ? t.dbh[i]^2 * t.tpa[i] :                     # live: backdated
                        (t.history[i] == 8 ? 0f0 : t.dbh[i]^2 * t.tpa[i]) for i in 1:ntot]  # dead: current/0
        ord = sortperm(rankd; rev = true)
        tot = sum(wk5); cum = 0f0
        if tot > 0f0
            @inbounds for k in ntot:-1:1
                ii = ord[k]
                cum += wk5[ii]
                ii <= nlive2 && (t.crown_ratio[ii] = cum / tot * 100f0)
            end
        end
    end
    if _fintr != 1f0                          # restore TRUE dead TPA (the FINT/FINTM inflation was calibration-
        @inbounds for (k, j) in enumerate((_livec + 1):(_livec + t.ndead)); t.tpa[j] = _saved_dead_tpa[k]; end
    end                                       # only — covers both density passes AND the PCTILE percentile)
    # GROWTH IDG=1/3 bark correction (dgdriv.f:330-333): the measured increment is now the
    # OUTSIDE-bark DBH difference (DBH−past), which DENSE used to backdast above. Convert it to
    # the INSIDE-bark increment with BRATIO at the CURRENT dbh (`saved_dbh`) before the calibration
    # term, so it matches the IDG=0 inside-bark basis. IDG=0/2 supply the inside increment directly.
    if s.control.growth_idg == 1 || s.control.growth_idg == 3
        @inbounds for i in 1:t.n
            t.diam_growth[i] > 0f0 &&
                (t.diam_growth[i] *= bark_ratio(bark_a, bark_b, t.species[i], saved_dbh[i]))
        end
    end

    # FORTYP is computed in the GROW path (per cycle), AFTER the LSTART calibration,
    # so the calibration prediction must NOT include the forest-type term (kuphd etc.
    # are all 0 at LSTART). Zero it for this dgf!, restore after. (dgf.f:453 reads
    # IFORTP, which is 0 until the first STKVAL/FORTYP call in the cycle loop.)
    saved_fortype = s.plot.forest_type
    s.plot.forest_type = 0
    dgf!(s, s.variant)                        # WK2 = DGF prediction at the PAST stand (variant dgf)
    s.plot.forest_type = saved_fortype
    wk2 = view(s.scratch.wk, 2, :)

    # calibration VMLT (autcor LSTART: new=old=YR, the measurement base period = 5 for SN —
    # NOT the cycle length FINT, which TIMEINT can change; the projection VMLT below uses FINT).
    # SERLCORR can override the ARMA(1,1) phi/theta; recompute BJRHO only when it does.
    bjrho = _stand_bjrho(s)
    _, vmlt = autcor(5, 5, bjrho); c.vmlt = vmlt

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
        # NOCALIB (LDGCAL=.FALSE.) suppresses the self-calibration for a species: skip the COR
        # fit so dg_cor / dg_cor_goal stay 0 (dgdriv.f:567 — correction terms not scaled).
        if s.control.dg_calib_sp[sp] && isct[sp, 1] != 0 && fn[sp] >= fnmin && snp[sp] > 0f0
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
            bound = s.control.dg_stddev_bound * c.sigma[sp]
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
        lim = s.control.dg_stddev_bound * c.sigma[t.species[i]]
        oldrn[i] > lim && (oldrn[i] = lim)
        oldrn[i] < -lim && (oldrn[i] = -lim)
    end

    # Small-tree height-growth calibration: HCOR_init (regent.f:411-516). For each
    # LHTCAL species (all by default) regress measured-vs-predicted small-tree (dbh<5)
    # height growth: HCOR_init = ln( Σ(HTG·P) / Σ(EDH·P) ) when ≥ NCALHT(5) observations,
    # where EDH = HTCALC-predicted increment (mode 0 age → mode 9 HTGR, ≥0.1) and the
    # measured HTG is scaled by SCALE3 = REGYR/FINTH (regent.f:406/462; REGYR=5, FINTH = the
    # GROWTH keyword's height measurement period, default 5 → SCALE3 = 1; RHCON = 1, no HCOR2). This
    # initial value seeds DIFH = HCOR_init − WCI; the per-cycle attenuation (diameter_growth!)
    # rides it like COR. Uses the CURRENT diameters (saved_dbh) for the dbh<5 filter, since
    # the regent calibration is independent of the large-tree diameter backdating.
    # SN-only: the small-tree height calibration uses the SN ht_curve (HTCALC) coefs + the REGENT
    # model. NE's height-dbh/growth model differs (htdbh + the NE htgf), so its small-tree height
    # calibration is a separate NE piece — skip this block for NE (htg_cor_init stays 0).
    s.variant isa Southern && let bc = (sd[:ht_curve_b1], sd[:ht_curve_b2], sd[:ht_curve_b3], sd[:ht_curve_b4], sd[:ht_curve_b5]),
        montane = !isempty(s.plot.eco_unit) && s.plot.eco_unit[1] == 'M'
        ncalht = 5
        scale3 = s.control.growth_finth > 0f0 ? 5f0 / s.control.growth_finth : 1f0  # REGYR/FINTH
        @inbounds for sp in 1:MAXSP
            i1 = isct[sp, 1]; i1 == 0 && continue
            i2 = isct[sp, 2]
            si = s.plot.sp_site_index[sp]
            snp = 0f0; snx = 0f0; sny = 0f0; nh = 0
            for k in i1:i2
                i = ind1[k]
                saved_dbh[i] >= 5f0 && continue            # large trees excluded
                hh = t.height[i] - t.ht_growth[i]          # start-of-period height (IHTG<2)
                hh < 0.01f0 && continue
                t.ht_growth[i] < 0.001f0 && continue       # no measured height growth
                aget = htcalc_age(bc, sp, si, hh, montane)
                htgr = htcalc_incr(bc, sp, si, aget, montane)
                htgr < 0.1f0 && (htgr = 0.1f0)
                edh = htgr                                  # ·RHCON(=1); EDH≥0.1 already
                p = t.tpa[i]
                term = t.ht_growth[i] * scale3            # TERM = HTG·SCALE3 (regent.f:462)
                snp += p; snx += edh * p; sny += term * p; nh += 1
            end
            nh < ncalht && continue
            cornew = sny / snx                              # (ΣHTG·P)/(ΣEDH·P)
            cornew <= 0f0 && (cornew = 1f-4)
            (cornew < 0.0821f0 || cornew > 12.1825f0) && (cornew = 1f0)
            c.htg_cor_init[sp] = log(cornew)
        end
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
const DG_FU =  1.271f0        # tripling upper-record factor (dgdriv.f FU)
const DG_FL = -1.549f0        # tripling lower-record factor (dgdriv.f FL)

function diameter_growth!(s::StandState, ::AbstractVariant; sfint::Float32 = 5f0,
                          tripling::Bool = true)
    t, c = s.trees, s.calib
    sd = s.coef.species
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    yr = htg_period(s.variant)   # DG model native period (gradd.f FINT/YR scale): 5 SN, 10 NE
    # DGBND DBH-range bounds are SN-only (NE's DGBND is just the SIZCAP cap, ne/dgbnd.f); `nothing`
    # ⇒ the per-tree bound skips the dlo/dhi adjustment and applies only the size cap.
    dlo_v = haskey(sd, :dg_bound_dbh_lo) ? sd[:dg_bound_dbh_lo] : nothing
    dhi_v = haskey(sd, :dg_bound_dbh_hi) ? sd[:dg_bound_dbh_hi] : nothing
    isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
    oldrn = t.old_random
    nlive = t.n
    # tripling scratch: per-tree upper/lower DG + their serial-correlation residual
    do_trip = tripling && 3 * nlive + t.ndead <= length(t.dbh)
    dgU = do_trip ? Vector{Float32}(undef, nlive) : Float32[]
    dgL = do_trip ? Vector{Float32}(undef, nlive) : Float32[]
    rnU = do_trip ? Vector{Float32}(undef, nlive) : Float32[]
    rnL = do_trip ? Vector{Float32}(undef, nlive) : Float32[]

    # Attenuate COR toward the calibration goal before predicting (dgdriv.f:76-79).
    # The attenuation clock is the cumulative elapsed time SINCE the inventory (FVS
    # SFINT = IY(icyc)−IY(1)), so the first projection cycle uses the FULL COR
    # (cormlt=1 at 0 elapsed years) and it decays thereafter. Use the actual elapsed
    # years from the IY schedule (= current_cycle_year − inventory), NOT `sfint·cycle`
    # — the two are equal only for UNIFORM cycles; a TIMEINT/CYCLEAT non-uniform
    # schedule (e.g. a 10-yr cycle 2) needs the true cumulative time (5, not 10).
    # NB (verified the hard way, -1823 tests): a debug-FVS dgdriv COR dump shows COR
    # ONE CYCLE AHEAD (1.0221 at cyc1) because that WRITE fires AFTER dgdriv updates COR
    # for the NEXT cycle, while dgf already baked the CURRENT (pre-update) COR into WK2.
    # jl's dg_cor at cycle N = the value FVS USES for cycle N's DG (the pre-update one) —
    # the START clock here is correct; do NOT "fix" it to elapsed+sfint.
    elapsed = Float32(current_cycle_year(s) - Int(s.control.cycle_year[1]))
    cormlt = exp(-0.02773f0 * elapsed)
    # The REGENT small-tree height calibration HCOR rides the SAME WCI attenuation as the
    # diameter COR (dgdriv.f:188-194) but on the elapsed-at-END-of-period clock (cycle+1):
    # HCOR = WCI + cormlt_h·DIFH, DIFH = HCOR_init − WCI (set at ICYC=1). This runs for LDGCAL
    # species, and LDGCAL defaults TRUE for ALL species (grinit.f:102; only the unported
    # NOCALIB keyword turns it off), so it always runs here. For a species with NO diameter
    # calibration WCI=0, so the height calibration DECAYS as cormlt_h·HCOR_init (it is NOT held
    # constant — that was the bug); for one with no height calibration HCOR_init=0, reducing to
    # the WCI·(1−cormlt_h) progression. HCOR is SEPARATE from the large-tree HTGF term HTCON
    # (`htg_cor`, from the HCOR2 keyword, 0 for snt01). HCOR_init is computed by the regent
    # regression in `calibrate_diameter_growth!`.
    cormlt_h = exp(-0.02773f0 * (elapsed + sfint))   # elapsed at END of this period (cumulative)
    @inbounds for sp in 1:MAXSP
        c.dg_cor[sp] = c.dg_cor_goal[sp] + cormlt * c.dg_cor_goal[sp]
        c.htg_cor_small[sp] = c.dg_cor_goal[sp] + cormlt_h * (c.htg_cor_init[sp] - c.dg_cor_goal[sp])
    end

    species_sort!(s)
    dgf!(s, s.variant)
    wk2 = view(s.scratch.wk, 2, :)

    # per-cycle ARMA multipliers: AUTCOR(new, old) where `new` = THIS cycle's period and
    # `old` = the PREVIOUS cycle's period (dgdriv.f). For uniform 5-yr cycles both are 5
    # (unchanged); a non-uniform TIMEINT/CYCLEAT schedule (e.g. a 10-yr cycle following a
    # 5-yr one) needs new=10, old=5 — using the base YR for both under-grows the long cycle.
    cyc = Int(s.control.cycle)
    newp = max(1, cycle_period_at(s.control, cyc))
    oldp = cyc == 0 ? 5 : max(1, cycle_period_at(s.control, cyc - 1))   # first cycle: 5-yr measurement base (dgdriv AUTCOR)
    covmlt, vmlt = autcor(newp, oldp, _stand_bjrho(s))
    pvmlt = c.vmlt > 0f0 ? c.vmlt : vmlt
    corr = covmlt / sqrt(vmlt * pvmlt)
    c.vmlt = vmlt   # FVS dgdriv.f:116 PVMLT=VMLT carry (uniform 5-yr unaffected)
    # BAIMULT (MULTS kind 1): per-species diameter-growth multiplier scaling DDS
    # (dgdriv.f XDGROW=ln(XDMULT) added to ln(DDS) ⇒ DDS·XDMULT).
    cur_year = current_cycle_year(s)   # IY schedule (TIMEINT/CYCLEAT-aware)

    @inbounds for sp in 1:MAXSP
        i1 = isct[sp, 1]; i1 == 0 && continue
        i2 = isct[sp, 2]
        xbai = active_multiplier(s.control, :bai, sp, cur_year)
        vardg = c.vardg[sp]
        evarp1 = (sqrt(1f0 + 4f0 * vardg * pvmlt) + 1f0) / 2f0
        sig1   = sqrt(log(max(evarp1, 1f0 + eps(Float32))))
        evarp2 = (sqrt(1f0 + 4f0 * vardg * vmlt) + 1f0) / 2f0
        ssigma = sqrt(log(max(evarp2, 1f0 + eps(Float32))))
        rho = (sig1 > 0f0 && ssigma > 0f0) ?
              log(1f0 + corr * sqrt((evarp1 - 1f0) * (evarp2 - 1f0))) / (sig1 * ssigma) : 0f0
        rhocp = sqrt(max(1f0 - rho * rho, 0f0))
        frmbase = DG_FM * ssigma * rhocp
        fru = DG_FU * ssigma * rhocp           # upper-triple FRM factor (dgdriv.f:91)
        frl = DG_FL * ssigma * rhocp           # lower-triple FRM factor (dgdriv.f:89)
        for k in i1:i2
            i = ind1[k]
            bark = bark_ratio(bark_a, bark_b, sp, t.dbh[i])
            d_ib = t.dbh[i] * bark
            # FVS bounds the 5-yr DG (DGBND, dgdriv.f:255-269) THEN scales to the cycle length
            # (gradd.f:79-90, DDS·(FINT/YR)) WITHOUT re-bounding. So DDS here is the 5-yr basis (BAIMULT
            # only); `bsc` bounds the 5-yr DG and then scales by sfint/5. FINT=5 ⇒ identity (no scale).
            dds5 = exp(wk2[i]) * xbai                       # YR-yr DDS (BAIMULT: DDS·XDMULT); YR=5 SN / 10 NE
            bsc(dg5) = _bound_scale(dlo_v, dhi_v, sp, t.dbh[i], d_ib, dg5, sfint, s.control.sp_size_cap, yr)
            if do_trip
                rnpar = oldrn[i]                            # original residual (dgdriv.f:116)
                frmt = frmbase + corr * rnpar; oldrn[i] = frmt
                t.diam_growth[i] = bsc(sqrt(d_ib * d_ib + dds5 * exp(frmt)) - d_ib)
                ru = fru + corr * rnpar; rnU[i] = ru
                dgU[i] = bsc(sqrt(d_ib * d_ib + dds5 * exp(ru)) - d_ib)
                rl = frl + corr * rnpar; rnL[i] = rl
                dgL[i] = bsc(sqrt(d_ib * d_ib + dds5 * exp(rl)) - d_ib)
            else
                if tripling
                    frmt = frmbase + corr * oldrn[i]       # deterministic (dgdriv.f:117)
                    oldrn[i] = frmt
                    frm = exp(frmt)
                else
                    frm = dgscor!(s.rng, oldrn, i, ssigma, rho, rhocp, wk2[i];
                                  dgsd = s.control.dg_stddev_bound)
                end
                t.diam_growth[i] = bsc(sqrt(d_ib * d_ib + dds5 * frm) - d_ib)
            end
        end
    end

    # Record tripling is performed LATER (after height growth + mortality), because
    # FVS runs MORTS before TRIPLE (grincr.f) — VARMRT must distribute mortality over
    # the original ITRN records, not the tripled set. Return the per-tree upper/lower
    # DG + serial-correlation residual so `triple_records!` can build the records.
    # htgU/htgL/is_small carry the small-tree (REGENT) tripled-record overrides that
    # `small_tree_growth!` fills after height growth (large-tree records keep is_small
    # = false and inherit the central HTG via copy_tree!).
    htgU = do_trip ? zeros(Float32, nlive) : Float32[]
    htgL = do_trip ? zeros(Float32, nlive) : Float32[]
    is_small = do_trip ? falses(nlive) : BitVector()
    return do_trip ? (nlive = nlive, dgU = dgU, dgL = dgL, rnU = rnU, rnL = rnL,
                      htgU = htgU, htgL = htgL, is_small = is_small) : nothing
end

"""
    triple_records!(state, stash)

TRIPLE (triple.f): split each of the `stash.nlive` original live records into 3
weighted records — central 0.60, upper 0.25 (FU growth), lower 0.15 (FL growth) —
using the DGs stashed by `diameter_growth!`. Run AFTER mortality so the split is of
the post-mortality TPA (matching FVS's MORTS-before-TRIPLE order). Dead records are
pushed to the end. No-op when `stash === nothing`.
"""
function triple_records!(s::StandState, stash)
    stash === nothing && return s
    t = s.trees; nlive = stash.nlive
    dgU = stash.dgU; dgL = stash.dgL; rnU = stash.rnU; rnL = stash.rnL
    htgU = stash.htgU; htgL = stash.htgL; is_small = stash.is_small
    @inbounds for k in t.ndead:-1:1
        copy_tree!(t, 3 * nlive + k, nlive + k)
    end
    @inbounds for i in 1:nlive
        # FVS TRIPLE (triple.f:18) appends BOTH new records per parent contiguously:
        # itfn = ITRN+2i-1 (weight .25), itfn+1 = ITRN+2i (weight .15). This interleaved
        # physical layout (not all-uppers-then-all-lowers) is what TREDEL's swap-from-end
        # walks after a thin, so it must match the oracle's append order exactly.
        u = nlive + 2 * i - 1; l = nlive + 2 * i
        copy_tree!(t, u, i); copy_tree!(t, l, i)
        t.tpa[u] = t.tpa[i] * 0.25f0; t.diam_growth[u] = dgU[i]; t.old_random[u] = rnU[i]
        t.tpa[l] = t.tpa[i] * 0.15f0; t.diam_growth[l] = dgL[i]; t.old_random[l] = rnL[i]
        # the record's period mortality (MortPA) splits with the surviving TPA (0.60/0.25/0.15)
        t.mort_pa[u] = t.mort_pa[i] * 0.25f0; t.mort_pa[l] = t.mort_pa[i] * 0.15f0
        # small-tree records carry per-record height increments (REGENT random effect)
        if is_small[i]
            t.ht_growth[u] = htgU[i]; t.ht_growth[l] = htgL[i]
        end
        t.tpa[i] *= 0.60f0; t.mort_pa[i] *= 0.60f0
        # lineage keys: upper=3K, central=3K+1, lower=3K+2 → species-sort then visits
        # records in the oracle's LNKCHN order (upper,central,lower depth-first), so the
        # untripled DGSCOR consumes the BACHLO stream in the same sequence (RNG-exact).
        kk = t.sort_key[i]
        t.sort_key[u] = 3 * kk; t.sort_key[i] = 3 * kk + 1; t.sort_key[l] = 3 * kk + 2
    end
    t.n = 3 * nlive
    return s
end
