# =============================================================================
# mortality.jl — Southern periodic mortality (MORTS)
#
# Ported from: sn/morts.f.
#
# Each cycle a fraction of each tree's trees-per-acre dies. SN combines:
#   * background mortality — a Hamilton logistic of DBH: ri = 1/(1+exp(b0+b1·DBH));
#   * density (Pretzsch) — when stand density `t` exceeds PMSDIL·(max density at the
#     stand QMD), a self-thinning-line rate `rn` applies instead.
# Per tree, deaths = PROB·(1−(1−rip)^FINT). Validated vs Oracle A on snt01 cycle-1
# (SDIMAX 348.4, t 589.65, dia0 4.70, d10 5.52, tn10 559, rn 0.0109, ~29 TPA).
# Tables (b0=PMSC, b1=PMD) copied verbatim.
# =============================================================================

const SDI_EXP = -1.605f0
const PRETZSCH_SDIK = 0.02483133f0
# Background-mortality coefficients (PMSC/PMD) and SDIMAX defaults live in
# data/southern/species_coefficients.csv (mort_bkgd_intercept/mort_bkgd_dbh).

"BA-weighted stand maximum SDI (SDICAL, pre-CLMAXDEN)."
function stand_sdimax(s::StandState)
    t = s.trees; p = s.plot
    num = 0f0; totba = 0f0
    @inbounds for i in 1:t.n
        tb = 0.0054542f0 * t.dbh[i]^2 * t.tpa[i]
        num   += p.sp_sdi_def[t.species[i]] * tb
        totba += tb
    end
    return totba <= 0f0 ? 1f0 : num / totba
end

# Pretzsch self-thinning target density tn10 (morts.f:200-343).
function _pretzsch_tn10(t, dia0, d10, const_v, pmsdil, pmsdiu)
    tmd0  = min(const_v * dia0^SDI_EXP, 35000f0)
    t85d0 = tmd0 * pmsdiu;  t55d0 = pmsdil * tmd0
    tmd10  = min(const_v * d10^SDI_EXP, 35000f0)
    t85d10 = tmd10 * pmsdiu; t55d10 = pmsdil * tmd10

    t > t85d0 && return min(t85d10, t)

    # solve the self-thinning line at a trial density → (slope, intercept)
    line(tem) = begin
        d55m = (log(tem) - log(pmsdil * const_v)) / SDI_EXP
        t55m = log(tem)
        d85m = d55m * 1.25f0
        local slp::Float32
        while true
            d85m = clamp(d85m, 0.125f0, 5f0)
            t85m = log(const_v * exp(d85m)^SDI_EXP * pmsdiu)
            slp = (t85m - t55m) / (d85m - d55m)
            (slp > -0.5f0 && d85m < 5f0) ? (d85m += 0.1f0) : break
        end
        (slp, t55m - slp * d55m)
    end

    if t > t55d0                                   # ipath 1: converge treeit
        abs(t85d0 - t) <= 5f0 && return min(t85d10, t)
        treeit = t + 0.1f0 * t; slp = 0f0; cept = 0f0
        for _ in 1:100
            slp, cept = line(treeit)
            diff = t - exp(cept + slp * log(dia0))
            (-5f0 <= diff <= 5f0) && break
            treeit += 0.5f0 * diff
        end
        return min(exp(cept + slp * log(d10)), t85d10)
    else                                           # t ≤ t55d0
        t <= t55d10 && return t
        slp, cept = line(t)                        # ipath 2
        return min(exp(cept + slp * log(d10)), t85d10)
    end
end

"""
    mortality!(state, ::Southern; fint=5f0)

Compute and apply periodic mortality, reducing `trees.tpa`. Combines background
(Hamilton) and density (Pretzsch self-thinning) rates. Runs after diameter growth
(uses `trees.diam_growth` for the projected end-of-cycle QMD).
"""
function mortality!(s::StandState, ::Southern; fint::Float32 = 5f0)
    p, t = s.plot, s.trees
    pmsdil = p.pct_sdimax_mort_lo > 0f0 ? p.pct_sdimax_mort_lo : 0.55f0
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0
    dbhstage = s.control.dbh_stage

    # SN uses the Zeide/Reineke mean diameter for density mortality (LZEIDE):
    #   dia0 = (Σ p·d^1.605 / Σ p)^(1/1.605), and d10 the same with grown diameters.
    zeide = s.control.zeide_sdi
    dthresh = zeide ? s.control.dbh_zeide : dbhstage
    bark_a = s.coef.species[:bark_intercept]; bark_b = s.coef.species[:bark_slope]
    mort_b0 = s.coef.species[:mort_bkgd_intercept]; mort_b1 = s.coef.species[:mort_bkgd_dbh]
    tt = 0f0; sdq0 = 0f0; sd2sq = 0f0; sumdr0 = 0f0; sumdr10 = 0f0
    @inbounds for i in 1:t.n
        d = t.dbh[i]
        d < dthresh && continue
        pr = t.tpa[i]
        bark = bark_ratio(bark_a, bark_b, t.species[i], d)
        g = (t.diam_growth[i] / bark) * (fint / 5f0)
        sd2sq += pr * (d * d + 2f0 * d * g + g * g)
        sdq0  += pr * d * d
        sumdr0  += pr * d^1.605f0
        sumdr10 += pr * (d + g)^1.605f0
        tt += pr
    end
    tt < 1f0 && return s
    tt > 35000f0 && (tt = 35000f0)
    dia0 = zeide ? (sumdr0  / tt)^(1f0 / 1.605f0) : sqrt(sdq0  / tt)
    d10  = zeide ? (sumdr10 / tt)^(1f0 / 1.605f0) : sqrt(sd2sq / tt)
    dia0 < 0.3f0 && (d10 = 0.3f0 + d10 - dia0; dia0 = 0.3f0)

    sdimax = stand_sdimax(s)
    density_on = false; rn = 0f0
    if sdimax >= 5f0
        const_v = sdimax / PRETZSCH_SDIK
        tn10 = _pretzsch_tn10(tt, dia0, d10, const_v, pmsdil, pmsdiu)
        tn10 = clamp(tn10, 0f0, tt); tn10 < 0.1f0 && (tn10 = 0f0)
        rn = 1f0 - (1f0 - (tt - tn10) / tt)^(1f0 / fint)
        tem_v2 = min(const_v * d10^SDI_EXP, 35000f0) * pmsdil
        density_on = !(tt <= tem_v2 || rn <= 0f0)
    end

    @inbounds for i in 1:t.n
        pr = t.tpa[i]; pr <= 0f0 && continue
        sp = t.species[i]; d = t.dbh[i]
        ri = 1f0 / (1f0 + exp(mort_b0[sp] + mort_b1[sp] * d))
        rip = density_on ? rn : ri
        rip > 1f0 && (rip = 1f0)
        deaths = min(pr * (1f0 - (1f0 - rip)^fint), pr)
        t.tpa[i] = pr - deaths
    end
    return s
end
