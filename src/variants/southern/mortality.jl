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

# Pretzsch self-thinning target density tn10 (morts.f:200-343). The self-thinning
# line (slope/intercept) is computed ONCE per stand and PERSISTED in `dens`
# (SLPMRT/CEPMRT, morts.f:317-322); subsequent cycles reuse it with the new d10.
function _pretzsch_tn10(dens::Density, t, dia0, d10, const_v, pmsdil, pmsdiu)
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

    local slp::Float32, cept::Float32
    if t > t55d0                                   # ipath 1: converge treeit
        abs(t85d0 - t) <= 5f0 && return min(t85d10, t)
        treeit = t + 0.1f0 * t; slp = 0f0; cept = 0f0
        for _ in 1:100
            slp, cept = line(treeit)
            diff = t - exp(cept + slp * log(dia0))
            (-5f0 <= diff <= 5f0) && break
            treeit += 0.5f0 * diff
        end
    else                                           # t ≤ t55d0
        t <= t55d10 && return t
        slp, cept = line(t)                        # ipath 2
    end
    # persist the line the first time it is solved; reuse it every later cycle
    if dens.mort_slope == 0f0
        dens.mort_slope = slp; dens.mort_intercept = cept
    end
    return min(exp(dens.mort_intercept + dens.mort_slope * log(d10)), t85d10)
end

# Per-species shade-tolerance scalar for the VARMRT mortality distribution (VARADJ,
# varmrt.f). Higher = more tolerant (survives suppression). (TODO: move to a CSV.)
const VARMRT_SHADE_ADJ = Float32[
    0.1, 0.7, 0.3, 0.7, 0.7, 0.7, 0.1, 0.7, 0.7, 0.7,
    0.7, 0.5, 0.7, 0.7, 0.5, 0.5, 0.1, 0.3, 0.3, 0.3,
    0.3, 0.1, 0.3, 0.7, 0.7, 0.1, 0.5, 0.7, 0.5, 0.3,
    0.1, 0.1, 0.1, 0.3, 0.7, 0.7, 0.3, 0.7, 0.3, 0.3,
    0.1, 0.7, 0.7, 0.7, 0.7, 0.3, 0.5, 0.3, 0.5, 0.3,
    0.7, 0.3, 0.7, 0.3, 0.7, 0.3, 0.3, 0.3, 0.5, 0.9,
    0.9, 0.7, 0.5, 0.9, 0.5, 0.7, 0.7, 0.3, 0.5, 0.7,
    0.7, 0.7, 0.7, 0.5, 0.5, 0.7, 0.7, 0.5, 0.5, 0.9,
    0.9, 0.7, 0.3, 0.5, 0.3, 0.5, 0.3, 0.5, 0.5, 0.5]

"""
    _varmrt!(killed, t, n, tokill) -> sumkil

VARMRT (varmrt.f): distribute `tokill` TPA of mortality across the `n` live records
by a geometric progression weighted toward suppressed trees. Per-tree efficiency
`efftr = peff(PCT)·shade_adj·0.1`, where `peff = 0.84525 − 0.01074·PCT +
2e-7·PCT³` (low percentile ⇒ high mortality). Fills `killed[i]`; returns the total.
"""
function _varmrt!(killed::Vector{Float32}, t::TreeList, n::Int, tokill::Float32)
    fill!(view(killed, 1:n), 0f0)
    tokill <= 0f0 && return 0f0
    pct = t.crown_ratio; tpa = t.tpa; sp = t.species
    efftr = Vector{Float32}(undef, n); temwk2 = zeros(Float32, n)
    pass1 = 0f0
    @inbounds for i in 1:n
        pe = clamp(0.84525f0 - 0.01074f0 * pct[i] + 0.0000002f0 * pct[i]^3f0, 0.01f0, 1f0)
        efftr[i] = pe * VARMRT_SHADE_ADJ[sp[i]] * 0.1f0
        pass1 += tpa[i] * efftr[i]
    end
    pass1 <= 0f0 && return 0f0
    npass = floor(Int, tokill / pass1) + 1
    sumkil = 0f0; temkil = tokill; short_v = 0f0; jpass = 0
    while true
        jpass += 1; jpass > 1 && (temkil = short_v)
        iswtch = 0; temsum = 0f0
        while true                                   # adjust npass into [0.8,1.2]
            temsum = 0f0
            @inbounds for i in 1:n
                tpalft = tpa[i] - killed[i]
                if tpalft > 0f0
                    temwk2[i] = -tpalft * ((1f0 - efftr[i])^npass - 1f0)
                    temsum += temwk2[i]
                end
            end
            minstp = npass > 50 ? 5 : (npass > 20 ? 2 : 1)
            adjust = temsum > 0f0 ? temkil / temsum : 1f0
            if adjust < 0.8f0 && iswtch != 2
                npass -= max(minstp, floor(Int, (temsum - temkil) / pass1)); iswtch = 1
                npass > 0 && continue
            elseif adjust > 1.2f0 && iswtch != 1
                npass += max(minstp, floor(Int, (temkil - temsum) / pass1)); iswtch = 2
                continue
            end
            break
        end
        short_v = 0f0
        adjust = temsum == 0f0 ? 1f0 : temkil / temsum
        @inbounds for i in 1:n
            tpalft = tpa[i] - killed[i]
            tpalft < 0.00001f0 && continue
            xkill = temwk2[i] * adjust
            if (tpa[i] - killed[i] - xkill) <= 0.00001f0
                xk = tpa[i] - killed[i]
                short_v += xkill - xk; pass1 -= efftr[i]
                killed[i] += xk; sumkil += xk
            else
                killed[i] += xkill; sumkil += xkill
            end
        end
        short_v <= 0f0 && break
        pass1 <= 0f0 && break
        npass = floor(Int, short_v / pass1) + 1
    end
    return sumkil
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
    n = t.n
    killed = zeros(Float32, n)

    # Background (Hamilton) mortality total — used when density mortality is off; it
    # depends only on the start-of-cycle TPA, so it is computed once.
    bg_tokill = 0f0
    @inbounds for i in 1:n
        pr = t.tpa[i]; pr <= 0f0 && continue
        sp = t.species[i]
        ri = 1f0 / (1f0 + exp(mort_b0[sp] + mort_b1[sp] * t.dbh[i]))
        ri > 1f0 && (ri = 1f0)
        bg_tokill += min(pr * (1f0 - (1f0 - ri)^fint), pr)
    end

    if sdimax < 5f0
        _varmrt!(killed, t, n, bg_tokill)
    else
        # MORTS QMD-convergence iteration (morts.f:184-481): solve tn10 for the
        # assumed end-of-cycle QMD d10, distribute the excess (t − tn10) by VARMRT,
        # recompute the post-mortality QMD d10n, and re-iterate with d10=d10n until it
        # converges (|d10−d10n|≤0.1) or the QMD would fall below dia0. The self-
        # thinning line is solved once (persisted in s.density) and reused each pass.
        const_v = sdimax / PRETZSCH_SDIK
        d10cur = d10
        @inbounds for _ in 1:10
            tn10 = _pretzsch_tn10(s.density, tt, dia0, d10cur, const_v, pmsdil, pmsdiu)
            tn10 = clamp(tn10, 0f0, tt); tn10 < 0.1f0 && (tn10 = 0f0)
            rn = 1f0 - (1f0 - (tt - tn10) / tt)^(1f0 / fint)
            tem_v2 = min(const_v * d10cur^SDI_EXP, 35000f0) * pmsdil
            density_on = !(tt <= tem_v2 || rn <= 0f0)
            tokill = density_on ? max(tt - tn10, 0f0) : bg_tokill
            _varmrt!(killed, t, n, tokill)
            density_on || break              # background ⇒ no d10 dependence, one pass
            # recompute the post-mortality QMD (d10n) from the surviving TPA
            ttn = 0f0; sdr = 0f0
            for i in 1:n
                d = t.dbh[i]; d < dthresh && continue
                pr = t.tpa[i] - killed[i]; pr <= 0f0 && continue
                bark = bark_ratio(bark_a, bark_b, t.species[i], d)
                g = (t.diam_growth[i] / bark) * (fint / 5f0)
                if zeide
                    sdr += pr * (d + g)^1.605f0
                else
                    sdr += pr * (d * d + 2f0 * d * g + g * g)
                end
                ttn += pr
            end
            d10n = ttn <= 0f0 ? 0f0 : (zeide ? (sdr / ttn)^(1f0 / 1.605f0) : sqrt(sdr / ttn))
            (abs(d10cur - d10n) <= 0.1f0 || d10n <= dia0) && break
            d10cur = d10n
        end
    end
    @inbounds for i in 1:n
        t.tpa[i] = max(0f0, t.tpa[i] - killed[i])
    end
    return s
end
