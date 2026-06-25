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

# Mortality self-thinning trajectory growth (morts.f:225/583/690). FVS runs MORTS
# *before* GRADD's concave sqrt rescaling, so every SDI/size-cap diameter trajectory
# uses the 5-yr DG *linearly* extrapolated to FINT years: G = DG_5/BARK · FINT/5.
# We only retain the sqrt-scaled fint-year growth, so recover DG_5 from it first.
# Identity at fint=5 ⇒ snt01 and all 5-yr scenarios stay bit-exact.
@inline function _mort_traj_g(dg_ib_fint::Float32, dbh::Float32, bark::Float32, fint::Float32)
    fint == 5f0 && return dg_ib_fint / bark   # exact identity — avoid the sqrt roundtrip's 1-ULP noise
    dib  = dbh * bark
    ddsf = (dg_ib_fint + dib)^2 - dib * dib   # fint-year inside-bark DDS
    dds5 = ddsf * (5f0 / fint)                # back to 5-yr DDS
    dg5  = sqrt(dib * dib + dds5) - dib       # 5-yr inside-bark DG
    return (dg5 / bark) * (fint / 5f0)        # outside-bark, linear FINT extrapolation
end
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

"""
    sdi_max_check!(state)

SDICHK (sdichk.f): if the INITIAL stand stocking exceeds 5% above the upper density
limit (TPROB > (PMSDIU+0.05)·temmax), the species SDI maxima are too low for this
stand — reset them all to a value `tem2` fitted to the observed (TPROB, QMD) so the
self-thinning mortality doesn't run away. Called once at setup. No-op for stands at or
below the limit (e.g. snt01). Without it, over-dense single-species stands over-kill.
"""
function sdi_max_check!(s::StandState)
    t = s.trees; n = t.n
    n == 0 && return s
    sdimax = stand_sdimax(s)
    sdimax < 5f0 && return s
    p = s.plot
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0
    zeide = s.control.zeide_sdi
    dthresh = zeide ? s.control.dbh_zeide : s.control.dbh_stage
    tprob = 0f0; sumdr = 0f0; sumd2 = 0f0
    @inbounds for i in 1:n
        d = t.dbh[i]; d < dthresh && continue
        pr = t.tpa[i]
        sumdr += pr * d^1.605f0; sumd2 += pr * d * d; tprob += pr
    end
    tprob < 1f0 && return s
    dq0 = zeide ? (sumdr / tprob)^(1f0 / 1.605f0) : sqrt(sumd2 / tprob)
    dq0 < 0.3f0 && (dq0 = 0.3f0)
    const_v = sdimax / PRETZSCH_SDIK
    upmax = min(pmsdiu + 0.05f0, 1f0)
    temmax = const_v * dq0^SDI_EXP
    tprob <= upmax * temmax && return s        # not over the upper limit ⇒ keep SDIDEF
    const_v2 = exp(log(tprob + 1f0) + 1.605f0 * log(dq0)) / pmsdiu
    tem2 = const_v2 * PRETZSCH_SDIK
    @inbounds for i in 1:MAXSP; p.sp_sdi_def[i] = tem2; end
    return s
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

"""
    _varmrt!(killed, efftr, temwk2, shade_adj, t, n, tokill) -> sumkil

VARMRT (varmrt.f): distribute `tokill` TPA of mortality across the `n` live records
by a geometric progression weighted toward suppressed trees. Per-tree efficiency
`efftr = peff(PCT)·shade_adj·0.1`, where `peff = 0.84525 − 0.01074·PCT +
2e-7·PCT³` (low percentile ⇒ high mortality). Fills `killed[i]`; returns the total.
"""
function _varmrt!(killed::AbstractVector{Float32}, efftr::AbstractVector{Float32},
                  temwk2::AbstractVector{Float32}, shade_adj::AbstractVector{Float32},
                  t::TreeList, n::Int, tokill::Float32)
    fill!(view(killed, 1:n), 0f0)
    tokill <= 0f0 && return 0f0
    pct = t.crown_ratio; tpa = t.tpa; sp = t.species
    pass1 = 0f0
    @inbounds for i in 1:n
        pe = clamp(0.84525f0 - 0.01074f0 * pct[i] + 0.0000002f0 * pct[i]^3f0, 0.01f0, 1f0)
        efftr[i] = pe * shade_adj[sp[i]] * 0.1f0
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
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    mort_b0 = s.coef.species[:mort_bkgd_intercept]; mort_b1 = s.coef.species[:mort_bkgd_dbh]
    shade_adj = s.coef.species[:varmrt_shade_adj]   # VARMRT per-species shade-tolerance scalar (CSV)
    tt = 0f0; sdq0 = 0f0; sd2sq = 0f0; sumdr0 = 0f0; sumdr10 = 0f0
    @inbounds for i in 1:t.n
        d = t.dbh[i]
        d < dthresh && continue
        pr = t.tpa[i]
        bark = bark_ratio(bark_a, bark_b, t.species[i], d)
        g = _mort_traj_g(t.diam_growth[i], d, bark, fint)   # morts.f:225 (linear FINT extrap)
        sd2sq += pr * (d * d + 2f0 * d * g + g * g)
        sdq0  += pr * d * d
        sumdr0  += pr * d^1.605f0
        sumdr10 += pr * (d + g)^1.605f0
        tt += pr
    end
    tt < 1f0 && return s
    # Reset the persisted self-thinning line when the stand TPA changed materially
    # since last cycle — i.e. after a thin or ingrowth (morts.f:160, |t−TPAMRT|>1).
    # For a closed stand t≈TPAMRT so the line persists (snt01 unaffected); after a
    # from-below thin it is recomputed for the new (large-tree) QMD.
    if s.density.tpa_mort > 0f0 && abs(tt - s.density.tpa_mort) > 1f0
        s.density.mort_slope = 0f0; s.density.mort_intercept = 0f0
    end
    tt > 35000f0 && (tt = 35000f0)
    dia0 = zeide ? (sumdr0  / tt)^(1f0 / 1.605f0) : sqrt(sdq0  / tt)
    d10  = zeide ? (sumdr10 / tt)^(1f0 / 1.605f0) : sqrt(sd2sq / tt)
    dia0 < 0.3f0 && (d10 = 0.3f0 + d10 - dia0; dia0 = 0.3f0)

    sdimax = stand_sdimax(s)
    n = t.n
    # preallocated VARMRT work buffers (sliced to the live count; no per-cycle allocation in the hot path)
    killed = @view s.scratch.mort_killed[1:n]; fill!(killed, 0f0)
    efftr  = @view s.scratch.mort_efftr[1:n]
    temwk2 = @view s.scratch.mort_temwk2[1:n]

    # Background (Hamilton) mortality total — used when density mortality is off; it
    # depends only on the start-of-cycle TPA, so it is computed once.
    # MORTMULT (MULTS kind 4): per-species multiplier on the BACKGROUND rate only
    # (morts.f:520-525: X=XMORT in [D1,D2], X=1 when the density rate is in effect).
    cur_year = current_cycle_year(s)   # IY schedule (TIMEINT/CYCLEAT-aware)
    bg_tokill = 0f0
    @inbounds for i in 1:n
        pr = t.tpa[i]; pr <= 0f0 && continue
        sp = t.species[i]
        ri = 1f0 / (1f0 + exp(mort_b0[sp] + mort_b1[sp] * t.dbh[i]))
        ri > 1f0 && (ri = 1f0)
        xmort = active_mort_mult(s.control, sp, cur_year, t.dbh[i])  # 1 outside the DBH window
        bg_tokill += min(pr * (1f0 - (1f0 - ri)^fint) * xmort, pr)
    end

    if sdimax < 5f0
        _varmrt!(killed, efftr, temwk2, shade_adj, t, n, bg_tokill)
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
            _varmrt!(killed, efftr, temwk2, shade_adj, t, n, tokill)
            density_on || break              # background ⇒ no d10 dependence, one pass
            # recompute the post-mortality QMD (d10n) from the surviving TPA
            ttn = 0f0; sdr = 0f0
            for i in 1:n
                d = t.dbh[i]; d < dthresh && continue
                pr = t.tpa[i] - killed[i]; pr <= 0f0 && continue
                bark = bark_ratio(bark_a, bark_b, t.species[i], d)
                # morts.f:583 — FVS uses the linear FINT-extrapolated 5-yr G here too,
                # but in this *iterative* d10n recompute the linear form amplifies the
                # DGBND-ordering gap (FVSjl bounds DG at the fint-year, FVS at 5-yr) and
                # over-thins fast cycles. Kept sqrt until DG is bounded at 5-yr in
                # diameter_growth.jl; identical at fint=5. See DIVERGENCES.md / line 223.
                g = t.diam_growth[i] / bark
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

    # Size-cap mortality (SIZCAP, morts.f:691-694): a tree whose grown DBH reaches the
    # per-species cap SIZCAP[1] (set by TREESZCP) gets a kill FLOOR of P·SIZCAP[2]·FINT/5
    # (≤ P), unless the no-mortality flag SIZCAP[3] truncates to 1. Applied after VARMRT,
    # before BAMAX. No-op unless TREESZCP set a cap (default SIZCAP[*,1]=999).
    let sc = s.control.sp_size_cap
        @inbounds for i in 1:n
            sp = t.species[i]
            (sc[sp, 1] >= 999f0 || trunc(Int, sc[sp, 3]) == 1) && continue
            d = t.dbh[i]
            # G is the OUTSIDE-bark, linearly FINT-extrapolated 5-yr increment
            # (sn/morts.f:690), same trajectory as the SDI self-thinning calc above.
            g = t.diam_growth[i] / bark_ratio(bark_a, bark_b, sp, d)
            if (d + g) >= sc[sp, 1]
                kc = min(t.tpa[i] * sc[sp, 2] * fint / 5f0, t.tpa[i])
                killed[i] < kc && (killed[i] = kc)
            end
        end
    end

    # BAMAX enforcement (morts.f:555-615): if the residual basal area exceeds the
    # stand maximum BAMAX = SDImax·0.5454154·PMSDIU, scale up every record's kill by
    # adjfac=(BA−BAMAX)/BAdead and re-test, iterating until BA ≤ BAMAX (max 100). This
    # caps the stand at its self-thinning BA — without it dense stands grow unbounded.
    bamax = sdimax * 0.5454154f0 * pmsdiu
    if bamax > 0f0
        @inbounds for _ in 1:100
            banew = 0f0; badead = 0f0
            for i in 1:n
                d = t.dbh[i]
                bark = bark_ratio(bark_a, bark_b, t.species[i], d)
                g = t.diam_growth[i] / bark   # morts.f:714
                de2 = 0.0054542f0 * (d + g)^2
                banew  += de2 * (t.tpa[i] - killed[i])
                badead += de2 * killed[i]
            end
            ((banew - bamax) > 1f0 && badead > 0f0) || break
            adjfac = (banew - bamax) / badead
            for i in 1:n
                killed[i] = min(t.tpa[i], killed[i] * (1f0 + adjfac))
            end
        end
    end

    # TPAMRT (morts.f:772): the surviving over-threshold TPA used for next cycle's self-thinning
    # line-reset test — locked HERE, from the BA-check survivors, BEFORE FIXMORT overrides the
    # kill (the forced FIXMORT mortality must NOT move the self-thinning line).
    surv = 0f0
    @inbounds for i in 1:n
        t.dbh[i] >= dthresh && (surv += t.tpa[i] - killed[i])
    end
    s.density.tpa_mort = surv

    # FIXMORT (morts.f:781): forced-mortality override, applied AFTER the BA-check and TPAMRT.
    apply_fixmort!(s, killed, n, fint)

    # FFE: trees killed by ordinary mortality become standing snags (FMSDIT, the FFE counterpart of
    # the fire-kill snags in fmburn!). The snags feed the Stand-Dead carbon pool and, as they fall
    # (update_snags!), the down-wood pool. No-op unless FFE is active.
    if s.fire !== nothing && s.fire.active
        yr = current_cycle_year(s); coef = s.coef
        v2t = coef_col(coef, :v2t); dkr = coef_col(coef, :dkr_cls)
        @inbounds for i in 1:n
            (killed[i] > 0f0 && t.dbh[i] > 0f0) || continue
            sp = Int(t.species[i])
            # FFE snag basis = death-time STEM-volume biomass (cuft·V2T), not whole-tree Jenkins
            # (fmdout.f SNVIS·V2T). V2T is the /2000-rescaled value, so multiply by P2T=1/2000.
            bolevol = t.cuft_vol[i] * v2t[sp] / 2000f0
            add_snag!(s.fire, sp, t.dbh[i], killed[i], yr; bolevol = bolevol)
            # the CROWN goes to CWD2B, scheduled to fall over TFALL years (FMSCRO)
            xv = crown_biomass(s, sp, t.dbh[i], t.height[i], Int(round(t.crown_pct[i])))
            fmscro!(s, sp, t.dbh[i], xv, killed[i], clamp(Int(dkr[sp]), 1, 4))
            # dead coarse roots accrue to the BIOROOT pool (FMSADD/FMCBIO RBIO, fmsadd.f:320)
            _, _, rbio = jenkins_biomass(coef, sp, t.dbh[i])
            s.fire.bioroot += rbio * killed[i]
        end
    end

    @inbounds for i in 1:n
        t.tpa[i] = max(0f0, t.tpa[i] - killed[i])
    end
    return s
end
