# =============================================================================
# crown_ratio.jl — Southern per-cycle crown-ratio update (CROWN)
#
# Ported from: sn/crown.f (CROWN). Runs once per cycle AFTER growth (gradd.f:172,
# i.e. after DG/HTG applied, the tree list re-sorted on dbh, and density recomputed)
# to update each live tree's crown ratio `trees.crown_pct` (ICR). The new crown
# ratio is a Weibull draw at the tree's diameter percentile, recentred on the
# species mean-crown-ratio curve `acrnew(relsdi)`, then limited to a ±1%/yr change
# from the prior crown and capped by the crown-length geometry. Coefficients in
# data/southern/crown_ratio_coeffs.csv (MCREQN mean-CR eqn + WEIBUL params).
#
# Crown ratio is bit-exact at cycle 1 only when this runs AFTER cycle 1's growth
# (so cycle-1 DGF still sees the input crown); without it, the frozen crown made
# every cycle ≥2 diverge (crown → log(ICR) in DGF wk2 → DG → mortality → volume).
# =============================================================================

"""
    init_crown_ratios!(s)

CRATET: estimate the INITIAL crown ratio for inventory trees that have no input crown
(`crown_pct == 0`), using the CROWN model on the inventory stand (FVS does this in
INITRE/CRATET so the first cycle's DGF and the FFE inventory crown are based on a real
crown, not 0). Trees WITH an input crown are left untouched. No-op if every tree already
has a crown. Idempotent; safe to call once before the first grow.
"""
function init_crown_ratios!(s::StandState)
    t = s.trees
    t.n == 0 && return s
    n = t.n
    any(@views t.crown_pct[1:n] .== 0) || return s   # all crowns already set
    compute_density!(s)
    # FVS's CRATET dubs missing crowns with CROWN using DENSE's CCF computed on the BACKDATED dbh
    # (DENSE/LBKDEN: past dbh = sqrt(d²·r), r from the measured DG; unmeasured trees use the stand-average
    # ratio). Replicate that backdated CCF (RELDEN) so the init crown matches — the percentile SCALE depends
    # on it, and using the current (denser) dbh gives a ~5% low crown.
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    bagr = 0f0; nb = 0f0
    @inbounds for i in 1:n
        g = t.diam_growth[i]; g <= 0f0 && continue
        d = t.dbh[i]; gadj = g / bark_ratio(bark_a, bark_b, t.species[i], d)
        gadj > d && continue
        bagr += 1f0 - (2f0 * d * gadj - gadj * gadj) / (d * d); nb += 1f0
    end
    nb > 0f0 && (bagr /= nb)
    saved_dbh = copy(@view t.dbh[1:n])
    @inbounds for i in 1:n
        d = t.dbh[i]; g = t.diam_growth[i]; r = bagr
        if g > 0f0
            gadj = min(g / bark_ratio(bark_a, bark_b, t.species[i], d), d)
            rr = 1f0 - (2f0 * d * gadj - gadj * gadj) / (d * d)
            rr > 0f0 && (r = rr)
        end
        r > 0f0 && (t.dbh[i] = sqrt(d * d * r))        # backdated dbh
    end
    compute_density!(s)
    bd_relden = stand_ccf(s)                            # CCF on the backdated stand
    @inbounds for i in 1:n; t.dbh[i] = saved_dbh[i]; end   # restore current dbh
    compute_density!(s)                                 # rank/SDI use current dbh (as FVS CROWN)
    saved = copy(@view t.crown_pct[1:n])
    crown_ratio_update!(s; fint = 5f0, relden_override = bd_relden)
    @inbounds for i in 1:n
        saved[i] != 0 && (t.crown_pct[i] = saved[i])   # restore input crowns; keep only the estimated 0s
    end
    return s
end

"""
    crown_ratio_update!(state; fint=5f0)

CROWN: update `trees.crown_pct` (ICR, %) for every live record from the post-growth
stand. No-op for an empty stand. Call once per cycle after growth + density.
"""
function crown_ratio_update!(s::StandState; fint::Float32 = 5f0, crown_sdi::Float32 = -1f0,
                             relden_override::Float32 = -1f0)
    t = s.trees; sd = s.coef.species; n = t.n
    n == 0 && return s
    # CRNMULT: cycle year for the persistent crown-ratio-change multiplier lookup.
    cur_year = current_cycle_year(s)   # IY schedule (TIMEINT/CYCLEAT-aware)
    sdiac = crown_sdi >= 0f0 ? crown_sdi : stand_sdi_reineke(s)  # pre-growth Reineke SDIBC (grincr.f:241)
    relden = relden_override >= 0f0 ? relden_override : stand_ccf(s)  # RELDEN — crown competition factor
                                         # (override = the DENSE-backdated CCF, used by CRATET init crown)
    sdidef = s.plot.sp_sdi_def
    # Ascending diameter rank: isort[i] = 1 (smallest) … n (largest); x = rank/n.
    ord = sortperm(view(t.dbh, 1:n))
    isort = Vector{Int32}(undef, n)
    @inbounds for r in 1:n; isort[ord[r]] = Int32(r); end
    scale = clamp(1f0 - 0.00167f0 * (relden - 100f0), 0.30f0, 1f0)

    eqn = sd[:mcr_eqn]; ma = sd[:mcr_a]; mc = sd[:mcr_c]; mb = sd[:mcr_b]
    mb2 = sd[:mcr_b2]; mb3 = sd[:mcr_b3]
    wa = sd[:wb_a]; wb0 = sd[:wb_b0]; wb1 = sd[:wb_b1]; wc = sd[:wb_c]
    # Per-species Weibull params (acrnew depends only on species via relsdi).
    Aw = zeros(Float32, MAXSP); Bw = zeros(Float32, MAXSP); Cw = zeros(Float32, MAXSP)
    seen = falses(MAXSP)
    @inbounds for i in 1:n
        sp = t.species[i]
        if !seen[sp]
            relsdi = sdidef[sp] > 0f0 ? sdiac / sdidef[sp] * 10f0 : 6f0
            relsdi = clamp(relsdi, 1f0, 12f0)
            ie = Int(eqn[sp])
            acrnew = ie == 1 ? exp(ma[sp] + mb[sp] * log(relsdi) + mc[sp] * relsdi) :
                     ie == 2 ? exp(ma[sp] + mb[sp] * log(relsdi)) :
                     ie == 3 ? ma[sp] + mc[sp] * relsdi :
                     ie == 4 ? ma[sp] + mb2[sp] * log10(relsdi) :
                     ie == 5 ? relsdi / (ma[sp] * relsdi + mb3[sp]) : 0f0
            bb = wb0[sp] + wb1[sp] * acrnew; bb < 3f0 && (bb = 3f0)
            cc = wc[sp]; cc < 2f0 && (cc = 2f0)
            Aw[sp] = wa[sp]; Bw[sp] = bb; Cw[sp] = cc; seen[sp] = true
        end
        icr_old = t.crown_pct[i]
        if icr_old < 0   # crown change already computed by the topkill/pest model
            t.crown_pct[i] = -icr_old   # (sn/crown.f:271): restore sign, bypass the recompute
            continue
        end
        d = t.dbh[i]
        x = d > 0f0 ? Float32(isort[i]) / Float32(n) * scale : 0.5f0 * scale
        x = clamp(x, 0.05f0, 0.95f0)
        crnew = Aw[sp] + Bw[sp] * ((-log(1f0 - x))^(1f0 / Cw[sp]))
        # Limit change to ±1%/yr of the prior crown (crown.f:442-459).
        if icr_old != 0
            chg = crnew - Float32(icr_old)
            pdifpy = chg / Float32(icr_old) / fint
            pdifpy >  0.01f0 && (chg = Float32(icr_old) *  0.01f0 * fint)
            pdifpy < -0.01f0 && (chg = Float32(icr_old) * -0.01f0 * fint)
            # CRNMULT (crown.f:319): scale the crown-ratio change for trees in the keyword's
            # DBH window (1.0 = no CRNMULT keyword, the common case).
            crnew = Float32(icr_old) + chg * active_crn_mult(s.control, sp, cur_year, d)
        end
        icri = trunc(Int32, crnew + 0.5f0)
        # Crown-length cap: the crown can't exceed (old length + HTG) over new height.
        if icr_old != 0
            crln = t.height[i] * Float32(icr_old) / 100f0
            crmax = (crln + t.ht_growth[i]) / (t.height[i] + t.ht_growth[i]) * 100f0
            (Float32(icri) > crmax || icri < 10) && (icri = trunc(Int32, crmax + 0.5f0))
        end
        icri > 95 && (icri = Int32(95))
        icri < 10 && (icri = Int32(10))
        icri < 1  && (icri = Int32(1))
        t.crown_pct[i] = icri
    end
    return s
end
