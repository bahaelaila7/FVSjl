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
    crown_ratio_update!(state; fint=5f0)

CROWN: update `trees.crown_pct` (ICR, %) for every live record from the post-growth
stand. No-op for an empty stand. Call once per cycle after growth + density.
"""
function crown_ratio_update!(s::StandState; fint::Float32 = 5f0)
    t = s.trees; sd = s.coef.species; n = t.n
    n == 0 && return s
    # CRNMULT: cycle year for the persistent crown-ratio-change multiplier lookup.
    cur_year = Int(s.control.cycle_year[1]) + Int(s.control.cycle) * round(Int, s.control.year)
    sdiac = stand_sdi(s)                 # SDIAC — post-growth stand SDI
    relden = stand_ccf(s)                # RELDEN — crown competition factor
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
