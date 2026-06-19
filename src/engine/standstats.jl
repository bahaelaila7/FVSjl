# =============================================================================
# standstats.jl — per-acre expansion (NOTRE) and stand summary statistics
#
# Ported from: base/notre.f (expansion) + base/stats.f / sumout (summary columns).
#
# NOTRE turns each record's tally (PROB) into trees-per-acre using the sampling
# design (variable-radius BAF plots vs fixed plots). The stand statistics
# (TPA, basal area, QMD, ...) are then simple weighted reductions over the trees —
# pure loops, autovectorizable, no allocation.
# =============================================================================

const BA_PER_TREE = 0.005454154f0     # ft² of basal area per inch² of DBH

"""
    finalize_design!(state)

INITRE end-of-keywords finalization (initre.f:360-414). The crucial, easy-to-miss
step: FVS overwrites the variable named `PI` with the plot count `IPTINV` — that is
the divisor NOTRE actually uses (NOT π). Also finalizes the stockable proportion
GROSPC into the reciprocal multiplier NOTRE applies.
"""
function finalize_design!(s::StandState)
    p = s.plot
    p.points_inv <= 0 && (p.points_inv = Int32(1))
    p.pi = Float32(p.points_inv)                      # PI := IPTINV (NOTRE divisor!)
    p.sample_weight < 0f0 && (p.sample_weight = Float32(p.points_inv))
    if p.gross_space < 0f0                             # default (set to -1 at stand start)
        g = (p.pi - Float32(p.nonstockable)) / p.pi
        g > 1f0 && (g = 1f0)
        (p.pi - Float32(p.nonstockable)) <= 0f0 && (g = 1f0)
        p.gross_space = g
    end
    p.gross_space = 1f0 / p.gross_space               # → reciprocal multiplier
    return s
end

"""
    notre!(state)

Compute trees-per-acre (`trees.tpa`) for every live record from the sampling
design. Mirrors notre.f: variable-plot trees (DBH ≥ BRK) scale by BAF/DBH²,
fixed-plot trees (DBH < BRK) by FPA, all times the stockable proportion GROSPC.
"""
function notre!(s::StandState)
    p, t = s.plot, s.trees
    fp  = p.fixed_plot_inv / p.pi
    p.total_fixed_plot > 0f0 && (fp = 1f0 / p.total_fixed_plot)
    vp  = p.baf * 183.3465f0 / p.pi
    fp2 = 0f0
    if p.baf <= 0f0
        vp = 0f0; fp2 = -p.baf / p.pi
    end
    brk = p.min_dbh_var_plot
    @inbounds for i in 1:t.n
        pr = t.tpa[i]; d = t.dbh[i]
        pr <= 0f0 && (pr = 1f0)
        pr = d < brk ? pr * fp : pr * vp / (d * d) + pr * fp2
        pr <= 0f0 && (pr = 9.0f-25)
        t.tpa[i] = pr * p.gross_space
    end
    return s
end

"""
    stand_tpa(state)   -> total trees per acre        (TPROB)
    stand_ba(state)    -> basal area ft²/acre         (BA)
    stand_qmd(state)   -> quadratic mean diameter, in (RMSQD)

Weighted reductions over the live trees. Pure (no mutation).
"""
function stand_tpa(s::StandState)
    t = s.trees; tot = 0f0
    @inbounds for i in 1:t.n; tot += t.tpa[i]; end
    return tot
end

function stand_ba(s::StandState)
    t = s.trees; ba = 0f0
    @inbounds for i in 1:t.n; ba += t.tpa[i] * BA_PER_TREE * t.dbh[i]^2; end
    return ba
end

function stand_qmd(s::StandState)
    t = s.trees; sd2 = 0f0; tpa = 0f0
    @inbounds for i in 1:t.n
        sd2 += t.tpa[i] * t.dbh[i]^2
        tpa += t.tpa[i]
    end
    return tpa > 0f0 ? sqrt(sd2 / tpa) : 0f0
end

"Reineke stand density index by summation: Σ tpa·(DBH/10)^1.605 (DR016/dense.f)."
function stand_sdi(s::StandState)
    t = s.trees; sdi = 0f0
    @inbounds for i in 1:t.n
        t.dbh[i] >= s.control.dbh_sdi && (sdi += t.tpa[i] * (t.dbh[i] / 10f0)^1.605f0)
    end
    return sdi
end
