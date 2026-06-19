# =============================================================================
# simulate.jl — the projection cycle loop (FVS / TREGRO orchestration)
#
# Ported from: base/fvs.f (cycle loop) + base/tregro.f (per-cycle driver).
#
# After initialization, FVS projects the stand cycle by cycle: recompute stand
# density, grow diameters and heights, apply mortality and regeneration, then
# recompute statistics. This is the skeleton; mortality (C4), regeneration (C4)
# and volume/output (C5) are wired in as those chunks land.
# =============================================================================

"""
    setup_growth!(state)

One-time growth setup after initialization (the FVS LSTART pass): compute the
site-dependent diameter-growth constants (DGCONS) and run the diameter-growth
calibration against the input measured growth (COR). Needs density set first.
"""
function setup_growth!(s::StandState)
    compute_density!(s)
    dgcons!(s)
    calibrate_diameter_growth!(s)        # LSTART — uses the input measured DG
    return s
end

"""
    compute_density!(state)

Recompute the per-cycle stand density quantities the growth models read:
basal area, average dominant height (AVH), and per-point basal area (PTBAA).
"""
function compute_density!(s::StandState)
    s.plot.basal_area = stand_ba(s)
    s.plot.avg_height = stand_top_height(s)
    point_basal_area!(s)
    stand_pct!(s)                      # PCT = stand BA percentile (for DGF competition)
    return s
end

"""
    grow_cycle!(state)

Advance the stand by one growth cycle: recompute density, compute diameter and
height growth, and apply them. (Mortality/regeneration/volume added in C4/C5.)
"""
function grow_cycle!(s::StandState)
    compute_density!(s)
    diameter_growth!(s, s.variant)
    height_growth!(s, s.variant)
    mortality!(s, s.variant)               # reduces tpa (uses the projected diameter)
    t = s.trees
    sd = s.coef.species
    bark_a = sd[:bark_intercept]; bark_b = sd[:bark_slope]
    @inbounds for i in 1:t.n
        # DG is the INSIDE-bark increment; outside-bark DBH grows by DG/bark, with
        # bark evaluated at the pre-growth DBH (update.f:115 / update.jl:75).
        bark = bark_ratio(bark_a, bark_b, t.species[i], t.dbh[i])
        t.dbh[i]    += t.diam_growth[i] / bark
        t.height[i] += t.ht_growth[i]
    end
    s.control.cycle += Int32(1)
    return s
end
