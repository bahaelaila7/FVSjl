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
    dub_missing_heights!(s)              # CRATET — dub HT=0 / resolve broken-top NORMHT
    setup_volume_equations!(s)           # VOLEQDEF — per-species NVEL equation ids
    compute_forest_type!(s)              # FORTYP — needed by dgf!'s forest-type term
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

"Number of early cycles that use deterministic record tripling (FVS ICL4)."
const TRIPLE_CYCLE_LIMIT = 2

"""
    grow_cycle!(state; fint=5f0) -> (; accretion, mortality)

Advance the stand by one growth cycle: recompute density, grow diameters/heights
(with record tripling), apply mortality, update dimensions, and recompute volumes.
Returns the period's per-acre cubic accretion and mortality (OACC/OMORT,
vols.f:190 / update.f:60): accretion = Σ(newCFV−oldCFV)·survivingTPA / fint,
mortality = Σ killedTPA·oldCFV / fint, both ÷ gross area. Requires the cycle-start
volumes to be present in `trees.cuft_vol` (run `compute_volumes!` once at setup).
"""
function grow_cycle!(s::StandState; fint::Float32 = 5f0)
    compute_density!(s)
    t = s.trees
    nlive = t.n                              # ORIGINAL live records (pre-tripling)
    # Cycle-start volume + TPA of the originals, for the period accounting.
    old_cfv = Float32[t.cuft_vol[i] for i in 1:nlive]
    old_tpa = Float32[t.tpa[i]      for i in 1:nlive]
    # Tripling is active only for the first few cycles (ICL4); afterwards growth is
    # the stochastic serial-correlation path.
    trip = Int(s.control.cycle) < TRIPLE_CYCLE_LIMIT
    stash = diameter_growth!(s, s.variant; tripling = trip)  # DGs only; no records yet
    height_growth!(s, s.variant)
    small_tree_growth!(s, stash; fint = fint)  # REGENT overrides DG/HTG for dbh < 3"
    mortality!(s, s.variant)               # MORTS on the ORIGINAL records (FVS order)
    g = s.plot.gross_space
    # Mortality volume (OMORT) is accounted on the originals, before tripling.
    mort = 0f0
    @inbounds for i in 1:nlive
        mort += (old_tpa[i] - t.tpa[i]) * old_cfv[i]
    end
    triple_records!(s, stash)              # TRIPLE after mortality (splits surviving TPA)
    # Per-record cycle-start CFV (tripled records inherit the originals' cycle-0 vol).
    n = t.n
    old_cfv2 = Float32[t.cuft_vol[i] for i in 1:n]
    sd = s.coef.species
    bark_a = sd[:bark_intercept]; bark_b = sd[:bark_slope]
    @inbounds for i in 1:n
        # DG is the INSIDE-bark increment; outside-bark DBH grows by DG/bark, with
        # bark evaluated at the pre-growth DBH (update.f:115 / update.jl:75).
        bark = bark_ratio(bark_a, bark_b, t.species[i], t.dbh[i])
        t.dbh[i]    += t.diam_growth[i] / bark
        t.height[i] += t.ht_growth[i]
    end
    compute_volumes!(s)                     # end-of-period volumes
    accr = 0f0
    @inbounds for i in 1:n
        accr += (t.cuft_vol[i] - old_cfv2[i]) * t.tpa[i]   # OACC over the tripled set
    end
    s.control.cycle += Int32(1)
    return (; accretion = accr / fint / g, mortality = mort / fint / g)
end
