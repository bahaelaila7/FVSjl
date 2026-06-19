# =============================================================================
# bark_and_bounds.jl — bark ratio (BRATIO) and diameter-growth bounds (DGBND)
#
# Ported from: sn/bratio.f, sn/dgbnd.f.
#
# `bark_ratio` converts outside-bark DBH to the inside-bark ratio used by the
# growth model (per-species linear, clamped [0.80, 0.99]). `dg_bound` clamps a
# predicted diameter increment by the species' large-tree taper and size cap.
# (The IFOR=20 Fort Bragg special-case bark equations are noted for later.)
# =============================================================================

# Bark + DG-bound coefficients now live in data/southern/species_coefficients.csv
# (columns bark_intercept/bark_slope, dg_bound_dbh_lo/dg_bound_dbh_hi). The kernels
# below take the loaded per-species columns so they stay pure and allocation-free.

"""
    bark_ratio(bark_a, bark_b, sp, d) -> Float32

Inside-bark / outside-bark DBH ratio for species `sp` at DBH `d` (bratio.f), where
`bark_a`/`bark_b` are the `:bark_intercept`/`:bark_slope` per-species columns.
"""
@inline function bark_ratio(bark_a, bark_b, sp::Integer, d::Real)::Float32
    df = Float32(d)
    ratio = df > 0f0 ? (bark_a[sp] + bark_b[sp] * df) / df : 0.99f0
    ratio > 0.99f0 && (ratio = 0.99f0)
    ratio < 0.80f0 && (ratio = 0.80f0)
    return ratio
end

# Convenience: pull the bark columns straight from a coefficient container.
@inline bark_ratio(c::SpeciesCoefficients, sp::Integer, d::Real) =
    bark_ratio(c.species[:bark_intercept], c.species[:bark_slope], sp, d)

"""
    dg_bound(dlo_v, dhi_v, sp, dbh, ddg, sizcap) -> Float32

Bound a predicted diameter increment `ddg` by the species large-tree taper
(`dlo_v`/`dhi_v` = `:dg_bound_dbh_lo`/`:dg_bound_dbh_hi` columns) and the size cap
(`sizcap` = MAXSP×4 SIZCAP matrix). dgbnd.f.
"""
@inline function dg_bound(dlo_v, dhi_v, sp::Integer, dbh::Real, ddg_in::Real, sizcap::AbstractMatrix)::Float32
    dlo = dlo_v[sp]; dhi = dhi_v[sp]
    ddg = Float32(ddg_in); d = Float32(dbh)
    if d <= dlo
        # no adjustment
    elseif d > dhi
        ddg = 0.048f0
    else
        ddg = ddg * (1f0 + 0.90f0 * ((d - dlo) / (dlo - dhi)))
        ddg < 0.048f0 && (ddg = 0.048f0)
    end
    if (d + ddg) > sizcap[sp, 1] && sizcap[sp, 3] < 1.5f0
        ddg = sizcap[sp, 1] - d
        ddg < 0.01f0 && (ddg = 0.01f0)
    end
    return ddg
end
