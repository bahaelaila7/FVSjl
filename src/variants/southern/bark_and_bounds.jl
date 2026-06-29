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
    ddg = Float32(ddg_in); d = Float32(dbh)
    if dlo_v !== nothing                          # SN DGBND: the dlo/dhi DBH-range adjustment (NE has none)
        dlo = dlo_v[sp]; dhi = dhi_v[sp]
        if d <= dlo
            # no adjustment
        elseif d > dhi
            ddg = 0.048f0
        else
            ddg = ddg * (1f0 + 0.90f0 * ((d - dlo) / (dlo - dhi)))
            ddg < 0.048f0 && (ddg = 0.048f0)
        end
    end
    if (d + ddg) > sizcap[sp, 1] && sizcap[sp, 3] < 1.5f0   # both variants: the SIZCAP cap (ne/dgbnd.f)
        ddg = sizcap[sp, 1] - d
        ddg < 0.01f0 && (ddg = 0.01f0)
    end
    return ddg
end

"""
    _bound_scale(dlo_v, dhi_v, sp, dbh, d_ib, dg5, sfint, sizcap) -> Float32

Apply `DGBND` to the **5-yr** diameter increment `dg5` (dgdriv.f:267-269), THEN scale the bounded
increment to the `sfint`-year cycle (gradd.f:79-90: `DDS=(DG·(2·d_ib+DG))·(FINT/5)`,
`DG=sqrt(d_ib²+DDS)−d_ib`) WITHOUT re-bounding. `sfint==5` ⇒ identity (just the bound).
"""
@inline function _bound_scale(dlo_v, dhi_v, sp::Integer, dbh::Real, d_ib::Real, dg5::Real,
                              sfint::Real, sizcap::AbstractMatrix, yr::Real = 5f0)::Float32
    dg = dg_bound(dlo_v, dhi_v, sp, dbh, dg5, sizcap)
    s = Float32(sfint); y = Float32(yr)             # YR = the DG model's native period (SN 5, NE 10)
    (s != y && dg > 0f0) || return dg               # gradd.f:79 IF(FINT.NE.YR .AND. DG.GT.0)
    dib = Float32(d_ib)
    dds = dg * (2f0 * dib + dg) * (s / y)
    return sqrt(dib * dib + dds) - dib
end
