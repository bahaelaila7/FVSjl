# =============================================================================
# fire/fuel_moisture.jl — fuel-moisture scenario + wind reduction (FFE chunk F5b)
#
# Ported from: bin/FVSsn_buildDir/fmmois.f (FMMOIS) + fmburn.f:390 / fmvinit.f
# (the canopy wind-reduction CANCLS/CORFAC).
#
# These supply the two environmental inputs the Rothermel model (F5-core) needs but the
# stand state doesn't carry: the dead/live fuel moistures (from a preset dryness model
# keyed to the fire's severity/season) and the midflame wind multiplier (the 20-ft wind
# reduced by canopy cover).
# =============================================================================

# FMMOIS preset fuel moistures by dryness model 1..4 (fmmois.f), fraction dry weight.
# Columns: dead 1hr(0–.25"), 10hr(.25–1"), 100hr(1–3"), 3+", duff; then live woody, herb.
const _FM_MOIS = Float32[
#   1hr   10hr  100hr 3+    duff  Lwoody Lherb
    0.05  0.07  0.12  0.17  0.40  0.55   0.55    # 1 very dry
    0.06  0.08  0.13  0.18  0.75  0.80   0.80    # 2 dry
    0.07  0.09  0.14  0.20  1.00  1.00   1.00    # 3 wet
    0.16  0.16  0.18  0.50  1.75  1.50   1.50]   # 4 very wet

"""
    fuel_moisture(fmois) -> Matrix{Float32}

Preset fuel-moisture values (fraction) for dryness model `fmois` (1 very dry … 4 very
wet) as a 2×5 matrix indexed `[1=dead/2=live, class]` (FMMOIS, fmmois.f). Dead classes
are 0–.25"/.25–1"/1–3"/3+"/duff; live are woody/herb in columns 1–2. The matrix is laid
out for direct use by `rothermel_surface_fire` (which reads dead 1–3 + dead-herb→1hr and
live 1–2).
"""
function fuel_moisture(fmois::Integer)::Matrix{Float32}
    @assert 1 <= fmois <= 4 "moisture model must be 1..4"
    r = @view _FM_MOIS[fmois, :]
    m = zeros(Float32, 2, 5)
    m[1, 1] = r[1]; m[1, 2] = r[2]; m[1, 3] = r[3]; m[1, 4] = r[4]; m[1, 5] = r[5]
    m[2, 1] = r[6]; m[2, 2] = r[7]
    return m
end

# Canopy-cover classes and their wind-reduction factors (fmvinit.f:36-44).
const _FM_CANCLS = (5.0f0, 17.5f0, 37.5f0, 75.0f0)
const _FM_CORFAC = (0.5f0, 0.3f0, 0.2f0, 0.1f0)

"""
    fire_wind_reduction(percov) -> Float32

Midflame wind multiplier (WMULT) for canopy cover `percov` (%): the 20-ft wind is
reduced toward midflame by interpolating the reduction factors CORFAC over the canopy-
cover classes CANCLS (fmburn.f:390 `WMULT=ALGSLP(PERCOV,CANCLS,CORFAC,4)`). Denser
canopy ⇒ smaller multiplier (more sheltering). `FWIND = wind · WMULT`.
"""
@inline function fire_wind_reduction(percov::Float32)::Float32
    percov < _FM_CANCLS[1]  && return _FM_CORFAC[1]
    percov >= _FM_CANCLS[4] && return _FM_CORFAC[4]
    @inbounds for i in 1:3
        if percov < _FM_CANCLS[i+1]
            return _FM_CORFAC[i] + (_FM_CORFAC[i+1] - _FM_CORFAC[i]) /
                   (_FM_CANCLS[i+1] - _FM_CANCLS[i]) * (percov - _FM_CANCLS[i])
        end
    end
    return _FM_CORFAC[4]
end
