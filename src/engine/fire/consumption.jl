# =============================================================================
# fire/consumption.jl — fire fuel consumption + carbon release (FFE F7/F8-rest)
#
# Ported from: bin/FVSsn_buildDir/fmcons.f (FMCONS, the natural-unpiled-fuels path).
#
# When a fire burns, it consumes a moisture-dependent fraction of each surface fuel size
# class; the consumed fuel is removed from the down-wood pools (`fire.cwd`) and, converted
# to carbon at 0.5, is the fire carbon *release*. Only the natural-fuels path is ported
# (the activity-fuels / piled-burn variants need the harvest-year context).
# =============================================================================

# PDIA: midpoints of the 6 large (>3") fuel size classes (fmcons.f:57).
const _FM_PDIA = (4.0f0, 8.0f0, 15.0f0, 15.0f0, 15.0f0, 15.0f0)

"""
    fire_consumption_fractions(mois) -> NTuple{11,Float32}

Consumed fraction of each of the 11 surface fuel classes for natural (unpiled) fuels
(FMCONS, fmcons.f:121-189). The <1" classes burn 90%, 1–3" 65%, litter 100%; the >3"
classes use a moisture-driven diameter-reduction `1 − ((PDIA−DIARED)/PDIA)²` and duff a
moisture-linear `(83.7 − 0.426·m_duff%)/100`. `mois` is the fuel-moisture matrix.
"""
function fire_consumption_fractions(mois::AbstractMatrix{Float32})::NTuple{11,Float32}
    m100 = mois[1, 4]                                  # 3+" (100-hr+) moisture drives the large classes
    diared = m100 > 1.25f0 ? 0f0 : max(0f0, 3.38f0 - 0.027f0 * m100 * 100f0)
    big = ntuple(i -> (pd = _FM_PDIA[i]; 1f0 - ((pd - diared) / pd)^2), 6)   # classes 4–9
    prduf = min(1f0, max(0f0, 83.7f0 - 0.426f0 * mois[1, 5] * 100f0) / 100f0)
    return (0.9f0, 0.9f0, 0.65f0, big[1], big[2], big[3], big[4], big[5], big[6], 1.0f0, prduf)
end

"""
    apply_fire_consumption!(fs, mois) -> Float32

Consume the surface fuel pools `fs.cwd` by the natural-fuels consumption fractions
(FMCONS) and return the fire carbon release (tons C/acre = consumed biomass × 0.5).
Mutates `fs.cwd` (the unburned remainder stays for the down-wood pools).
"""
function apply_fire_consumption!(fs::FireState, mois::AbstractMatrix{Float32})::Float32
    pr = fire_consumption_fractions(mois)
    released = 0f0
    @inbounds for isz in 1:11
        f = pr[isz]
        f > 0f0 || continue
        for k in 1:2, l in 1:4
            consumed = fs.cwd[isz, k, l] * f
            fs.cwd[isz, k, l] -= consumed
            released += consumed
        end
    end
    return released * 0.5f0
end
