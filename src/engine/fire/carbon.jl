# =============================================================================
# fire/carbon.jl — standing live-tree carbon pools (FFE chunk F8 — carbon)
#
# Ported from: bin/FVSsn_buildDir/fmcrbout.f (the live-tree pools of the FFE Stand
# Carbon Report), which sums the Jenkins biomass (FMCBIO, ported as `jenkins_biomass`)
# over the tree list and converts biomass → carbon at 0.5 (fmcrbout.f:89/158).
#
# These are the live aboveground / merchantable / belowground carbon pools (tons C/ac).
# The dead pools (snags, down wood, forest floor) build on the F7 snag/CWD model; this
# is the live-tree foundation. Carbon appears in the DBS Carbon report, not the `.sum`.
# =============================================================================

"""
    stand_live_carbon(s) -> (; aboveground, merch, belowground)

Standing live-tree carbon pools in tons C/acre (FFE Stand Carbon Report, fmcrbout.f):
the per-tree Jenkins aboveground / merchantable / belowground (root) biomass summed over
the tree list (weighted by TPA) and converted to carbon at the 0.5 biomass→carbon ratio.
"""
function stand_live_carbon(s::StandState)
    t = s.trees; coef = s.coef
    above = 0f0; merch = 0f0; root = 0f0
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        a, m, r = jenkins_biomass(coef, t.species[i], t.dbh[i])
        above += a * t.tpa[i]
        merch += m * t.tpa[i]
        root  += r * t.tpa[i]
    end
    return (aboveground = above * 0.5f0, merch = merch * 0.5f0, belowground = root * 0.5f0)
end

"""
    standing_dead_carbon(s) -> Float32

Standing-dead (snag) carbon pool in tons C/acre: the Jenkins aboveground biomass of each
snag cohort (`fire.snags`, F7) weighted by its still-standing density (hard + soft),
converted to carbon at 0.5 (fmcrbout.f). Zero when FFE is off / no snags.
"""
function standing_dead_carbon(s::StandState)::Float32
    fs = s.fire; fs === nothing && return 0f0
    sn = fs.snags; coef = s.coef
    c = 0f0
    @inbounds for i in eachindex(sn.sp)
        den = sn.den_hard[i] + sn.den_soft[i]
        den > 0f0 || continue
        a, _, _ = jenkins_biomass(coef, sn.sp[i], sn.dbh[i])
        c += a * den
    end
    return c * 0.5f0
end

"""
    down_wood_carbon(s) -> Float32

Down dead wood + forest-floor carbon pool in tons C/acre: the FFE surface fuel pools
(`fire.cwd`, the dead down-wood/litter/duff loadings in tons/ac, F3) converted to carbon
at 0.5. Zero when FFE is off.
"""
function down_wood_carbon(s::StandState)::Float32
    fs = s.fire; fs === nothing && return 0f0
    return sum(fs.cwd) * 0.5f0
end

"""
    stand_carbon(s) -> (; live_above, live_below, standing_dead, down_wood, total)

All the main stand carbon pools (tons C/acre): live aboveground + belowground (trees),
standing dead (snags), and down dead wood / forest floor (FFE Stand Carbon Report).
"""
function stand_carbon(s::StandState)
    lc = stand_live_carbon(s)
    sd = standing_dead_carbon(s)
    dw = down_wood_carbon(s)
    return (; live_above = lc.aboveground, live_below = lc.belowground,
            standing_dead = sd, down_wood = dw,
            total = lc.aboveground + lc.belowground + sd + dw)
end
