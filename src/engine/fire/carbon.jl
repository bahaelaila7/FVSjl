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
