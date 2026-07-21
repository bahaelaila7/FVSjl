# =============================================================================
# mortality.jl (centralrockies) — CR VARMRT efficiency (cr/varmrt.f)
#
# CR distributes the SDI/background mortality across trees by a per-tree
# "efficiency" EFFTR = PEFF·((100−CRI)/100)·VARADJ·0.01, where PEFF is a cubic
# in the BA percentile PCT (NOT relative height, unlike NE/CS/LS), CRI is the
# crown %, and VARADJ is species shade-tolerance. Oaks (23-27) cap CRI at 50.
# The rest of the mortality driver is shared (southern/mortality.jl).
# =============================================================================

function _varmrt_efftr!(efftr, s, ::CentralRockies, t::TreeList, n::Int)
    varadj = s.coef.species[:varmrt_varadj]
    pass1 = 0.0f0
    @inbounds for i in 1:n
        sp = Int(t.species[i])
        cri = Float32(t.crown_pct[i])
        (23 <= sp <= 27 && cri > 50.0f0) && (cri = 50.0f0)   # oak crown-ratio ramp-down (varmrt.f:113)
        pct = t.crown_ratio[i]                                # PCT = BA percentile
        peff = 0.84525f0 - 0.01074f0 * pct + 0.0000002f0 * fpow(pct, 3.0f0)
        peff > 1.0f0 && (peff = 1.0f0); peff < 0.01f0 && (peff = 0.01f0)
        efftr[i] = peff * ((100.0f0 - cri) / 100.0f0) * varadj[sp] * 0.01f0
        pass1 += t.tpa[i] * efftr[i]
    end
    return pass1
end
