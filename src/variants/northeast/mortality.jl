# =============================================================================
# mortality.jl (northeast) — NE-specific pieces of MORTS.
#
# The MORTS driver (SDI/Pretzsch density mortality, BAMAX cap, background rate,
# size-cap floor, MSBMRT) is shared from variants/southern/mortality.jl, made
# variant-generic via three hooks (htg_period / mort_ri_scale / mort_dbh_threshold)
# plus this variant-dispatched VARMRT efficiency fill.
#
# NE VARMRT (ne/varmrt.f) weights mortality toward suppressed trees by RELATIVE
# HEIGHT (HT/AVH) and species shade tolerance (1−VARADJ), where SN uses crown
# percentile and a shade-adjustment scalar. The geometric-progression distribution
# loop itself is identical and stays in the shared `_varmrt!`.
# =============================================================================

# CS shares NE's VARMRT efficiency (cs/varmrt.f ≡ ne/varmrt.f modulo the VARADJ DATA, which CS
# supplies via its own varmrt_varadj column). One method serves both eastern variants.
function _varmrt_efftr!(efftr, s, ::Union{Northeast,CentralStates,LakeStates}, t::TreeList, n::Int)
    tpa = t.tpa; sp = t.species; avh = s.plot.avg_height
    varadj = s.coef.species[:varmrt_varadj]
    pass1 = 0f0
    @inbounds for i in 1:n
        relht = avh > 0f0 ? min(t.height[i] / avh, 1f0) * 100f0 : 100f0
        peff = clamp(0.84525f0 - 0.01074f0 * relht + 0.0000002f0 * relht^3f0, 0.01f0, 1f0)
        efftr[i] = peff * (1f0 - varadj[sp[i]]) * 0.1f0
        pass1 += tpa[i] * efftr[i]
    end
    return pass1
end
