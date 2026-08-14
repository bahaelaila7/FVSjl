# =============================================================================
# mortality.jl (westsierra) — WS mortality VARMRT efficiency (ws/varmrt.f). Chunk 7.
#
# WS's morts.f is the STANDARD Prognosis SDI mortality — identical in structure to the shared driver
# `mortality!(::AbstractVariant)`: background RI = 0.5/(1+exp(PMSC + PMD·D)) (mort_ri_scale(::WestSierra)=0.5,
# ws/morts.f:558+569), merged with the SDI density self-thin (WS is Zeide-SDI), excess distributed by VARMRT
# (shade-tolerance efficiency). So WS reuses the shared driver + (1) mort_ri_scale=0.5 (chunk 0), (2) the
# PMSC/PMD/VARADJ CSV columns (mort_bkgd_intercept/mort_bkgd_dbh/varmrt_varadj), and (3) this WS VARMRT.
#
# ★★ MEASURE CATCH (ws/varmrt.f:171-187): PEFF is byte-identical to SO (0.84525 − 0.01074·PCT + 2e-7·PCT³),
#   but the EFFTR scale FACTOR differs — WS CASE DEFAULT uses ×0.1 (NOT SO's ×0.01; SO's note: "SO uses ×0.01
#   unlike CA/SN ×0.1" — WS is the ×0.1 form), and GB(21) uses ×((100−CRI)/100)·VARADJ·0.01 (CRI=crown %).
# =============================================================================

# ws/varmrt.f — EFFTR = PEFF·VARADJ·0.1 (DEFAULT) ; GB(21) = PEFF·((100−CRI)/100)·VARADJ·0.01.
# PEFF = 0.84525 − 0.01074·PCT + 2e-7·PCT³. PCT = the BA percentile (t.crown_ratio after stand_pct!).
function _varmrt_efftr!(efftr, s, ::WestSierra, t::TreeList, n::Int)
    varadj = s.coef.species[:varmrt_varadj]
    pass1 = 0.0f0
    @inbounds for i in 1:n
        pct = t.crown_ratio[i]
        peff = 0.84525f0 - 0.01074f0 * pct + 0.0000002f0 * fpow(pct, 3f0)
        peff > 1.0f0 && (peff = 1.0f0); peff < 0.01f0 && (peff = 0.01f0)
        sp = Int(t.species[i])
        if sp == 21                                            # GB special: ·((100−CRI)/100)·VARADJ·0.01
            cri = Float32(t.crown_pct[i])
            efftr[i] = peff * ((100f0 - cri) / 100f0) * varadj[sp] * 0.01f0
        else                                                   # DEFAULT ·VARADJ·0.1 (WS ×0.1, unlike SO ×0.01)
            efftr[i] = peff * varadj[sp] * 0.1f0
        end
        pass1 += t.tpa[i] * efftr[i]
    end
    return pass1
end
