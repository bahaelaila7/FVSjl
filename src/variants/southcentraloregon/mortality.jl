# =============================================================================
# mortality.jl (southcentraloregon) — SO mortality VARMRT efficiency (so/scomrt.f). Chunk 7.
#
# SO's morts.f is the STANDARD Prognosis/Reineke SDI mortality (LZEIDE=.FALSE. for R6 ⇒ Reineke QMD),
# identical in structure to the shared driver `mortality!(::AbstractVariant)` (southern/mortality.jl):
# background RI = 0.5/(1+exp(PMSC + PMD·D)) (mort_ri_scale(::SouthCentralOregon)=0.5, so/morts.f:490)
# merged with the SDI density self-thin, the excess distributed by VARMRT (shade-tolerance efficiency).
# So SO reuses the shared driver — it needs only (1) mort_ri_scale=0.5 (chunk 0), (2) the so_bratio bark
# in the DG/self-thin trajectory (the 4 southern/diameter_growth.jl bark dispatches, wired as `_so_*`),
# (3) the PMSC/PMD/VARADJ CSV columns (mort_bkgd_intercept/mort_bkgd_dbh/varmrt_varadj), and (4) this
# SO VARMRT efficiency. so/scomrt.f: EFFTR = PEFF·VARADJ·0.01 (SO uses ×0.01, unlike CA/SN ×0.1); PEFF
# is byte-identical (0.84525 − 0.01074·PCT + 2e-7·PCT³). PCT = the BA percentile (t.crown_ratio after stand_pct!).
# =============================================================================

# so/scomrt.f — EFFTR = PEFF·VARADJ[JSPC=ISPC]·0.01, PEFF = 0.84525 − 0.01074·PCT + 2e-7·PCT³.
function _varmrt_efftr!(efftr, s, ::SouthCentralOregon, t::TreeList, n::Int)
    varadj = s.coef.species[:varmrt_varadj]
    pass1 = 0.0f0
    @inbounds for i in 1:n
        pct = t.crown_ratio[i]
        peff = 0.84525f0 - 0.01074f0 * pct + 0.0000002f0 * fpow(pct, 3f0)
        peff > 1.0f0 && (peff = 1.0f0); peff < 0.01f0 && (peff = 0.01f0)
        efftr[i] = peff * varadj[Int(t.species[i])] * 0.01f0
        pass1 += t.tpa[i] * efftr[i]
    end
    return pass1
end
