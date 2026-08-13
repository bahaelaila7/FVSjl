# =============================================================================
# mortality.jl (eastcascades) — EC mortality (ec/morts.f + ec/varmrt.f). Chunk 7.
#
# EC's morts.f is the STANDARD Stage/Reineke SDI mortality (LZEIDE=.FALSE. ⇒ Reineke QMD), IDENTICAL in
# structure to the shared N-Rockies driver `mortality!(::AbstractVariant)` (southern/mortality.jl): a
# background rate RI = 0.5/(1+exp(PMSC + PMD·D)) merged with the SDI density rate (55%/85% Reineke self-
# thin), the excess distributed by VARMRT (shade-tolerance-weighted geometric progression). So EC reuses
# the shared driver — it needs only (1) mort_ri_scale = 0.5 (ec/morts.f RI=0.5·RI), (2) the POWER bark in
# the DG/BARK trajectory (wired in southern/mortality.jl `_mbark`), (3) PMSC/PMD/VARADJ in the species CSV
# (mort_bkgd_intercept/mort_bkgd_dbh/varmrt_varadj), and (4) this EC VARMRT efficiency method. ec/varmrt.f
# PEFF is byte-identical to SN's; only the per-species VARADJ shade-tolerance differs.
# =============================================================================

# ec/varmrt.f — EFFTR = PEFF·VARADJ·0.1, PEFF = 0.84525 − 0.01074·PCT + 2e-7·PCT³ (PCT = BA percentile,
# t.crown_ratio after stand_pct!). Same PEFF as SN; EC's own VARADJ (varmrt_varadj column).
function _varmrt_efftr!(efftr, s, ::EastCascades, t::TreeList, n::Int)
    varadj = s.coef.species[:varmrt_varadj]
    pass1 = 0.0f0
    @inbounds for i in 1:n
        pct = t.crown_ratio[i]
        peff = 0.84525f0 - 0.01074f0 * pct + 0.0000002f0 * fpow(pct, 3f0)
        peff > 1.0f0 && (peff = 1.0f0); peff < 0.01f0 && (peff = 0.01f0)
        efftr[i] = peff * varadj[Int(t.species[i])] * 0.1f0
        pass1 += t.tpa[i] * efftr[i]
    end
    return pass1
end
