# =============================================================================
# mortality.jl (bluemountains) — BM VARMRT self-thin kill-distribution efficiency (bm/bmtmrt.f)
#
# BM routes through the SHARED MORTS driver (southern/mortality.jl): background Hamilton RI (×0.5)
# + SDI self-thinning via `_pretzsch_tn10` (the 55%/85%-line target WITH the CEPMRT/SLPMRT
# persistence) + the IPASS QMD-convergence loop + BAMAX/SIZCAP caps — exactly like CR. The only
# BM-specific piece is the per-tree kill-distribution efficiency EFFTR (bm/bmtmrt.f): original BM
# species (1-5,7-10,17) use a DBH cubic; the added SO/UT/TT/WC species (6,11-16,18) use the SO
# percentile form. EFFTR = PEFF·VARADJ·0.01 (★ ·0.01, cf. the SN ·0.1). PCT = t.crown_ratio (the BA
# percentile filled by stand_pct!). Background PMSC/PMD + VARADJ are in the species CSV
# (mort_bkgd_intercept/mort_bkgd_dbh/varmrt_varadj). Bark is bm_bratio (POWER) via the shared _mbark.
# =============================================================================

mort_ri_scale(::BlueMountains) = 0.5f0     # bm/morts.f:471 RI = 0.5·RI (Hamilton half-rate, like CR/NE)

# bm/bmtmrt.f efficiency: original species DBH cubic; added species (6,11-16,18) SO percentile form.
function _varmrt_efftr!(efftr, s, ::BlueMountains, t::TreeList, n::Int)
    varadj = s.coef.species[:varmrt_varadj]
    pass1 = 0.0f0
    @inbounds for i in 1:n
        sp = Int(t.species[i]); d = t.dbh[i]
        if sp == 6 || (11 <= sp <= 16) || sp == 18
            pct = t.crown_ratio[i]                                     # BA percentile (stand_pct!)
            peff = 0.84525f0 - 0.01074f0 * pct + 0.0000002f0 * (pct * pct * pct)  # PCT**3 → X·X·X
        else
            peff = (14.94435f0 - 0.69929f0 * d + 0.00868f0 * d * d) * 0.1f0
            d > 40.0f0 && (peff = 0.086f0)
        end
        peff > 1.0f0 && (peff = 1.0f0); peff < 0.01f0 && (peff = 0.01f0)
        efftr[i] = peff * varadj[sp] * 0.01f0                          # ★ ·0.01 (bmtmrt.f), not the SN ·0.1
        pass1 += t.tpa[i] * efftr[i]
    end
    return pass1
end
