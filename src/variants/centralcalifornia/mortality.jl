# =============================================================================
# mortality.jl (centralcalifornia) — CA mortality (ca/morts.f + ca/varmrt.f). Chunk 7.
#
# CA's morts.f is the STANDARD Prognosis/Reineke SDI mortality (LZEIDE=.FALSE. ⇒ Reineke QMD; CA's
# sitset sets LZEIDE off when CALCSDI blank), IDENTICAL in structure to the shared driver
# `mortality!(::AbstractVariant)` (southern/mortality.jl): background RI = 0.5/(1+exp(PMSC + PMD·D))
# merged with the SDI density rate (55%/85% Reineke self-thin, CEPMRT/SLPMRT persistence, IPASS
# QMD-convergence), the excess distributed by VARMRT (shade-tolerance geometric progression). CA's
# PMSC/PMD are the SO-variant 7-group background constants mapped to the 50 CA species by IBGMAP
# (folded into the per-species mort_bkgd_intercept/mort_bkgd_dbh CSV columns). So CA reuses the shared
# driver — it needs only (1) mort_ri_scale = 0.5 (ca/morts.f:495 "RI = 0.5·RI"), (2) the POWER bark in
# the DG/BARK self-thin trajectory (wired in southern/mortality.jl `_mbark`, #140/EC class), (3) the CSV
# columns, and (4) this CA VARMRT efficiency. ca/varmrt.f PEFF is byte-identical to SN/EC's; only the
# per-species VARADJ shade-tolerance differs (varmrt_varadj column, from ca/varmrt.f DATA VARADJ).
#
# NOTE (cornered): ca/morts.f floors RW/GS (sp 23/50) RI at 0.0001 before the ×0.5 halving. Their
# group-7 background (B0=2.5968, B1=0.51261) already gives RI≈0 for any D>0, so the floor is inert on
# any real stand and is not separately wired; cat01 (Sierra mixed-conifer) carries no coast redwood.
# =============================================================================

# mort_ri_scale(::CentralCalifornia)=0.5 is defined in centralcalifornia.jl (ca/morts.f:495 RI=0.5·RI).

# ca/varmrt.f — EFFTR = PEFF·VARADJ·0.1, PEFF = 0.84525 − 0.01074·PCT + 2e-7·PCT³ (PCT = BA percentile in
# t.crown_ratio after stand_pct!; no crown term, unlike CR). Same PEFF as SN/EC; CA's own VARADJ.
function _varmrt_efftr!(efftr, s, ::CentralCalifornia, t::TreeList, n::Int)
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
