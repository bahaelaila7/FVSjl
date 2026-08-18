# =============================================================================
# ontario/crown.jl — Ontario crown-ratio dub/update (canada/on/crown.f).
#
# ON's crown.f IS the shared TWIGS NC-125 model (GTR NC-125) used by NE/CS/LS:
#     CRNEW = 10·( BCR1/(1 + BCR2·BA) + BCR3·(1 − exp(BCR4·D)) )        (crown.f:180)
# with the same ±1%/yr change limit, CRNMLT window, crown-length cap, top-kill reduction and
# 1/10/95 bounds. BCR1..BCR4 are per species (data/ontario/species_coefficients.csv crown_bcr1..4,
# transcribed from canada/on/crown.f DATA via a dump of FVSon_g16). BCR4's sign is baked into the
# constant (ON uses exp(BCR4·D) directly), matching the shared kernel's `fexp(bcr4·d)`.
#
# VALIDATED: the LSTART crown dub reproduces FVSon_g16's cyc0 ICR 8/8 given the oracle's crown-time
# stand BA (scratchpad/on/check_crown_formula.jl). End-to-end cyc0 crown bit-exactness is gated on
# the metric plot/TPA expansion (stand BA/QMD), the documented next density blocker.
# =============================================================================

crown_ratio_update!(s::StandState, ::Ontario; kwargs...) = _twigs_crown_update!(s; kwargs...)
