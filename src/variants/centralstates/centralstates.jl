# =============================================================================
# centralstates.jl — the CS (Central States) variant singleton + registration
#
# Ported from: cs/*.f (the FVS Central States variant). ~94% of the Fortran is
# shared with SN/NE and already implemented variant-agnostically; this directory
# holds only the CS-specific equations + data. See docs/CS_VARIANT_PORT_SCOPE.md
# and docs/CS_GOAL.md.
#
# Validation oracle: the live Fortran relinked from `bin/FVScs_buildDir/*.o`
# (`test/harness/cs_oracle.sh`); canonical stand `tests/FVScs/cst01.key`. There is
# NO FVSjulia ("Oracle A") for CS — the live binary + Fortran source are the SOLE
# ground truth, bit-exact differential per chunk (same doctrine as SN/NE).
#
# CS infra (cs/blkdat.f / cs/grinit.f): MAXSP=96, RNG seed S0/SS=55329 (same as
# SN/NE), YR=10 (10-yr cycle, like NE), LZEIDE=.TRUE. (Zeide SDI, like NE).
# Most variant routines are NEAR-IDENTICAL to NE (htgf 96%, cratet 98%, regent 96%,
# varmrt 94%, crown 88%, dgdriv shared, balmod framework); the one genuinely-new
# model is cs/dgf.f — an SN-family ln(DDS) regression (BA-percentile/QMD), NOT NE's
# BAL-potential iteration. Equations + data (`data/centralstates/`) land chunk by
# chunk; until a hook is implemented, dispatching it on `CentralStates` errors loudly.
# =============================================================================

"""
    CentralStates <: AbstractVariant

The FVS Central States variant (VARACD `"CS"`, MAXSP = 96). Pass `CentralStates()`
as the `variant` to `run_keyfile`/`each_stand`.
"""
struct CentralStates <: AbstractVariant end

variant_code(::CentralStates) = "CS"
nspecies(::CentralStates) = 96
htg_period(::CentralStates) = 10f0   # /CONTRL/ YR = 10 (cs/blkdat.f:68), like NE
mort_ri_scale(::CentralStates) = 0.5f0                   # background-rate halving (morts.f), as NE — verify vs cs
mort_dbh_threshold(s, ::CentralStates) = s.control.dbh_sdi

const CS_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "centralstates"))

"Cached Central States coefficient load (first call reads `data/centralstates/*.csv`)."
coefficients(::CentralStates) = get!(() -> load_species_coefficients(CS_DATADIR), _COEF_CACHE, "CS")
