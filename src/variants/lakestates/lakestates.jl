# =============================================================================
# lakestates.jl — the LS (Lake States) variant singleton + registration
#
# Ported from: ls/*.f (the FVS Lake States variant). ~90%+ of the Fortran is shared
# with SN/NE/CS and already implemented variant-agnostically; this directory holds
# only the LS-specific equations + data. See docs/LS_VARIANT_PORT_SCOPE.md and
# docs/LS_GOAL.md.
#
# Validation oracle: the live Fortran relinked from `bin/FVSls_buildDir/*.o`
# (`test/harness/ls_oracle.sh`); canonical stand `tests/FVSls/lst01.key`. There is
# NO FVSjulia ("Oracle A") for LS — the live binary + Fortran source are the SOLE
# ground truth, bit-exact differential per chunk (same doctrine as SN/NE/CS).
#
# LS infra (ls/blkdat.f / ls/grinit.f): MAXSP=68, RNG seed S0/SS=55329 (same as
# SN/NE/CS), YR=10 (10-yr cycle, like NE/CS), LZEIDE=.TRUE. (Zeide SDI, like NE/CS).
# Mortality (ls/morts.f): RI=0.5·RI background halving (line 504, like NE/CS) and the
# SDI DBH gate reads DBHSDI (line 203, like CS). Most variant routines are near-NE/CS
# (htcalc.f IDENTICAL to NE); the genuinely-LS-specific models are ls/dgf.f (an
# SN-family ln(DDS) per-species regression, ~508 lines — extend the Southern dgf!
# framework, NOT NE's BAL-potential iteration), ls/crown.f, and the volume/merch path.
# Equations + data (`data/lakestates/`) land chunk by chunk; until a hook is
# implemented, dispatching it on `LakeStates` errors loudly.
# =============================================================================

"""
    LakeStates <: AbstractVariant

The FVS Lake States variant (VARACD `"LS"`, MAXSP = 68). Pass `LakeStates()`
as the `variant` to `run_keyfile`/`each_stand`.
"""
struct LakeStates <: AbstractVariant end

variant_code(::LakeStates) = "LS"
nspecies(::LakeStates) = 68
htg_period(::LakeStates) = 10f0   # /CONTRL/ YR = 10 (ls/blkdat.f:71), like NE/CS
mort_ri_scale(::LakeStates) = 0.5f0                      # background-rate halving (ls/morts.f:504), as NE/CS
mort_dbh_threshold(s, ::LakeStates) = s.control.dbh_sdi  # ls/morts.f:203 IF(D.LT.DBHSDI), as CS

const LS_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "lakestates"))

"Cached Lake States coefficient load (first call reads `data/lakestates/*.csv`)."
coefficients(::LakeStates) = cached_coefficients(() -> load_species_coefficients(LS_DATADIR), "LS")
