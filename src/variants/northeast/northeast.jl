# =============================================================================
# northeast.jl — the NE (Northeast) variant singleton + registration
#
# Ported from: ne/*.f (the FVS Northeast variant). ~93% of the Fortran is shared
# with SN and already implemented variant-agnostically; this directory holds only
# the NE-specific equations + data. See docs/NE_VARIANT_PORT_SCOPE.md.
#
# Validation oracle: the live Fortran `bin/FVSne_buildDir/FVSne` + the committed
# baseline `tests/FVSne/net01.sum.save` (the NE analog of snt01). NB FVSjulia
# ("Oracle A") does NOT have NE, so the live binary + Fortran source are the SOLE
# ground truth — bit-exact differential per chunk, same discipline as SN.
# =============================================================================

"""
    Northeast <: AbstractVariant

The FVS Northeast variant (VARACD `"NE"`, MAXSP = 108). Pass `Northeast()` as the
`variant` to `run_keyfile`/`each_stand`. The variant-specific equations (diameter
growth incl. the BAL competition term, height growth, mortality, crown, volume,
site) and data (`data/northeast/`) are ported chunk by chunk; until a given
interface method is implemented, dispatching it on `Northeast` errors loudly.
"""
struct Northeast <: AbstractVariant end

variant_code(::Northeast) = "NE"
nspecies(::Northeast) = 108
htg_period(::Northeast) = 10f0    # /CONTRL/ YR = 10 for NE (blkdat.f:71); HTCALC increment is 10-yr
mort_ri_scale(::Northeast) = 0.5f0                       # morts.f:504 — NE halves the background rate
mort_dbh_threshold(s, ::Northeast) = s.control.dbh_sdi   # morts.f:203 — NE gates on DBHSDI (default 0)

# NE coefficient/data directory (mirrors data/southern/). Built up chunk by chunk:
# species_translation.csv (the 108-species roster) is the first, most-upstream piece.
const NE_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "northeast"))

"Cached Northeast coefficient load (first call reads `data/northeast/*.csv`)."
coefficients(::Northeast) = cached_coefficients(() -> load_species_coefficients(NE_DATADIR), "NE")
