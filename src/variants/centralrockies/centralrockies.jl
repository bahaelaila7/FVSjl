# =============================================================================
# centralrockies.jl — the CR (Central Rockies) variant singleton + registration
#
# Ported from: cr/*.f (the FVS Central Rockies variant). This is the FIRST WESTERN
# variant (SN/NE/CS/LS are eastern). ~most of the Fortran is the shared, already-
# implemented variant-agnostic engine; this directory holds only CR-specific
# equations + data. Charter/doctrine: docs/CR_VARIANT_PORT_GOAL.md.
#
# Validation oracle: live Fortran relinked from bin/FVScr_buildDir/*.o
# (tmp/oracles/FVScr_new). NO FVSjulia oracle — the live binary + Fortran source are
# the SOLE ground truth, bit-exact differential per chunk (same doctrine as LS).
#
# CR infra (cr/blkdat.f / cr/grinit.f): VARACD='CR', MAXSP=38, YR=10 (10-yr cycle,
# like NE/CS/LS), RNG seed S0/SS=55329 (SAME as eastern). The genuinely-CR-specific
# models are the WESTERN Prognosis/Wykoff large-tree DDS (cr/dgf.f: DDS = f(ln D, BAL,
# CCF, CR, CR², site, habitat-type)), cr/htgf.f, cr/crown.f/cratet.f/ccfcal.f,
# cr/regent.f (small tree), and cr/sitset.f/habtyp.f (site index + HABITAT-TYPE groups,
# a coefficient dimension the eastern variants lack). Volume routes through the shared
# NVEL/NATCRS driver. Equations + data (data/centralrockies/) land chunk by chunk;
# until a hook is implemented, dispatching it on CentralRockies() errors loudly.
# =============================================================================

"""
    CentralRockies <: AbstractVariant

The FVS Central Rockies variant (VARACD "CR", MAXSP = 38). Pass `CentralRockies()`
as the `variant` to `run_keyfile`/`each_stand`. PORT IN PROGRESS (see CR_VARIANT_PORT_GOAL.md).
"""
struct CentralRockies <: AbstractVariant end

variant_code(::CentralRockies) = "CR"
nspecies(::CentralRockies) = 38
htg_period(::CentralRockies) = 10f0          # /CONTRL/ YR = 10 (cr/blkdat.f:126), like NE/CS/LS
# mort_ri_scale / mort_dbh_threshold: TODO — reconcile vs cr/morts.f in the mortality chunk.

const CR_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "centralrockies"))
