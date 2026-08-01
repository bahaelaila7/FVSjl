# =============================================================================
# easternmontana.jl — the EM (Eastern Montana) variant singleton + registration
#
# Ported from: em/*.f (the FVS "EASTERN MONTANA PROGNOSIS" variant). EM is the
# NEXT western Rockies CLUSTER variant after CR/KT/IE, and — like KT/IE — a STANDARD
# WESTERN WYKOFF DDS (em/dgf.f reads the classic coefficient arrays DGLD/DGBAL/DGCR/
# DGCRSQ/DGDBAL/DGHAB/DGFOR/DGCCF with NO `CALL GEMDG`), NOT CR's GENGYM model. So EM
# is a DISCOUNT of the KT/IE engine: same diameter_growth!/height/crown/regent hooks,
# differing mainly in (a) the 19-species list + coefficients and (b) the habitat/forest
# array DIMS — DGHAB(8,MAXSP) and DGFOR(6,MAXSP) (vs KT's 9/7). Measured, not inferred.
#
# EM grinit (em/grinit.f) is IDENTICAL to KT: LZEIDE=.FALSE. (Stage SDI, not CR's Zeide),
# DGSD=2.0 (DG serial-corr + ZZRAN), IFINT=10 (10-yr cycle), IFINTH=5, LHTDRG=.TRUE.,
# RNG seed 55329 (em/blkdat.f S0/SS — SAME as all variants).
#
# Validation oracle: live Fortran relinked from bin/FVSem_buildDir/*.o
# (/workspace/.emwork/FVSem_clean; relink_em.sh mirrors the KT recipe), VERIFIED vs
# tests/FVSem/emt01.sum.save. NO FVSjulia oracle — the live binary + Fortran source are
# the SOLE ground truth, bit-exact differential per chunk (same doctrine as CR/KT).
#
# EM infra (em/blkdat.f / PRGPRM.F77 / grinit.f): VARACD='EM', MAXSP=19, YR=10. Species
# JSP (19): WB WL DF LM LL RM LP ES AF PP GA AS CW BA PW NC PB OS OH. Equations + data
# (data/easternmontana/) land chunk by chunk; until a hook is implemented, dispatching it
# on EasternMontana() errors loudly (like CR/KT/LS did during their ports).
# =============================================================================

"""
    EasternMontana <: AbstractVariant

The FVS Eastern Montana variant (VARACD "EM", MAXSP = 19) — the next western Rockies
cluster variant after CR/KT/IE. Pass `EasternMontana()` as the `variant` to
`run_keyfile`/`each_stand`. PORT IN PROGRESS (chunk 0 scaffold; western Wykoff DDS,
KT-engine discount).
"""
struct EasternMontana <: AbstractVariant end

variant_code(::EasternMontana) = "EM"
nspecies(::EasternMontana) = 19
htg_period(::EasternMontana) = 10f0     # /CONTRL/ YR = 10 (em IFINT=10), like CR/KT/NE/CS/LS
# mort_ri_scale / mort_dbh_threshold: TODO — reconcile vs em/morts.f in the mortality chunk.

const EM_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "easternmontana"))

coefficients(::EasternMontana) = cached_coefficients(() -> load_species_coefficients(EM_DATADIR), "EM")
