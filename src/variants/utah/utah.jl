# =============================================================================
# utah.jl — the UT (Utah) variant singleton + registration
#
# Ported from: ut/*.f (the FVS "UTAH PROGNOSIS" variant, Region-4). UT is a STANDARD
# WESTERN WYKOFF DDS (ut/dgf.f reads the classic coefficient arrays DGLD/DGBAL/DGCR/
# DGCRSQ/DGDBAL/DGCCF/DGFOR/DGDS with NO `CALL GEMDG`), NOT CR's GENGYM model — so, like
# KT/IE/EM, it is a DISCOUNT of that shared Wykoff engine: same diameter_growth!/height/
# crown/regent hooks, differing mainly in (a) the 24-species list + coefficients and
# (b) array DIMS — DGFOR(5,MAXSP), DGDS(4,MAXSP).
#
# UT grinit (ut/grinit.f): VARACD='UT', MAXSP=24, LZEIDE=.TRUE. (Zeide SDI — like CR,
# UNLIKE KT/EM's Stage SDI), DGSD=2.0 (DG serial-corr + ZZRAN), FINT=10 (10-yr cycle),
# FINTH=5, FINTM=5, LHTDRG=.TRUE. except MC(20)/BI(21)=.FALSE., RNG seed 55329
# (ut/blkdat.f S0/SS — SAME as all variants). Heavy pinyon-juniper/woodland content
# (Utah specialty): PI/WJ/PM/RM/UJ/GB junipers+pinyons, GO gambel oak, MC mahogany.
#
# Species JSP (24): WB LM DF WF BS AS LP ES AF PP PI WJ GO PM RM UJ GB NC FC MC BI BE OS OH
# FIA:  101 113 202 015 096 746 108 093 019 122 106 064 814 133 066 065 142 749 748 475 322 313 299 998
# SPCTRN: UT maps to shared ASPT column 4 (ut/spctrn.f `SPCOUT=ASPT(I,4)`).
#
# Validation oracle: live Fortran relinked from bin/FVSut_buildDir/*.o
# (/workspace/.utwork/FVSut_clean via relink_ut.sh, gfortran-16 + isoc23 shim). VERIFIED
# BIT-EXACT ON GROWTH vs tests/FVSut/utt01.sum.save (TPA/BA/SDI/CCF/TopHt/QMD all identical);
# VOLUME cols differ ~6% = the gfortran-16-vs-original Float32 precision delta in the NVEL
# library (growth code is not precision-sensitive; volume is). So FVSut_clean is a faithful
# GROWTH oracle. NO FVSjulia oracle — the live binary + Fortran source are the SOLE ground
# truth, bit-exact differential per chunk (same doctrine as CR/KT/IE/EM/TT).
#
# PORT IN PROGRESS (chunk 0 scaffold). Until a hook is implemented, dispatching it on
# Utah() errors loudly (like CR/KT/EM did during their ports).
# =============================================================================

"""
    Utah <: AbstractVariant

The FVS Utah variant (VARACD "UT", MAXSP = 24, Region-4) — a western Wykoff-DDS variant
(Zeide SDI) with heavy pinyon-juniper/woodland content. Pass `Utah()` as the `variant`
to `run_keyfile`/`each_stand`. PORT IN PROGRESS (chunk 0 scaffold).
"""
struct Utah <: AbstractVariant end

variant_code(::Utah) = "UT"
nspecies(::Utah) = 24
htg_period(::Utah) = 10f0     # /CONTRL/ YR = 10 (ut FINT=10), like CR/KT/EM

const UT_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "utah"))

coefficients(::Utah) = cached_coefficients(() -> load_species_coefficients(UT_DATADIR), "UT")
