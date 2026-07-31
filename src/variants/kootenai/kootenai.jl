# =============================================================================
# kootenai.jl — the KT (Kootenai / Kaniksu / Tally Lake) variant singleton + registration
#
# Ported from: kt/*.f (the FVS "KOOTENAI, KANIKSU, TALLY LAKE PROGNOSIS" variant).
# KT is the FIRST of the western Rockies CLUSTER (KT/IE/EM/BM/TT/UT) that follow CR
# at a discount. Unlike CR (which is a GENGYM growth-and-yield model), KT uses the
# STANDARD WESTERN WYKOFF DDS — measured, not inferred (kt/dgf.f reads the classic
# coefficient arrays DGLD/DGCR/DGCRSQ/DGDBAL/DGHAB(9,MAXSP)/DGFOR(7,MAXSP)/DGDS/DGEL/
# DGEL2/DGSASP, with NO `CALL GEMDG`). So each western-cluster variant is a discount of
# BOTH (a) the eastern SN/NE/CS/LS Wykoff-DDS engine (structurally identical dgf) AND
# (b) CR's western framework (crown/ccfcal/regent/htgf), plus western terms the eastern
# variants lack (elevation DGEL/DGEL2, slope-aspect DGSASP, 9-group habitat DGHAB).
#
# Validation oracle: live Fortran relinked from bin/FVSkt_buildDir/*.o
# (/workspace/.ktwork/FVSkt_clean; relink_kt.sh mirrors the CR recipe). NO FVSjulia
# oracle — the live binary + Fortran source are the SOLE ground truth, bit-exact
# differential per chunk (same doctrine as CR/LS).
#
# KT infra (kt/blkdat.f / kt/grinit.f): VARACD='KT', MAXSP=11, IFINT=10 (10-yr cycle,
# like CR/NE/CS/LS), IFINTH=5, RNG seed 55329 (SAME as the others). Species JSP (11):
# WP WL DF GF WH RC LP ES AF PP OT (western Montana conifers). Equations + data
# (data/kootenai/) land chunk by chunk; until a hook is implemented, dispatching it on
# Kootenai() errors loudly (like CR/LS did during their ports).
# =============================================================================

"""
    Kootenai <: AbstractVariant

The FVS Kootenai / Kaniksu / Tally Lake variant (VARACD "KT", MAXSP = 11) — the first
of the western Rockies cluster (KT/IE/EM/BM/TT/UT). Pass `Kootenai()` as the `variant`
to `run_keyfile`/`each_stand`. PORT IN PROGRESS (chunk 0 scaffold; western Wykoff DDS).
"""
struct Kootenai <: AbstractVariant end

variant_code(::Kootenai) = "KT"
nspecies(::Kootenai) = 11
htg_period(::Kootenai) = 10f0          # /CONTRL/ YR = 10 (kt IFINT=10), like CR/NE/CS/LS
# mort_ri_scale / mort_dbh_threshold: TODO — reconcile vs kt/morts.f in the mortality chunk.

const KT_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "kootenai"))

coefficients(::Kootenai) = cached_coefficients(() -> load_species_coefficients(KT_DATADIR), "KT")
