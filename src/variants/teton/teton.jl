# =============================================================================
# teton.jl — the TT (Teton) variant singleton + registration
#
# Ported from: tt/*.f (the FVS "TETON" variant, Region-4 Intermountain). TT is the
# NEXT western Rockies CLUSTER variant after CR/KT/IE/EM, and the cleanest remaining
# DISCOUNT — it reuses BOTH already-built engines (measured, not inferred; see
# fvsjl-tt-variant-port memory):
#   • DG = EM's Wykoff-DDS structure. tt/dgf.f is inline DDS (NO `CALL GEMDG`) with the
#     same habitat+forest+DSQ dimensions as EM: DGHAB(7)/DGFOR(5,MAXSP)/DGDS(4,MAXSP)
#     (EM is DGHAB(8)/DGFOR(6)). ⇒ clone em_dgcons!/dgf!, resize the arrays, swap coeffs.
#   • Density = CR's Zeide SDI. tt/grinit.f sets LZEIDE=.TRUE. (EM/KT are .FALSE.=Stage) —
#     the ONE cross-engine difference. TT = EM-DG + CR-Zeide-density.
#
# TT grinit (tt/grinit.f): DGSD=2.0, IFINT=10 (10-yr cycle), IFINTH=5, LZEIDE=.TRUE.,
# LHTDRG=.TRUE. (except sp13 BI/sp16 MC .FALSE.), RNG seed 55329 (shared).
#
# Validation oracle: live Fortran relinked from bin/FVStt_buildDir/*.o (670 .o) →
# /workspace/.ttwork/FVStt_clean (relink_tt.sh mirrors the EM recipe), VERIFIED bit-exact
# vs tests/FVStt/ttt01.sum.save (all 11 cycles). NO FVSjulia oracle — the live binary +
# Fortran source are the SOLE ground truth, bit-exact differential per chunk (CR/KT/EM doctrine).
# ttt01 = the SAME 248112 conifer stand as EM's emt01 ⇒ chunk-3 DG reuses EM's main-conifer path.
#
# TT infra (tt/blkdat.f / PRGPRM.F77 / grinit.f): VARACD='TT', MAXSP=18, YR=10. Species JSP
# (18): WB LM DF PM BS AS LP ES AF PP UJ RM BI MM NC MC OS OH; FIA 101 113 202 133 096 746
# 108 093 019 122 065 066 322 321 749 475 299 998. Equations + data (data/teton/) land chunk
# by chunk; until a hook is implemented, dispatching it on Teton() errors loudly (like EM did).
# =============================================================================

"""
    Teton <: AbstractVariant

The FVS Teton variant (VARACD "TT", MAXSP = 18) — the next western Rockies cluster
variant after CR/KT/IE/EM. Pass `Teton()` as the `variant` to `run_keyfile`/`each_stand`.
PORT IN PROGRESS (chunk 0 scaffold; western Wykoff DDS = EM-engine discount + CR Zeide SDI).
"""
struct Teton <: AbstractVariant end

variant_code(::Teton) = "TT"
nspecies(::Teton) = 18
htg_period(::Teton) = 10f0     # /CONTRL/ YR = 10 (tt IFINT=10), like CR/KT/EM/NE/CS/LS

const TT_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "teton"))

coefficients(::Teton) = cached_coefficients(() -> load_species_coefficients(TT_DATADIR), "TT")
