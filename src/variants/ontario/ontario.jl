# =============================================================================
# ontario.jl — the ON (Ontario) variant singleton + registration.
# DRAFT chunk-0 scaffold (UNVALIDATED). Mirrors src/variants/britishcolumbia/britishcolumbia.jl.
# Place at: src/variants/ontario/ontario.jl  (create the dir).
#
# Ported from: canada/on/*.f (FVS "ONTARIO"). ON is NOT a Wykoff-DDS discount — it uses the
# **Penner** growth family (own model), so diameter_growth! is a fresh implementation, not a
# template discount. Distinguishing infra (MEASURED, canada/on/PRGPRM.F77 + grinit.f + blkdat.f):
#   • VARACD='ON', MAXSP=72, MAXTRE=6000, DGSD=2.0 (grinit.f:219), metric variant (METRIC.F77).
#   • Density: Reineke SDI exponent 1.605 (morts.f:45,263) with optional LZEIDE Zeide branch.
#   • Large-tree DG: canada/on/dgf.f — Penner (2006) ANNUAL diameter-increment model, 35 equations
#     (MAXEQ=35), OSPMAP maps the 72 FVS species → 35 growth species; 8 species carry an AGS/UGS
#     hardwood-quality code (2 eqn variants). Metric: DBHM=DBH*INtoCM, HTM=HT*FTtoM, SIM, BAM, QMDM.
#   • Height: canada/on/htont.f — Penner diameter-height models (LHYBRID / broken-down fallback).
#   • Site index: canada/on/sitset.f — Ontario SI equations (sitset.f:392+).
#   • Bark: canada/on/bratio.f (metric, BRATIO), IMAP-style.
# Species JSP (72): PJ PS RN RP PW SW SN BF SB TA CE HE SO CR AB AR CW MV MR CB EW ES ER BY BD MH
#   MB BE AW OW OP OB OC OR BO PN HB HP HU PG PT PB BW CH BT WB IW LB NC MM MS MT BB CA BH DF HT ML
#   GB SY CP CC PL WI WI WI SS AM JP WP SP BP  (blkdat.f:104-171)
# Oracle: /workspace/.onwork/FVSon_g16 (build_g16_on.sh). RUNS inline/keyword input to normal STOP
#   + emits .sum; the SQLite DATABASE tree-read path (dbstreesin) SEGFAULTS under gcc-16 — use the
#   CR-style inline TREEDATA (.tre) beachhead, or rebuild the oracle with gcc-15.
#
# STATUS: DGF go/no-go PASSED — the Penner large-tree DG core (diameter_growth.jl on_penner_dds)
#   reproduces FVSon_g16's per-tree DBHM/DIAGR/DDS 8/8 BIT-EXACT (Float32-hex) on an inline stand
#   incl. both AGS/UGS quality species. Remaining hooks (htg/crown/mort/vol/site/regen) + the
#   dgf! stand-context wiring + maple/beech HYBRID adjustment are un-ported (error loudly).
# =============================================================================

"""
    Ontario <: AbstractVariant

The FVS Ontario variant (VARACD "ON", MAXSP = 72) — own **Penner** annual diameter-increment
growth model (NOT Wykoff-DDS). DRAFT chunk-0 scaffold (UNVALIDATED).
"""
struct Ontario <: AbstractVariant end

variant_code(::Ontario) = "ON"
nspecies(::Ontario) = 72
htg_period(::Ontario) = 10f0   # CONFIRM vs canada/on/grinit.f IFINT/YR before relying on this

const ON_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "ontario"))
coefficients(::Ontario) = cached_coefficients(() -> load_species_coefficients(ON_DATADIR), "ON")
