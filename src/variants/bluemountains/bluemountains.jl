# =============================================================================
# bluemountains.jl — the BM (Blue Mountains) variant singleton + registration
#
# Ported from: bm/*.f (the FVS "BLUE MOUNTAINS PROGNOSIS" variant, Region-6). BM is a STANDARD
# WESTERN WYKOFF DDS (bm/dgf.f reads DGLD/DGBAL/DGCR/DGCRSQ/DGSIC/... with NO `CALL GEMDG`,
# DDS=CONSPP+DGLD·lnD+DGBAL·BAL+… at bm/dgf.f:512), a discount of the KT/IE/EM/UT engine —
# same diameter_growth!/height/crown/regent/mortality hooks, differing in the 18-species list +
# coefficients. BM grinit (bm/grinit.f): VARACD='BM', MAXSP=18, LZEIDE=.FALSE. (STAGE SDI, like
# EM/KT — UNLIKE UT/CR's Zeide), DGSD=1.5 (NOT 2.0!), FINT=10 (10-yr cycle), FINTH=5, seed 55329.
# BM bark is a POWER model DIB=BARK1·D^BARK2 (bm/bratio.f, NOT the IMAP-dispatched linear form) ⇒
# needs a bm_bratio in the DG chunk.
#
# Species JSP (18): WP WL DF GF MH WJ LP ES AF PP WB LM PY YC AS CW OS OH
# FIA:  119 073 202 017 264 064 108 093 019 122 101 113 231 042 746 747 299 998
# SPCTRN: BM uses a FOREST-DEPENDENT ASPT column (bm/spctrn.f SPCOUT=ASPT(I,4..13)); default col 4.
#
# Oracle: /workspace/.bmwork/FVSbm_clean (relink_bm.sh, gfortran-16 + isoc23 shim). VERIFIED
# BIT-EXACT ON GROWTH vs tests/FVSbm/bmt01.sum.save (volume may differ by the gfortran-16 delta).
# PORT IN PROGRESS (chunk 0 scaffold).
# =============================================================================

"""
    BlueMountains <: AbstractVariant

The FVS Blue Mountains variant (VARACD "BM", MAXSP = 18, Region-6) — a western Wykoff-DDS
variant (Stage SDI, DGSD=1.5). PORT IN PROGRESS (chunk 0 scaffold).
"""
struct BlueMountains <: AbstractVariant end

variant_code(::BlueMountains) = "BM"
nspecies(::BlueMountains) = 18
htg_period(::BlueMountains) = 10f0     # /CONTRL/ YR = 10 (bm FINT=10)

const BM_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "bluemountains"))

coefficients(::BlueMountains) = cached_coefficients(() -> load_species_coefficients(BM_DATADIR), "BM")
