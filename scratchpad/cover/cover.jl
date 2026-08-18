# cover.jl — FVS COVER understory-vegetation / canopy-cover extension (covr/ vcovr/ pg/).
#
# SCOPE (see scratchpad/cover/HANDOFF.md for full evidence):
#   COVER is a REPORT-ONLY extension. It reads the tree list + stand density and
#   emits a text report (CANOPY COVER STATISTICS / SHRUB STATISTICS / SUMMARY /
#   SHRUB-SMALL CONIFER COMPETITION) to output unit JOSHRB at end of projection.
#   It does NOT modify tree DBH/HT/PROB/growth/mortality — grep of covr/*.f
#   vcovr/*.f shows ZERO writes to any tree array (verified). Therefore:
#     - INERT  = a run without the COVER keyword is byte-identical (LCVATV=.FALSE.
#                → CVBROW/CVCNOP/CVOUT all early-return). Proven at the oracle:
#                FVSem_g16 cover_off.out ≡ stock (no COVER section emitted).
#     - VALID  = the emitted report text matches the live relinked oracle.
#   COVER draws NO random numbers (grep RANN/RANDOM/ran( over covr/vcovr/pg = none;
#   the only "random" refs are pg/{get,put}std serializing the ENGINE seed, not a
#   COVER LCG). So there is NO RNG to port.
#
# LINKAGE: real cvin.o (active) is compiled into ALL 18 western variant binaries
#   (ws ak bm ca cr ec op pn so tt wc ut ci em kt ie oc nc); the base/excov.f stub
#   is linked instead only in the 6 eastern variants (ls ne on sn bc cs). So the
#   oracle needs NO relink — the shipped western binary already fires COVER when the
#   base keyword COVER (keywds.f opt 46) is present.
#
# THIS FILE STAGES THE VALIDATED BEACHHEAD SLICE ONLY:
#   * CoverState (subset of CVCOM needed for the canopy-cover half)
#   * cover_init!         (CVINIT, covr/cvinit.f)
#   * cover_cvcw!         (CVCW, covr/cvcw.f)  — DUMP-REPLAY BIT-EXACT vs FVSem_g16
#                          (scratchpad/cover/cvcw_g16_dump.txt, 3 cycles, hex-exact).
#   The CVIN keyword reader, CVSHAP/CVCBMS/CVSUM/CVOUT report path, and the shrub
#   half (CVBROW/CVSCON/CVBCAL/CVCLAS) are NOT ported here — see HANDOFF "NEXT".

abstract type AbstractCoverState end

"""
    CoverState

Subset of the Fortran CVCOM common (common/CVCOM.F77) needed for the canopy-cover
computation beachhead. `active` mirrors LCVATV (set when the COVER activity 900 is
scheduled). Full port will extend this with the per-cycle/per-height-class arrays
(TXHT/VOLXHT/CFBXHT/PROXHT) and the shrub arrays (SCOV/TRSH/SHT/SPB).
"""
mutable struct CoverState <: AbstractCoverState
    active::Bool          # LCVATV: COVER activity scheduled → routinely called
    lcov::Bool            # LCOV:   activity has fired this/earlier cycle
    lcnop::Bool           # LCNOP:  CANOPY sub-keyword (crown shape+biomass) on
    lbrow::Bool           # LBROW:  SHRUBS sub-keyword on
    lcover::Bool          # LCOVER: emit CANOPY COVER STATISTICS table  (NOCOVOUT→false)
    lshrub::Bool          # LSHRUB: emit SHRUB STATISTICS table         (NOSHBOUT→false)
    lcvsum::Bool          # LCVSUM: emit SUMMARY table                  (NOSUMOUT→false)
    covopt::Int           # COVOPT: foliage-biomass eqn option (default 2)
    idt::Int              # cycle/date the COVER activity turns on
    joshrb::Int           # output unit for the report (default = main out)
    crarea::Float32       # CRAREA: sum of crown projection areas (sq ft/ac), CVCW out
    trecw::Vector{Float32}# TRECW(MAXTRE): per-tree crown width used by cover (=CRWDTH)
end

CoverState() = CoverState(false, false, false, false, true, true, true,
                          2, 1, 0, 0.0f0, Float32[])

# CVINIT (covr/cvinit.f): initialise per-stand cover flags/arrays. Called from INITRE.
function cover_init!(cv::CoverState)
    cv.lcov   = false
    cv.lcnop  = false
    cv.lbrow  = false
    cv.lcover = true
    cv.lshrub = true
    cv.lcvsum = true
    cv.crarea = 0.0f0
    empty!(cv.trecw)
    return cv
end

"""
    cover_cvcw!(cv, crown_width, prob, itrn) -> Float32

CVCW (covr/cvcw.f): per-tree crown width for cover = the base CRWDTH array (the
open-grown crown width FVS already computed for CCF), and the stand crown-projection
area CRAREA = (Σ TRECW(i)² · PROB(i)) · π/4.  Constant is FVS's literal 0.785398f0
(NOT Float32(pi)/4), and the accumulation is a Float32 running sum in tree order —
both are load-bearing for bit-exactness.

DUMP-REPLAY VALIDATED bit-exact (Float32 hex) vs FVSem_g16 over 3 cycles:
  icyc0 26605.12 (46cfda3d), icyc1 30998.373 (46f22cbf), icyc2 34740.656 (4707b4a8).
"""
function cover_cvcw!(cv::CoverState, crown_width::AbstractVector{Float32},
                     prob::AbstractVector{Float32}, itrn::Integer)
    resize!(cv.trecw, itrn)
    acc = 0.0f0
    @inbounds for i in 1:itrn
        cw = crown_width[i]            # TRECW(I) = CRWDTH(I)
        cv.trecw[i] = cw
        acc = acc + cw * cw * prob[i]  # Float32 running sum, tree order
    end
    cv.crarea = acc * 0.785398f0
    return cv.crarea
end
