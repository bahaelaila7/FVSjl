# =============================================================================
# mortality.jl (pacificnorthwest) — PN mortality (pn/morts.f). Chunk 7.
#
# FVSpn compiles WC's morts.f (FVSpn_buildDir/morts.f == FVSwc_buildDir/morts.f): the ORGANON logistic RIP, the
# integer-PASS density self-thin with PASS 1's carried TA, the SIZCAP floor and MISMRT. The ONLY PN difference is
# morts.f:337 — XSITE1 is the RAW DF SITEAR (VARACD≠'WC' ⇒ no Curtis→King conversion). Shared body:
# _wcpn_mortality! (westcascades/mortality.jl).
# =============================================================================

_wcpn_mort_xsite1(::PacificNorthwest, si16::Float32) = si16

mortality!(s::StandState, v::PacificNorthwest; fint::Float32 = 10.0f0, book_snags::Bool = true) =
    _wcpn_mortality!(s, v; fint = fint, book_snags = book_snags)
