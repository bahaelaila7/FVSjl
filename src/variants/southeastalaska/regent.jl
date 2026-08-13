# =============================================================================
# regent.jl (southeastalaska) — AK small-tree height/diameter growth (ak/regent.f). Chunk 6 — STUB.
#
# ak/regent.f overrides growth for small trees (site-curve HTCALC-based height increment + a diameter
# back-out) and inserts establishment records. That model (HTCALC Hegyi/Payandeh site curves + the
# regen HD Curtis-Arney) is a separate chunk. Until it lands, this is a NO-OP: small trees keep the
# large-tree DGF/HTGF increments. akt01 is a mature stand (few sub-breast-height records), so the
# effect on its .sum is bounded; a seedling-heavy stand would need the real REGENT.
# =============================================================================

small_tree_growth!(s::StandState, stash, ::SoutheastAlaska; fint::Float32 = 10.0f0) = s
