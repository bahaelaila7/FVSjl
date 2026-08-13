# =============================================================================
# mortality.jl (westcascades) — WC mortality (vwc/morts.f). Chunk 7 — NOT YET PORTED.
#
# ⚠ PLACEHOLDER. WC mortality is NOT the shared MORTS self-thin driver: vwc/morts.f is a distinct
# ORGANON-style model (per-species BM0..BM5 logistic annual RIP + MCLASS/MVALUES size classes +
# SDIMAX self-thin credit). Porting/validating it is chunk 7 (separate from this task's chunks 6+8).
#
# This no-op method exists ONLY so the end-to-end wct01 run COMPLETES and the CYCLE-0 (year-1990) .sum
# row — the START-of-simulation inventory stand, which is pre-growth AND pre-mortality — can be
# validated bit-exact for TPA/BA/SDI/TopHt/QMD + the chunk-8 volume columns. Cycles 1+ (year ≥2000)
# will NOT match the oracle until chunk 7 lands (TPA will not decline), so multi-cycle rows are
# EXPECTED to diverge. Do not treat this as a faithful mortality port.
# =============================================================================

function mortality!(s::StandState, ::WestCascades; fint::Float32 = 10.0f0, book_snags::Bool = true)
    return s   # chunk 7 pending (vwc/morts.f ORGANON RIP) — no-op so cyc0 .sum is producible
end
