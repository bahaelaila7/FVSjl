# =============================================================================
# crown.jl (bluemountains) — CCF / crown width (bm/ccfcal.f). Chunk 5 (partial: the CCF
# piece is a prerequisite for chunk-3 DG, since RELDEN = stand CCF feeds CONSPP).
#
# bm/ccfcal.f MODE=1 per-tree CCFT (Paine & Hann, Oregon RP46):
#   CASE(1:12,15,17):  D>=1  → RD1 + D·RD2 + D²·RD3 ; 0.1<D<1 → RDA·D^RDB ; D<=0.1 → 0.001
#   CASE(13,14,16,18): D<1   → D·(RD1+RD2+RD3)      ; D>=1    → RD1 + RD2·D + RD3·D²
# Stand CCF = Σ CCFT·TPA = RELDEN.  Coeffs from data/bluemountains/ccf_coeffs_bm.csv (verified vs source).
# =============================================================================

let
    path = joinpath(BM_DATADIR, "ccf_coeffs_bm.csv")
    rows = [split(strip(l), ',') for l in readlines(path)[2:end]]
    col(j) = Float32[parse(Float32, rows[sp][j]) for sp in 1:18]
    global const BM_RD1 = col(2); global const BM_RD2 = col(3); global const BM_RD3 = col(4)
    global const BM_RDA = col(5); global const BM_RDB = col(6)
end

@inline function bm_tree_ccf(sp::Integer, d::Real)::Float32
    d <= 0f0 && return 0f0
    dd = Float32(d)
    if sp == 13 || sp == 14 || sp == 16 || sp == 18
        return dd < 1f0 ? dd * (BM_RD1[sp] + BM_RD2[sp] + BM_RD3[sp]) :
                          BM_RD1[sp] + BM_RD2[sp] * dd + BM_RD3[sp] * dd * dd
    else
        return dd >= 1f0 ? BM_RD1[sp] + dd * BM_RD2[sp] + dd * dd * BM_RD3[sp] :
               dd > 0.1f0 ? BM_RDA[sp] * dd ^ BM_RDB[sp] : 0.001f0
    end
end
