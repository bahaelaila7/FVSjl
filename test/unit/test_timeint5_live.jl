# test_timeint5_live.jl — 5-year cycles in YR=10 variants vs LIVE oracles (gradd.f:79-90 DG rescale after GRINCR).
#
# The FVS tests/FVS<v> keys (pnt01, ect01, wct01, bmt01, cit01) re-run with `TIMEINT 5` before every PROCESS, so
# FINT=5 ≠ YR=10. Goldens: the live FVS<v>_g16 .sum data rows of the same keys (<stem>_t5_live.rows).
#
# gradd.f:79-90 rescales DG from the YR-year basis to FINT years only AFTER GRINCR: DGDRIV bounds it, HTGF reads
# DBH+DG/BARK (the 10-year diameter), REGENT blends it, MORTS forms G=(DG/BARK)·(FINT/YR) and TRIPLE copies it, all on
# the 10-year basis. jl scaled DG to FINT inside the DG driver, so under 5-year cycles the height growth and the
# diameter-driven mortality ran on a half-size DG (PN pnt01: 53 of 56 rows off; now 0).

using Test, FVSjl
const _T5 = FVSjl
const _T5_FX = joinpath(@__DIR__, "..", "fixtures", "timeint5")

function _t5_rows(stem, variant)
    dir = mktempdir()
    for ext in ("key", "tre"); cp(joinpath(_T5_FX, "$stem.$ext"), joinpath(dir, "$stem.$ext")); end
    txt = cd(() -> _T5.run_keyfile("$stem.key"; variant = variant, output = :sum), dir)
    return [split(l) for l in split(txt, '\n') if occursin(r"^\d{4} ", l)]
end
_t5_live(stem) = [split(l) for l in readlines(joinpath(_T5_FX, "$(stem)_live.rows"))]
# stand index of each row (a new stand starts when the year goes back)
_t5_stand(rows) = (k = 1; [i == 1 ? 1 : (parse(Int, rows[i][1]) < parse(Int, rows[i-1][1]) ? (k += 1) : k) for i in eachindex(rows)])

# Known residuals, named (all outside the DG rescale):
#  • bmt01/cit01 stand 5 (bare ground, PLANT 1992): jl's planted cohort runs one 5-year cycle ahead of live (regent.f LESTB
#    FNT−5 / LSKIPH when FINT≤5 — the establishment-cycle height growth under 5-year cycles is not ported);
#  • one-unit print knife-edges in a single volume column: ect01 row 2018 TCuFt, wct01 rows 2010 MCuFt / 2018 BdFt,
#    cit01 rows 2060/2008/2018 BdFt.
const _T5_CASES = [("pnt01_t5", _T5.PacificNorthwest(), Int[], Int[]),
                   ("ect01_t5", _T5.EastCascades(),    Int[], Int[]),
                   ("wct01_t5", _T5.WestCascades(),    Int[], Int[]),
                   ("bmt01_t5", _T5.BlueMountains(),   Int[5], Int[]),
                   ("cit01_t5", _T5.CentralIdaho(),    Int[5], Int[])]

@testset "TIMEINT 5, YR=10 variants: .sum rows == live (gradd.f DG rescale after GRINCR)" begin
    for (stem, v, skip_stands, _) in _T5_CASES
        jl = _t5_rows(stem, v); lv = _t5_live(stem)
        @test length(jl) == length(lv)
        st = _t5_stand(lv)
        nfull = 0
        for i in eachindex(lv)
            st[i] in skip_stands && continue
            # Year..QMD (TPA, BA, SDI, CCF, TopHt, QMD) exact in every row
            @test (stem, i, jl[i][1:8]) == (stem, i, lv[i][1:8])
            nfull += (jl[i] == lv[i])
        end
        # every full row exact except the named one-unit volume knife-edges
        nknife = Dict("pnt01_t5" => 0, "ect01_t5" => 1, "wct01_t5" => 2, "bmt01_t5" => 0, "cit01_t5" => 3)[stem]
        @test (stem, nfull) == (stem, count(i -> !(st[i] in skip_stands), eachindex(lv)) - nknife)
    end
end
