# test_west_fortyp_live.jl — FORTYP's California mixed-conifer test and the R5 variants' grinit TLAT/TLONG, vs LIVE.
#
# fortyp.f:1117-1139: in California (ISTATE 6) or Region 5 a Douglas-fir type away from the north coast, sugar
# pine/incense cedar, and ponderosa/Jeffrey or white/red fir types under 80% of the stocking are typed 371
# (California mixed conifer); with no county the north-coast test is the TLAT/TLONG box, whose inputs for a
# keyword stand are grinit.f's defaults (WS 39/120, CA 42/124, NC 42/123, OC 42/124). Goldens: the stock
# tests/FVSws wst01 and tests/FVSnc nct01 keys' live FVS{ws,nc}_g16 .sum rows. Their 4th stand ("FFE TEST") has an
# open FFE fire-outcome residual and is not compared here.

using Test, FVSjl
const _FT = FVSjl
const _FT_FX = joinpath(@__DIR__, "..", "fixtures", "west_fortyp")

function _ft_rows(stem, v)
    dir = mktempdir()
    for ext in ("key", "tre"); cp(joinpath(_FT_FX, "$stem.$ext"), joinpath(dir, "$stem.$ext")); end
    txt = cd(() -> _FT.run_keyfile("$stem.key"; variant = v, output = :sum), dir)
    return [split(l) for l in split(txt, '\n') if occursin(r"^\d{4} ", l) && length(split(l)) >= 20]
end
_ft_stand(rows) = (k = 1; [i == 1 ? 1 : (parse(Int, rows[i][1]) < parse(Int, rows[i-1][1]) ? (k += 1) : k) for i in eachindex(rows)])

@testset "R5 forest type (FORTYP California test) + grinit TLAT/TLONG: stock keys == live" begin
    for (stem, v) in (("wst01", _FT.WestSierra()), ("nct01", _FT.Klamath()))
        jl = _ft_rows(stem, v)
        lv = [split(l) for l in readlines(joinpath(_FT_FX, "$(stem)_live.rows"))]
        @test (stem, length(jl)) == (stem, length(lv))
        st = _ft_stand(lv)
        for i in eachindex(lv)
            (st[i] == 4 || i > length(jl)) && continue          # the FFE TEST stand
            @test (stem, i, jl[i]) == (stem, i, lv[i])
        end
    end
end
