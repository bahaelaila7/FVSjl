# test_west_kcwa_tiered.jl — KT / WS / CA / AK fixes vs the LIVE western tiered goldens (test/fixtures/tiered/<v>: stands.db +
# <stand>_<regime>.key + the live .sum / DBS tables of FVSkt_clean / FVS<v>_g16) and live keyfile goldens. Each testset names
# the Fortran it follows and the case it was measured on; a testset only asserts what its mechanism moved.
module WestKcwaTieredTest
using Test
using FVSjl
include(joinpath(@__DIR__, "..", "harness", "tiered", "tiered_runner.jl"))

function _case(v, cn, r)
    d = mktempdir()
    txt, db, crashed, _ = run_case(v, cn, r; dir = d)
    return (txt = txt, db = db, crashed = crashed, ms = compare_case(v, cn, r, txt, db))
end
_yr(m) = tryparse(Int, m.year)
_before(ms, y) = [m for m in ms if (k = _yr(m); k !== nothing && k < y)]

# kt/htgf.f, dgf.f, crown.f, morts.f, regent.f, ktfctr.f: REAL*4 EXP/ALOG/** are glibc expf/logf/powf and 0.985**K is
# libgcc __powisf2 (doctrine §4). Julia's own Float32 exp/log/^ put the 1999 HTG of four records 1-2 ULP off live
# (3021216010690: HtG 2.92662287 live / 2.92662311 jl), compounding through every later cycle.
@testset "KT growth kernels in glibc expf/logf/powf (kt/htgf.f:107-114, dgf.f, crown.f, morts.f)" begin
    c = _case("KT", "3021216010690", "none")
    @test !c.crashed
    @test isempty(_before(c.ms, 2000))
end

# kt/regent.f:560-562: a sub-4.5 ft record gets DBH(K)=0.1+DIAM·0.01+HK·0.001 with DG(K)=0 in its OWN slot — the tripled
# copies too (kt/dgdriv.f:253/261 put the pre-REGENT D in their slots; TRIPLE does not copy DBH). jl gave the copies the
# central DBH plus an equivalent increment (TreeList DG 7.9E-4 vs live 0) and a direct-set central's DBH for the others.
@testset "KT REGENT sub-4.5 ft tripled copies keep their own DBH(K) (kt/regent.f:560-562)" begin
    c = _case("KT", "4718785010690", "none")
    @test !c.crashed
    @test count(m -> m.file == "FVS_TreeList" && m.col == "DG" && m.year == "2014", c.ms) == 0
    c = _case("KT", "1627650292290487", "none")
    @test isempty(_before(c.ms, 2050))          # 8325 cells from 2030 before
end
end # module
