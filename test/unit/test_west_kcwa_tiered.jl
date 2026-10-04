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

# dense.f:168-229 sums CCFT·P species-major in IND1 order into RELDSP(ISPC), then RELDT=RELDT+RELDSP(ISPC). The flat
# record-order sum put KT's RELDEN 1 ULP off (3021216010690 cycle 2: DGF CCF2 11904.7773 live / 11904.774 jl), and the
# tripled copies' DG with it (TreeIndex 114 DG 2.35536575 live / 2.3553638 jl).
@testset "KT stand CCF summed per species (dense.f:168-229 RELDSP)" begin
    c = _case("KT", "3021216010690", "none")
    @test count(m -> m.file == "FVS_TreeList" && m.col == "DG" && m.year == "2009", c.ms) == 0
end

# kt/morts.f:201-203 CIOBDS=(2.0*D*G+G*G); SD2SQ=SD2SQ+P*(D*D+CIOBDS) — D²+(2DG+G²), not (D²+2DG)+G².
@testset "KT MORTS SD2SQ association (kt/morts.f:201-203)" begin
    c = _case("KT", "4718785010690", "none")
    @test isempty(_before(c.ms, 2034))          # first divergence 2024 before
end

# kt/morts.f:288-307: X = XMORT inside the MORTMULT DBH window, then the ESTAB best-tree IESTAT guard. jl ran X=1 for every
# record. ktt01 with MORTMULT 1990 all-species ×2.5 on 5-15" and 2010 DF ×0.4: every .sum row equals live FVSkt_clean
# (all 7 rows off without the window).
@testset "KT MORTS MORTMULT window (kt/morts.f:288-292) vs live FVSkt" begin
    dir = joinpath(@__DIR__, "..", "fixtures", "kootenai", "ktt01_mortmult")
    out = mktempdir()
    for f in ("ktmm.key", "ktmm.tre"); cp(joinpath(dir, f), joinpath(out, f)); end
    rows(t) = [l for l in split(t, '\n') if occursin(r"^\s?\d{4} ", l)]
    jl = rows(cd(() -> FVSjl.run_keyfile("ktmm.key"; variant = FVSjl.Kootenai(), output = :sum), out))
    lv = rows(read(joinpath(dir, "ktmm.live.sum"), String))
    @test length(jl) == length(lv)
    for (g, j) in zip(lv, jl); @test rstrip(j) == rstrip(g); end
end
end # module
