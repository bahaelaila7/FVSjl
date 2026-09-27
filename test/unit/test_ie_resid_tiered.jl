# test_ie_resid_tiered.jl — IE per-record fixes vs the LIVE FVSie_g16 tiered goldens (test/fixtures/tiered/ie: stands.db +
# <stand>_<regime>.key + the live .sum / DBS tables). Each testset names the Fortran it follows and the measured case.
module IEResidTieredTest
using Test
using FVSjl
include(joinpath(@__DIR__, "..", "harness", "tiered", "tiered_runner.jl"))

_case(cn, r) = (d = mktempdir(); txt, db, crashed, _ = run_case("IE", cn, r; dir = d);
                (txt = txt, db = db, crashed = crashed, ms = compare_case("IE", cn, r, txt, db)))

# fvsvol.f:90-96 hands VOLINIT IREGN=KODFOR/100; mrules.f then gives region 6 (the Colville, KODFOR 621) COR='N' (raw
# Scribner) and OPT 23, region 1 COR='Y' and OPT 22. jl hard-wired region 1 for IE (MEASURED FVSie_g16 374547584489998
# 2015: RC 11.8" BdFt live 76 / jl 80, LP 12.6" MCuFt 36.9 / 36.6 — 32 inventory-year volume cells).
@testset "IE FW2 volumes follow the KODFOR region's NVEL merch rules (fvsvol.f/mrules.f) vs FVSie_g16" begin
    c = _case("374547584489998", "none")
    @test !c.crashed
    @test count(m -> m.col in ("MCuFt", "BdFt", "TCuFt"), c.ms) == 0
end

end # module
