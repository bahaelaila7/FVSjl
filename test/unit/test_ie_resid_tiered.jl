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

# estab.f:545-549 sets PNN(NCOUNT)=ESA on every plot of a calibrated fresh tally — the cycle-1 ingrowth tally too — and
# estab.f:583 floors each plot's PROB1 at PNN+0.0001, which the ingrowth NSTORE (:587-589) then divides by. jl left PNN=0
# on ingrowth (MEASURED FVSie_g16 3285544010690 2011: point 1 logistic 0.3833 < PNN 0.429568 ⇒ live PROB1 0.429668 and
# NSTORE 3, jl 0.3833 and NSTORE 4 ⇒ 2022 TPA live 1028 / jl 949, 20,187 tiered cells).
@testset "IE AUTOES ingrowth PROB1 floored at PNN=ESA (estab.f:545-589) vs FVSie_g16" begin
    for cn in ("3285544010690", "51032748020004")
        c = _case(cn, "none")
        @test !c.crashed
        @test count(m -> m.file in ("sum", "FVS_TreeList") || m.col in ("Tpa", "BA", "SDI", "CCF"), c.ms) == 0
    end
end
# dense.f:179-188 sums TPROB/TSUMD2 over IND1 species-major with WK5=D*(D*P) (dense.f is byte-identical in the IE build);
# RMSQD=SQRT(TSUMD2/TPROB) is FVS_Summary QMD/ATQMD. jl summed IE in record order with P·D² (MEASURED FVSie_g16
# 3285544010690 2012 QMD live 7.00555182, jl 7.00555038 — 1,280 QMD/ATQMD cells across the IE suite).
@testset "IE QMD (RMSQD) over IND1 with dense.f's WK5 (dense.f:179-188) vs FVSie_g16" begin
    for cn in ("3285544010690", "51032748020004", "374547584489998")
        @test count(m -> m.col in ("QMD", "ATQMD"), _case(cn, "none").ms) == 0
    end
end
end # module
