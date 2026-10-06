# test_west_tail2.jl — western-tail round 2, checked against the LIVE tiered goldens (test/fixtures/tiered/<v>).
module WestTail2Test
using Test
using FVSjl
include(joinpath(@__DIR__, "..", "harness", "tiered", "tiered_runner.jl"))

function _wt2_case(v, cn, r)
    d = mktempdir()
    txt, db, crashed, err = run_case(v, cn, r; dir = d)
    return (crashed = crashed, ms = crashed ? Mismatch[] : compare_case(v, cn, r, txt, db))
end
_wt2_in(ms, files) = [m for m in ms if m.file in files]

# CA CWHR canopy cover: FMCFMD's `CALL CWHR(ISP,DBH,HT,FMICR,FMPROB,CRWDTH,...)` (ca/fmcfmd.f:230) reads the common-block
# CRWDTH(I) — the very stored value FMCBA used (R5CRWD on IFOR<=5, load/end-of-cycle dims, ca/cwcalc.f:385) — not a
# ca_cwcalc recomputed from the current tree dims. Before: 15320267010497 salvage PotFire 14 cells off (live fuel model
# 3P→FM6 vs jl 2M), 22959652010497 / 22960323010497 / 23676707010900 likewise (CA 234 → 23 tiered cells).
@testset "CA/NC/WS CWHR reads the FMCBA CRWDTH (ca/fmcfmd.f:230)" begin
    for (cn, r, files) in (("15320267010497", "salvage", ["FVS_PotFire"]),
                           ("15320267010497", "simfire", ["FVS_BurnReport", "FVS_PotFire"]),
                           ("22959652010497", "salvage", ["FVS_PotFire"]),
                           ("22959652010497", "simfire", ["FVS_BurnReport", "FVS_PotFire"]),
                           ("23676707010900", "simfire", ["FVS_BurnReport", "FVS_PotFire", "FVS_Carbon"]))
        c = _wt2_case("CA", cn, r)
        @test !c.crashed
        @test isempty(_wt2_in(c.ms, files))
    end
end

# CI habitat-moisture decay multiplier: ci/fmcba.f:407-445 scales every DKR(I,J) by DCYMLT (MAPDRY(ICINDX): 0 -> 0.66,
# 1 -> 1, 2 -> 1.33) at the first FFE year. Before: 5385337010690 / 5388215010690 litter 0.508 -> 0.254 vs live 0.3404 after
# year 1 (DKR(10) 0.33), Forest_Floor/Forest_Down_Dead_Wood/fuel model weights off (CI 179 tiered cells).
@testset "CI DCYMLT decay-rate habitat multiplier (ci/fmcba.f:407-445)" begin
    for (cn, r) in (("5385337010690", "salvage"), ("5385337010690", "simfire"),
                    ("5388215010690", "salvage"), ("5388215010690", "simfire"))
        c = _wt2_case("CI", cn, r)
        @test !c.crashed
        @test isempty(_wt2_in(c.ms, ["FVS_Carbon", "FVS_PotFire", "FVS_BurnReport"]))
    end
    @test FVSjl.ci_adjusted_dkr(5)[10, 1] == 0.5f0 * 0.66f0     # MAPDRY(5) = 0 dry
    @test FVSjl.ci_adjusted_dkr(15)[10, 1] == 0.5f0             # MAPDRY(15) = 1 mesic
    @test FVSjl.ci_adjusted_dkr(57)[11, 4] == 0.002f0 * 1.33f0  # MAPDRY(57) = 2 moist
end

# KT/IE FMKOD = KODTYP = MTYPE(ITYPE) (kt/habtyp.f:121, ie/habtyp.f:147; fmsdit.f:70) feeds NIFMHAB's MAPDRY. jl passed the KT
# KKTYPE index (41) instead of 310, so IDRY stayed 0 (model 8) where live blends 9 (PSME/SYAL, dry shrubby).
@testset "KT/IE FMKOD = MTYPE(ITYPE) for the MAPDRY fuel-model class (kt/habtyp.f:121)" begin
    for r in ("salvage", "simfire")
        c = _wt2_case("KT", "4718785010690", r)
        @test !c.crashed
        @test isempty(_wt2_in(c.ms, ["FVS_PotFire", "FVS_BurnReport"]))
    end
end
end
