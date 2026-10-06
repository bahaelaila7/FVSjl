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
end
