# test_west_tail2.jl — western-tail round 2, checked against the LIVE tiered goldens (test/fixtures/tiered/<v>).
module WestTail2Test
using Test
using FVSjl
using SQLite, DBInterface
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

# A cycle whose only due activities are FFE ones (SALVAGE/SALVSP/PILEBURN/FUELMOVE) is not a cut: FVS compacts the zero-PROB
# records at CUTS entry only AFTER the start-of-cycle DENSE (grincr.f), and FMSALV runs at cuts.f:386. jl compacted early inside
# cuts!, re-ordering RDPSRT ties: AK 10706662010497 salvage Acc 63 vs live 62, QMD/BA drifting from 2017 (the same run without
# SALVAGE matched; live salvage .sum == live none .sum).
@testset "SALVAGE-only cycle does not compact the tree list early (cuts.f:255 vs grincr.f DENSE)" begin
    c = _wt2_case("AK", "10706662010497", "salvage")
    @test !c.crashed
    @test isempty(_wt2_in(c.ms, ["FVS_Summary", "sum", "FVS_Carbon", "FVS_Fuels", "FVS_SnagSum"]))
end

# CR MISTPINF on a stand with no infected tree: mistoe.f:193 clears DMFLAG each MISTOE call, misinf.f:182 sets it, and MISPRT
# (misprt.f:265-393) prints the cycle's zero-infection row (one FVS_DM_Stnd_Sum row, an empty FVS_DM_Spp_Sum table).
# jl's CR path never set DMFLAG and the report gate excluded CR: both tables were missing (5 CR stands, 10 PRESENCE cells).
@testset "CR MISTPINF zero-infection DM report (mistoe.f:193/267, misinf.f:182, misprt.f)" begin
    for cn in ("2463020010690", "3026069010690")
        c = _wt2_case("CR", cn, "mistletoe")
        @test !c.crashed
        @test isempty(_wt2_in(c.ms, ["FVS_DM_Stnd_Sum", "FVS_DM_Spp_Sum"]))
    end
end

# CR DM spread: the overstory-intensification branch's second CALL RANN(XNUM) overwrites XNUM (mistoe.f:415), and the DMR
# decrease test PMINUS.GT.XNUM (mistoe.f:462) reads that second draw. jl's cr_mistoe! kept the first draw (ie_mistoe! already
# had this): FIA 46279527020004 cycle 2 record 447 PMINUS 0.00328 vs draw 0.0029 (first) / 0.746 (second): DMR 3 live, 2 jl.
@testset "CR mistoe DMR decrease reads the overwritten XNUM (mistoe.f:415/462)" begin
    c = _wt2_case("CR", "46279527020004", "mistletoe")
    @test !c.crashed
    @test isempty(c.ms)
end

# CR/UT/TT FMCFMD LPPDOM reads FMTBA — the FMCOM array FMCBA filled — not a per-call recomputation (cr/fmcfmd.f:805-820). After a
# SIMFIRE kill FMPROB is 0 (CTBA rebuilt empty) but FMTBA still holds the pre-fire stand, so a non-ponderosa mixed-conifer stand is
# "not LPPDOM" (IFMST 0 -> model 2); jl's recomputed zeros made every species a tie -> LPPDOM -> model 9 (CR 31673493010690 2013).
@testset "CR LPPDOM uses the stale FMCBA FMTBA after the fire kill (fmcfmd.f:805-820)" begin
    c = _wt2_case("CR", "31673493010690", "simfire")
    @test !c.crashed
    @test isempty(_wt2_in(c.ms, ["FVS_PotFire", "FVS_BurnReport", "FVS_Carbon"]))
end

# AK ESTAB's BAAA/PTPAA/PRDA are the arrays of the gradd.f:192 DENSE, which precedes ESNUTR: ESUCKR's fire-kill sprout records
# (appended before ESTAB) are not in them. jl recomputed them over the list WITH the sprouts: AK 10709344010497 simfire ESTOCK
# XTPA 210.06 (432 rec) vs live 200.54 (396 rec, DENSE) => post-fire regen TPA +3.5%, Tpa/BA/TCuFt/PotFire/Carbon all drifting
# (live DENSE/EST trace: private FVSak relink of dense.f/estab.f, /workspace/.postswap/wt2/ak_relink).
@testset "AK ESTAB density stats exclude ESUCKR sprouts (gradd.f:192 DENSE vs ESNUTR)" begin
    c = _wt2_case("AK", "10709344010497", "simfire")
    @test !c.crashed
    @test isempty([m for m in c.ms if !(m.col in ("QMD", "ATQMD"))])      # QMD/ATQMD: 1-2 ULP (3e-6 rel), open
end

# crown.f:118/143 IF((ITRN.LE.0).AND.(IREC2.LT.MAXTP1)) GO TO 74 — a stand whose inventory is ALL dead records still dubs their
# missing crowns (DUBSCR BACHLO draws, DO 79). jl's CA/WS crown_ratio_update! returned at n == 0: PctCr 0 on every dead record
# and the RNG short of those draws (WS 7689156010901: 52+52+3 tiered cells, dead trees PctCr 29/10/53 live).
# ca_alldead: FIA CA stand 15320267010497 with every record HISTORY 6 and CR blank, goldens = live FVSca_g16 TreeList.
@testset "all-dead inventory still dubs dead crowns (WS/CA crown.f GO TO 74)" begin
    c = _wt2_case("WS", "7689156010901", "plant_cal")
    @test !c.crashed
    @test isempty(c.ms)
    fx = joinpath(@__DIR__, "..", "fixtures", "west_tail")
    dir = mktempdir(); cp(joinpath(fx, "ca_alldead.stands.db"), joinpath(dir, "stands.db"))
    key = [l == "out.db" ? joinpath(dir, "out.db") : l == "stands.db" ? joinpath(dir, "stands.db") : l
           for l in readlines(joinpath(fx, "ca_alldead.key"))]
    write(joinpath(dir, "s.key"), join(key, '\n'))
    FVSjl.run_keyfile(joinpath(dir, "s.key"); variant = FVSjl.CentralCalifornia())
    gl = readlines(joinpath(fx, "ca_alldead.FVS_TreeList.csv")); hdr = split(gl[1], ',')
    gold = sort([(split(l, ',')[findfirst(==("TreeId"), hdr)], split(l, ',')[findfirst(==("SpeciesFVS"), hdr)],
                  parse(Int, split(l, ',')[findfirst(==("PctCr"), hdr)])) for l in gl[2:end]])
    db = SQLite.DB(joinpath(dir, "out.db"))
    jl = sort([(string(r.TreeId), string(r.SpeciesFVS), Int(r.PctCr)) for r in DBInterface.execute(db, "SELECT TreeId,SpeciesFVS,PctCr FROM FVS_TreeList")])
    SQLite.close(db)
    @test gold == jl
end
end
