# test_west_shared_tiered.jl — shared western fixes vs the LIVE western tiered goldens (test/fixtures/tiered/<v>: stands.db +
# <stand>_<regime>.key + the live .sum / DBS tables written by FVS<v>_g16 / FVScr_clean). Each testset names the Fortran
# it follows and the measured case; cells are compared at a 1e-5 relative tolerance (Float32 drift from other open
# residuals is not what these tests pin — the named mechanism is).
module WestSharedTieredTest
using Test
using FVSjl
include(joinpath(@__DIR__, "..", "harness", "tiered", "tiered_runner.jl"))

function _case(v, cn, r)
    d = mktempdir()
    txt, db, crashed, _ = run_case(v, cn, r; dir = d)
    return (txt = txt, db = db, crashed = crashed, ms = compare_case(v, cn, r, txt, db))
end
_rel(m) = (a = tryparse(Float64, m.gold); b = tryparse(Float64, m.got);
           (a === nothing || b === nothing) ? Inf : abs(a - b) / max(abs(a), 1e-9))
_material(ms) = [m for m in ms if _rel(m) > 1e-5]

# tt/fmcrow.f and ut/fmcrow.f call FMCROWW(SPIE=ISPMAP(SPIW)) for the conifers and FMCROWE only for the hardwoods
# (TT CASE(6,15,16,18), UT CASE(6,18:20,22)); both fmcroww.f/fmcrowe.f are byte-identical to CR's. jl had no TT/UT route,
# so every species took the eastern FMCROWE with `ls_spi` (MEASURED FVStt_g16 2780339010690 2000, the inventory year:
# Aboveground_Total_Live 13.17 live / 14.30 jl, Canopy_Density 0.102 / 0.061 — 196 inventory-year FFE cells on TT).
@testset "TT/UT crown weights route FMCROWW by ISPMAP (tt|ut/fmcrow.f) vs live at the inventory year" begin
    for (v, cn, y) in (("TT", "2780339010690", "2000"), ("TT", "2839796010690", "1999"),
                       ("UT", "42642675010690", "2010"), ("UT", "2390216010690", "2002"))
        c = _case(v, cn, "salvage")
        @test !c.crashed
        @test isempty(_material([m for m in c.ms if m.year == y && m.file in ("FVS_Carbon", "FVS_PotFire")]))
    end
end

# fmcba.f bare-stand cover type: with no basal area the first FFE year takes COVTYP = COVINI(ITYPE) (the habitat's seral
# cover) else the variant's "NO VALID HABITAT" default (tt/ut 7 LP, wc/pn COVINI6→16, ec →3); the CA-FFE top-2 variants
# then load the one cover at weight 1 (nc:359-360 COVCA(1)=COVTYP, COVCAWT(1)=1). jl used DF 3 for TT/UT, 16 for every
# WC/PN habitat, and left NC's COVCAWT at 0 (MEASURED: TT 3333677010690 1992 live COVTYP 7 — DDW 4.05 / jl 2.40; NC
# 450603388489998 2016 DDW 1.40 / jl 0.0; WC 374601116489998 COVTYP 19; PN 26379272010900 19; EC 450507010497 10).
@testset "bare-stand FFE cover type COVINI(ITYPE) + top-2 COVCA (fmcba.f) vs live at the first FFE year" begin
    for (v, cn, y) in (("TT", "3333677010690", "1992"), ("UT", "471768556489998", "2016"), ("NC", "450603388489998", "2016"),
                       ("WC", "374601116489998", "2015"), ("PN", "26379272010900", "2002"), ("EC", "450507010497", "2009"))
        c = _case(v, cn, "salvage")
        @test !c.crashed
        @test isempty(_material([m for m in c.ms if m.year == y && m.file == "FVS_Carbon"]))
    end
end

# ci/habtyp.f:44-90: PVREF4 crosswalks (PV_CODE, PV_REF_CODE) → KODTYP; an unknown pair (FVS34/33/32) or a KODTYP outside
# 10..999 (FVS14) keeps ci/grinit.f ICINDX=21 / ITYPE=4 (habitat 260). jl used PV_CODE mod 1000 and ICINDX 1 (MEASURED
# FVSci_g16: 3369538010690 9999999/491 and 12276084010690 45101/494 "MAPPED TO 260"; jl 999 / 101 ⇒ DGHAB/ITYPE wrong,
# e.g. 3369538010690 2006 HtG 11.67 live / 8.64 jl, BdFt 1246 / 1201).
@testset "CI habitat via PVREF4 + ICINDX 21 default (ci/habtyp.f, ci/pvref4.f) vs FVSci_g16" begin
    for cn in ("3369538010690", "12276084010690")
        c = _case("CI", cn, "none")
        @test !c.crashed
        @test count(m -> m.file == "sum", c.ms) == 0
        @test count(m -> m.file == "FVS_Error", c.ms) == 0          # FVS33 + FVS32 (unknown pair)
    end
    c = _case("CI", "3159852010690", "salvage")                        # bare, no PV code: FVS14 + CIPVG(ICINDX 21)
    @test count(m -> m.file in ("FVS_Error", "FVS_PotFire", "FVS_Carbon"), c.ms) == 0
end

# fortyp.f:1117-1140 California mixed-conifer test (ISTATE 6 or region 5): DF off the north coast (needs ICNTY, read
# from the FIA COUNTY column — dbsstandin.f:392-395), SP/IC, PP/JP with PP < 80%, WF/RF with true fir < 80% → 371
# (MEASURED: NC 44777283020004 2011 live 371 / jl 221, SO 7690240010901 371/261, CA 23742358010900 371/201, WS
# 850400255290487 371/222).
@testset "FORTYP California mixed conifer 371 (fortyp.f:1117-1140) vs live" begin
    for (v, cn) in (("NC", "44777283020004"), ("SO", "7690240010901"), ("CA", "23742358010900"), ("WS", "850400255290487"))
        c = _case(v, cn, "none")
        @test !c.crashed
        @test count(m -> m.col == "ForTyp", c.ms) == 0
    end
end

end # module
