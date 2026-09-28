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

end # module
