# test_kwcov_east_ffe.jl — eastern (CS/LS/NE) FFE keyword-isolation coverage vs the LIVE oracle.
#
# test/fixtures/kwcov_ffe/<v>/<v>_<case>.{key,tre,live.sum,<Table>.csv,tables}: one stand × one isolated FFE keyword
# (FLAMEADJ, MOISTURE, FUELMOVE, PILEBURN, FUELTRET, POTF*, SNAG*, SALVAGE, DROUGHT, FUELINIT, FUELFOTO, CARBCALC,
# FMODLIST, FUELMODL, DEFULMOD, FIRECALC, SOILHEAT, FMORTMLT, CANCALC, …) with every FFE report + DBS table requested.
# Goldens come from the stock-flag LIVE oracle (test/harness/kwcov_ffe/make_cases.jl, PROVENANCE.toml). Every .sum field
# and every DBS table cell must match at print precision — no tolerance.
#   KWFFE_VARIANTS=cs,ls   KWFFE_ONLY=regex(on stem)   restrict the run (default: all)
# A case that is not yet bit-exact goes in _KWFFE_OPEN with the measured diagnosis (reported as @test_broken so an
# unexpected pass surfaces); never silently skipped.
using Test, FVSjl
include(joinpath(@__DIR__, "..", "harness", "kwcov_ffe", "runner.jl"))

# OPEN causes (2026-10-06): removed-carbon cells (Total_Removed_Carbon/Biomass_Removed; ported on west-tail 5956ebda, not duplicated
# here), 1-ULP DWD 305to508 cells after FUELMOVE/PILEBURN, and 1-ULP PotFire flame cells (cs/ne firecalc_40/_off/_sav, ls_simfire_nomort).
const _KWFFE_OPEN = Dict{String,String}(
    "cs_carbcut" => "OPEN (measured 2026-10-06, stock-flag live oracle): 24 cells in FVS_Carbon,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_PotFire_East",
    "cs_firecalc_40" => "OPEN (measured 2026-10-06, stock-flag live oracle): 1 cells in FVS_PotFire_East",
    "cs_firecalc_off" => "OPEN (measured 2026-10-06, stock-flag live oracle): 1 cells in FVS_PotFire_East",
    "cs_fuelmove_imp" => "OPEN (measured 2026-10-06, stock-flag live oracle): 2 cells in FVS_Carbon,FVS_Fuels",
    "cs_fueltret" => "OPEN (measured 2026-10-06, stock-flag live oracle): 11 cells in FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Vol,FVS_Fuels",
    "cs_fueltret2" => "OPEN (measured 2026-10-06, stock-flag live oracle): 55 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_SnagDet",
    "cs_pileburn" => "OPEN (measured 2026-10-06, stock-flag live oracle): 14 cells in FVS_Carbon,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels",
    "cs_pileburn_jack" => "OPEN (measured 2026-10-06, stock-flag live oracle): 13 cells in FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels",
    "cs_salvage" => "OPEN (measured 2026-10-06, stock-flag live oracle): 2 cells in FVS_Carbon,FVS_Down_Wood_Cov",
    "cs_salvage_sp" => "OPEN (measured 2026-10-06, stock-flag live oracle): 1 cells in FVS_Down_Wood_Cov",
    "ls_carbcut" => "OPEN (measured 2026-10-06, stock-flag live oracle): 28 cells in FVS_Carbon,FVS_Down_Wood_Vol,FVS_Fuels,FVS_PotFire",
    "ls_fuelmove_imp" => "OPEN (measured 2026-10-06, stock-flag live oracle): 2 cells in FVS_Carbon,FVS_Fuels",
    "ls_fueltret" => "OPEN (measured 2026-10-06, stock-flag live oracle): 12 cells in FVS_Carbon,FVS_Consumption,FVS_Fuels,FVS_PotFire",
    "ls_fueltret2" => "OPEN (measured 2026-10-06, stock-flag live oracle): 13 cells in FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Fuels,FVS_PotFire",
    "ls_pileburn" => "OPEN (measured 2026-10-06, stock-flag live oracle): 27 cells in FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Vol,FVS_Fuels,FVS_PotFire",
    "ls_pileburn_jack" => "OPEN (measured 2026-10-06, stock-flag live oracle): 19 cells in FVS_Carbon,FVS_Down_Wood_Vol,FVS_Fuels,FVS_PotFire",
    "ls_salvage" => "OPEN (measured 2026-10-06, stock-flag live oracle): 1 cells in FVS_Carbon",
    "ls_simfire_nomort" => "OPEN (measured 2026-10-06, stock-flag live oracle): 1 cells in FVS_Down_Wood_Cov",
    "ne_carbcut" => "OPEN (measured 2026-10-06, stock-flag live oracle): 9 cells in FVS_Carbon,FVS_Down_Wood_Vol,FVS_Fuels",
    "ne_firecalc_40" => "OPEN (measured 2026-10-06, stock-flag live oracle): 5 cells in FVS_PotFire",
    "ne_firecalc_sav" => "OPEN (measured 2026-10-06, stock-flag live oracle): 1 cells in FVS_PotFire",
    "ne_fuelmove_imp" => "OPEN (measured 2026-10-06, stock-flag live oracle): 1 cells in FVS_Carbon",
    "ne_fueltret" => "OPEN (measured 2026-10-06, stock-flag live oracle): 12 cells in FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_PotFire",
    "ne_fueltret2" => "OPEN (measured 2026-10-06, stock-flag live oracle): 12 cells in FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_PotFire",
    "ne_pileburn" => "OPEN (measured 2026-10-06, stock-flag live oracle): 11 cells in FVS_Carbon,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels",
    "ne_pileburn_jack" => "OPEN (measured 2026-10-06, stock-flag live oracle): 22 cells in FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_PotFire",
    "ne_salvage" => "OPEN (measured 2026-10-06, stock-flag live oracle): 1 cells in FVS_Carbon",
)

@testset "eastern FFE keyword-isolation vs live" begin
    root = KCASE_ROOTS["ffe"]
    vs = split(get(ENV, "KWFFE_VARIANTS", "cs,ls,ne"), ",")
    pat = Regex(get(ENV, "KWFFE_ONLY", "."))
    for v in vs, stem in kcase_stems(root, v)
        occursin(pat, stem) || continue
        @testset "$stem" begin
            ms, crashed, err = run_kcase(root, v, stem)
            if haskey(_KWFFE_OPEN, stem)
                @test_broken isempty(ms)
            else
                isempty(ms) || foreach(m -> println("  MISMATCH ", key_str(m), " gold=", m.gold, " jl=", m.got), first(ms, 12))
                @test isempty(ms)
            end
        end
    end
end
