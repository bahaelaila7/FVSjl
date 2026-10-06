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

const _KWFFE_OPEN = Dict{String,String}(
    "cs_carbcut" => "OPEN (measured 2026-10-06, stock-flag live oracle): 47 cells in FVS_Carbon,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire_East",
    "cs_firecalc" => "OPEN (measured 2026-10-06, stock-flag live oracle): 510 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire_East,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "cs_firecalc_ml" => "OPEN (measured 2026-10-06, stock-flag live oracle): 510 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire_East,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "cs_fuelfoto" => "OPEN (measured 2026-10-06, stock-flag live oracle): 581 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire_East,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "cs_fuelmove_imp" => "OPEN (measured 2026-10-06, stock-flag live oracle): 2 cells in FVS_Carbon,FVS_Fuels",
    "cs_fueltret" => "OPEN (measured 2026-10-06, stock-flag live oracle): 34 cells in FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon",
    "cs_fueltret2" => "OPEN (measured 2026-10-06, stock-flag live oracle): 78 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_Mortality,FVS_SnagDet",
    "cs_pileburn" => "OPEN (measured 2026-10-06, stock-flag live oracle): 344 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_Mortality,FVS_PotFire_East,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "cs_pileburn_jack" => "OPEN (measured 2026-10-06, stock-flag live oracle): 331 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_Mortality,FVS_PotFire_East,FVS_SnagDet,FVS_SnagSum",
    "cs_salvage" => "OPEN (measured 2026-10-06, stock-flag live oracle): 2 cells in FVS_Carbon,FVS_Down_Wood_Cov",
    "cs_salvage_sp" => "OPEN (measured 2026-10-06, stock-flag live oracle): 1 cells in FVS_Down_Wood_Cov",
    "ls_carbcut" => "OPEN (measured 2026-10-06, stock-flag live oracle): 62 cells in FVS_Carbon,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire",
    "ls_firecalc" => "OPEN (measured 2026-10-06, stock-flag live oracle): 278 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ls_firecalc_ml" => "OPEN (measured 2026-10-06, stock-flag live oracle): 318 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ls_fuelfoto" => "OPEN (measured 2026-10-06, stock-flag live oracle): 346 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary",
    "ls_fuelmove_imp" => "OPEN (measured 2026-10-06, stock-flag live oracle): 2 cells in FVS_Carbon,FVS_Fuels",
    "ls_fueltret" => "OPEN (measured 2026-10-06, stock-flag live oracle): 281 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ls_fueltret2" => "OPEN (measured 2026-10-06, stock-flag live oracle): 48 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire",
    "ls_pileburn" => "OPEN (measured 2026-10-06, stock-flag live oracle): 206 cells in FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire",
    "ls_pileburn_jack" => "OPEN (measured 2026-10-06, stock-flag live oracle): 209 cells in FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire",
    "ls_salvage" => "OPEN (measured 2026-10-06, stock-flag live oracle): 1 cells in FVS_Carbon",
    "ls_simfire_nomort" => "OPEN (measured 2026-10-06, stock-flag live oracle): 1 cells in FVS_Down_Wood_Cov",
    "ne_carbcut" => "OPEN (measured 2026-10-06, stock-flag live oracle): 45 cells in FVS_Carbon,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire",
    "ne_firecalc" => "OPEN (measured 2026-10-06, stock-flag live oracle): 437 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ne_firecalc_ml" => "OPEN (measured 2026-10-06, stock-flag live oracle): 291 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ne_fuelfoto" => "OPEN (measured 2026-10-06, stock-flag live oracle): 370 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ne_fuelmove_imp" => "OPEN (measured 2026-10-06, stock-flag live oracle): 1 cells in FVS_Carbon",
    "ne_fueltret" => "OPEN (measured 2026-10-06, stock-flag live oracle): 48 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire",
    "ne_fueltret2" => "OPEN (measured 2026-10-06, stock-flag live oracle): 381 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ne_pileburn" => "OPEN (measured 2026-10-06, stock-flag live oracle): 204 cells in FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire",
    "ne_pileburn_jack" => "OPEN (measured 2026-10-06, stock-flag live oracle): 205 cells in FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire",
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
