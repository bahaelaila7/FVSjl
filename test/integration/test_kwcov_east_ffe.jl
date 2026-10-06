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
    "cs_carbcalc" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 10 cells in FVS_Carbon",
    "cs_carbcut" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 47 cells in FVS_Carbon,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire_East",
    "cs_defulmod" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 31 cells in FVS_PotFire_East",
    "cs_firecalc" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 510 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire_East,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "cs_firecalc_ml" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 510 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire_East,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "cs_flameadj_len" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 414 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire_East,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "cs_flameadj_mult" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 2 cells in FVS_BurnReport",
    "cs_fmortmlt" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 8 cells in FVS_PotFire_East",
    "cs_fmortmlt_dbh" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 481 cells in FVS_Carbon,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire_East,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "cs_fmortmlt_early" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 495 cells in FVS_Carbon,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire_East,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "cs_fuelfoto" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 581 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire_East,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "cs_fuelmodl" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 31 cells in FVS_PotFire_East",
    "cs_fuelmove_imp" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 329 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire_East,FVS_SnagDet,FVS_SnagSum,FVS_Summary",
    "cs_fuelpool" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 327 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire_East,FVS_SnagDet,FVS_SnagSum,FVS_Summary",
    "cs_fueltret" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 34 cells in FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon",
    "cs_fueltret2" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 78 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_Mortality,FVS_SnagDet",
    "cs_mortclas" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 42 cells in FVS_Mortality",
    "cs_pileburn" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 344 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_Mortality,FVS_PotFire_East,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "cs_pileburn_jack" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 331 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_Mortality,FVS_PotFire_East,FVS_SnagDet,FVS_SnagSum",
    "cs_salvage" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 228 cells in FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire_East,FVS_SnagDet,FVS_SnagSum",
    "cs_salvage_sp" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 1 cells in FVS_Down_Wood_Cov",
    "cs_simfire_nomort" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 173 cells in FVS_Carbon,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_PotFire_East,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "cs_snagfall_sp" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 1 cells in FVS_SnagDet",
    "cs_snagpsft" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 765 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire_East,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "cs_statfuel" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 500 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire_East,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ls_cancalc" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 8 cells in FVS_PotFire",
    "ls_canfprof" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 792 cells in FVS_CanProfile",
    "ls_carbcalc" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 10 cells in FVS_Carbon",
    "ls_carbcut" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 62 cells in FVS_Carbon,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire",
    "ls_defulmod" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 29 cells in FVS_PotFire",
    "ls_firecalc" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 278 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ls_firecalc_ml" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 318 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ls_flameadj_len" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 305 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ls_fmortmlt" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 8 cells in FVS_PotFire",
    "ls_fmortmlt_dbh" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 486 cells in FVS_Carbon,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ls_fmortmlt_early" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 486 cells in FVS_Carbon,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ls_fueldcay" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 390 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ls_fuelfoto" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 346 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary",
    "ls_fuelmodl" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 27 cells in FVS_PotFire",
    "ls_fuelmove_imp" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 200 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary",
    "ls_fuelmult" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 390 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ls_fueltret" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 281 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ls_fueltret2" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 48 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire",
    "ls_mortclas" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 46 cells in FVS_Mortality",
    "ls_pileburn" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 206 cells in FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire",
    "ls_pileburn_jack" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 209 cells in FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire",
    "ls_salvage" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 244 cells in FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire,FVS_SnagDet,FVS_SnagSum",
    "ls_simfire_nomort" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 60 cells in FVS_Carbon,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_PotFire",
    "ls_snagpsft" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 672 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary",
    "ls_statfuel" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 257 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ne_cancalc" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 10 cells in FVS_PotFire",
    "ne_canfprof" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 841 cells in FVS_CanProfile",
    "ne_carbcalc" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 10 cells in FVS_Carbon",
    "ne_carbcut" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 45 cells in FVS_Carbon,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire",
    "ne_defulmod" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 24 cells in FVS_PotFire",
    "ne_firecalc" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 437 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ne_firecalc_ml" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 291 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ne_flameadj_len" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 251 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ne_fmortmlt" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 8 cells in FVS_PotFire",
    "ne_fmortmlt_dbh" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 465 cells in FVS_Carbon,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ne_fmortmlt_early" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 464 cells in FVS_Carbon,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ne_fueldcay" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 312 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary",
    "ne_fuelfoto" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 370 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ne_fuelmodl" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 24 cells in FVS_PotFire",
    "ne_fuelmove_imp" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 174 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary",
    "ne_fuelmult" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 312 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary",
    "ne_fueltret" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 48 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire",
    "ne_fueltret2" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 381 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary,sum",
    "ne_mortclas" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 28 cells in FVS_Mortality",
    "ne_pileburn" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 204 cells in FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire",
    "ne_pileburn_jack" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 205 cells in FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_PotFire",
    "ne_salvage" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 482 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Hrv_Carbon,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary",
    "ne_snagbrk" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 377 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary",
    "ne_snagbrk_sp" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 259 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary",
    "ne_snagpsft" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 674 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary",
    "ne_statfuel" => "OPEN (first measured 2026-10-06, stock-flag live oracle): 201 cells in FVS_BurnReport,FVS_Carbon,FVS_Consumption,FVS_Down_Wood_Cov,FVS_Down_Wood_Vol,FVS_Fuels,FVS_Mortality,FVS_PotFire,FVS_SnagDet,FVS_SnagSum,FVS_Summary",
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
