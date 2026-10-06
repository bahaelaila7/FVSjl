# test_kwcov_on.jl — Ontario (ON, metric) keyword-isolation coverage vs the LIVE oracle.
#
# test/fixtures/kwcov_on/on/on_<stand>_[ffe_]<case>.{key,live.sum,<Table>.csv,tables} (+ stands.db): two FIA stands from
# test/fixtures/tiered/on/stands.db (DB-read path) × one isolated keyword — base-model keywords (thinnings, multipliers,
# volume/defect, SETSITE, ESTAB, ...) and the FFE keywords the ON build accepts (every FFE report + *_Metric DBS table
# requested). Goldens = /workspace/.onwork/tiered/FVSon_tiered (test/harness/kwcov_on/make_cases.jl, PROVENANCE.toml). Every
# .sum field and every DBS table cell must match at print precision — no tolerance.
#   KWON_ONLY=regex(on stem)   restrict the run (default: all)
# A case that is not yet bit-exact goes in _KWON_OPEN with the measured diagnosis (reported as @test_broken so an unexpected
# pass surfaces); never silently skipped.
using Test, FVSjl
include(joinpath(@__DIR__, "..", "harness", "kwcov_ffe", "runner.jl"))

const _KWON_OPEN = Dict{String,String}(
    "on_21769405010661_bamax" => "OPEN (measured 2026-10-06, FVSon_tiered): 23609 cells in FVS_StrClass_Metric,FVS_TreeList_East_Metric,sum",
    "on_21769405010661_estab" => "OPEN (measured 2026-10-06, FVSon_tiered): 143 cells in FVS_TreeList_East_Metric",
    "on_21769405010661_eventmon" => "OPEN (measured 2026-10-06, FVSon_tiered): 3313 cells in FVS_Error,FVS_StrClass_Metric,FVS_TreeList_East_Metric,sum",
    "on_21769405010661_fertiliz" => "OPEN (measured 2026-10-06, FVSon_tiered): 21825 cells in FVS_Error,FVS_StrClass_Metric,FVS_TreeList_East_Metric,sum",
    "on_21769405010661_ffe_carbcut" => "OPEN (measured 2026-10-06, FVSon_tiered): 39 cells in FVS_Carbon_Metric,FVS_Down_Wood_Cov_Metric,FVS_Down_Wood_Vol_Metric,FVS_Hrv_Carbon",
    "on_21769405010661_ffe_fmortmlt" => "OPEN (measured 2026-10-06, FVSon_tiered): 542 cells in FVS_Carbon_Metric,FVS_Down_Wood_Cov_Metric,FVS_Down_Wood_Vol_Metric,FVS_Mortality_Metric,FVS_SnagDet_Metric,FVS_SnagSum_Metric,sum",
    "on_21769405010661_ffe_fuelinit" => "OPEN (measured 2026-10-06, FVSon_tiered): 150 cells in FVS_Carbon_Metric,FVS_Consumption_Metric,FVS_Down_Wood_Cov_Metric,FVS_Down_Wood_Vol_Metric",
    "on_21769405010661_ffe_fuelmove" => "OPEN (measured 2026-10-06, FVSon_tiered): 52 cells in FVS_Carbon_Metric,FVS_Consumption_Metric,FVS_Down_Wood_Cov_Metric,FVS_Down_Wood_Vol_Metric",
    "on_21769405010661_ffe_fuelmove_imp" => "OPEN (measured 2026-10-06, FVSon_tiered): 15 cells in FVS_Carbon_Metric,FVS_Consumption_Metric",
    "on_21769405010661_ffe_fuelsoft" => "OPEN (measured 2026-10-06, FVSon_tiered): 116 cells in FVS_Carbon_Metric,FVS_Consumption_Metric,FVS_Down_Wood_Cov_Metric,FVS_Down_Wood_Vol_Metric",
    "on_21769405010661_ffe_fueltret" => "OPEN (measured 2026-10-06, FVSon_tiered): 45 cells in FVS_Carbon_Metric,FVS_Consumption_Metric,FVS_Down_Wood_Cov_Metric,FVS_Down_Wood_Vol_Metric,FVS_Hrv_Carbon",
    "on_21769405010661_ffe_mortclas" => "OPEN (measured 2026-10-06, FVSon_tiered): 62 cells in FVS_Mortality_Metric",
    "on_21769405010661_ffe_pileburn" => "OPEN (measured 2026-10-06, FVSon_tiered): 32 cells in FVS_Carbon_Metric,FVS_Down_Wood_Cov_Metric,FVS_Down_Wood_Vol_Metric,FVS_Hrv_Carbon",
    "on_21769405010661_ffe_pileburn_jack" => "OPEN (measured 2026-10-06, FVSon_tiered): 33 cells in FVS_Carbon_Metric,FVS_Down_Wood_Cov_Metric,FVS_Down_Wood_Vol_Metric,FVS_Hrv_Carbon",
    "on_21769405010661_ffe_simfire_nomort" => "OPEN (measured 2026-10-06, FVSon_tiered): 116 cells in FVS_Carbon_Metric,FVS_SnagDet_Metric,FVS_SnagSum_Metric,sum",
    "on_21769405010661_ffe_snagpbn" => "OPEN (measured 2026-10-06, FVSon_tiered): 62 cells in FVS_Carbon_Metric,FVS_Down_Wood_Cov_Metric,FVS_Down_Wood_Vol_Metric,FVS_SnagDet_Metric,FVS_SnagSum_Metric",
    "on_21769405010661_htgstop" => "OPEN (measured 2026-10-06, FVSon_tiered): 672 cells in FVS_TreeList_East_Metric,sum",
    "on_21769405010661_natural" => "OPEN (measured 2026-10-06, FVSon_tiered): 84 cells in FVS_TreeList_East_Metric",
    "on_21769405010661_thinbba_estab" => "OPEN (measured 2026-10-06, FVSon_tiered): 56 cells in FVS_TreeList_East_Metric",
    "on_21769405010661_thincc" => "OPEN (measured 2026-10-06, FVSon_tiered): 18769 cells in FVS_StrClass_Metric,FVS_TreeList_East_Metric,sum",
    "on_21769405010661_thinpt" => "OPEN (measured 2026-10-06, FVSon_tiered): 25834 cells in FVS_StrClass_Metric,FVS_TreeList_East_Metric,sum",
    "on_21769405010661_thinqfa" => "OPEN (measured 2026-10-06, FVSon_tiered): 15672 cells in FVS_StrClass_Metric,FVS_TreeList_East_Metric",
    "on_21769405010661_topkill" => "OPEN (measured 2026-10-06, FVSon_tiered): 11945 cells in FVS_StrClass_Metric,FVS_TreeList_East_Metric,sum",
    "on_21769405010661_treeszcp" => "OPEN (measured 2026-10-06, FVSon_tiered): 18904 cells in FVS_StrClass_Metric,FVS_TreeList_East_Metric,sum",
    "on_21769405010661_yardloss" => "OPEN (measured 2026-10-06, FVSon_tiered): 9738 cells in FVS_StrClass_Metric,FVS_TreeList_East_Metric,sum",
    "on_64767919010661_ffe_carbcut" => "OPEN (measured 2026-10-06, FVSon_tiered): 31 cells in FVS_Carbon_Metric,FVS_Down_Wood_Vol_Metric,FVS_Hrv_Carbon",
    "on_64767919010661_ffe_flameadj_len" => "OPEN (measured 2026-10-06, FVSon_tiered): 1 cells in FVS_Down_Wood_Cov_Metric",
    "on_64767919010661_ffe_fmortmlt" => "OPEN (measured 2026-10-06, FVSon_tiered): 365 cells in FVS_Carbon_Metric,FVS_Down_Wood_Cov_Metric,FVS_Down_Wood_Vol_Metric,FVS_Mortality_Metric,FVS_SnagDet_Metric,FVS_SnagSum_Metric,sum",
    "on_64767919010661_ffe_fuelinit" => "OPEN (measured 2026-10-06, FVSon_tiered): 147 cells in FVS_Carbon_Metric,FVS_Consumption_Metric,FVS_Down_Wood_Cov_Metric,FVS_Down_Wood_Vol_Metric",
    "on_64767919010661_ffe_fuelmove" => "OPEN (measured 2026-10-06, FVSon_tiered): 58 cells in FVS_Carbon_Metric,FVS_Consumption_Metric,FVS_Down_Wood_Cov_Metric,FVS_Down_Wood_Vol_Metric",
    "on_64767919010661_ffe_fuelmove_imp" => "OPEN (measured 2026-10-06, FVSon_tiered): 27 cells in FVS_Carbon_Metric,FVS_Consumption_Metric,FVS_Down_Wood_Vol_Metric",
    "on_64767919010661_ffe_fuelsoft" => "OPEN (measured 2026-10-06, FVSon_tiered): 106 cells in FVS_Carbon_Metric,FVS_Consumption_Metric,FVS_Down_Wood_Cov_Metric,FVS_Down_Wood_Vol_Metric",
    "on_64767919010661_ffe_fueltret" => "OPEN (measured 2026-10-06, FVSon_tiered): 30 cells in FVS_Carbon_Metric,FVS_Consumption_Metric,FVS_Hrv_Carbon",
    "on_64767919010661_ffe_mortclas" => "OPEN (measured 2026-10-06, FVSon_tiered): 42 cells in FVS_Mortality_Metric",
    "on_64767919010661_ffe_pileburn" => "OPEN (measured 2026-10-06, FVSon_tiered): 31 cells in FVS_Carbon_Metric,FVS_Down_Wood_Vol_Metric,FVS_Hrv_Carbon",
    "on_64767919010661_ffe_pileburn_jack" => "OPEN (measured 2026-10-06, FVSon_tiered): 28 cells in FVS_Carbon_Metric,FVS_Hrv_Carbon",
    "on_64767919010661_ffe_salvage" => "OPEN (measured 2026-10-06, FVSon_tiered): 24 cells in FVS_Carbon_Metric,FVS_Hrv_Carbon",
    "on_64767919010661_ffe_simfire_nomort" => "OPEN (measured 2026-10-06, FVSon_tiered): 6 cells in FVS_Carbon_Metric,FVS_SnagDet_Metric",
    "on_64767919010661_ffe_snagpbn" => "OPEN (measured 2026-10-06, FVSon_tiered): 61 cells in FVS_Carbon_Metric,FVS_Down_Wood_Cov_Metric,FVS_Down_Wood_Vol_Metric,FVS_SnagDet_Metric,FVS_SnagSum_Metric",
)

@testset "ON keyword-isolation vs live" begin
    root = KCASE_ROOTS["on"]
    pat = Regex(get(ENV, "KWON_ONLY", "."))
    for stem in kcase_stems(root, "on")
        occursin(pat, stem) || continue
        @testset "$stem" begin
            ms, crashed, err = run_kcase(root, "on", stem)
            if haskey(_KWON_OPEN, stem)
                @test_broken isempty(ms)
            else
                isempty(ms) || foreach(m -> println("  MISMATCH ", key_str(m), " gold=", m.gold, " jl=", m.got), first(ms, 12))
                @test isempty(ms)
            end
        end
    end
end
