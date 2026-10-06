# test_kwcov_ext.jl — extension-STUB parity for CS/LS/NE/ON vs the LIVE oracles.
#
# test/fixtures/kwcov_ext/<v>/<v>_<case>.{key,tre|stands.db,live.sum,<Table>.csv,tables}: the canonical stand with ONE extension keyword block
# (MPB DFTM WSBW COVER BRUST MISTOE DFB RDIN BGCIN ECON CLIMATE ORGANON, + FERTILIZ), opener alone / with an unknown sub-keyword + END. The eastern
# builds link the ex*.f stubs (FVS11 for MPB DFTM COVER BRUST MISTOE DFB RDIN CLIMATE; ON also FERTILIZ and ECON; WSBW is the silent exbudl.f stub), so the
# run must equal the no-extension run and the FVS_Error rows (FVS11 / FVS01 per record) must match live. Goldens: test/harness/kwcov_ext/make_cases.jl
# (PROVENANCE.toml). Every .sum field and every DBS cell, no tolerance.   KWEXT_ONLY=regex(stem)   restrict the run.
using Test, FVSjl
include(joinpath(@__DIR__, "..", "harness", "kwcov_ffe", "runner.jl"))

const _KWEXT_OPEN = Dict{String,String}(
    "cs_brust_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "cs_brust_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "cs_dfb_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "cs_dfb_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "cs_dftm_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "cs_dftm_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "cs_econ_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "cs_fertiliz" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 4 cells in FVS_Summary",
    "cs_fertiliz_none" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 148 cells in FVS_Summary,sum",
    "cs_mpb_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "cs_mpb_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "cs_wsbw_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "cs_wsbw_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ls_brust_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 155 cells in FVS_Error,FVS_Summary,sum",
    "ls_brust_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 155 cells in FVS_Error,FVS_Summary,sum",
    "ls_dfb_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ls_dfb_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ls_dftm_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ls_dftm_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ls_econ_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ls_fertiliz" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 4 cells in FVS_Summary",
    "ls_fertiliz_none" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 162 cells in FVS_Summary,sum",
    "ls_mpb_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ls_mpb_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ls_wsbw_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ls_wsbw_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ne_brust_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ne_brust_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ne_dfb_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ne_dfb_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ne_dftm_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ne_dftm_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ne_econ_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ne_fertiliz_none" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 166 cells in FVS_Summary,sum",
    "ne_mpb_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ne_mpb_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ne_wsbw_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "ne_wsbw_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "on_brust_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "on_brust_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "on_dfb_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "on_dfb_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "on_dftm_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "on_dftm_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "on_econ_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "on_econ_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "on_fertiliz" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 50 cells in FVS_Error,sum",
    "on_fertiliz_none" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "on_mpb_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "on_mpb_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "on_wsbw_end" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
    "on_wsbw_sub" => "OPEN (measured 2026-10-06, stock-flag live oracle, before the stub fixes): 1 cells in FVS_Error",
)

@testset "extension stubs vs live (CS/LS/NE/ON)" begin
    root = KCASE_ROOTS["ext"]
    pat = Regex(get(ENV, "KWEXT_ONLY", "."))
    for v in ("cs", "ls", "ne", "on"), stem in kcase_stems(root, v)
        occursin(pat, stem) || continue
        @testset "$stem" begin
            ms, crashed, err = run_kcase(root, v, stem)
            if haskey(_KWEXT_OPEN, stem)
                @test_broken isempty(ms)
            else
                isempty(ms) || foreach(m -> println("  MISMATCH ", key_str(m), " gold=", m.gold, " jl=", m.got), first(ms, 12))
                @test isempty(ms)
            end
        end
    end
end
