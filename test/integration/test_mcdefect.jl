# test_mcdefect.jl — MCDEFECT cubic-volume defect (sdefet.f + FVSsn vols.f:294-332) vs Fortran.
#
# MCDEFECT sets a per-species CUBIC defect curve (CFDEFT) — defect fractions at DBH 5/10/15/20/25"
# (the 25" value extends flat to 40"). In the SN volume path the pulpwood/topwood part of merch
# cubic (MCFV−SCFV) is reduced by ICDF% = NINT(ALGSLP(DBH, CFDEFT)·100) clamped to [0,99]; the
# sawtimber cubic is untouched. The override is undated here, so (like sdefet.f's "date not
# defined" branch) it applies from cycle 0. Two checks:
#   1. FIRES — vs the same stand with the MCDEFECT line removed, merch cubic (col 10) drops by
#      >100 cuft/ac in every cycle, the 1990 inventory included (immediate, not deferred).
#   2. FORTRAN — TPA/BA and the defect-reduced merch-cubic column match live Fortran bit-exact
#      (bar the documented ±Float32 total-cuft/board-foot/SDI noise). NOTE: per-tree DEFECT input
#      and the CFLA0/CFLA1 log-linear form model (default no-op) are deferred (see DIVERGENCES.md).

using Test, FVSjl

const _MD_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_md_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && tryparse(Int, first(split(l))) !== nothing]
_md_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_mdcol(r, c) = parse(Float64, r[c])

@testset "MCDEFECT cubic defect vs Fortran" begin
    key = joinpath(_MD_DIR, "mcdefect_override.key")
    sav = joinpath(_MD_DIR, "mcdefect_override.sum.save")
    tre = joinpath(_MD_DIR, "mcdefect_override.tre")
    if !isfile(key) || !isfile(sav)
        @test_skip "mcdefect_override scenario not available"
    else
        jl = _md_rows(FVSjl.run_keyfile(key; faithful = true))

        # 1. the defect FIRES from cycle 0: build a no-MCDEFECT twin, merch cubic must drop everywhere.
        novkey = tempname() * ".key"
        write(novkey, join(filter(l -> !startswith(l, "MCDEFECT"), readlines(key)), "\n") * "\n")
        cp(tre, replace(novkey, ".key" => ".tre"); force = true)
        nd = _md_rows(FVSjl.run_keyfile(novkey; faithful = true))
        @test length(jl) == length(nd)
        @test all(_mdcol(nd[i], 10) - _mdcol(jl[i], 10) > 100 for i in 1:length(jl))  # incl. row 1 (1990)

        # 2. matches live Fortran on the structural + defect-reduced merch-cubic columns.
        ft = _md_base(sav)
        @test length(jl) == length(ft)
        if length(jl) == length(ft)
            for i in 1:length(jl)
                @test abs(_mdcol(jl[i], 3) - _mdcol(ft[i], 3)) <= 1      # TPA
                @test abs(_mdcol(jl[i], 4) - _mdcol(ft[i], 4)) <= 1      # BA
                @test abs(_mdcol(jl[i], 10) - _mdcol(ft[i], 10)) <= 2    # merch cubic (defect target)
            end
        end
    end
end
