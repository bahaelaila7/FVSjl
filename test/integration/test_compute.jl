# test_compute.jl — COMPUTE event-monitor user variables (vbase/initre.f:1266) vs live Fortran.
#
# COMPUTE defines named variables from arithmetic expressions, re-evaluated each cycle and
# readable by IF/THEN conditions. compute_cycle.key is snt01_alpha (a multi-stand IF/THEN scenario
# whose `(FRAC(CYCLE/3)==0)` condition fires a THINDBH every 3rd cycle) rewritten to compute
# `MYCYC = CYCLE` and drive the condition off MYCYC instead. Two checks:
#   1. EQUIVALENCE — the COMPUTE-fed run is byte-identical to the original direct-`CYCLE` run, so
#      the user variable is computed and consumed correctly (and the firing thins fire identically).
#   2. FORTRAN — the COMPUTE-fed run matches live Fortran on the bit-exact lead stand.

using Test, FVSjl

const _CP_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_cp_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 11 && tryparse(Int, first(split(l))) !== nothing]
_cp_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 11 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]

@testset "COMPUTE event-monitor variables vs Fortran" begin
    ckey = joinpath(_CP_DIR, "compute_cycle.key")
    dkey = joinpath(_CP_DIR, "snt01_alpha.key")
    if !isfile(ckey) || !isfile(dkey)
        @test_skip "compute_cycle scenario not available"
    else
        comp = _cp_rows(FVSjl.run_keyfile(ckey; faithful = true))
        # 1. COMPUTE(MYCYC=CYCLE) feeding the IF condition ≡ the direct CYCLE reference, exactly.
        dir = _cp_rows(FVSjl.run_keyfile(dkey; faithful = true))
        @test length(comp) == length(dir)
        @test all(comp[i] == dir[i] for i in 1:min(length(comp), length(dir)))
        # the thins must actually fire (TPA pulled down at the ÷3 cycles) — not a silent no-op.
        @test any(parse(Float64, r[3]) < 410 for r in comp)
        # 2. lead stand matches live Fortran (the COMPUTE-fed firing thins reproduce the .sum).
        sav = joinpath(_CP_DIR, "compute_cycle.sum.save")
        if isfile(sav)
            ft = _cp_base(sav)
            nlead = 11                                   # first stand's 11 rows
            for i in 1:min(nlead, length(comp), length(ft)), c in (3, 4, 8)
                @test abs(parse(Float64, comp[i][c]) - parse(Float64, ft[i][c])) <= 1
            end
        end
    end
end
