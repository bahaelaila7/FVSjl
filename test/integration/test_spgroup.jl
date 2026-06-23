# test_spgroup.jl — SPGROUP species groups (vbase/initre.f:4726) vs live Fortran.
#
# SPGROUP defines a named species group from a next-record species list (alpha or numeric
# codes). Groups are numbered in definition order and referenced from any keyword's species
# field by the negative index −N (the ISPCC<0 branch). A group reference must be byte-identical
# to naming the member species directly. Each .sum.save is the live-Fortran group run:
#   * spgroup_fixmort — group {SM} + FIXMORT −1 0.9 (DBH≥10): exercises the FIXMORT group branch;
#     verified equal to FIXMORT 22 (SM) in Fortran.
#   * spgroup_fixdg   — group {SM, AB} + FIXDG −1 0.3: exercises a MULTI-species group parse and
#     the FIXDG (apply_fix_scalers!) group branch.
#   * spgroup_thindbh — group {SM} + THINDBH −1: the thin-method species filter (_cut_eligible
#     threads s.control.sp_groups); verified equal to THINDBH 22 in Fortran.

using Test, FVSjl

const _SG_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_sg_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 11 && tryparse(Int, first(split(l))) !== nothing]
_sg_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 11 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_sgcol(r, c) = parse(Float64, r[c])

@testset "SPGROUP species groups vs Fortran" begin
    have(nm) = isfile(joinpath(_SG_DIR, nm * ".key")) && isfile(joinpath(_SG_DIR, nm * ".sum.save"))
    for nm in ("spgroup_fixmort", "spgroup_fixdg", "spgroup_thindbh")
        if !have(nm); @test_skip "$nm scenario not available"; continue; end
        @testset "$nm" begin
            jl = _sg_rows(FVSjl.run_keyfile(joinpath(_SG_DIR, nm * ".key"); faithful = true))
            ft = _sg_base(joinpath(_SG_DIR, nm * ".sum.save"))
            @test length(jl) == length(ft)
            if length(jl) == length(ft)
                for i in 1:length(jl), c in (3, 4, 7, 8)   # TPA / BA / TopHt / QMD
                    @test abs(_sgcol(jl[i], c) - _sgcol(ft[i], c)) <= 1
                end
            end
        end
    end
end
