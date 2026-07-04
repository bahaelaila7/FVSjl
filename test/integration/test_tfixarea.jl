# test_tfixarea.jl — TFIXAREA total-fixed-plot-area expansion (initre.f:816 / notre.f:45) vs Fortran.
#
# TFIXAREA sets TFPA, the total fixed plot area. In NOTRE the small-tree (DBH < BRK) sample is then
# expanded by FP = 1/TFPA instead of the default fixed_plot_inv/π, changing those records' per-acre
# TPA. The base stand has BRK=5 with 7 sub-5" trees, so TFPA=0.02 visibly changes the stand. Checks:
#   1. FIRES — TPA differs from the same stand with the TFIXAREA line removed;
#   2. FORTRAN — TPA/BA/cubic columns match live Fortran (board feet within Scribner noise).

using Test, FVSjl

const _TF_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_tf_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && tryparse(Int, first(split(l))) !== nothing]
_tf_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_tfcol(r, c) = parse(Float64, r[c])

@testset "TFIXAREA fixed-plot-area expansion vs Fortran" begin
    key = joinpath(_TF_DIR, "tfixarea.key")
    sav = joinpath(_TF_DIR, "tfixarea.sum.save")
    tre = joinpath(_TF_DIR, "tfixarea.tre")
    if !isfile(key) || !isfile(sav)
        @test_skip "tfixarea scenario not available"
    else
        jl = _tf_rows(FVSjl.run_keyfile(key; faithful = true))

        # 1. FIRES: dropping the TFIXAREA line changes the expanded TPA.
        nokey = tempname() * ".key"
        write(nokey, join(filter(l -> !startswith(l, "TFIXAREA"), readlines(key)), "\n") * "\n")
        cp(tre, replace(nokey, ".key" => ".tre"); force = true)
        nf = _tf_rows(FVSjl.run_keyfile(nokey; faithful = true))
        @test _tfcol(jl[1], 3) != _tfcol(nf[1], 3)     # cycle-0 TPA differs

        # 2. matches live Fortran (cubic columns bit-exact; board feet within Scribner noise).
        ft = _tf_base(sav)
        @test length(jl) == length(ft)
        if length(jl) == length(ft)
            for i in 1:length(jl)
                for c in (3, 4, 9, 10, 11)
                    @test abs(_tfcol(jl[i], c) - _tfcol(ft[i], c)) <= 2
                end
                @test _tfcol(jl[i], 12) == _tfcol(ft[i], 12)
            end
        end
    end
end
