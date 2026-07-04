# test_fertiliz.jl — FERTILIZE / FFERT fertilizer growth response (ffert.f) vs live Fortran.
#
# FERTILIZE applies a 200-lb-N fertilizer treatment whose effect boosts each tree's diameter (DDS)
# and height growth for up to 10 years after application (scaled by the application efficacy). The
# scenario fertilizes the base stand in 1995. Two checks:
#   1. FIRES — basal area is visibly higher than the same stand without the FERTILIZE line;
#   2. FORTRAN — TPA/BA/QMD + cubic columns match live Fortran (board feet within Scribner noise,
#      larger here because the growth boost compounds).
# (SN is outside the model's calibrated DF/GF range — Fortran warns but still applies the
# species-agnostic factor, so the bit-exact comparison is valid.)

using Test, FVSjl

const _FZ_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_fz_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && tryparse(Int, first(split(l))) !== nothing]
_fz_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_fzcol(r, c) = parse(Float64, r[c])

@testset "FERTILIZE growth response vs Fortran" begin
    key = joinpath(_FZ_DIR, "fertiliz.key")
    sav = joinpath(_FZ_DIR, "fertiliz.sum.save")
    tre = joinpath(_FZ_DIR, "fertiliz.tre")
    if !isfile(key) || !isfile(sav)
        @test_skip "fertiliz scenario not available"
    else
        jl = _fz_rows(FVSjl.run_keyfile(key; faithful = true))

        # 1. FIRES: the fertilizer boost raises basal area above the no-FERTILIZE twin.
        nokey = tempname() * ".key"
        write(nokey, join(filter(l -> !startswith(l, "FERTILIZ"), readlines(key)), "\n") * "\n")
        cp(tre, replace(nokey, ".key" => ".tre"); force = true)
        nf = _fz_rows(FVSjl.run_keyfile(nokey; faithful = true))
        @test length(jl) == length(nf)
        @test any(_fzcol(jl[i], 4) - _fzcol(nf[i], 4) > 5 for i in 1:length(jl))   # BA boosted

        # 2. matches live Fortran on the structural + cubic columns; board feet within Scribner noise.
        ft = _fz_base(sav)
        @test length(jl) == length(ft)
        if length(jl) == length(ft)
            for i in 1:length(jl)
                for c in (3, 4, 7, 8, 9, 10, 11)   # TPA/BA/TopHt/QMD/total/merch/sawtimber cubic
                    @test abs(_fzcol(jl[i], c) - _fzcol(ft[i], c)) <= 2
                end
                # Board feet: BIT-EXACT every cycle bar a single print-boundary ULP (the per-acre Scribner
                # sum lands within one ULP of the +0.5 integer-render knife-edge). Bound = exactly 1 (one step).
                @test abs(_fzcol(jl[i], 12) - _fzcol(ft[i], 12)) <= 1
            end
        end
    end
end
