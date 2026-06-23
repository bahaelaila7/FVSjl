# test_volume_override.jl — VOLUME merch-standard override (volkey.f:9915) vs live Fortran.
#
# VOLUME overrides the per-stand cubic merch standards (DBHMIN/TOPD/STMP/SCFMIND/SCFTOPD/
# SCFSTMP) for a species (0=all, <0=SPGROUP). The scenario raises DBHMIN 4"→8" for all
# species (every present species shares the same other defaults, so ONLY the merch-cubic
# DBH gate changes — no ripple into the shared R8 Clark taper call). Two checks:
#   1. FIRES — vs the same stand with the VOLUME line removed, the merch-cubic column (MCF)
#      drops by >100 cuft/ac once the override is in effect (it gates out 4-8" trees). The
#      1990 inventory + the 1995 cycle keep the defaults (volkey.f ICYC.EQ.0 skip + the
#      activity taking effect from the cycle starting at its 1995 date).
#   2. FORTRAN — TPA/BA and the gated MCF column match live Fortran (bit-exact bar the
#      documented ±Float32 volume/SDI noise). NOTE: VOLUME's merch-top/stump params and the
#      whole BFVOLUME keyword are wired but only partially reproduce Fortran for board feet —
#      FVSjl derives board feet from the shared sawtimber R8 Clark call rather than a separate
#      board-foot taper (see docs/DIVERGENCES.md); the DBHMIN gate is the clean, exact part.

using Test, FVSjl

const _VO_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_vo_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && tryparse(Int, first(split(l))) !== nothing]
_vo_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_vocol(r, c) = parse(Float64, r[c])

@testset "VOLUME merch-standard override vs Fortran" begin
    key = joinpath(_VO_DIR, "volume_override.key")
    sav = joinpath(_VO_DIR, "volume_override.sum.save")
    tre = joinpath(_VO_DIR, "volume_override.tre")
    if !isfile(key) || !isfile(sav)
        @test_skip "volume_override scenario not available"
    else
        jl = _vo_rows(FVSjl.run_keyfile(key; faithful = true))

        # 1. the DBHMIN gate FIRES: build a no-VOLUME twin and compare merch cubic (col 10).
        novkey = tempname() * ".key"
        write(novkey, join(filter(l -> !startswith(l, "VOLUME"), readlines(key)), "\n") * "\n")
        cp(tre, replace(novkey, ".key" => ".tre"); force = true)
        nv = _vo_rows(FVSjl.run_keyfile(novkey; faithful = true))
        @test length(jl) == length(nv)
        # last row: the override has zeroed merch cubic for the 4-8" cohort → a real drop.
        @test _vocol(nv[end], 10) - _vocol(jl[end], 10) > 100
        # the cycle-0 inventory (1990) is untouched by VOLUME (volkey.f ICYC.EQ.0 skip).
        @test _vocol(jl[1], 10) == _vocol(nv[1], 10)

        # 2. matches live Fortran on the structural + gated columns (±Float32 volume noise).
        ft = _vo_base(sav)
        @test length(jl) == length(ft)
        if length(jl) == length(ft)
            for i in 1:length(jl)
                @test abs(_vocol(jl[i], 3) - _vocol(ft[i], 3)) <= 1     # TPA
                @test abs(_vocol(jl[i], 4) - _vocol(ft[i], 4)) <= 1     # BA
                @test abs(_vocol(jl[i], 10) - _vocol(ft[i], 10)) <= 2   # merch cubic (the gated column)
            end
        end
    end
end
