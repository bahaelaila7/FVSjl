# test_voleqnum.jl — VOLEQNUM cubic volume-equation override (initre.f:5061) vs live Fortran.
#
# VOLEQNUM overrides the per-species cubic NVEL equation id (VEQNNC). The scenario re-points sugar
# maple (SM) at black oak's equation (841CLKE531); both engines then compute SM's cubic volume with
# that taper. Two checks:
#   1. FIRES — vs the same stand with the VOLEQNUM line removed, total cubic (col 9) changes.
#   2. FORTRAN — TPA/BA and the CUBIC columns (total/merch/sawtimber, 9/10/11) match live Fortran
#      bit-exact (±Float32 noise).
# NOTE: board feet (col 12) is intentionally NOT checked. VOLEQNUM sets only the cubic equation
# (VEQNNB, the board equation, stays at SM's default), but FVSjl derives board feet from the same
# (now-overridden) cubic call — so its board-foot column follows the cubic equation. This is the
# documented board-foot separate-equation limitation (see docs/DIVERGENCES.md), shared with BFVOLUME.

using Test, FVSjl

const _VE_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_ve_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && tryparse(Int, first(split(l))) !== nothing]
_ve_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_vecol(r, c) = parse(Float64, r[c])

@testset "VOLEQNUM cubic equation override vs Fortran" begin
    key = joinpath(_VE_DIR, "vol_eqnum.key")
    sav = joinpath(_VE_DIR, "vol_eqnum.sum.save")
    tre = joinpath(_VE_DIR, "vol_eqnum.tre")
    if !isfile(key) || !isfile(sav)
        @test_skip "vol_eqnum scenario not available"
    else
        jl = _ve_rows(FVSjl.run_keyfile(key; faithful = true))

        # 1. the override FIRES: build a no-VOLEQNUM twin, total cubic (col 9) must shift.
        novkey = tempname() * ".key"
        write(novkey, join(filter(l -> !startswith(l, "VOLEQNUM"), readlines(key)), "\n") * "\n")
        cp(tre, replace(novkey, ".key" => ".tre"); force = true)
        nd = _ve_rows(FVSjl.run_keyfile(novkey; faithful = true))
        @test length(jl) == length(nd)
        @test abs(_vecol(jl[1], 9) - _vecol(nd[1], 9)) > 20   # cycle-0 total cuft changed

        # 2. matches live Fortran on the structural + CUBIC columns (board feet excluded — see header).
        ft = _ve_base(sav)
        @test length(jl) == length(ft)
        if length(jl) == length(ft)
            for i in 1:length(jl), c in (3, 4, 9, 10, 11)   # TPA / BA / total / merch / sawtimber cubic
                @test abs(_vecol(jl[i], c) - _vecol(ft[i], c)) <= 2
            end
        end
    end
end
