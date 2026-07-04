# test_voleqnum.jl — VOLEQNUM cubic volume-equation override (initre.f:5061) vs live Fortran.
#
# VOLEQNUM overrides the per-species cubic NVEL equation id (VEQNNC). The scenario re-points sugar
# maple (SM) at black oak's equation (841CLKE531); both engines then compute SM's cubic volume with
# that taper. Two checks:
#   1. FIRES — vs the same stand with the VOLEQNUM line removed, total cubic (col 9) changes.
#   2. FORTRAN — TPA/BA and the CUBIC columns (total/merch/sawtimber, 9/10/11) match live Fortran
#      bit-exact (±Float32 noise).
# Board feet (col 12) is ALSO checked: VOLEQNUM sets only the cubic equation (VEQNNB, the board
# equation, stays at SM's default), and FVSjl now keeps a board-equation snapshot (sp_bf_vol_eq) and
# recomputes board feet from a separate board-equation call (BFPFLG=0 path, fvsvol.f:362) — so board
# feet stays on the default equation, matching Fortran (only the documented ±Float32 Scribner noise).

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

        # 2. matches live Fortran on the structural + ALL volume columns (board feet now included).
        ft = _ve_base(sav)
        @test length(jl) == length(ft)
        if length(jl) == length(ft)
            for i in 1:length(jl)
                @test _vecol(jl[i], 3) == _vecol(ft[i], 3)    # TPA — BIT-EXACT
                @test _vecol(jl[i], 4) == _vecol(ft[i], 4)    # BA  — BIT-EXACT
                # Cubic (9/10/11) is BIT-EXACT every cycle except a single print-boundary ULP at 2020 (col 9):
                # the per-acre cubic sum lands within one ULP of the `+0.5` integer-render knife-edge, so the
                # rendered integer rounds to the opposite side of live's. Bound = exactly 1 (one integer step).
                for c in (9, 10, 11)   # total / merch / sawtimber cubic
                    @test abs(_vecol(jl[i], c) - _vecol(ft[i], c)) <= 1
                end
                # Board feet: BIT-EXACT bar a single print-boundary ULP at 2030. Previously carried a systematic
                # −16→−23 residual at the largest cycles; root-caused (BFDUMP per-tree trace) to the BROKEN-TOP
                # board top-kill (BFTOPK) being fit to the WRONG equation's total cubic. VOLEQNUM splits the board
                # equation (SM, VEQNNB) from the cubic (black oak, VEQNNC) → BFPFLG=0; FVS vols.f:391 passes BFMAX
                # (the BOARD call's total, from BFVOL) to BFTOPK, but jl passed v[1] (the cubic call's total). On a
                # broken-top tree that scaled the SM board by black-oak's taper. Fixed to use the board call's vb[1].
                # Now every cycle is bit-exact except the 2030 render knife-edge → bound = exactly 1.
                @test abs(_vecol(jl[i], 12) - _vecol(ft[i], 12)) <= 1
            end
        end
    end
end
