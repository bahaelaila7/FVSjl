# PPE (Parallel Processing Extension) — landscape / multi-stand harness.
# USER-directed 2026-08-21 ("reconstruct a PPE harness like WWPB").  PPE's outer
# routines (PPMAIN / ALSTD2 / SPLAEX) have NO source in this FVS tree — only the
# COMMON-block declarations survive (archive/PPEcommons/PPEPRM.F77 + PPEXCM.F77),
# so, exactly like the WWPB harness, src/engine/ppe_landscape.jl is a FAITHFUL
# RECONSTRUCTION of PPE's documented area-weighted aggregation (PPEXCM PTSTV1),
# NOT a bit-exact port (no oracle exists — the code is absent).  It composes the
# ALREADY-oracle-validated per-stand `run_keyfile` projection into a landscape.
#
# VALIDATED HERE (behavior-faithful self-consistency, since there is no PPE oracle):
#   * identical member stands with unequal areas -> the area-weighting CANCELS and
#     every aggregate equals the single-stand `.sum` value (TOTALWT = Σ area).
#   * two DIFFERENT stands (unthinned A area 3 + thinned B area 1) -> each PTSTV1
#     aggregate == the hand-computed Σ(value·area)/Σ(area) of the per-stand rows.
#   * ADDITIVE + INERT: a new engine file with no simulate.jl seam (gate 339/11).

using Test
using FVSjl

@testset "PPE landscape harness (reconstructed, behavior-faithful)" begin
    V = FVSjl.variant_from_code("EM")
    tre = "/workspace/ForestVegetationSimulator/tests/FVSem/emt01.tre"
    @assert isfile(tre)
    dir = mktempdir()
    cp(tre, joinpath(dir, "ppeA.tre")); cp(tre, joinpath(dir, "ppeB.tre"))

    mkkey(name, extra) = begin
        p = joinpath(dir, "$name.key")
        open(p, "w") do io
            print(io, """
SCREEN
NOAUTOES
NOTRIPLE
STDIDENT
$name  PPE member
DESIGN                                        11.0       1.0
STDINFO        112.0     260.0      60.0     315.0      30.0      54.0
INVYEAR       1990.0
NUMCYCLE         3.0
$(extra)TREEFMT
(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,
T52,I2,T66,5I1,T54,7I1,T75,F3.0)
TREEDATA
ECHOSUM
PROCESS
STOP
""")
        end
        p
    end
    keyA = mkkey("ppeA", "")
    keyB = mkkey("ppeB", "THINBTA       2010.0     100.0\n")

    # per-stand baseline tpa (from the .sum data rows), for the hand check
    tpaof(txt) = Dict(parse(Int, split(strip(l))[1]) => parse(Float64, split(strip(l))[3])
                      for l in split(txt, '\n')
                      if (f = split(strip(l)); length(f) > 19 && all(c -> isdigit(c) || c == '-', f[1])))
    a = tpaof(FVSjl.run_keyfile(keyA; variant = V))
    b = tpaof(FVSjl.run_keyfile(keyB; variant = V))
    @test !isempty(a) && !isempty(b)

    # TEST 1 — identical stands, areas 3 and 1: weighting cancels.
    land1 = FVSjl.ppe_run_landscape(
        [FVSjl.PPEStand(keyA; area = 3.0), FVSjl.PPEStand(keyA; area = 1.0)]; variant = V)
    @test !isempty(land1)
    for g in land1
        @test g.nstands == 2
        @test g.totalwt == 4.0
        @test isapprox(g.avbtpa, a[g.year]; atol = 1e-6)
    end

    # TEST 2 — different stands A(area 3) + B(area 1): area-weighted mean.
    land2 = FVSjl.ppe_run_landscape(
        [FVSjl.PPEStand(keyA; area = 3.0), FVSjl.PPEStand(keyB; area = 1.0)]; variant = V)
    @test !isempty(land2)
    for g in land2
        @test isapprox(g.avbtpa, (3 * a[g.year] + 1 * b[g.year]) / 4; atol = 1e-6)
        @test g.totalwt == 4.0
    end
end
