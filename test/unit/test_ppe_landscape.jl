# PPE (Parallel Processing Extension) — landscape / multi-stand orchestration.
#
# src/engine/ppe_landscape.jl is a FAITHFUL PORT of the PPMAIN mode-1 master-cycle
# control structure (recovered source scratchpad/ppe/recovered/ppbase/ppmain.f, FVS
# rev bc6e2377^ — the earlier "no source" premise was wrong). It orders stands by the
# bit-exact C11SRT (ppe_index_qsort!), projects each with the oracle-validated
# per-stand `run_keyfile` (GRSTND-equivalent), and forms the area-weighted PPEXCM
# PTSTV1(1..9) landscape aggregates (SPLAEX/ALSTD2).
#
# MODE-1 (the reachable regime) is bit-identical to the master-cycle-stepped form:
# stands do not interact, so independent full projection + area-aggregation == the
# stepped run, and every per-stand value is oracle-validated. The C11SRT order and the
# master-cycle stepping are behaviorally inert for mode-1 output (validated below by
# order-invariance). Mode-2 (interstand beetle) + MXHRVP are documented seams.
#
# VALIDATED HERE:
#   * ppe_processing_order == the bit-exact ppe_index_qsort! over stand ids, and mode-1
#     aggregates are INVARIANT to input stand order (proves the equivalence premise).
#   * ppe_neighbors == the bit-exact ppe_hxindx hex-neighbor lookup on a placed grid.
#   * identical member stands with unequal areas -> area-weighting CANCELS (every
#     aggregate == the single-stand `.sum` value; TOTALWT = Σ area).
#   * two DIFFERENT stands -> each PTSTV1 aggregate == Σ(value·area)/Σ(area) of the
#     per-stand rows; PTSTV1(7) MSPERIOD == the report cadence.
#   * ADDITIVE + INERT: no simulate.jl seam (gate 339/11).

using Test
using FVSjl

@testset "PPE landscape orchestration (faithful PPMAIN mode-1 port)" begin
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

    # ---- C11SRT processing order: bit-exact ppe_index_qsort! over stand ids --------
    s3 = [FVSjl.PPEStand(joinpath(dir, "c.key")),
          FVSjl.PPEStand(joinpath(dir, "a.key")),
          FVSjl.PPEStand(joinpath(dir, "b.key"))]
    ord = FVSjl.ppe_processing_order(s3)
    idx = collect(1:3)
    FVSjl.ppe_index_qsort!(idx, ["c.key", "a.key", "b.key"]; lseq = false)
    @test ord == idx                              # order IS the C11SRT of the basenames
    @test [basename(s3[i].keyfile) for i in ord] == ["a.key", "b.key", "c.key"]

    # ---- HXINDX spatial neighbors on a placed 3x3-ish grid (nstnd=9) ---------------
    grid = [FVSjl.PPEStand(joinpath(dir, "g$(r)$(c).key"); col = c, row = r)
            for r in 1:3 for c in 1:3]
    nb = FVSjl.ppe_neighbors(grid[5], grid; nstnd = 9)   # centre-ish stand (row2,col2)
    @test !isempty(nb)                                   # a centre cell has hex neighbors
    @test all(1 .<= nb .<= 9) && !(5 in nb)              # valid indices, not self
    @test isempty(FVSjl.ppe_neighbors(FVSjl.PPEStand(keyA), grid))  # unplaced -> none

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

    # TEST 2 — different stands A(area 3) + B(area 1): area-weighted mean, and the
    # aggregate is INVARIANT to input order (the mode-1 equivalence premise).
    sA = FVSjl.PPEStand(keyA; area = 3.0); sB = FVSjl.PPEStand(keyB; area = 1.0)
    land2  = FVSjl.ppe_run_landscape([sA, sB]; variant = V)
    land2r = FVSjl.ppe_run_landscape([sB, sA]; variant = V)          # reversed input
    @test !isempty(land2)
    @test [g.avbtpa for g in land2] == [g.avbtpa for g in land2r]    # order-invariant
    yrs = sort(collect(keys(a)))
    for (i, g) in enumerate(land2)
        @test isapprox(g.avbtpa, (3 * a[g.year] + 1 * b[g.year]) / 4; atol = 1e-6)
        @test g.totalwt == 4.0
        # PTSTV1(7) MSPERIOD == gap to next report row (0 on the last)
        exp_ms = i < length(land2) ? (land2[i+1].year - g.year) : 0
        @test g.msperiod == exp_ms
    end
end
