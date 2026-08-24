# =============================================================================
# PPE mode-2 LIVE interstand-beetle coupling — the master-cycle LOCKSTEP BARRIER.
#
# ppe_run_landscape_live! (src/engine/ppe_landscape.jl) steps N stands in lockstep:
# every member is projected by the ordinary per-stand engine (each_stand → setup_growth!
# → write_sum_file → grow_cycle!) but PAUSES at the WWPB seam (simulate.jl, post-MORTS/
# pre-GRADD, un-tripled) via the `wwpb_barrier` hook. When all stands reach that seam for
# a master cycle, ONE landscape bmdrv_multi! dispersal runs across the placed neighbors
# (SPLALO geometry), and each stand's beetle kill is handed back (bmkill! → WK2 → t.tpa)
# before it resumes. The barrier is a cooperative Julia-Task + Channel rendezvous; the
# BMRANN RNG is advanced ONLY inside the coordinator (bmdrv_multi!/bmistd!) in fixed stand
# order, so the run is deterministic. `wwpb_barrier === nothing` on ordinary runs keeps the
# single-stand path byte-identical (guard test_multicycle stays 339/11).
#
# The dispersal KERNELS (bmsdit!/bmdrv_multi!/bmatct_multi!/bmkill!) are ALREADY bit-exact
# vs pristine Fortran (test_wwpb_bmatct_multi.jl). This suite proves the NEW ORCHESTRATION
# PLUMBING (the lockstep barrier + state extraction at the seam + kill handback) adds ZERO
# divergence, oracle-free but rigorously:
#
#   (A) EQUIVALENCE — the LIVE in-flight-coupled run == a PREMADE-decisions run bit-exact:
#       * ppe_standalone_decisions recomputes every per-cycle dispersal decision OUTSIDE the
#         engine, feeding the ALREADY-VALIDATED standalone bmdrv_multi! cascade the exact
#         per-cycle barrier treelists the LIVE run recorded → bit-exact to the LIVE kills.
#       * ppe_landscape_replay re-projects each stand INDEPENDENTLY (no Tasks, no channel, no
#         coupling) injecting those fixed decisions at the seam → byte-identical .sum to LIVE.
#       Two paths, differing in MECHANISM (in-flight lockstep vs pre-computed injection),
#       sharing the bit-exact kernels ⇒ equality proves the barrier/extraction/handback.
#   (B) COLLAPSE — MXSTND=1 live-coupled == the single-stand DISPERSE path (run_keyfile /
#       wwpb_apply!) BYTE-IDENTICAL (bmdrv_multi! degenerates to bmatct_single!).
#   (C) USE-CASE — an outbreak seeded in ONE stand disperses to its adjacent neighbor over
#       master cycles (neighbor beetle mortality > 0, and > the neighbor projected in
#       isolation; monotone with the seed), while a stand beyond the dispersal radius is
#       unaffected (zero beetle mortality, identical to isolation).
# =============================================================================
using Test
using FVSjl
const F = FVSjl

@testset "PPE mode-2 LIVE interstand coupling (lockstep barrier)" begin
    V = F.variant_from_code("EM")
    tre = "/workspace/ForestVegetationSimulator/tests/FVSem/emt01.tre"
    @assert isfile(tre)
    dir = mktempdir()

    # a member keyfile; TREEDATA reads <stem>.tre, so give each member its own copy.
    mkkey(name, disperse) = begin
        cp(tre, joinpath(dir, "$name.tre"); force = true)
        p = joinpath(dir, "$name.key")
        open(p, "w") do io
            print(io, """
SCREEN
NOAUTOES
NOTRIPLE
STDIDENT
$name
DESIGN                                        11.0       1.0
STDINFO        112.0     260.0      60.0     315.0      30.0      54.0
INVYEAR       1990.0
NUMCYCLE         3.0
BMIN
DISPERSE$disperse
END
TREEFMT
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
    # DISPERSE start=1 dur=3 class=6; seeded stand tpa=40, neighbors tpa=0 (receive only).
    seeded = "         1.0       3.0       6.0      40.0"
    noseed = "         1.0       3.0       6.0       0.0"
    kA = mkkey("standA", seeded)   # outbreak seeded here
    kB = mkkey("standB", noseed)   # adjacent neighbor (500 m — within URMAX=15 mi)
    kC = mkkey("standC", noseed)   # far stand (40 km — beyond the dispersal radius)

    mk(k, x, tpa) = F.PPELiveMember(k; area = 100.0, xloc = x, yloc = 0.0,
                                    seed_class = 6, seed_tpa = tpa)
    members = [mk(kA, 0.0, 40.0), mk(kB, 500.0, 0.0), mk(kC, 40000.0, 0.0)]

    res = F.ppe_run_landscape_live!(members; variant = V, iyr1 = 1, iyr2 = 3)
    @test res.ncyc == 3
    @test length(res.sums) == 3
    @test all(!isempty, res.sums)

    hxv(v) = UInt32[reinterpret(UInt32, x) for x in v]        # Float32-exact comparison

    # ---- (A) EQUIVALENCE — standalone-decisions == in-flight kills (bit-exact) ----------
    sd = F.ppe_standalone_decisions(res; variant = V)
    @test all(length(sd[i]) == res.ncyc for i in 1:3)
    @test all(all(hxv(sd[i][c]) == hxv(res.kills[i][c]) for c in 1:res.ncyc) for i in 1:3)

    # ---- (A) EQUIVALENCE — independent replay of the fixed decisions == LIVE (byte-id) --
    rep = F.ppe_landscape_replay(members, sd; variant = V)
    @test rep == res.sums                                    # every stand's .sum byte-identical
    # and replay off the LIVE-recorded kills agrees too (the handback is the recorded WK2).
    @test F.ppe_landscape_replay(members, res.kills; variant = V) == res.sums

    # ---- (B) COLLAPSE — MXSTND=1 live == single-stand run_keyfile DISPERSE (byte-id) ----
    D = "01-01-2026"; Tm = "12:00:00"
    ref1 = F.run_keyfile(kA; variant = V, date = D, time = Tm)
    res1 = F.ppe_run_landscape_live!([mk(kA, 0.0, 40.0)]; variant = V, iyr1 = 1, iyr2 = 3,
                                     date = D, time = Tm)
    @test res1.sums[1] == ref1

    # ---- (C) USE-CASE — spread to the neighbor, none to the far stand -------------------
    # beetle-attributable mortality = the WK2 excess over the FVS density/background baseline.
    bexc(r, i) = sum(sum(max.(r.kills[i][c] .- r.fvsmort[i][c], 0f0)) for c in 1:r.ncyc)

    excA = bexc(res, 1); excB = bexc(res, 2); excC = bexc(res, 3)
    resB0 = F.ppe_run_landscape_live!([mk(kB, 0.0, 0.0)]; variant = V, iyr1 = 1, iyr2 = 3)
    excB_iso = bexc(resB0, 1)

    @test excA > 0f0                       # the seeded stand has an outbreak
    @test excB_iso == 0f0                  # the neighbor, in ISOLATION + unseeded, has no beetle kill
    @test excB > 0f0                       # ...but in the landscape it DOES — dispersed from A
    @test excB > excB_iso                  # the spread is the coupling (not a local artifact)
    @test excC == 0f0                      # the far stand (beyond URMAX) is UNAFFECTED

    # monotone: a stronger seed at A raises the neighbor's dispersed mortality.
    res2 = F.ppe_run_landscape_live!([mk(kA, 0.0, 80.0), mk(kB, 500.0, 0.0),
                                      mk(kC, 40000.0, 0.0)]; variant = V, iyr1 = 1, iyr2 = 3)
    @test bexc(res2, 2) > excB             # neighbor mortality grows with the seed
    @test bexc(res2, 3) == 0f0             # far stand still unaffected
end
