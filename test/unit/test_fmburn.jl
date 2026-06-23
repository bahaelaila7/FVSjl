# Tests for the FFE fire event driver (F5b, FMBURN/FMEFF) — the composed fire that
# applies the .sum-affecting kill to tree TPA. End-to-end bit-exactness vs Fortran needs
# the SIMFIRE keyword wiring + the live projection; these pin the driver's behavior:
# it kills trees size-dependently, honors %-burned and the mortality switch, is
# deterministic, and is a no-op when FFE is off.
using FVSjl: fmburn!, FireResult, StandState, init_blockdata!, init_merch_standards!,
             FireState, Southern

function _burn_stand()
    s = StandState(Southern()); init_blockdata!(s, s.variant); init_merch_standards!(s)
    s.plot.forest_type = Int32(520)
    s.plot.latitude = 35f0; s.plot.longitude = -80f0; s.plot.elevation = 10f0
    t = s.trees; t.n = 5
    for (i, (sp, d, h, ic)) in enumerate(((65,14f0,72f0,40),(33,8f0,55f0,45),
                                          (22,4f0,18f0,50),(65,2f0,8f0,55),(33,1f0,6f0,60)))
        t.species[i]=Int32(sp); t.dbh[i]=d; t.height[i]=h; t.tpa[i]=30f0; t.crown_pct[i]=Int32(ic)
    end
    s.fire = FireState(); s.fire.active = true
    s
end

@testset "FFE fire driver (FMBURN/FMEFF)" begin
    s = _burn_stand(); tpa0 = copy(s.trees.tpa[1:5])
    res = fmburn!(s; wind = 20f0, fmois = 1, psburn = 100f0, burnseas = 3)

    @testset "fire runs and kills trees" begin
        @test res isa FireResult
        @test res.flame > 0f0 && res.byram > 0f0 && res.scorch > 0f0
        @test res.killed > 0f0
        # every record loses some TPA (whole stand burned), none goes negative
        @test all(s.trees.tpa[i] < tpa0[i] for i in 1:5)
        @test all(s.trees.tpa[i] >= 0f0 for i in 1:5)
        # killed total equals the summed TPA reduction
        @test res.killed ≈ sum(tpa0[i] - s.trees.tpa[i] for i in 1:5)
    end

    @testset "size-dependent mortality" begin
        # surviving fraction: the big thick-barked oak survives far better than the saplings
        surv(i) = s.trees.tpa[i] / tpa0[i]
        @test surv(1) > surv(3)     # d14 oak vs d4 sapling
        @test surv(1) > surv(5)     # d14 oak vs d1 stem (≤1" → fully killed)
        @test s.trees.tpa[5] ≈ 0f0  # the ≤1" hardwood is killed outright
    end

    @testset "controls" begin
        # %-burned = 0 ⇒ no record is in the burned portion ⇒ no kill
        s0 = _burn_stand(); r0 = fmburn!(s0; psburn = 0f0)
        @test r0.killed == 0f0 && all(s0.trees.tpa[i] == 30f0 for i in 1:5)
        # mortality switch off ⇒ no kill (behavior still computed)
        sm = _burn_stand(); rm = fmburn!(sm; mortcode = 0)
        @test rm.killed == 0f0 && rm.flame > 0f0
        # FFE inactive ⇒ no-op
        si = _burn_stand(); si.fire.active = false
        @test fmburn!(si).killed == 0f0
    end

    @testset "determinism" begin
        a = _burn_stand(); ra = fmburn!(a; wind = 20f0, fmois = 1, burnseas = 3)
        b = _burn_stand(); rb = fmburn!(b; wind = 20f0, fmois = 1, burnseas = 3)
        @test ra.killed == rb.killed && a.trees.tpa[1:5] == b.trees.tpa[1:5]
        # a wetter fuel model kills fewer trees than a very-dry one
        d = _burn_stand(); rd = fmburn!(d; fmois = 1)
        w = _burn_stand(); rw = fmburn!(w; fmois = 4)
        @test rw.killed <= rd.killed
    end
end
