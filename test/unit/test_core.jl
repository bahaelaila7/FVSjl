# C0 unit tests — foundation: state construction, RNG bit-exactness, variant, units.

using Test
using FVSjl
using FVSjl: Control, PlotData, FVSRng, rann!, esrann!, ranseed!, bachlo,
             nint, MAXTRE, MAXSP, MAXCYC

@testset "parameters" begin
    @test MAXTRE == 3000
    @test MAXSP  == 90
    @test MAXCYC == 40
    @test nint(2.5) == 3      # Fortran NINT rounds half away from zero (not banker's)
    @test nint(-2.5) == -3
    @test nint(2.4) == 2
end

@testset "StandState construction & faithful defaults" begin
    s = StandState(Southern())
    @test variant_code(s.variant) == "SN"
    @test s.control.variant_code == "SN"
    @test s.control.faithful == true
    @test s.trees.n == 0
    @test length(s.trees.dbh) == MAXTRE
    @test length(s.trees.height) == MAXTRE
    # faithful BLOCK DATA / initializer defaults
    @test s.control.unit_stand == 6
    @test s.control.unit_summary == 8
    @test s.control.auto_max == 60f0
    @test s.control.auto_min == 45f0
    @test s.control.cut_eff == 0.98f0
    @test s.control.mort_period == 5f0
    @test s.plot.baf == 40f0
    @test s.plot.fixed_plot_inv == 300f0
    @test s.plot.min_dbh_var_plot == 5f0
    @test s.plot.gross_space == 1f0
    @test isapprox(s.plot.pi, 3.14159265f0; atol=1f-6)
    # extensions inactive until their keywords fire
    @test s.fire === nothing
    @test s.econ === nothing
    # faithful=false path
    s2 = StandState(Southern(); faithful=false)
    @test s2.control.faithful == false
end

@testset "RNG bit-exactness vs Oracle A" begin
    r = FVSRng()
    # establishment stream seeded 55329 (ESBLKD)
    @test r.es0 == 55329.0
    # main stream defaults to 0 → first draw is 0 (16807*0 mod m)
    @test rann!(r) == 0f0
    # GOLDEN: seed 12345 then two draws — verified identical to FVSjulia.RANN
    ranseed!(r, true, 12345f0)
    @test rann!(r) == 0.09661653f0
    @test rann!(r) == 0.8339946f0
    # determinism: same seed → same sequence
    r2 = FVSRng(); ranseed!(r2, true, 12345f0)
    @test rann!(r2) == 0.09661653f0
    @test rann!(r2) == 0.8339946f0
    # establishment stream is independent of the main stream
    r3 = FVSRng()
    e1 = esrann!(r3); e2 = esrann!(r3)
    @test 0f0 < e1 < 1f0 && 0f0 < e2 < 1f0
    @test e1 != e2
    # bachlo respects stdev<=0 shortcut
    @test bachlo(FVSRng(), 7.0, 0.0) == 7f0
end

@testset "thread independence (no shared state)" begin
    # two states must not share any mutable array
    a = StandState(Southern()); b = StandState(Southern())
    @test a.trees.dbh !== b.trees.dbh
    @test a.rng !== b.rng
    a.trees.dbh[1] = 9.9f0
    @test b.trees.dbh[1] == 0f0
end
