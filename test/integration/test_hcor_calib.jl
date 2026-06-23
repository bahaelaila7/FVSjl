# test_hcor_calib.jl — small-tree height-growth calibration (REGENT HCOR_init) vs Fortran.
#
# The REGENT small-tree height calibration regression (regent.f:411-516) computes
# HCOR_init = ln( Σ(measured HTG·P) / Σ(predicted EDH·P) ) for species with ≥ NCALHT(5)
# small (dbh<5) trees that have measured height growth. No standard SN scenario carries
# such data, so this is a purpose-built stand: s29's 30 records plus 6 small SM (sp 22)
# trees with measured HTG, all on an existing plot (hcor_smalltree.key/.tre).
#
# Species 22 ends up height-calibrated but NOT diameter-calibrated (LDGCAL=false). The
# resulting HCOR_init = -0.893823 is BIT-EXACT to live Fortran (verified by an instrumented
# regent.f HCORF dump on this exact stand). Switching FVSjl from the old HCOR_init≡0
# assumption to this regression moves cycle-1 TPA from 4986 to 5564 against the Fortran
# 5511 — i.e. a −525 error becomes +53, confirming the port's direction and magnitude.
# (A small residual remains in the broader small-tree height-growth blend; tracked
# separately — it is NOT in this regression, whose value matches to all printed digits.)

using Test, FVSjl

const _HCOR_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
const _HCOR_KEY = joinpath(_HCOR_DIR, "hcor_smalltree.key")
const _HCOR_TRE = joinpath(_HCOR_DIR, "hcor_smalltree.tre")

@testset "REGENT small-tree HCOR calibration vs Fortran" begin
    if !isfile(_HCOR_KEY) || !isfile(_HCOR_TRE)
        @test_skip "hcor_smalltree scenario not available"
    else
        s, _ = FVSjl.initialize(_HCOR_KEY)
        FVSjl.notre!(s); FVSjl.setup_growth!(s)
        # sp 22 (SM) — HCOR_init from the regent regression, bit-exact to Fortran.
        @test isapprox(s.calib.htg_cor_init[22], -0.893823f0; atol = 1f-4)
        # height-calibrated but not diameter-calibrated → constant-HCOR (non-attenuating) path
        @test s.calib.ldgcal[22] == false
        # only SM clears the NCALHT(5) threshold; every other species stays HCOR_init = 0
        @test count(!=(0f0), s.calib.htg_cor_init) == 1
    end
end
