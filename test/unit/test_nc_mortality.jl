# NC (Klamath) density-mortality distribution = NWCMRT (nc/nwcmrt.f), regression 2026-08-24.
#
# The western full-population FIA sweep flagged NC with a ~7.5% needs_dig rate whose worst column was TPA:
# on MATURE over-max-SDI stands the port UNDER-killed TPA ~3× (both-sides-traced CN 504553467126144:
# live 345→49, old jl 345→175; CN 25028875010900: live 101→18, old jl 101→52). Root cause: the port applied
# a UNIFORM per-tree self-thinning rate, but NC's morts.f distributes the density (or background) mortality
# TOTAL across records by NWCMRT — a BA-percentile/shade-tolerance geometric progression that CONCENTRATES
# kills on suppressed (low-percentile = small) trees — and then runs the QMD-convergence loop. Uniform matches
# the TOTAL first-pass kill but not the distribution; feeding the BAMAX residual-BA cap a uniform (BA-heavy)
# kill removes BA too fast per tree, leaving far too much TPA. Signed ≥10-stand tally: 0 over / 4 eq / 11 under
# (one-directional bias) BEFORE → 2 over / 5 eq / 8 under (mixed, few-%) AFTER, magnitudes 3×→~few%. The
# residual is downstream of the separately-cornered NC diameter growth (#207): the pre-growth Reineke mean
# (dia0) is bit-exact vs FVSnc_clean (23.3376 vs 23.338), only the grown dr10 drifts (growth, not mortality).
#
# This is a self-contained unit test of the two invariants the fix restores (no DB / no oracle):
#   (1) `_varmrt_efftr!(::Klamath)` = PEFF·VARADJ·0.1 with PEFF DECREASING in the BA percentile (nc/nwcmrt.f:99-102);
#   (2) `_varmrt!` distributes a fixed TOKILL CONCENTRATED on the low-percentile (small) trees, NOT uniformly.
using Test
using FVSjl

@testset "NC (Klamath) NWCMRT mortality distribution" begin
    # NC VARADJ (nc/nwcmrt.f:49-51): DF(sp3)=0.65, RF(sp9)=0.50, BO(sp7)=1.00.
    peff(pct) = clamp(0.84525f0 - 0.01074f0 * pct + 0.0000002f0 * Float32(pct)^3, 0.01f0, 1.00f0)

    @testset "efftr = PEFF·VARADJ·0.1, decreasing in BA percentile (nwcmrt.f:99-102)" begin
        # 3 DF (sp3) records: small (PCT 10), mid (PCT 50), dominant (PCT 100).
        t = FVSjl.TreeList(3)
        for (i, pct) in enumerate((10f0, 50f0, 100f0))
            t.species[i] = Int32(3); t.dbh[i] = Float32(2 + i); t.tpa[i] = 100f0
            t.crown_ratio[i] = pct                       # crown_ratio holds FVS's PCT (stand_pct!)
        end
        t.n = 3
        efftr = zeros(Float32, 3)
        pass1 = FVSjl._varmrt_efftr!(efftr, nothing, FVSjl.Klamath(), t, 3)
        # exact NWCMRT efficiency for DF (VARADJ 0.65)
        @test efftr[1] ≈ peff(10f0)  * 0.65f0 * 0.1f0
        @test efftr[2] ≈ peff(50f0)  * 0.65f0 * 0.1f0
        @test efftr[3] ≈ peff(100f0) * 0.65f0 * 0.1f0
        # concentrated on the suppressed (small) tree; dominant hits the 0.01 PEFF floor
        @test efftr[1] > efftr[2] > efftr[3]
        @test efftr[3] ≈ 0.01f0 * 0.65f0 * 0.1f0        # PEFF floored at PCT=100
        @test pass1 ≈ sum(t.tpa[i] * efftr[i] for i in 1:3)
        # species tolerance scales it: an intolerant BO (VARADJ 1.0) is killed harder than a tolerant RF (0.50)
        t.species[1] = Int32(7); efb = zeros(Float32, 3); FVSjl._varmrt_efftr!(efb, nothing, FVSjl.Klamath(), t, 3)
        t.species[1] = Int32(9); efr = zeros(Float32, 3); FVSjl._varmrt_efftr!(efr, nothing, FVSjl.Klamath(), t, 3)
        @test efb[1] > efr[1]
        @test efb[1] ≈ peff(10f0) * 1.00f0 * 0.1f0
        @test efr[1] ≈ peff(10f0) * 0.50f0 * 0.1f0
    end

    @testset "_varmrt! distributes TOKILL concentrated on small trees, total preserved" begin
        # equal-TPA trees spanning the percentile range; NWCMRT must kill more of the low-percentile records.
        t = FVSjl.TreeList(4)
        for (i, pct) in enumerate((8f0, 30f0, 60f0, 100f0))
            t.species[i] = Int32(3); t.dbh[i] = Float32(i); t.tpa[i] = 100f0; t.crown_ratio[i] = pct
        end
        t.n = 4
        killed = zeros(Float32, 4); efftr = zeros(Float32, 4); temwk2 = zeros(Float32, 4)
        tokill = 120f0
        sumkil = FVSjl._varmrt!(killed, efftr, temwk2, nothing, FVSjl.Klamath(), t, 4, tokill)
        @test sumkil ≈ tokill rtol = 1e-3               # NWCMRT iterates to the requested total
        @test sum(killed) ≈ tokill rtol = 1e-3
        @test all(killed[i] <= t.tpa[i] + 1f-4 for i in 1:4)   # never kill more than a record holds
        # STRICTLY monotone-decreasing kill down the size gradient — the whole point of the fix
        # (a uniform rate would give killed[1]==killed[4]).
        @test killed[1] > killed[2] > killed[3] > killed[4]
        @test killed[1] > 3 * killed[4]                 # strongly concentrated on the suppressed tree
    end
end
