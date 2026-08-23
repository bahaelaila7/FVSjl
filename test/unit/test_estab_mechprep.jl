# MECHPREP / BURNPREP site preparation (esin.f opt 5/4 → estb/esetpr.f + estab.f:246-399) ESTAB-packet keywords
# (2026-08-23). A MECHPREP/BURNPREP keyword schedules mechanical / broadcast-burn site prep on the AUTAL
# DISTURBANCE tally (NTALLY==1): a fraction of the DUPNPT replicate plots is assigned IPPREP∈{1 NONE,2 MECH,
# 3 BURN} by a sample-without-replacement draw off the WK6 site-prep RNG vector the tally already consumes
# (estab.f:333). Each prepped plot's regen uses its IPREP in the advance/excess species mix (ie_espadv/ie_espxcs
# CPRE prep term) — a species-COMPOSITION shift (count near-invariant, PROB1's SPRE term is 0 for this habitat).
#
# LIVE bit-exact vs FVSie_g16 on the under-stocked IE fixture 753189105290487 with THINBTA 2029 (near-clearcut,
# drives the AUTAL disturbance tally). Once AUTOES #143 aligned the establishment seed stream, the disturbance
# tally draws the SAME seed (61997) and hence the SAME WK6 vector as the oracle, so the per-plot IPPREP — and
# the whole .sum — is byte-identical for MECHPREP 100 / BURNPREP 100 / MECHPREP 50 (the WK6-sampled mix):
#   jl h_{base,mech100,burn100,mech50}.sum == oracle, full data rows byte-identical; and jl h_mech100 differs
#   from jl h_base exactly as the oracle does (2069 BA 190→197) — the prep effect is real and bit-for-bit.

using Test
using FVSjl
using FVSjl: ie_esetpr, ie_esetpr_normalize, ie_esetpr_sample, ie_esrann!, IEEstabRNG

@testset "MECHPREP / BURNPREP site preparation (esetpr.f + estab.f sampler)" begin
    @testset "ie_esetpr keyword parse (esetpr.f)" begin
        m100 = ie_esetpr(100f0, nothing)          # MECHPREP 100
        @test m100.pmech == 1f0 && m100.pburn == 0f0 && m100.ialn2 == 1 && m100.ialn3 == 0 && m100.any_kw
        b100 = ie_esetpr(nothing, 100f0)          # BURNPREP 100
        @test b100.pmech == 0f0 && b100.pburn == 1f0 && b100.ialn2 == 0 && b100.ialn3 == 1
        m50 = ie_esetpr(50f0, nothing)            # MECHPREP 50
        @test m50.pmech == 0.5f0 && m50.ialn2 == 1
        bare = ie_esetpr(nothing, nothing)        # no keyword
        @test !bare.any_kw
    end

    @testset "ie_esetpr_normalize → SUMUP bucket weights (estab.f:249-375)" begin
        @test ie_esetpr_normalize(0f0, 1f0, 0f0, 1, 0) == (0f0, 1f0, 0f0)      # MECH100 → all-mech
        @test ie_esetpr_normalize(0f0, 0f0, 1f0, 0, 1) == (0f0, 0f0, 1f0)      # BURN100 → all-burn
        @test ie_esetpr_normalize(0f0, 0.5f0, 0f0, 1, 0) == (0.5f0, 0.5f0, 0f0) # MECH50 → half none/half mech
        # PMECH+PBURN>1 renormalizes (estab.f:353-356)
        s1, s2, s3 = ie_esetpr_normalize(0f0, 0.8f0, 0.8f0, 1, 1)
        @test isapprox(s2, 0.5f0) && isapprox(s3, 0.5f0) && isapprox(s1, 0f0; atol = 1f-6)
    end

    @testset "ie_esetpr_sample per-plot IPPREP (estab.f:382-399)" begin
        wk6 = Float32[i / 51 for i in 1:50]
        @test all(==(2), ie_esetpr_sample((0f0, 1f0, 0f0), wk6, 50, 1))   # MECH100 → every plot MECH
        @test all(==(3), ie_esetpr_sample((0f0, 0f0, 1f0), wk6, 50, 1))   # BURN100 → every plot BURN
    end

    @testset "LIVE IPPREP off the aligned disturbance seed == oracle FVSie_estabdump (h_mech50)" begin
        # The AUTOES-aligned disturbance tally draws seed0 = 61997 (bit-exact vs the live oracle, post-#143).
        # WK6 = the first 50 draws off that seed's odd-adjusted stream (the plot-1 site-prep prefix, estab.f:333).
        seed0 = 61997; s = seed0; iseven(s) && (s += 1)
        rng = IEEstabRNG(Float64(s))
        wk6 = Float32[ie_esrann!(rng) for _ in 1:50]
        sumup = ie_esetpr_normalize(0f0, 0.5f0, 0f0, 1, 0)                 # MECHPREP 50%
        ipprep = ie_esetpr_sample(sumup, wk6, 50, 1)
        # oracle ZIPPREP dump (FVSie_estabdump, h_mech50 disturbance tally, SUMUP 0.5/0.5/0):
        oracle = [1,2,1,2,2,2,2,2,1,1,2,1,1,1,1,2,1,1,2,1,2,2,2,1,2,1,2,1,2,2,
                  1,2,1,2,2,1,2,1,1,2,2,1,2,2,1,1,1,1,1,2]
        @test ipprep == oracle
    end
end
