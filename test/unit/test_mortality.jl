# C4 — mortality model (MORTS), snt01 cycle 1. RE-GROUNDED 2026-07-02 (forget Oracle-A):
# these pin jl's DETERMINISTIC internal density-model intermediates (SDIMAX, Zeide dia0,
# Pretzsch t, isolated mortality! TPA loss). They are NOT .sum columns; their LIVE grounding
# is that snt01.sum cyc1 = 507/103/202 BIT-EXACT vs live FVSsn (verified via sn_oracle.sh),
# which these intermediates produce. Tolerances below are DISPLAY-ROUNDING of the pinned
# Float32 value (e.g. SDIMAX 348.43875 → shown 348.44), not live-divergence slack.

using Test
using FVSjl
using FVSjl: notre!, setup_growth!, diameter_growth!, mortality!, stand_sdimax,
             stand_tpa, _pretzsch_tn10, PRETZSCH_SDIK, bark_ratio

const KEYM = "/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key"

@testset "SDIMAX (SDICAL) — exact internal pin (live via .sum)" begin
    s, _ = initialize(KEYM); notre!(s)
    @test stand_sdimax(s) == 348.43875f0                  # BA-weighted SDIDEF — exact Float32 (was atol=0.1 display-round)
end

@testset "density mortality (Pretzsch) — Zeide dia0 + activation" begin
    s, _ = initialize(KEYM); notre!(s)
    setup_growth!(s); diameter_growth!(s, s.variant)
    t = s.trees; fint = 5f0
    @test s.control.zeide_sdi                              # SN uses Zeide/Reineke SDI
    # Zeide/Reineke start diameter (no growth) = exact Oracle A value
    tt = 0f0; sumdr0 = 0f0
    for i in 1:t.n
        tt += t.tpa[i]; sumdr0 += t.tpa[i] * t.dbh[i]^1.605f0
    end
    dia0 = (sumdr0 / tt)^(1f0/1.605f0)
    # both are deterministic sums over the (bit-exact) tree data; they RENDER to the reference exactly
    # (measured dia0 4.7009993→4.701 Δ9.5e-7 Float32 ULP; tt 589.6528→589.65 Δ0.0027). Compare rounded
    # value `==` (was atol 0.005 / 0.5 — ~5000×/180× the print-resolution floor).
    @test round(Float64(dia0); digits=3) == 4.701         # vs quadratic QMD 5.14
    @test round(Float64(tt);   digits=2) == 589.65        # density t (per stockable)
end

@testset "mortality reduces TPA in the right regime" begin
    s, _ = initialize(KEYM); notre!(s)
    setup_growth!(s); diameter_growth!(s, s.variant)
    g = s.plot.gross_space
    before = stand_tpa(s) / g
    mortality!(s, s.variant)
    after = stand_tpa(s) / g
    @test after < before                                  # trees died
    # This isolated mortality! step loses 25.98935 TPA (536.048→510.05862), FULLY DETERMINISTIC. The FULL
    # cycle loses 29 (536→507) — and jl's snt01.sum cyc1 is BIT-EXACT vs live FVSsn (507/103/202, verified
    # via sn_oracle.sh), so jl is NOT behind live; the 3-TPA difference is the full cycle's tripling/regen
    # ordering, not a DG gap. before-after == 25.98935f0 BIT-FOR-BIT (Float32, measured Δ=0.0 — both are
    # stand_tpa/g Float32). Was atol 1f-3 (~300× a real ULP at 26) then ±1.5 — both padding a bit-exact
    # reproducible quantity that masked regressions.
    @test before - after == 25.98935f0
end
