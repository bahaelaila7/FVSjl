# NC (Klamath) small-tree growth cascade — regression 2026-08-25 (extreme-dense >10k-TPA seedling under-kill).
#
# The western full-population sweep flagged 782 NC stands that CONSISTENTLY, one-directionally UNDER-killed in
# the extreme-dense (>10k TPA) SEEDLING regime (worst_col SDI/TPA): the stand held ~16k TPA forever instead of
# self-thinning. Both-sides trace vs FVSnc_clean (CN 13191881010497, 16,750 TPA): live 16,514→3,771 by 2015
# while old jl held 16,281; QMD frozen at ~1.3 vs live 3.7→6.9. The primary hypothesis (a gap in the NWCMRT/
# QMD-loop/BAMAX mortality fix) was REFUTED by instrumenting FVSnc_g16: given the SAME 2010 stand the mortality
# math was correct — it was STARVED by an upstream diameter growth of ~0. The Reineke grown diameter DR10 was
# 0.59 in jl vs 3.24 in FVS, keeping TMD10 pinned at its 35000 cap so the density self-thin never engaged.
#
# THREE cascading root causes (all upstream of mortality; all FVS-source-backed, none cornered):
#   (1) nc/htgr5.f CR argument: regent.f:183 passes CR=ICR(I)/10.0 (0–10 scale). The port passed crown_pct/100
#       (0–1). Since IMETH=1 height growth has a 0.0566·CR² term, the 10× CR error is a 100× term error: +5.11
#       (FVS) vs +0.0005 (jl) ⇒ seedling HTGR flipped −2.70 (jl) vs +2.41 (FVS).
#   (2) nc/htgr5.f floor `IF(HTGR.LE.0.0)HTGR=0.01` was missing.
#   (3) nc/crown.f LSTART crown dub was NEVER CALLED for Klamath (simulate.jl init skipped compute_density! +
#       crown_ratio_update!(lstart=true)), and the d<1 branch had no DUBSCR ⇒ sub-1" inventory seedlings kept
#       crown_pct=0 ⇒ CR²=0 ⇒ negative HTGR ⇒ never crossed 4.5' breast height ⇒ no DBH growth ⇒ QMD frozen ⇒
#       self-thin never fires. This is the IE #137 / EM sibling bug. nc/dubscr.f dubs those to CR 95%.
#
# Signed extreme-dense (>10k TPA) 20-stand tally vs FVSnc_clean: BEFORE = one-directional under-kill (16 over-
# retain / 4 under, worst 891%); AFTER = 10 over / 10 under (balanced straddle, worst 222% on the single
# FVS-flagged >1000-TPA mega-record). The systematic under-kill is eliminated; the residual is the cornered
# mega-record instability FVS itself warns about. Gate test_multicycle stays 339/11 byte-identical.
#
# Self-contained (no DB / no oracle): asserts the two htgr5 invariants + the DUBSCR seedling crown, vs values
# dumped from the instrumented FVSnc_g16 (ZZHTGR / ZZDUB) for CN 13191881010497 cycle 1.
using Test
using FVSjl

@testset "NC (Klamath) small-tree growth cascade" begin
    @testset "nc_htgr5 CR is on the ICR/10 (0–10) scale, not /100 (regent.f:183)" begin
        # FVSnc_g16 ZZHTGR (CN 13191881010497, cyc1, sp4 WF seedling): SSITE=49, BAA=133.541, RELHT=0.0120,
        # CR=9.5, H=1.010 ⇒ HTGR=2.4087 (POSITIVE). The 0.0566·CR² term (+5.11 at CR=9.5) is what lifts it.
        htgr = FVSjl.nc_htgr5(4, 49.0f0, 133.541f0, 0.0120f0, 9.5f0, 1.010f0)
        @test htgr ≈ 2.4087f0 rtol = 1e-3
        # The old /100 bug passed CR≈0.095 for the same tree ⇒ HTGR goes NEGATIVE (0.0566·CR²≈0).
        wrong = FVSjl.nc_htgr5(4, 49.0f0, 133.541f0, 0.0120f0, 0.095f0, 1.010f0)
        @test wrong ≈ 0.01f0                              # negative raw ⇒ floored (see next)
        @test htgr > 2.0f0                                # correct scale ⇒ real, positive height growth
    end

    @testset "nc_htgr5 floors HTGR≤0 to 0.01 (htgr5.f: IF(HTGR.LE.0.0)HTGR=0.01)" begin
        # CR=0 (a still-crownless seedling) ⇒ deterministic raw HTGR≈−2.70 ⇒ must floor to exactly 0.01, never
        # a negative height increment.
        @test FVSjl.nc_htgr5(4, 49.0f0, 133.541f0, 0.0120f0, 0.0f0, 1.010f0) == 0.01f0
    end

    @testset "nc_dubscr dubs sub-1\" seedlings to CR 0.95 (nc/dubscr.f, RNG-independent clamp)" begin
        # FVSnc_g16 ZZDUB (same stand, cyc1 LSTART): sp4 & sp9 seedlings D=0.1, H=1.01, BA=133.541,
        # TPCCF=79.916, AVH=123.12, RMAI=50 ⇒ deterministic logit ≈ −6 ⇒ CR clamps to 0.9500 for ANY FCR draw
        # (|FCR|≤CRSD), so the crown VALUE is bit-exact without reproducing FVS's RANN stream.
        cr4 = FVSjl.nc_dubscr(4, 0.1f0, 1.01f0, 133.541f0, 79.916f0, 123.12f0, 0.0f0, 1.0f0)
        cr9 = FVSjl.nc_dubscr(9, 0.1f0, 1.01f0, 133.541f0, 79.916f0, 123.12f0, 0.0f0, 1.0f0)
        @test cr4 == 0.95f0
        @test cr9 == 0.95f0
        # Clamp bounds hold for both tails.
        @test 0.05f0 <= FVSjl.nc_dubscr(3, 0.9f0, 8.0f0, 5.0f0, 0.0f0, 10.0f0, 0.0f0, 1.0f0) <= 0.95f0
    end
end
