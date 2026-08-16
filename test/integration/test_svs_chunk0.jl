# test_svs_chunk0.jl — SVS (Stand Visualization System) data path, chunk 0.
#
# Validates the cycle-0 inventory SVS picture against the LIVE relinked FVSkt oracle
# (stand S248112, kt0.key/kt0.tre; SVS 0 = IPLGEM=0 square acre, integer TPA). The oracle
# `kt0_001.svs` was produced by /workspace/.ktwork/FVSkt_clean and is committed as a fixture.
#
# Exercises the whole chunk-0 path: SVKEY flag (kw_svs!), the 2nd Park–Miller stream
# (svrann!, divisor 2^31) seeded from the main stream at the SVSTART seam, SVGTPL(IPLGEM=0),
# SVESTB integer placement, SVGTPT rectangle draws, SVOBOL/SVCROL overlap (inert here),
# KT western crown width (CWCALC Crookston R1), and the SVOUT header + live-tree object loop.
#
# The RNG seam parity is exact and unforced: the main stream's s0 at the SVSTART seam already
# equals the oracle's SVS0 (1211791048) with no adjustment, because the Julia setup consumes
# the identical main-stream draws as FVS up to fvs.f:333.

using Test, FVSjl

const _SVS_DIR = joinpath(@__DIR__, "..", "fixtures", "svs")

@testset "SVS chunk 0 — cyc0 object list bit-exact vs live FVSkt" begin
    key = joinpath(_SVS_DIR, "kt0.key")
    s = FVSjl.each_stand(key; variant = FVSjl.Kootenai())[1]
    FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)

    # SVKEY parsed the flag.
    @test s.control.svs_on
    @test Int(s.control.svs_iplgem) == 0
    @test Int(s.control.svs_imetric) == 0

    # RNG seam parity: main stream s0 at the SVSTART seam == oracle SVS0 (measured from FVSkt_svsdbg).
    @test s.rng.s0 == 1211791048.0

    # 2nd stream matches svrann.f exactly (Park–Miller 16807/2^31-1, divisor 2^31) from that seed.
    FVSjl.svs_seed!(s.rng)
    draw1 = FVSjl.svrann!(s.rng)
    @test draw1 ≈ 0.924207f0 rtol=0 atol=1f-6      # first SVGTPT x-draw ⇒ xloc 192.89 on the 208.71 acre

    # Re-setup (the draw above advanced the SVS stream) and render the picture end-to-end.
    s = FVSjl.each_stand(key; variant = FVSjl.Kootenai())[1]
    FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
    body   = FVSjl.svs_render_cycle0(s)
    oracle = read(joinpath(_SVS_DIR, "kt0_001.svs.oracle"), String)

    @test body == oracle                            # BYTE-IDENTICAL cyc0 object block
end
