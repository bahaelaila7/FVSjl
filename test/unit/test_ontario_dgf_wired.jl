# test_ontario_dgf_wired.jl — ON wired large-tree DGF end-to-end vs production FVSon_g16.
#
# Complements test_ontario.jl (which replays the isolated on_penner_dds against a per-tree
# dump). Here we drive the ENGINE-WIRED path: load the oracle's exact internal stand state
# (DIAM/HT/PCT/IMC per tree + stand BA/RMSQD + per-species SITEAR, dumped as Float32-hex from
# an instrumented canada/on/dgf.f), run species_sort!(s) + dgf!(s, Ontario()), and assert each
# tree's WK2 (scratch.wk[2,:]) is bit-exact to the PRODUCTION WK2. This exercises the stand-
# context assembly glue (BAM=BA·FT2pACRtoM2pHA, QMDM=RMSQD·INtoCM, SIM=SITEAR·FTtoM,
# BALM=(1−PCT/100)·BAM), OSPMAP/AGS selection (sugar-maple/beech get the +1 AGS eqn), the
# on_bratio grown-diameter bark, and the faithful dgf.f DO-303 quirk where the DDS term reuses
# the leftover BARK scalar from the last species processed (beech, 0.94) for every tree.
#
# The StandState is built with the Southern() container (a fresh Ontario() has no coefficients/
# species table yet — those are later chunks); dgf! dispatches on the Ontario() argument, not
# s.variant, so the container is irrelevant to the code path under test.

using Test

@testset "ON — wired dgf! WK2 bit-exact vs FVSon_g16 (ont01 cyc0)" begin
    _h2f(h) = reinterpret(Float32, parse(UInt32, h; base = 16))
    golden = joinpath(@__DIR__, "ontario_dgf_wired_dump.txt")

    rows = NamedTuple[]
    for ln in eachline(golden)
        startswith(strip(ln), "#") && continue
        f = split(strip(ln)); isempty(f) && continue
        I = parse(Int, f[1]); isp = parse(Int, f[2]); imc = parse(Int, f[3]); hx = f[4]
        push!(rows, (; I, isp, imc,
            diam = _h2f(hx[1:8]), ht = _h2f(hx[9:16]), pct = _h2f(hx[17:24]),
            sitear = _h2f(hx[25:32]), ba = _h2f(hx[33:40]), rmsqd = _h2f(hx[41:48]),
            wk2 = _h2f(f[5])))
    end
    @test !isempty(rows)

    s = FVSjl.StandState(FVSjl.Southern())      # container; dgf! dispatches on the Ontario() arg
    t = s.trees; p = s.plot; c = s.calib
    t.n = length(rows)
    p.basal_area = rows[1].ba                    # BA (ft²/ac) — stand-level, identical across trees
    p.qmd = rows[1].rmsqd                        # RMSQD (in)
    for r in rows
        t.species[r.I] = Int32(r.isp)
        t.dbh[r.I] = r.diam
        t.height[r.I] = r.ht
        t.crown_ratio[r.I] = r.pct               # FVS PCT
        t.mort_code[r.I] = Int32(r.imc)          # IMC
        t.tpa[r.I] = 1f0
        t.sort_key[r.I] = Int32(r.I)
        p.sp_site_index[r.isp] = r.sitear        # SITEAR(sp)
        c.dg_cor[r.isp] = 0f0                     # ground-truth COR = 0 at cyc0
    end

    FVSjl.species_sort!(s)
    FVSjl.dgf!(s, FVSjl.Ontario())
    wk2 = view(s.scratch.wk, 2, :)

    nok = 0
    for r in rows
        @test reinterpret(UInt32, wk2[r.I]) == reinterpret(UInt32, r.wk2)
        nok += reinterpret(UInt32, wk2[r.I]) == reinterpret(UInt32, r.wk2)
    end
    @test nok == length(rows)                    # all trees bit-exact
end
