# test_ontario_growth_wired.jl — ON runnable-growth wiring end-to-end vs FVSon_g16.
#
# Three things this asserts, all against the LIVE oracle (FVSon_g16 / FVSon_wkidump):
#   (1) SPECIES TABLE + TRANSLATION: coefficients(::Ontario) loads data/ontario, and the shared
#       reader resolve_species() maps the 8 ont01 alpha codes to the SAME ON species indices the
#       oracle used (PW→5 SW→6 MH→26 BE→28 SB→9 CE→11 PJ→1 BF→8, from the g16 dgf dump ISPC).
#   (2) coefficients(::Ontario): StandState(Ontario()) constructs STANDALONE (no Southern container).
#   (3) DDS→DG on_bratio BRANCH: on the oracle's exact internal stand state, the DBH-increment path
#       d_ib = DBH·on_bratio(sp,DBH,HT) (dgdriv.f:201, ORIGINAL DBH), DDS = exp(WK2), and the
#       deterministic expected increment WKI = √(d_ib²+DDS)−d_ib are each bit-exact (Float32-hex) to
#       FVSon_wkidump (an instrumented canada/on/dgdriv.f). This is the DDS→DG round-trip the shared
#       diameter_growth! driver performs; the VARDG/COR serial-correlation spread + full-cycle .sum
#       (site/crown/mort/htg/vol) are downstream chunks.

using Test

@testset "ON — growth wiring: species table+translation + coefficients + DDS→DG on_bratio (ont01 cyc0)" begin
    _h2f(h) = reinterpret(Float32, parse(UInt32, h; base = 16))
    _u(x::Float32) = reinterpret(UInt32, x)

    # (2) StandState(Ontario()) constructs standalone via coefficients(::Ontario)
    s = FVSjl.StandState(FVSjl.Ontario())
    @test FVSjl.nspecies(FVSjl.Ontario()) == 72
    FVSjl.load_species_coefficients!(s, s.variant)

    # (1) species translation through the shared reader matches the oracle's ISPC mapping
    want = Dict("PW"=>5, "SW"=>6, "MH"=>26, "BE"=>28, "SB"=>9, "CE"=>11, "PJ"=>1, "BF"=>8)
    for (code, idx) in want
        got, _ = FVSjl.resolve_species(code, s.variant, s.species, s.coef)
        @test Int(got) == idx
    end

    # (3) DDS→DG branch bit-exact vs FVSon_wkidump.
    # Load the oracle internal stand state (DIAM/HT/PCT/IMC + BA/RMSQD/SITEAR) from the wired golden.
    stategold = joinpath(@__DIR__, "ontario_dgf_wired_dump.txt")
    rows = NamedTuple[]
    for ln in eachline(stategold)
        startswith(strip(ln), "#") && continue
        f = split(strip(ln)); isempty(f) && continue
        I=parse(Int,f[1]); isp=parse(Int,f[2]); imc=parse(Int,f[3]); hx=f[4]
        push!(rows, (; I, isp, imc,
            diam=_h2f(hx[1:8]), ht=_h2f(hx[9:16]), pct=_h2f(hx[17:24]),
            sitear=_h2f(hx[25:32]), ba=_h2f(hx[33:40]), rmsqd=_h2f(hx[41:48])))
    end
    # WKI golden: I ISPC D DDS WKI
    wki = Dict{Int,NamedTuple}()
    for ln in eachline(joinpath(@__DIR__, "ontario_dgdriv_wki_dump.txt"))
        startswith(strip(ln), "#") && continue
        f = split(strip(ln)); isempty(f) && continue
        wki[parse(Int,f[1])] = (; D=_h2f(f[3]), dds=_h2f(f[4]), wki=_h2f(f[5]))
    end

    sd = FVSjl.StandState(FVSjl.Southern())   # container; dgf! + on_bratio dispatch on the value/arg
    t=sd.trees; p=sd.plot; c=sd.calib
    t.n = length(rows)
    p.basal_area = rows[1].ba; p.qmd = rows[1].rmsqd
    for r in rows
        t.species[r.I]=Int32(r.isp); t.dbh[r.I]=r.diam; t.height[r.I]=r.ht
        t.crown_ratio[r.I]=r.pct; t.mort_code[r.I]=Int32(r.imc); t.tpa[r.I]=1f0
        t.sort_key[r.I]=Int32(r.I); p.sp_site_index[r.isp]=r.sitear; c.dg_cor[r.isp]=0f0
    end
    FVSjl.species_sort!(sd)
    FVSjl.dgf!(sd, FVSjl.Ontario())
    wk2 = view(sd.scratch.wk, 2, :)

    nD=nDDS=nW=0
    for r in rows
        i=r.I
        bark = FVSjl.on_bratio(Int(r.isp), t.dbh[i], t.height[i])
        d_ib = t.dbh[i]*bark
        dds  = FVSjl.on_expf(wk2[i])
        wkival = sqrt(d_ib*d_ib + dds) - d_ib
        g = wki[i]
        @test _u(d_ib) == _u(g.D)
        @test _u(dds) == _u(g.dds)
        @test _u(wkival) == _u(g.wki)
        nD += _u(d_ib)==_u(g.D); nDDS += _u(dds)==_u(g.dds); nW += _u(wkival)==_u(g.wki)
    end
    @test nD == length(rows) && nDDS == length(rows) && nW == length(rows)
end
