# Golden-replay regression for the LPMPB LPOPDY chain (src/engine/lpmpb_lpopdy.jl) vs the
# oracle dumps in test/fixtures/lpopdy/ (FVSie_lpmpb, lp_popdy stand). Validates BETIN, GARBEL,
# SURFCE, and the MPBMOD epidemic (total kill) — the engine-ported functions reproduce the oracle.

using Test
using FVSjl: mpb_betin, mpb_garbel, mpb_surflp, mpb_mpbmod

const _FX = joinpath(@__DIR__, "..", "fixtures", "lpopdy")
_h2d(s) = reinterpret(Float64, parse(UInt64, s, base=16))
_h2f(s) = reinterpret(Float32, parse(UInt32, s, base=16))
_fhex(x::Float32) = uppercase(string(reinterpret(UInt32, x), base=16, pad=8))

@testset "LPOPDY chain (golden replay vs FVSie_lpmpb)" begin
    # --- BETIN (betin.txt: a b x result, double hex) ---
    @testset "BETIN incomplete-beta" begin
        n = 0; worst = 0.0
        for ln in eachline(joinpath(_FX, "betin.txt"))
            p = split(ln); p[1] == "BETIN" || continue
            a = _h2d(p[2]); b = _h2d(p[3]); x = _h2d(p[4]); r = _h2d(p[5])
            got = mpb_betin(a, b, x)
            worst = max(worst, r != 0 ? abs(got-r)/abs(r) : abs(got-r)); n += 1
        end
        @test n == 30
        @test worst == 0.0            # BIT-EXACT (Fortran single-precision Stirling constant matched)
    end

    # --- load per-tree inputs (inputs.txt: IN i DBH XPT WK3 PROB) ---
    DBH = Float32[]; XPT = Float32[]; PROB = Float32[]
    for ln in eachline(joinpath(_FX, "inputs.txt"))
        t = split(ln); t[1] == "IN" || continue
        push!(DBH, _h2f(t[3])); push!(XPT, _h2f(t[4])); push!(PROB, _h2f(t[6]))
    end
    @test length(DBH) == 18

    mp1, mp2, ipt, clsprob, avgdbh, avgxpt = mpb_garbel(DBH, XPT, PROB, 10)

    # --- GARBEL (garbel.txt: CLS c mp1 mp2 probhex ; MEM c ipt) ---
    @testset "GARBEL classifier" begin
        orc_mem = Dict{Int,Vector{Int}}(); orc_hex = Dict{Int,String}()
        for ln in eachline(joinpath(_FX, "garbel.txt"))
            p = split(ln)
            if p[1] == "CLS"; orc_hex[parse(Int,p[2])] = uppercase(p[5])
            elseif p[1] == "MEM"; push!(get!(orc_mem, parse(Int,p[2]), Int[]), parse(Int,p[3])); end
        end
        @test length(mp1) == 10
        for c in 1:length(mp1)
            @test [ipt[jj] for jj in mp1[c]:mp2[c]] == orc_mem[c]
            @test _fhex(clsprob[c]) == orc_hex[c]
        end
    end

    # --- SURFCE (surfce.txt: SURF c avgDBHhex surfhex) ---
    @testset "SURFCE surface areas" begin
        for ln in eachline(joinpath(_FX, "surfce.txt"))
            p = split(ln); p[1] == "SURF" || continue
            c = parse(Int, p[2])
            @test _fhex(avgdbh[c]) == uppercase(p[3])
            @test _fhex(mpb_surflp(avgdbh[c])) == uppercase(p[4])
        end
    end

    # --- MPBMOD epidemic → total kill (mortality.txt CLPOP survivors << cap threshold) ---
    @testset "MPBMOD epidemic (total kill)" begin
        # params from params.txt
        P = Dict{String,Float32}()
        for ln in eachline(joinpath(_FX, "params.txt"))
            t = split(ln); t[1] == "P" && (P[t[2]] = _h2f(t[3]))
        end
        surf = Float32[mpb_surflp(avgdbh[c]) for c in 1:length(mp1)]
        trees0 = copy(clsprob)
        surv = mpb_mpbmod(trees0, surf, avgdbh, avgxpt, P["TA"], P["EFELEV"], P["EFLAT"])
        # every class driven below the 3.6e-7 mortality-cap threshold (survival ratio ~0)
        maxratio = maximum(abs.(surv ./ clsprob))
        @test maxratio < 3.6e-7
    end
end
