# =============================================================================
# test_ontario_allspecies_cw.jl — open-grown crown width (cwcalc.f IWHO=1) bit-exact across ALL 72 ON species.
#
# Open-grown crown width feeds the .sum CCF (ccfcal.f 0.001803·CW²·P). ON's cwcalc.f remaps each ISPC to
# a US 2-char code (ON_JSP2) and selects an Ek/Bechtold/Krajicek equation. Before this chunk only the 8
# species ont01 exercises had equations wired (ON_CW_OPEN/ON_CW_EQS); the other 64 returned CW=0 → the
# stand CCF was ~8× too low (a real bug the full-species stand exposed: the .sum CCF read 1322 where the
# oracle overflowed its I4 field). This locks the whole set: an instrumented FVSon_cwdump (unit 772; .sum
# data rows byte-identical to clean) dumps per tree D / CW / HILAT / HILONG / HIELEV / HI; the SHIPPED
# on_open_crown_width reproduces CW bit-exact for every species.
#
# Golden = test/fixtures/ontario/ont_all_cwdump.txt. Validated 2026-08-20: 70/70 distinct species (the 3
# willow WI codes collapse to one ISPC). HIELEV in the dump is full feet; on_open_crown_width takes
# elevation in HUNDREDS of feet (FVS ELEV storage), hence the /100.
# =============================================================================

using Test, FVSjl
const F = FVSjl

@testset "ON — open-grown crown width bit-exact across all 72 species (on_open_crown_width vs FVSon_cwdump)" begin
    dump = joinpath(@__DIR__, "..", "fixtures", "ontario", "ont_all_cwdump.txt")
    key  = joinpath(@__DIR__, "..", "fixtures", "ontario", "ont_all.key")
    if !isfile(dump) || !isfile(key)
        @test_skip "ont_all crown-width fixture not present"
    else
        h2f(h) = reinterpret(Float32, parse(UInt32, h; base = 16))
        u(x) = reinterpret(UInt32, Float32(x))
        # map tree index → species from the loaded ont_all stand (one tree per species)
        sp_of = Int[]
        for s in F.each_stand(key; variant = F.Ontario())
            F.notre!(s)
            for i in 1:s.trees.n; push!(sp_of, Int(s.trees.species[i])); end
            break
        end
        seen = Set{Int}(); n = ok = 0
        for ln in eachline(dump)
            startswith(strip(ln), "#") && continue
            f = split(strip(ln))
            length(f) >= 8 || continue
            tree = parse(Int, f[1])
            tree <= length(sp_of) || continue
            d, cw_o, lat, long, elev, _hi = h2f.(f[3:8])
            sp = sp_of[tree]
            sp in seen && continue
            push!(seen, sp)
            cw_j = F.on_open_crown_width(sp, d, lat, long, elev / 100f0)  # HIELEV full-ft → hundreds-of-ft
            n += 1
            ok += u(cw_j) == u(cw_o)
        end
        @test n >= 60             # most of the 72 species (duplicate alpha codes collapse a few)
        @test ok == n             # open-grown crown width bit-exact for every species
    end
end
