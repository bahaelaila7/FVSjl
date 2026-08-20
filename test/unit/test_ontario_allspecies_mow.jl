# =============================================================================
# test_ontario_allspecies_mow.jl — Mowraski net-merch cull (volont.f) bit-exact across ALL 72 ON species.
#
# ON's .sum board column is the Mowraski net-merch volume (NMV) — GMV culled by an age-based factor
# CV = VM·(1 − (1 − exp(−AM·ABIRTH))^BM), with AM/BM selected per species via ON_ISPMPM. ont01 exercised
# 8 species; this locks the whole set on ont_big (72 trees, 50 cm/25 m so every species passes the board
# gate). An instrumented FVSon_mowdump (volont.f MOWRASKI entry, unit 775; .sum data rows byte-identical
# to clean) dumps per tree VM / ABIRTH / CV; the SHIPPED on_mowraski, fed the oracle's ABIRTH, reproduces
# CV bit-exact (Float32-hex, glibc expf/powf) for every species.
#
# Golden = test/fixtures/ontario/ont_big_mowdump.txt. Validated 2026-08-20: 70/70 distinct species (the 3
# willow WI codes collapse to one ISPC). This isolates the cull equation + ON_ISPMPM/AM/BM tables; age
# dubbing (on_tree_age) is validated separately (ont01 29/29 + ON_MAPLS structurally complete 72/72).
# =============================================================================

using Test, FVSjl
const F = FVSjl

@testset "ON — Mowraski net-merch cull bit-exact across all 72 species (on_mowraski vs FVSon_mowdump)" begin
    dump = joinpath(@__DIR__, "..", "fixtures", "ontario", "ont_big_mowdump.txt")
    if !isfile(dump)
        @test_skip "ont_big Mowraski dump fixture not present"
    else
        h2f(h) = reinterpret(Float32, parse(UInt32, h; base = 16))
        u(x) = reinterpret(UInt32, Float32(x))
        seen = Set{Int}(); n = ok = 0
        for ln in eachline(dump)
            startswith(strip(ln), "#") && continue
            f = split(strip(ln))
            length(f) >= 4 || continue
            sp = parse(Int, f[1])
            vm, ab, cv_o = h2f.(f[2:4])
            sp in seen && continue
            push!(seen, sp)
            cv_j = F.on_mowraski(sp, vm, ab)
            n += 1
            ok += u(cv_j) == u(cv_o)
        end
        @test n >= 60           # most of the 72 species (duplicate willow codes collapse a few)
        @test ok == n           # net-merch cull bit-exact for every species
    end
end
