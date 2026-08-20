# =============================================================================
# test_ontario_allspecies_htg.jl — Penner diameter-height (htont.f) bit-exact across ALL 72 ON species.
#
# Companion to test_ontario_allspecies_dgf.jl. The ont01 stand exercised only 8 species; this locks
# the height model on the full 72-species ont_all stand. An instrumented FVSon_htdump (canada/on/
# htgf.f, unit 773; .sum data rows byte-identical to clean) dumps, per tree, the two HTONT calls'
# inputs (DBH/RMSQD/BA for the current height, DBH10/QMD10/BA10 for the projected) plus the constant
# SIM (= SITEAR(ISISP)·FTtoM) and the resulting HTNOW / HT10. The SHIPPED on_htont must reproduce
# both bit-exact (Float32-hex) for every species — validating the OSPMAP_SP/OSPMAP_P height-group
# branch and the α/β/δ/φ/γ + B0..BBA coefficient tables across the whole variant.
#
# Golden = test/fixtures/ontario/ont_all_htdump.txt. Validated 2026-08-20: HTNOW 72/72, HT10 72/72.
# =============================================================================

using Test, FVSjl
const F = FVSjl

@testset "ON — Penner htont bit-exact across all 72 species (on_htont vs FVSon_htdump)" begin
    dump = joinpath(@__DIR__, "..", "fixtures", "ontario", "ont_all_htdump.txt")
    if !isfile(dump)
        @test_skip "ont_all htont dump fixture not present"
    else
        h2f(h) = reinterpret(Float32, parse(UInt32, h; base = 16))
        n = okN = ok10 = 0
        for ln in eachline(dump)
            startswith(strip(ln), "#") && continue
            f = split(strip(ln))
            length(f) >= 12 || continue
            _, ispc = parse.(Int, f[1:2])
            dbh, dbh10, rmsqd, ba, qmd10, ba10, sim, htnow_o, ht10_o, _ = h2f.(f[3:12])
            htnow_j = F.on_htont(ispc, dbh, rmsqd, ba, sim)
            ht10_j  = F.on_htont(ispc, dbh10, qmd10, ba10, sim)
            n += 1
            okN  += reinterpret(UInt32, htnow_j) == reinterpret(UInt32, htnow_o)
            ok10 += reinterpret(UInt32, ht10_j)  == reinterpret(UInt32, ht10_o)
        end
        @test n >= 72
        @test okN == n        # current height HTONT bit-exact for every species
        @test ok10 == n       # projected height HTONT bit-exact
    end
end
