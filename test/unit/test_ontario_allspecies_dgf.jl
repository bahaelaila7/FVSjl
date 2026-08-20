# =============================================================================
# test_ontario_allspecies_dgf.jl — Penner large-tree DGF bit-exact across ALL 72 ON species.
#
# The ont01 reference stand exercises only 8 of the 72 Ontario species (PW/SW/MH/BE/SB/CE/PJ/BF),
# hitting 8 of the 27 distinct Penner DGF equations. This locks the FULL species set: a 72-tree
# stand (one tree per species, uniform 25 cm DBH / 18 m HT so per-species coefficient differences
# isolate) is dumped from an instrumented FVSon_g16 (canada/on/dgf.f, unit 771; the instrumented
# .sum is byte-identical to clean). The SHIPPED on_penner_dds is fed the oracle's per-tree inputs
# and must reproduce DBHM_final / DIAGR / DDS bit-exact (Float32-hex) for every record.
#
# Golden = test/fixtures/ontario/ont_all_dgfdump.txt. Validated 2026-08-20: 216/216 (72 sp × 3
# subcycle records) all three columns bit-exact — the Penner coefficient tables (OSPMAP 72→35,
# B0..BSI + AGS/UGS quality) are complete and correct across the whole variant.
# =============================================================================

using Test, FVSjl
const F = FVSjl

@testset "ON — Penner DGF bit-exact across all 72 species (on_penner_dds vs FVSon_g16)" begin
    dump = joinpath(@__DIR__, "..", "fixtures", "ontario", "ont_all_dgfdump.txt")
    if !isfile(dump)
        @test_skip "ont_all DGF dump fixture not present"
    else
        h2f(h) = reinterpret(Float32, parse(UInt32, h; base = 16))
        n = okB = okG = okD = 0
        for ln in eachline(dump)
            startswith(strip(ln), "#") && continue
            f = split(strip(ln))
            length(f) >= 15 || continue
            _, ispc, ksp, ags, _icyc = parse.(Int, f[1:5])
            diam, dbhm_o, sim, bam, qmdm, balm, htm, bark, diagr_o, dds_o = h2f.(f[6:15])
            dbhm_j, diagr_j, dds_j = F.on_penner_dds(ksp, ags, diam, sim, bam, qmdm, balm, htm, bark)
            n += 1
            okB += reinterpret(UInt32, dbhm_j)  == reinterpret(UInt32, dbhm_o)
            okG += reinterpret(UInt32, diagr_j) == reinterpret(UInt32, diagr_o)
            okD += reinterpret(UInt32, dds_j)   == reinterpret(UInt32, dds_o)
        end
        @test n >= 72                 # at least one record per species
        @test okB == n                # DBHM_final bit-exact for every record
        @test okG == n                # DIAGR bit-exact
        @test okD == n                # DDS bit-exact
    end
end
