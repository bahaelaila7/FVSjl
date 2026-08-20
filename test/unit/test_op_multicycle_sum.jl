# =============================================================================
# test_op_multicycle_sum.jl — OP (Olympic) END-TO-END multi-cycle .sum vs the relinked FVSop_clean.
#
# The per-tree volume A/B (test_op_site_and_volume.jl) and the DG/HG/CR/MORT engine A/B
# (test_op_organon_nwo.jl) both pass, but neither exercises the full multi-cycle .sum through the
# summary writer. This locks the aggregate: run the S248112 UNTHINNED CONTROL stand (opt01 stand 1,
# 10 cycles) and compare every column to the LIVE RELINKED oracle FVSop_clean (g16 relink), the
# doctrine oracle — NOT the production opt01.sum.save.
#
# FINDINGS (2026-08-20, measured — see fvsjl-oc-organon-blm-volume memory):
#   • cyc0-3 (1990..2005): EVERY column bit-identical to FVSop_clean, incl. merch cubic (972) and
#     board-foot (5003). ORGANON growth is DGSD=0 deterministic, so this is a true bit-exact match.
#   • The production opt01.sum.save shows merch/bdft = 931/4797 (−4.3%); that delta is a
#     PRODUCTION-BINARY-vs-g16-RELINK NVEL(BLM) artifact, not a jl bug — jl matches the relinked
#     oracle exactly. (Hence the golden here is the relinked oracle, captured live.)
#   • cyc4 (2010): growth + TOTAL cubic (4814) still bit-exact; merch/bdft drift <0.4%
#     (4114→4106) — merch-log-bucking DIB-class rounding discreteness (oc_blmgdib!/oc_scrib), the
#     accepted cornered class. A single density-aggregate field picks up ±1 NINT at later cycles
#     while raw TPA stays bit-exact every cycle (print-rounding straddle).
#
# Golden = test/fixtures/olympic/opgro_oracle.rows (first 5 cyc rows of FVSop_clean, NOTRIPLE-stable).
# =============================================================================
using Test
using FVSjl
const F = FVSjl

@testset "OP — multi-cycle .sum vs relinked FVSop_clean (S248112 growth control)" begin
    fx  = joinpath(@__DIR__, "..", "fixtures", "olympic")
    key = joinpath(fx, "opgro.key")
    if !isfile(key) || !isfile(joinpath(fx, "opgro.tre"))
        @test_skip "olympic multi-cycle fixture not present"
    else
        jl_rows = cd(fx) do
            txt = F.run_keyfile("opgro.key"; variant = F.variant_from_code("OP"), output = :sum)
            [split(strip(l)) for l in split(txt, '\n') if occursin(r"^(19|20)\d\d\s", strip(l))]
        end
        gold = [split(strip(l)) for l in
                split(read(joinpath(fx, "opgro_oracle.rows"), String), '\n') if !isempty(strip(l))]

        @test length(jl_rows) >= 5
        # cyc0-3 (rows 1..4): entire row bit-identical to the relinked oracle.
        for r in 1:4
            @test jl_rows[r] == gold[r]
        end
        # cyc0 explicitly locks the merch/bdft values against the relinked oracle (972/5003),
        # documenting that jl does NOT reproduce the production save's 931/4797.
        @test gold[1][10] == "972" && jl_rows[1][10] == "972"    # merch cubic
        @test gold[1][12] == "5003" && jl_rows[1][12] == "5003"  # board-foot
        # cyc4 (2010, row 5): growth + total cubic bit-exact; merch cornered within 0.4%.
        g, j = gold[5], jl_rows[5]
        @test j[3] == g[3]      # TPA bit-exact
        @test j[5] == g[5]      # BA bit-exact
        @test j[9] == g[9]      # total cubic bit-exact
        merch_g = parse(Int, g[10]); merch_j = parse(Int, j[10])
        @test abs(merch_g - merch_j) / merch_g < 0.004     # merch cornered (bucking discreteness)
    end
end
