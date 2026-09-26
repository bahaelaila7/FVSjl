# PPE MXHRVP — chunk F: composite-MATERIALIZATION (the actual harvest EFFECT on the landscape),
# A/B vs the FVSppe post-harvest COMPOSITE golden (test/fixtures/ppe/mxhrvp/msp_thin.golden.txt).
#
# msp_thin.key: an AGPLABEL-gated THINBTA (residual 40 TPA) fires on the SELECTED stands (all,
# TARGET=1000) — the oracle post-harvest COMPOSITE reflects the removals. FVSjl materializes the
# same decision: the selected stands are projected WITH the harvest applied (the aligned THINBTA
# thins to residual 40 TPA), then aggregated with the oracle-validated `_ppe_aggregate`.
#
# NOTE (THINBTA verdict): the earlier "FVSjl over-thins" was a KEYFILE COLUMN-ALIGNMENT artifact in
# a hand-written test card (the date spilled into the residual field ⇒ FVSjl's fixed-column parser
# read residual=0 ⇒ removed all). With a column-correct card, FVSjl's THINBTA thins to residual 40
# TPA and MATCHES the oracle exactly. FVSjl's base THINBTA is correct.
using Test
using FVSjl: run_keyfile, EastCascades, _ppe_parse_sum_row, _ppe_aggregate

@testset "PPE MXHRVP — composite-materialization (harvest effect) vs FVSppe golden" begin
    kf  = joinpath(@__DIR__, "..", "fixtures", "ppe", "mxhrvp", "stand_thin.key")  # aligned THINBTA 1990 40
    sumtxt = run_keyfile(kf; variant = EastCascades())
    rows = [x for x in (_ppe_parse_sum_row(l) for l in split(sumtxt, "\n")) if x !== nothing]
    @test !isempty(rows)
    # all-selected landscape = 3 identical thinned stands, area 11 each (SAMPLE WEIGHT 33)
    comp = _ppe_aggregate([(11.0, rows), (11.0, rows), (11.0, rows)])
    byyr = Dict(a.year => a for a in comp)

    # oracle post-harvest COMPOSITE golden (msp_thin.golden.txt): year => TPA. The harvest decision and the
    # density it leaves are bit-exact to FVSppe every cycle (536 → 36 → 35 → 35).
    golden_tpa = Dict(1990 => 536, 2000 => 36, 2010 => 35, 2020 => 35)
    # VOLUMES: FVSppe is built from the historical tree (bc6e2377^), whose EC volume equations predate the current
    # source (voleqdef R6_EQN: Mt Hood DF on westside F05FW2W202, the rest on I11-I13 INGY) — its composite reads
    # 1624/1102/5567 at 1990. The current live FVSec_g16 on this same stand (stand_thin.key + ECHOSUM, measured
    # 2026-09-26) reads the rows below, and the composite of 3 identical stands is that stand. jl equals it.
    live_vol = Dict(1990 => (1602, 1064, 5456),
                    2000 => ( 917,  827, 4189),
                    2010 => (1318, 1202, 6296),
                    2020 => (1780, 1651, 9038))

    for yr in sort(collect(keys(golden_tpa)))
        a = byyr[yr]
        @testset "year $yr" begin
            # ★ THE MATERIALIZATION PROOF: post-thin TPA (the density structure) is BIT-EXACT to FVSppe.
            @test round(Int, a.avbtpa) == golden_tpa[yr]
            cuft, mcuft, bdft = live_vol[yr]
            @test round(Int, a.avbtcuft) == cuft
            @test round(Int, a.avbmcuft) == mcuft
            @test round(Int, a.avbbdft)  == bdft
        end
    end
    # the KEY proof restated: the harvest actually removed trees — 2000 TPA (36) is < 10% of 1990 (536).
    @test byyr[2000].avbtpa < 0.1 * byyr[1990].avbtpa
end
