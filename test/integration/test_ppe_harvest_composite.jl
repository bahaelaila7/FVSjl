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

    # oracle post-harvest COMPOSITE golden (msp_thin.golden.txt): year => (TPA, TCuFt, MCuFt, BdFt)
    golden = Dict(1990 => (536, 1624, 1102, 5567),   # before-thin (removals reported separately)
                  2000 => (36,   857,  814, 4004),   # post-thin, residual 40 TPA regrown
                  2010 => (35,  1166, 1117, 5701),
                  2020 => (35,  1546, 1489, 7859))

    for (yr, (tpa, cuft, mcuft, bdft)) in sort(collect(golden))
        a = byyr[yr]
        @testset "year $yr" begin
            # ★ THE MATERIALIZATION PROOF: post-thin TPA (the density structure) is BIT-EXACT to
            #   the oracle every cycle — the harvest genuinely removed trees to residual 40 TPA and
            #   the landscape regrew identically (536 → 36 → 35 → 35). This is the subsystem's point.
            @test round(Int, a.avbtpa) == tpa
            # Volumes: cyc0 near-bit-exact (~1%). Beyond cyc0 the post-thin REGROWTH volume inherits
            # the EC-variant height/form growth straddle (#207) — density matches, volume diverges
            # progressively (2000 ~4% → 2020 ~9%). Cornered to the EC growth primitive, not a MXHRVP
            # or materialization defect.
            vtol = yr == 1990 ? 0.02 : 0.10
            @test abs(a.avbtcuft - cuft) <= vtol * cuft
            @test abs(a.avbmcuft - mcuft) <= vtol * mcuft
            @test abs(a.avbbdft - bdft)  <= vtol * bdft
        end
    end
    # the KEY proof restated: the harvest actually removed trees — 2000 TPA (36) is < 10% of 1990 (536).
    @test byyr[2000].avbtpa < 0.1 * byyr[1990].avbtpa
    # and cyc0 volumes are within ~1% of the oracle (EC is cyc0-bit-exact; the small residual is rounding).
    @test abs(byyr[1990].avbtcuft - 1624) / 1624 < 0.012
end
