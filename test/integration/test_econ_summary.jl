# ECON FVS_EconSummary / FVS_EconHarvestValue end-to-end vs live FVS (eccalc.f ECCALC → dbsecsum.f / dbsecharv.f).
# Golden rows were measured from the live oracles (FVSsn_g16 / FVSbm_g16 with main.f built WITHOUT -std=legacy — the
# legacy-runtime oracle zeroes the backward-T-tab TREEFMT fields these .tre fixtures carry; == FVS*_clean):
#   * econ_strtecon / econ_u5 (SN): STRTECON 5% + ANNUCST + HRVRVN log-graded board (unit 4) / cubic (unit 5).
#   * econ_bm_rich (BM bmt01 stand 3): STRTECON 4% + computed SEV, ANNUCST w/ '&' appreciation, ANNURVN, PCTSPEC/
#     PCTFXCST/PCTVRCST (the VAR_PCT-reads-varHrv quirk), HRVFXCST/HRVVRCST, whole-tree HRVRVN, SPECCST/SPECRVN.
#     Periods 1990–2010 (before the stand's own growth diverges): the 1990 THINPRSC removal is valued as a PCT.
#   * econ_bm_logs (BM bmt01 stand 3): HRVRVN unit-4 log grades on the R6 FW2 per-log ECVOL arrays; 1990 period.
# The FVS table-writer quirks are part of the contract: HarvestValue's prepared INSERT keeps the previous row's
# bindings for columns FVS does not re-bind (<0), and values are REAL→REAL*8 widened.
using SQLite, DBInterface

const _ECON_SCEN = joinpath(@__DIR__, "..", "harness", "scenarios")

function _econ_run(key::AbstractString, variant)
    tre = replace(key, r"\.key$" => ".tre")
    (isfile(key) && isfile(tre)) || return nothing
    d = mktempdir()
    cp(key, joinpath(d, basename(key))); cp(tre, joinpath(d, basename(tre)))
    cd(() -> FVSjl.run_keyfile(basename(key); variant = variant), d)
    return joinpath(d, "FVSOut.db")
end

_econ_rows(db, t, w = "") = [Tuple(r) for r in DBInterface.execute(SQLite.DB(db), "SELECT * FROM $t $w ORDER BY rowid")]

_eqv(a, b) = (ismissing(a) || ismissing(b)) ? (ismissing(a) && ismissing(b)) :
             (a isa AbstractFloat || b isa AbstractFloat) ? isapprox(Float64(a), Float64(b); rtol = 1e-6, atol = 1e-6) :
             a == b

function _check_rows(got, want, skip)
    @test length(got) == length(want)
    for (g, w) in zip(got, want)
        gs = g[skip+1:end]
        @test length(gs) == length(w)
        @test all(_eqv(x, y) for (x, y) in zip(gs, w))
    end
end

const M = missing
@testset "ECON FVS_EconSummary / HarvestValue vs live FVS" begin
    @testset "SN econ_strtecon (log-graded board)" begin
        db = _econ_run(joinpath(_ECON_SCEN, "econ_strtecon.key"), FVSjl.Southern())
        if db === nothing
            @test_skip "econ_strtecon fixture not available"
        else
            _check_rows(_econ_rows(db, "FVS_EconSummary"), [
                (1990, 5, "NO ", 15.0, 0.0, 13.637852668762207, 0.0, -13.637852668762207, M, 0.0, M, M, M, M, 0, 0, 5.0, M),
                (1995, 10, "NO ", 30.0, 0.0, 24.323471069335938, 0.0, -24.323471069335938, M, 0.0, M, M, M, M, 0, 0, 5.0, M),
                (2000, 15, "NO ", 45.0, 0.0, 32.6959342956543, 0.0, -32.6959342956543, M, 0.0, M, M, M, M, 0, 71, 5.0, M)], 2)
            _check_rows(_econ_rows(db, "FVS_EconHarvestValue"), [
                (2000, "SM", "ACSA3", "318", 10.0, 999.9000244140625, M, M, M, M, M, M, M, 16, 5, 5),
                (2000, "HI", "CARYA", "400", 10.0, 999.9000244140625, M, M, M, M, M, M, M, 9, 3, 3),
                (2000, "AB", "FAGR", "531", 10.0, 999.9000244140625, M, M, M, M, M, M, M, 41, 12, 12),
                (2000, "SK", "QUFA", "812", 10.0, 999.9000244140625, M, M, M, M, M, M, M, 5, 1, 1)], 1)
        end
    end
    @testset "SN econ_u5 (log-graded cubic, SCFV(I))" begin
        db = _econ_run(joinpath(_ECON_SCEN, "econ_u5.key"), FVSjl.Southern())
        if db === nothing
            @test_skip "econ_u5 fixture not available"
        else
            _check_rows(_econ_rows(db, "FVS_EconSummary", "WHERE Year = 2000"), [
                (2000, 15, "NO ", 45.0, 0.0, 32.6959342956543, 0.0, -32.6959342956543, M, 0.0, M, M, M, M, 14, 0, 5.0, M)], 2)
            _check_rows(_econ_rows(db, "FVS_EconHarvestValue"), [
                (2000, "SM", "ACSA3", "318", 10.0, 999.9000244140625, M, M, M, M, M, 3, 10, M, M, 10),
                (2000, "HI", "CARYA", "400", 10.0, 999.9000244140625, M, M, M, M, M, 2, 5, M, M, 5),
                (2000, "AB", "FAGR", "531", 10.0, 999.9000244140625, M, M, M, M, M, 8, 23, M, M, 23),
                (2000, "SK", "QUFA", "812", 10.0, 999.9000244140625, M, M, M, M, M, 1, 3, M, M, 3)], 1)
        end
    end
    @testset "BM econ_bm_rich (PCT valuation, SEV, IRR/RRR, specials, appreciation)" begin
        db = _econ_run(joinpath(_ECON_SCEN, "econ_bm_rich.key"), FVSjl.BlueMountains())
        if db === nothing
            @test_skip "econ_bm_rich fixture not available"
        else
            _check_rows(_econ_rows(db, "FVS_EconSummary", "WHERE Year <= 2010"), [
                (1990, 10, "NO ", 82.84915161132812, 15.0, 77.5325927734375, 12.166346549987793, -65.36624908447266,
                 -1.025473466143012e-05, 0.15691912174224854, -13.582658767700195, -267.1967468261719, M, M, 0, 0, 4.0, M),
                (2000, 20, "NO ", 182.7735595703125, 30.0, 136.45059204101562, 20.385494232177734, -116.06509399414062,
                 -1.025473466143012e-05, 0.14939835667610168, -5.430608749389648, -237.3234100341797, M, M, 0, 0, 4.0, M),
                (2010, 30, "NO ", 279.8040771484375, 45.0, 173.80516052246094, 25.938058853149414, -147.86709594726562,
                 -1.025473466143012e-05, 0.14923641085624695, -2.3896634578704834, -211.9961395263672, M, M, 0, 0, 4.0, M)], 2)
        end
    end
    @testset "BM econ_bm_logs (R6 FW2 per-log ECVOL grades)" begin
        db = _econ_run(joinpath(_ECON_SCEN, "econ_bm_logs.key"), FVSjl.BlueMountains())
        if db === nothing
            @test_skip "econ_bm_logs fixture not available"
        else
            _check_rows(_econ_rows(db, "FVS_EconSummary", "WHERE Year = 1990"), [
                (1990, 10, "NO ", 30.0, 0.0, 30.0, 0.0, -30.0, M, 0.0, M, M, M, M, 0, 1149, 0.0, M)], 2)
            _check_rows(_econ_rows(db, "FVS_EconHarvestValue", "WHERE Year = 1990"), [
                (1990, "WL", "LAOC", "073", 4.0, 10.0, M, M, M, M, M, M, M, 458, 9, 9),
                (1990, "DF", "PSME", "202", 6.0, 10.0, M, M, M, M, M, M, M, 308, 15, 15),
                (1990, "LP", "PICO", "108", 4.0, 10.0, M, M, M, M, M, M, M, 383, 8, 8)], 1)
        end
    end
end
