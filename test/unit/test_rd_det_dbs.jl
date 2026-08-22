# FVS_RD_Det — Western Root Disease per-species patch detail DBS table (dbs/dbsrd.f DBSRD2,
# via the rd/rddout.f report ported in engine/root_disease.jl `rd_det_report`), gated by the
# RDDETAIL DATABASE keyword + the RRDOUT report keyword. (2026-08-22.)
#
# Per active disease type, one row per variant species with a tree record in the patch: the DBH
# at the 10/30/50/70/90/100 percentile points of the killed-tree and live(in-patch)-tree
# distributions (RDPSRT→PCTILE→RDDST), plus per-species Mort/UnInf/Inf TPA and mean %-roots-
# infected (PRINF·100). Validated vs FVSkt_clean with RDDETAIL on the S248112 Armillaria stand
# (oracle at /workspace/.ktwork/rdrun/rddet_oracle.db, 66 rows = 6 species × 11 years). KT is
# imperial (LMTRIC only for BC/ON) so DBH values are raw inches (no ·INTOCM).
#
# The **1990 inventory row is BIT-EXACT** across all 17 numeric columns for every species (RDPR#1,
# fvs.f:347, pre-projection) — this validates the percentile pipeline + PRINF + per-species
# accumulators. The 2000 cycle is bit-exact-or-cornered; later cycles are cornered by the WRD
# DGSCOR/OLDRN #206 growth+tree-list straddle (WRD growth is cornered on KT). The DBH-percentile
# columns amplify that straddle (a percentile point is a discrete tree DBH, so a small growth
# difference can jump to an adjacent tree), so only the aggregate TPA columns are checked past 1990.

using Test
using FVSjl
using FVSjl: run_keyfile, Kootenai
using SQLite, DBInterface

# 1990 inventory rows (BIT-EXACT): (SpeciesFVS, RD_Area, [Mort_10..100, Mort_TPA, Live_10..100, UnInf, Inf, Root%])
const _RDDET_1990 = [
    ("DF", 10.0267, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 1.2, 1.9, 9.4, 12.7, 149.9104, 57.0489, 9.0211]),
    ("ES", 10.0267, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 3.2, 3.2, 5.0, 5.8, 5.8, 63.5674, 47.5679, 7.4601]),
    ("GF", 10.0267, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 6.1, 6.6, 10.9, 120.3287, 57.575, 11.9515]),
    ("LP", 10.0267, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5, 8.5, 9.5, 9.6, 11.5, 11.5, 7.7439, 24.0358, 12.8451]),
    ("WH", 10.0267, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 15.0633, 2.2948, 8.0]),
    ("WL", 10.0267, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.9, 8.0, 8.0, 8.2, 8.4, 8.4, 36.4841, 8.0265, 9.6205]),
]

# 2000 first-projected-cycle aggregates (cornered): (SpeciesFVS, Mort_TPA, UnInf_TPA, Inf_TPA, Pct_Roots_Inf)
const _RDDET_2000_AGG = [
    ("DF", 68.906, 117.237, 20.212, 24.758),
    ("ES", 52.056, 47.075, 12.004, 40.693),
    ("GF", 74.088, 79.198, 24.618, 68.597),
    ("LP", 6.6, 6.573, 16.82, 19.174),
    ("WH", 2.216, 12.985, 2.156, 52.235),
    ("WL", 0.0, 31.346, 10.447, 39.569),
]

const _RDDET_COLS = ["Mort_10Pctile_DBH", "Mort_30Pctile_DBH", "Mort_50Pctile_DBH",
    "Mort_70Pctile_DBH", "Mort_90Pctile_DBH", "Mort_100Pctile_DBH", "Mort_TPA_Total",
    "Live_10Pctile_DBH", "Live_30Pctile_DBH", "Live_50Pctile_DBH", "Live_70Pctile_DBH",
    "Live_90Pctile_DBH", "Live_100Pctile_DBH", "UnInf_TPA_Total", "Inf_TPA_Total", "Pct_Roots_Inf"]

function _run_rddet(dir)
    cp(joinpath(@__DIR__, "..", "fixtures", "rd", "rdsum.tre"), joinpath(dir, "rd.tre"))
    dbp = joinpath(dir, "rd.db"); key = joinpath(dir, "rd.key")
    open(key, "w") do io
        print(io, """
STDIDENT
S248112  RD DET
DATABASE
DSNOUT
$dbp
RDDETAIL
END
SCREEN
NOAUTOES
NOTRIPLE
STATS
DESIGN                                        11.0       1.0
STDINFO     11406001     570.0      60.0     315.0      30.0      34.0
INVYEAR       1990.0
NUMCYCLE        10.0
TREEFMT
(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,
T52,I2,T66,5I1,T54,7I1,T75,F3.0)
TREEDATA
RDIN
RRTYPE             3
RRINIT             0        10        10        20       0.1        10         3
SAREA            100
RRDOUT
BBCLEAR
END
ECHOSUM
PROCESS
STOP
""")
    end
    run_keyfile(key; variant = Kootenai())
    return dbp
end

@testset "FVS_RD_Det WRD per-species patch detail (RDDETAIL) vs FVSkt_clean" begin
    mktempdir() do dir
        dbp = _run_rddet(dir)
        @test isfile(dbp)
        db = SQLite.DB(dbp)
        try
            @test "FVS_RD_Det" in [r.name for r in DBInterface.execute(db,
                "SELECT name FROM sqlite_master WHERE type='table' AND name='FVS_RD_Det'")]
            @test first(DBInterface.execute(db, "SELECT COUNT(*) c FROM FVS_RD_Det")).c == 66

            # 1990 inventory: BIT-EXACT across all 17 numeric columns for every species.
            for (sp, area, vals) in _RDDET_1990
                r = first(DBInterface.execute(db,
                    "SELECT * FROM FVS_RD_Det WHERE Year=1990 AND SpeciesFVS='$sp'"))
                @test isapprox(Float32(r.RD_Area), Float32(area); atol = 5f-3)
                for (c, want) in zip(_RDDET_COLS, vals)
                    got = Float32(getproperty(r, Symbol(c)))
                    @test isapprox(got, Float32(want); atol = 5f-3)
                end
            end

            # 2000 first projected cycle: aggregate TPA + root-% columns bit-exact-or-cornered.
            for (sp, mort, un, inf, root) in _RDDET_2000_AGG
                r = first(DBInterface.execute(db,
                    "SELECT * FROM FVS_RD_Det WHERE Year=2000 AND SpeciesFVS='$sp'"))
                @test isapprox(Float32(r.Mort_TPA_Total), Float32(mort); atol = 0.05f0, rtol = 0.03f0)
                @test isapprox(Float32(r.UnInf_TPA_Total), Float32(un); rtol = 0.03f0)
                @test isapprox(Float32(r.Inf_TPA_Total), Float32(inf); rtol = 0.05f0)
                @test isapprox(Float32(r.Pct_Roots_Inf), Float32(root); rtol = 0.03f0)
            end

            # every year emits exactly the 6 host species (structure preserved through the run).
            for yr in 1990:10:2090
                @test first(DBInterface.execute(db,
                    "SELECT COUNT(*) c FROM FVS_RD_Det WHERE Year=$yr")).c == 6
            end
        finally
            SQLite.close(db)
        end
    end
end
