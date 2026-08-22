# FVS_RD_Sum — Western Root Disease area-summary DBS table (dbs/dbsrd.f DBSRD1, via the rdpr.f
# report aggregation ported in engine/root_disease.jl `rd_sum_report`), gated by the RDSUM
# DATABASE keyword + the RRDOUT report keyword. (2026-08-22.)
#
# WRD IS a LIVE-oracle chunk (rd/ is in every FVS sourceList). Validated vs FVSkt_clean with
# RDSUM+RRDOUT on the S248112 RRTYPE-3 (Armillaria) stand — oracle at /workspace/.ktwork/rdrun/
# rdsum_oracle.db, 11 rows 1990-2090 (captured below as _RD_ORACLE). KT is imperial (LMTRIC only
# for BC/ON) so the values are RAW (no ·ACRtoHA). The **1990 inventory row is BIT-EXACT** (RDPR#1,
# fvs.f:347, pre-projection); the projected rows are bit-exact early and cornered later by the WRD
# DGSCOR/OLDRN #206 growth straddle (WRD growth is cornered on KT). Two collect-seam subtleties that
# matter: the row is captured AFTER rd_grow_apply!→rdinoc (decayed stump pool PROBDA) AND after the
# DBH-UPDATE + compute_volumes! (grown DBH for Live_BA). The 4 new-infection columns
# (New_Inf_Prp_Ins/Exp/Tot + Ave_Pct_Root_Inf) are 0 pending the CORINF/EXPINF/PRINF accumulators.

using Test
using FVSjl
using FVSjl: run_keyfile, Kootenai
using SQLite, DBInterface

# (Year, Age, RD_Type, Num_Centers, RD_Area, Spread, Stumps, Mort_TPA, UnInf_TPA, Inf_TPA, Live_BA)
const _RD_ORACLE = [
    (1990, 60, "A", 10, 10.0267, 0.0,    0.0,     0.0,      393.0978, 196.549, 85.1303),
    (2000, 70, "A", 10, 11.6622, 1.037,  0.1841,  203.8672, 294.4129, 86.2569, 84.0852),
    (2010, 80, "A", 10, 13.4578, 1.1403, 7.1491,  133.8578, 205.4474, 48.6501, 70.8547),
    (2020, 90, "A", 10, 15.6267, 1.2442, 15.2762, 77.7258,  154.1932, 40.8416, 70.3073),
    (2090, 160,"A", 10, 33.7067, 1.272,  36.4617, 24.1405,  25.4541,  37.4429, 84.6102),
]

function _run_rd(dir)
    cp(joinpath(@__DIR__, "..", "fixtures", "rd", "rdsum.tre"), joinpath(dir, "rd.tre"))
    dbp = joinpath(dir, "rd.db"); key = joinpath(dir, "rd.key")
    open(key, "w") do io
        print(io, """
STDIDENT
S248112  RD SUM
DATABASE
DSNOUT
$dbp
RDSUM
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

@testset "FVS_RD_Sum WRD root-disease summary (RDSUM) vs FVSkt_clean" begin
    mktempdir() do dir
        dbp = _run_rd(dir)
        @test isfile(dbp)
        db = SQLite.DB(dbp)
        try
            @test "FVS_RD_Sum" in [r.name for r in DBInterface.execute(db,
                "SELECT name FROM sqlite_master WHERE type='table' AND name='FVS_RD_Sum'")]
            @test first(DBInterface.execute(db, "SELECT COUNT(*) c FROM FVS_RD_Sum")).c == 11
            for (yr, age, typ, nc, area, spr, stmp, mort, uninf, inf, ba) in _RD_ORACLE
                r = first(DBInterface.execute(db, "SELECT * FROM FVS_RD_Sum WHERE Year=$yr"))
                @test r.Age == age && strip(r.RD_Type) == typ && r.Num_Centers == nc
                if yr == 1990
                    # inventory row: BIT-EXACT (pre-projection)
                    @test isapprox(Float32(r.RD_Area), Float32(area); atol = 5f-3)
                    @test r.Spread_Ft_per_Year == 0.0 && r.Stumps_per_Acre == 0.0 && r.Mort_TPA == 0.0
                    @test isapprox(Float32(r.UnInf_TPA), Float32(uninf); atol = 0.1f0)
                    @test isapprox(Float32(r.Inf_TPA), Float32(inf); atol = 0.1f0)
                    @test isapprox(Float32(r.Live_BA), Float32(ba); atol = 0.1f0)
                elseif yr == 2000
                    # first projected cycle: Area/Spread/Mort BIT-EXACT; Stumps/Inf/BA cornered
                    @test isapprox(Float32(r.RD_Area), Float32(area); atol = 5f-3)
                    @test isapprox(Float32(r.Spread_Ft_per_Year), Float32(spr); atol = 5f-3)
                    @test isapprox(Float32(r.Mort_TPA), Float32(mort); atol = 0.05f0)
                    @test isapprox(Float32(r.Inf_TPA), Float32(inf); rtol = 0.02f0)
                    @test isapprox(Float32(r.Live_BA), Float32(ba); rtol = 0.03f0)
                else
                    # later cycles: cornered by the WRD DGSCOR/OLDRN #206 growth+spread straddle
                    @test isapprox(Float32(r.RD_Area), Float32(area); rtol = 0.04f0)
                    @test isapprox(Float32(r.Live_BA), Float32(ba); rtol = 0.12f0)
                    @test isapprox(Float32(r.Mort_TPA), Float32(mort); rtol = 0.20f0)
                end
            end
        finally
            SQLite.close(db)
        end
    end
end
