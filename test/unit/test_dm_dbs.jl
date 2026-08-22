# FVS_DM_Stnd_Sum + FVS_DM_Spp_Sum dwarf-mistletoe summary DBS tables (dbs/dbsmis.f DBSMIS2/DBSMIS1,
# via the misprt.f MISPRT aggregation ported in engine/mistletoe_report.jl), gated by the MISRPTS
# database keyword. (2026-08-22.)
#
# VALIDATED BIT-EXACT vs the relinked FVSie_clean .out mistletoe tables on a controlled 6-tree IE stand
# (3 DF infected DMR 3/4/2, 2 DF uninfected, 1 PP uninfected; damage code 33 = DF mistletoe). The cyc0
# (1990, inventory) row — no growth projection, so a clean deterministic A/B — matches every column of
# BOTH oracle tables:
#   Species DF: MEAN DMR 2.4, MEAN DMI 3.1, INF-TPA 1110, MORT-TPA 56, %INF 75, %MORT 4, COMPOSITION 85
#   Stand:      AGE 60, TREES 1729, BA 660, VOL 14666, INF 1110/353/7854, MORT 56/17/382,
#               %TPA-inf 64, %VOL-inf 54, %TPA-mort 3, %VOL-mort 3, MEAN DMR 2.0, MEAN DMI 3.1
# The DM MORTALITY columns exercise ie_dm_mortality_rate (mismrt.f, fint=10). Multi-cycle rows inherit
# the IE establishment (ESRANN)/OLDRN projection straddle (the oracle books AUTOES DF regen cohorts) —
# the report faithfully mirrors jl's own tree list, cornered like every IE multi-cycle validation.

using Test
using FVSjl
using FVSjl: run_keyfile, InlandEmpire
using SQLite, DBInterface

@testset "FVS_DM_* dwarf-mistletoe summary DBS tables (MISRPTS)" begin
    tre_src = joinpath(@__DIR__, "..", "fixtures", "dm", "dm6.tre")
    mktempdir() do dir
        cp(tre_src, joinpath(dir, "dm6.tre"))       # tree file read by keyfile basename (dm6.key → dm6.tre)
        dbp = joinpath(dir, "dm6.db")
        key = joinpath(dir, "dm6.key")
        open(key, "w") do io
            print(io, """
STDIDENT
DM6 IE DWARF MISTLETOE TEST
DESIGN                                        11.0       1.0
STDINFO          303    001010      60.0     315.0      30.0      88.0
INVYEAR         1990
NUMCYCLE           4
TREEFMT
(I4,T1,I7,F6.0,I1,A3,F4.1,F3.1,2F3.0,F4.1,I1,3(I2,I2),2I1,I2,2I3,2I1,F3.0)
OPEN            55
$(joinpath(dir, "dm6.tre"))
TREEDATA          55
DATABASE
DSNOUT
$dbp
SUMMARY
MISRPTS
END
PROCESS
STOP
""")
        end
        run_keyfile(key; variant = InlandEmpire())

        @test isfile(dbp)
        db = SQLite.DB(dbp)
        try
            tabs = [r.name for r in DBInterface.execute(db,
                "SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'FVS_DM%'")]
            @test "FVS_DM_Stnd_Sum" in tabs
            @test "FVS_DM_Spp_Sum" in tabs

            # --- Stand composite table (DBSMIS2), cyc0 = 1990, bit-exact vs FVSie_clean ---
            st = first(DBInterface.execute(db, "SELECT * FROM FVS_DM_Stnd_Sum WHERE Year=1990"))
            @test st.Age          == 60
            @test st.Stnd_TPA     == 1729
            @test st.Stnd_BA      == 660
            @test st.Stnd_Vol     == 14666
            @test st.Inf_TPA      == 1110
            @test st.Inf_BA       == 353
            @test st.Inf_Vol      == 7854
            @test st.Mort_TPA     == 56
            @test st.Mort_BA      == 17
            @test st.Mort_Vol     == 382
            @test st.Inf_TPA_Pct  == 64
            @test st.Inf_Vol_Pct  == 54
            @test st.Mort_TPA_Pct == 3
            @test st.Mort_Vol_Pct == 3
            @test round(st.Mean_DMR, digits = 1) == 2.0
            @test round(st.Mean_DMI, digits = 1) == 3.1

            # --- Top-4 species table (DBSMIS1), cyc0 = 1990; only DF is infected ---
            nsp = first(DBInterface.execute(db, "SELECT COUNT(*) c FROM FVS_DM_Spp_Sum WHERE Year=1990")).c
            @test nsp == 1
            df = first(DBInterface.execute(db, "SELECT * FROM FVS_DM_Spp_Sum WHERE Year=1990"))
            @test strip(df.Spp)               == "DF"
            @test round(df.Mean_DMR, digits=1) == 2.4
            @test round(df.Mean_DMI, digits=1) == 3.1
            @test df.Inf_TPA       == 1110
            @test df.Mort_TPA      == 56
            @test df.Inf_TPA_Pct   == 75
            @test df.Mort_TPA_Pct  == 4
            @test df.Stnd_TPA_Pct  == 85

            # DM report present each reported cycle (inventory + 4 cycles = 5 rows)
            n = first(DBInterface.execute(db, "SELECT COUNT(*) c FROM FVS_DM_Stnd_Sum")).c
            @test n == 5
        finally
            SQLite.close(db)
        end
    end
end
