# FVS_DM_Stnd_Sum + FVS_DM_Spp_Sum dwarf-mistletoe summary DBS tables (dbs/dbsmis.f DBSMIS2/DBSMIS1,
# via the misprt.f MISPRT aggregation ported in engine/mistletoe_report.jl), gated by the MISRPTS
# database keyword. (2026-08-22.)
#
# VALIDATED BIT-EXACT vs the relinked FVSie_clean / FVScr_clean .out mistletoe tables on a controlled
# 6-tree stand (3 DF infected DMR 3/4/2, 2 DF uninfected, 1 PP uninfected; damage code 33 = DF
# mistletoe). The cyc0 (1990, inventory) row — no growth projection, so a clean deterministic A/B —
# matches every column of BOTH oracle tables. CR carries its own volume equations + cr_dm_mortality_rate
# (VOL 14715 / INF-VOL 7988 / MORT-VOL 387 vs IE's 14666 / 7854 / 382), so the two variants independently
# exercise the report's volume + per-tree-DM-mortality dispatch. Multi-cycle rows inherit the projection
# straddle (the oracle books AUTOES regen cohorts) — cornered, the report mirrors jl's own tree list.

using Test
using FVSjl
using FVSjl: run_keyfile, InlandEmpire, CentralRockies
using SQLite, DBInterface

# Oracle cyc0 (1990) values per variant: (stand NamedTuple, species NamedTuple). From FVS{ie,cr}_clean .out.
const _DM_ORACLE = Dict(
    "IE" => (variant = InlandEmpire(),
        st = (Age=60, Stnd_TPA=1729, Stnd_BA=660, Stnd_Vol=14666, Inf_TPA=1110, Inf_BA=353,
              Inf_Vol=7854, Mort_TPA=56, Mort_BA=17, Mort_Vol=382, Inf_TPA_Pct=64, Inf_Vol_Pct=54,
              Mort_TPA_Pct=3, Mort_Vol_Pct=3, Mean_DMR=2.0, Mean_DMI=3.1),
        sp = (Spp="DF", Mean_DMR=2.4, Mean_DMI=3.1, Inf_TPA=1110, Mort_TPA=56,
              Inf_TPA_Pct=75, Mort_TPA_Pct=4, Stnd_TPA_Pct=85)),
    "CR" => (variant = CentralRockies(),
        st = (Age=60, Stnd_TPA=1729, Stnd_BA=660, Stnd_Vol=14715, Inf_TPA=1110, Inf_BA=358,
              Inf_Vol=7988, Mort_TPA=56, Mort_BA=17, Mort_Vol=387, Inf_TPA_Pct=64, Inf_Vol_Pct=54,
              Mort_TPA_Pct=3, Mort_Vol_Pct=3, Mean_DMR=2.0, Mean_DMI=3.1),
        sp = (Spp="DF", Mean_DMR=2.4, Mean_DMI=3.1, Inf_TPA=1110, Mort_TPA=56,
              Inf_TPA_Pct=75, Mort_TPA_Pct=4, Stnd_TPA_Pct=85)),
)

function _run_dm(variant, dir)
    tre_src = joinpath(@__DIR__, "..", "fixtures", "dm", "dm6.tre")
    cp(tre_src, joinpath(dir, "dm6.tre"))          # tree file read by keyfile basename (dm6.key → dm6.tre)
    dbp = joinpath(dir, "dm6.db"); key = joinpath(dir, "dm6.key")
    open(key, "w") do io
        print(io, """
STDIDENT
DM6 DWARF MISTLETOE TEST
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
    run_keyfile(key; variant = variant)
    return dbp
end

@testset "FVS_DM_* dwarf-mistletoe summary DBS tables (MISRPTS)" begin
    for (vc, spec) in _DM_ORACLE
        @testset "$vc cyc0 bit-exact vs FVS$(lowercase(vc))_clean" begin
            mktempdir() do dir
                dbp = _run_dm(spec.variant, dir)
                @test isfile(dbp)
                db = SQLite.DB(dbp)
                try
                    tabs = [r.name for r in DBInterface.execute(db,
                        "SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'FVS_DM%'")]
                    @test "FVS_DM_Stnd_Sum" in tabs
                    @test "FVS_DM_Spp_Sum" in tabs

                    st = first(DBInterface.execute(db, "SELECT * FROM FVS_DM_Stnd_Sum WHERE Year=1990"))
                    for (col, val) in pairs(spec.st)
                        got = st[col]
                        if col in (:Mean_DMR, :Mean_DMI)
                            @test round(got, digits = 1) == val
                        else
                            @test got == val
                        end
                    end

                    nsp = first(DBInterface.execute(db,
                        "SELECT COUNT(*) c FROM FVS_DM_Spp_Sum WHERE Year=1990")).c
                    @test nsp == 1
                    sp = first(DBInterface.execute(db, "SELECT * FROM FVS_DM_Spp_Sum WHERE Year=1990"))
                    @test strip(sp.Spp) == spec.sp.Spp
                    for col in (:Mean_DMR, :Mean_DMI)
                        @test round(sp[col], digits = 1) == spec.sp[col]
                    end
                    for col in (:Inf_TPA, :Mort_TPA, :Inf_TPA_Pct, :Mort_TPA_Pct, :Stnd_TPA_Pct)
                        @test sp[col] == spec.sp[col]
                    end

                    # inventory + 4 cycles = 5 reported rows
                    @test first(DBInterface.execute(db, "SELECT COUNT(*) c FROM FVS_DM_Stnd_Sum")).c == 5
                finally
                    SQLite.close(db)
                end
            end
        end
    end
end
