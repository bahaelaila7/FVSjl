# =============================================================================
# test_ontario_database_reader.jl — ON DATABASE (SQLite FVS_TreeInit) reader metric conversion.
#
# ON is a metric variant: its FVS_TreeInit DB stores DBH in cm, Ht in m, Tree_Count per hectare (exactly
# like BC). The DATABASE reader (fia_database.jl) originally gated its metric cm→in / m→ft conversion and
# per-hectare→per-acre expansion on `s.variant isa BritishColumbia` ONLY — so an Ontario DB was ingested
# as if cm were inches (~2.5× too large), the same class of bug as the BC DB-reader fix (42eb555). The
# gate is now `BritishColumbia || Ontario`.
#
# This builds a tiny in-memory-style Ontario DB (3 trees at 25 cm / 18 m) and asserts jl ingests DBH =
# 9.842 in (= 25 cm · 0.3937) and Ht = 59.06 ft (= 18 m · 3.28084) — the exact internal values the
# oracle-validated inline `.tre` path produces (the DGF dump-replay hex 411D7AE1 = 9.842). Validated
# 2026-08-20. The inline path was already correct; this covers the DB path.
# =============================================================================

using Test, FVSjl
const F = FVSjl
using .FVSjl.SQLite

@testset "ON — DATABASE reader metric conversion (cm→in / m→ft)" begin
    dir = mktempdir()
    dbf = joinpath(dir, "on_test.db")
    db = SQLite.DB(dbf)
    SQLite.execute(db, """CREATE TABLE FVS_StandInit (Stand_CN TEXT, Stand_ID TEXT, Variant TEXT,
        Inv_Year INT, Latitude REAL, Longitude REAL, Region INT, Forest INT, Location INT, Age INT,
        Aspect REAL, Slope REAL, Elevation REAL, Basal_Area_Factor REAL, Inv_Plot_Size REAL, Brk_DBH REAL,
        Num_Plots INT, Sam_Wt REAL, Site_Species TEXT, Site_Index REAL)""")
    SQLite.execute(db, """INSERT INTO FVS_StandInit VALUES ('ONDB1','ONDB1','ON',2004, 48.0,-89.0,10,10,
        1010, 15,180.0,30.0,300.0, 1.0,0.0,0.0, 11, 1.0,'PW',15.0)""")
    SQLite.execute(db, """CREATE TABLE FVS_TreeInit (Stand_CN TEXT, Stand_ID TEXT, Plot_ID INT, Tree_ID INT,
        Tree_Count REAL, History INT, Species TEXT, DBH REAL, Ht REAL, CrRatio INT)""")
    for (i, sp) in enumerate(("PW", "SW", "MH"))
        SQLite.execute(db, "INSERT INTO FVS_TreeInit VALUES ('ONDB1','ONDB1',1,$i,50.0,0,'$sp',25.0,18.0,40)")
    end
    close(db)

    key = joinpath(dir, "ondb1.key")
    write(key, """STDIDENT
ONDB1
STDINFO          301       1010      15.0     180.0     30.0      300.0
INVYEAR       2004.0
NUMCYCLE           1.0
DESIGN                                        11.0       1.0
DATABASE
DSNIN
$dbf
StandSQL
SELECT * FROM FVS_StandInit WHERE Stand_ID = 'ONDB1'
EndSQL
TreeSQL
SELECT * FROM FVS_TreeInit WHERE Stand_ID = 'ONDB1'
EndSQL
END
PROCESS
STOP
""")

    got = nothing
    for s in F.each_stand(key; variant = F.Ontario())
        F.notre!(s)
        got = (n = s.trees.n, dbh = copy(s.trees.dbh[1:s.trees.n]),
               ht = copy(s.trees.height[1:s.trees.n]), sp = copy(s.trees.species[1:s.trees.n]))
        break
    end
    @test got !== nothing && got.n == 3
    # DBH 25 cm → 9.842 in (metric-converted, not 25 in); Ht 18 m → 59.06 ft
    @test all(reinterpret(UInt32, d) == reinterpret(UInt32, 25f0 * 0.3937f0) for d in got.dbh)
    @test all(isapprox(h, 18f0 * 3.28084f0; atol = 0.01f0) for h in got.ht)
    @test got.sp == Int32[5, 6, 26]      # PW/SW/MH resolve to ON species indices
end
