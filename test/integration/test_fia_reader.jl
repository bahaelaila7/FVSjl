# test_fia_reader.jl — native FIA "FVS-ready" database INPUT path (src/io/fia_database.jl).
#
# Self-contained: a tiny extracted SQLite fixture (test/fixtures/fia/ls_sample.db, two real
# LS stands) + the golden live-FVS .sum captured from that same fixture. Proves FVSjl's
# DATABASE/DSNIN reader ingests an FIA stand and reproduces live's cycle-0 inventory
# BIT-EXACT on the core columns (TPA/BA/SDI/CCF/TopHt/QMD/CuFt) AND BdFt across all four
# variants — including the R9 per-national-forest board-type gate (Scribner vs International ¼").

using Test
using FVSjl
using .FVSjl.SQLite
using .FVSjl.DBInterface

const FIA_DIR = joinpath(@__DIR__, "..", "fixtures", "fia")
const FIA_DB  = abspath(joinpath(FIA_DIR, "ls_sample.db"))

# Parse the first (cycle-0) data row of a .sum into its whitespace fields.
_sum_cycle0(text) = for ln in split(text, '\n')
    f = split(strip(ln)); length(f) >= 8 || continue
    y = tryparse(Int, f[1]); (y === nothing || y < 1000) && continue
    return f
end

_fia_keyfile(cn) = """
STDIDENT
$cn
DATABASE
DSNin
$FIA_DB
StandSQL
SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
TreeSQL
SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
END
NUMCYCLE         1.0
ECHOSUM
PROCESS
STOP
"""

@testset "FIA-ready database reader (DATABASE/DSNIN)" begin
    # Cross-variant: one+ real FIA stand per variant, cycle-0 vs the golden live .sum.
    # Density/size (TPA/BA/SDI/TopHt/QMD) + cubic volume + CCF are BIT-EXACT on every
    # stand/variant — CCF once species resolve (the _fia_spcode 3-digit-FIA fix) AND the
    # stand's lat/long feed the Hopkins index (HI-dependent hardwood crowns). BdFt is
    # bit-exact except the LS hardwood stand, which keeps a small hardwood board-vol residual.
    density = [(3, "TPA"), (4, "BA"), (5, "SDI"), (6, "CCF"), (7, "TopHt"), (8, "QMD"),
               (9, "TotCuFt"), (10, "MerchCuFt")]
    # (stand_cn, variant) — BdFt is BIT-EXACT on every stand once the R9 per-national-forest
    # board-type gate (volinit.f:434-451; _R9_INTL_BDFT_FORESTS) is honored: LS conifer
    # IFORST=10 → Scribner, LS hardwood IFORST=24 → International, NE/CS/SN → their native board.
    cases = [("100180735010661", FVSjl.LakeStates()),    # LS conifer (IFORST=10 → Scribner)
             ("255262523010854", FVSjl.Southern()),      # SN
             ("14173137020004",  FVSjl.CentralStates()), # CS
             ("657546100126144", FVSjl.Northeast()),     # NE (CCF fixed by lat/long)
             ("55482390010661",  FVSjl.LakeStates())]    # LS hardwood (IFORST=24 → R9 International board)
    for (cn, var) in cases
        golden = read(joinpath(FIA_DIR, "ls_$(cn).live.sum"), String)
        lv = _sum_cycle0(golden)
        key = joinpath(mktempdir(), "fia_$(cn).key")
        write(key, _fia_keyfile(cn))
        out = FVSjl.run_keyfile(key; variant = var)
        jl = _sum_cycle0(out)
        @test jl !== nothing && lv !== nothing
        @test jl[1] == lv[1]                       # inventory year
        for (i, name) in density
            @test jl[i] == lv[i]                   # BIT-EXACT rendered inventory stat (incl CCF)
        end
        @test jl[12] == lv[12]                     # BdFt BIT-EXACT (per-forest board type)
    end
end

@testset "FIA SITE_SPECIES 2-digit FIA code → SITSET interpolates (not fill-all)" begin
    # Regression: a 2-digit numeric FIA SITE_SPECIES ("93" = Engelmann spruce) was zero-padded
    # to "093" by _fia_spcode but the variant coef stores the FIA code UNPADDED ("93"), so the
    # STRICT site-species match at apply_fia_stand! missed it ⇒ ISISP fell to 0 ⇒ SITE_INDEX
    # was filled to ALL species, freezing SITSET's per-species interpolation. On a dense TT
    # (Teton) woodland aspen stand this raised aspen SITEAR to the site index (86) vs live's
    # interpolated 60.67, driving RSIMOD 0.9 vs 0.719 ⇒ TT small-tree aspen height/DBH over-grew
    # from cycle 1 (worst dig CN 1856529552290487: BA 150 vs live 91 by 2072). The match now
    # normalizes leading zeros like resolve_species does for trees.
    dir = mktempdir(); dbf = joinpath(dir, "tt_site.db")
    db = SQLite.DB(dbf)
    DBInterface.execute(db, """CREATE TABLE FVS_STANDINIT_COND (STAND_CN TEXT, VARIANT TEXT,
        INV_YEAR INT, SITE_SPECIES TEXT, SITE_INDEX REAL, AGE INT, SLOPE REAL, ASPECT REAL,
        ELEVATION REAL, LATITUDE REAL, LONGITUDE REAL, ECOREGION TEXT, LOCATION INT)""")
    DBInterface.execute(db, "INSERT INTO FVS_STANDINIT_COND VALUES " *
        "('T1','TT',2022,'93',86.0,50,20,180,65,44.0,-110.0,'M331Df',415)")
    DBInterface.execute(db, """CREATE TABLE FVS_TREEINIT_COND (STAND_CN TEXT, PLOT_ID INT,
        SPECIES TEXT, DIAMETER REAL, HT REAL, TREE_COUNT REAL, CRRATIO REAL)""")
    DBInterface.execute(db, "INSERT INTO FVS_TREEINIT_COND VALUES ('T1',1,'746',0.1,NULL,50.0,NULL)")
    DBInterface.execute(db, "INSERT INTO FVS_TREEINIT_COND VALUES ('T1',1,'93',10.0,60.0,5.0,55.0)")
    SQLite.close(db)
    key = joinpath(dir, "t.key")
    write(key, "STDIDENT\nT1\nDATABASE\nDSNin\n$dbf\nStandSQL\n" *
        "SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\n" *
        "SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\n" *
        "NUMCYCLE 1.0\nECHOSUM\nPROCESS\nSTOP\n")
    s, _ = FVSjl.initialize(key; variant = FVSjl.Teton())
    @test Int(s.plot.site_species) == 8                      # FIA 93 → TT Engelmann spruce (idx 8), not 0/fill-all
    # aspen (idx 6) SITEAR is INTERPOLATED from the site species SI (SITSET), not the raw 86:
    #   SITELO[6]=30 + (86−40)/(100−40)·(70−30) = 60.667 (bit-matches live FVStt SITEAR(6)).
    @test isapprox(Float64(s.plot.sp_site_index[6]), 60.667; atol = 0.05)
    @test s.plot.sp_site_index[6] < 85f0                     # decisively NOT the fill-all value 86
end
