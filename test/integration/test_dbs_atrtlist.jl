# test_dbs_atrtlist.jl — DBS FVS_ATRTList (dbsatrtls.f) + DSNOUT redefinition (dbsin.f:116-122).
#
# FVS_ATRTList: written only when DATABASE ATRTLIDB set IATRTLIST>0 AND an ATRTLIST activity (initre.f opt 135,
# code 198) is due in the cycle of an applied cut — cuts.f:1740 CALL PRTRLS(3) right after the CUTLIST's PRTRLS(2).
# Rows = every record with PROB>0 AFTER the thin (dbsatrtls.f `IF (P.LE.0.0) CYCLE`), species-major IND1 order,
# TPA = PROB/GROSPC, MortPA = 0. Live A/B (FVSbm_g16.new, 7 BM FIA stands × 7 gating scenarios, 2026-09-19): table
# presence identical in every case; rows, order and every column identical except Ht2TDCF/Ht2TDBF (merch-top heights
# not yet computed by jl, 0 in every list table).
#
# DSNOUT redefinition: DSNOUT (dbsin.f:127-173) and a stand read (dbsstandin.f:216) both call DBSCASE(1), whose
# IFORSURE=1 forces the case ⇒ CASEID assigned (opening DSNOUT, default 'FVSOut.db'). A later DSNOUT then draws FVS16
# "DSNOUT DATA BASE CAN NOT BE REDEFINED" and output stays where it is. Live FVSbm (2026-09-19) on 6 block orders:
#   input block with SUMMARY, later DSNOUT b.db      ⇒ FVSOut.db      input block, later DSNOUT b.db + SUMMARY ⇒ FVSOut.db
#   DSNOUT a.db + SUMMARY + input, later DSNOUT b.db ⇒ a.db           DSNOUT a.db + input, later DSNOUT b.db+SUMMARY ⇒ a.db
#   DSNOUT b.db + SUMMARY block first, then input    ⇒ b.db           input block, later SUMMARY + DSNOUT b.db ⇒ FVSOut.db

using Test, FVSjl, SQLite, DBInterface

_tables(db) = isfile(db) ? (d = SQLite.DB(db); try [r.name for r in SQLite.tables(d)] finally SQLite.close(d) end) : String[]

@testset "DBS FVS_ATRTList + DSNOUT redefinition" begin
    tre = joinpath(@__DIR__, "..", "harness", "scenarios", "dbs_compute.tre")
    if !isfile(tre)
        @test_skip "dbs scenario not available"
    else
        thin = rpad("THINBTA", 10) * lpad("1995", 10) * lpad("80", 10)   # residual BA 80 @ 1995
        function key(dir, extra, dbkw)
            db = joinpath(dir, "out.db")
            cp(tre, joinpath(dir, "a.tre"); force = true)
            k = joinpath(dir, "a.key")
            write(k, "STDIDENT\nATRTDB\nSTDINFO        80106   231Dd        60.0     315.0      30.0       7.0\n" *
                  "INVYEAR       1990.0\nNUMCYCLE         3.0\nSITECODE          63      60.\n" *
                  "DESIGN                                        11.0       1.0\n$thin\n" *
                  join(extra, "\n") * (isempty(extra) ? "" : "\n") *
                  "TREEFMT\n(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,\n" *
                  "T52,I2,T66,5I1,T54,7I1,T75,F3.0)\nDATABASE\nDSNOUT\n$db\nSUMMARY\n$(dbkw)END\nTREEDATA\nPROCESS\nSTOP\n")
            return k, db
        end
        # positive: ATRTLIST + CUTLIST together ⇒ after-treatment + removed lists partition the pre-thin stand
        dir = mktempdir()
        k, db = key(dir, [rpad("ATRTLIST", 10) * lpad("1995", 10), rpad("CUTLIST", 10) * lpad("1995", 10)], "ATRTLIDB\nCUTLIDB\n")
        sumtxt = FVSjl.run_keyfile(k; faithful = true)
        r = nothing
        for ln in split(sumtxt, "\n"); f = split(ln); length(f) >= 18 && f[1] == "1995" && (r = f); end
        @test r !== nothing
        tpa = parse(Int, r[3]); atba = parse(Int, r[18])
        d = SQLite.DB(db)
        try
            @test "FVS_ATRTList" in [t.name for t in SQLite.tables(d)]
            at = [NamedTuple(x) for x in DBInterface.execute(d, "SELECT TPA,DBH,MortPA,SpeciesFVS FROM FVS_ATRTList WHERE Year=1995 ORDER BY rowid")]
            ct = [NamedTuple(x) for x in DBInterface.execute(d, "SELECT TPA FROM FVS_CutList WHERE Year=1995")]
            @test !isempty(at) && all(x -> x.TPA > 0 && x.MortPA == 0, at)
            @test round(Int, sum(0.005454154 * x.TPA * x.DBH^2 for x in at)) == atba          # Σ residual BA = .sum after-thin BA
            @test round(Int, sum(x.TPA for x in at) + sum(x.TPA for x in ct)) == tpa           # residual + removed = pre-thin TPA
            sp = [x.SpeciesFVS for x in at]; @test issorted(sp; by = s -> findfirst(==(s), unique(sp)))   # species-major blocks
        finally
            SQLite.close(d)
        end
        # negatives (live-verified gating): no ATRTLIDB / no ATRTLIST / ATRTLIST in a cut-free cycle ⇒ no table
        for (extra, dbkw) in (([rpad("ATRTLIST", 10) * lpad("0", 10)], ""), (String[], "ATRTLIDB\n"),
                              ([rpad("ATRTLIST", 10) * lpad("1990", 10)], "ATRTLIDB\n"))
            gd = mktempdir(); gk, gdb = key(gd, extra, dbkw)
            FVSjl.run_keyfile(gk; faithful = true)
            @test !("FVS_ATRTList" in _tables(gdb))
        end
    end

    # DSNOUT redefinition (live-measured, see header). Input = the tracked BM seedling fixture DB.
    fdb = joinpath(@__DIR__, "..", "fixtures", "bm_seedling", "bm_seedling.db")
    if !isfile(fdb)
        @test_skip "bm_seedling fixture db missing"
    else
        inblk(pre) = "DATABASE\n" * pre * "DSNin\n$(abspath(fdb))\nStandSQL\n" *
                     "SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\n" *
                     "SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\n"
        cases = (("sum-then-dsnout",   inblk("SUMMARY\n") * "DATABASE\nDSNOUT\nb.db\nSUMMARY\nEND\n",       "FVSOut.db"),
                 ("a-then-b",          inblk("DSNOUT\na.db\nSUMMARY\n") * "DATABASE\nDSNOUT\nb.db\nEND\n",  "a.db"),
                 ("noflag-then-b",     inblk("") * "DATABASE\nDSNOUT\nb.db\nSUMMARY\nEND\n",                "FVSOut.db"),
                 ("a-noflag-then-b",   inblk("DSNOUT\na.db\n") * "DATABASE\nDSNOUT\nb.db\nSUMMARY\nEND\n",  "a.db"),
                 ("b-first",           "DATABASE\nDSNOUT\nb.db\nSUMMARY\nEND\n" * inblk(""),                "b.db"),
                 ("sum-in-after",      inblk("") * "DATABASE\nSUMMARY\nDSNOUT\nb.db\nEND\n",                "FVSOut.db"))
        for (nm, blk, want) in cases
            dir = mktempdir(); kf = joinpath(dir, "k.key")
            write(kf, "STDIDENT\n449747082489998\n" * blk * "NUMCYCLE         2.0\nECHOSUM\nPROCESS\nSTOP\n")
            cd(() -> FVSjl.run_keyfile(kf; variant = FVSjl.BlueMountains()), dir)
            got = [f for f in ("a.db", "b.db", "FVSOut.db") if "FVS_Summary" in _tables(joinpath(dir, f))]
            @test got == [want]
        end
    end
end
