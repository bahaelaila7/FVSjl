# =============================================================================
# test_ontario_db_path_equiv.jl — ON DATABASE-read path validated against the inline `.tre` path
# on a REAL Ontario FVS input DB (FVSDataHardwood.db, stand LD3001: 94 metric trees, ON variant).
#
# WHY inline-equivalence and not a live DB oracle: the ON oracle (FVSon_g16, gfortran-16/gcc-16
# rebuild of canada/on) SEGFAULTS on the SQLite tree-read — the crash is __memset_avx2 in
# dbstreesin_'s prologue, reproducible with fresh gfortran-16 objects, the shipped Jun-4
# dbstreesin.o, and a fully-static link; the stand-read (dbsstandsin, 86 columns) succeeds, so it
# is specifically the tree-read that faults. gcc-15/gfortran-15 (the known workaround for this
# toolchain-skew class) is unavailable in this environment and not installable (no apt candidate),
# and the isoc23 sscanf shim is unrelated to the memset fault. The archived Hardwood.sum.save is
# from a DIFFERENT (non-FVSon_g16) FVS build — it reports MCuFt=178 where FVSon_g16's faithful
# metric-quirk behaviour (and FVSjl) give MCuFt=0 — so it is NOT a valid FVSon_g16 oracle. The
# ON port is therefore validated via the inline path; this test closes the "DB reader validated
# only through the synthetic 3-tree case" gap on the real 94-tree stand.
#
# The DB reader (fia_database.jl apply_fia_trees!) and the inline `.tre` loader (treeinput.jl
# load_trees!) both route records through the SAME ingest_tree_records!(...; metric=…), so they
# MUST produce byte-identical tree state. The ONLY DB-specific tree transform is the per-hectare→
# per-acre TREE_COUNT rescale (×ACRtoHA = 0.40468564). This test builds an inline stand from the
# SAME DB rows and asserts:
#   • species / DBH / height / tree-id / plot-index arrays are bit-exact between the two paths;
#   • the DB path's ingested (pre-notre) TPA == inline TPA × 0.40468564, bit-exact (Float32);
#   • the stand attributes the DB reader sets (age/slope/aspect/elevation/latitude) match the DB;
#   • the DB path runs end-to-end and emits a cyc0 `.sum` row.
# =============================================================================

using Test, FVSjl
const F = FVSjl
using .FVSjl.SQLite
using .FVSjl.DBInterface

const _ONDB_HA2AC = 0.40468564f0
_u32(x) = reinterpret(UInt32, x)

@testset "ON — DATABASE-read path == inline `.tre` path (real DB, stand LD3001)" begin
    dbf = joinpath(@__DIR__, "..", "fixtures", "ontario", "FVSDataHardwood.db")
    if !isfile(dbf)
        @warn "FVSDataHardwood.db fixture missing — skipping ON DB-path equivalence (fresh-worktree provisioning)"
        @test_skip false
    else
        dir = mktempdir()

        # --- DB-read keyfile (the path under test) ---
        dbkey = joinpath(dir, "ld3001_db.key")
        write(dbkey, """STDIDENT
Hardwood
DATABASE
DSNin
$(abspath(dbf))
StandSQL
SELECT * FROM FVS_StandInit WHERE Stand_ID = 'LD3001'
EndSQL
TreeSQL
SELECT * FROM FVS_TreeInit WHERE Stand_ID = 'LD3001'
EndSQL
End
TIMEINT                      5
NUMCYCLE          10
PROCESS
STOP
""")

        # --- Build an INLINE `.tre`(.csv sibling) from the SAME raw DB rows (independent of the
        #     DB reader): named-column CSV carries raw metric values (cm DBH / m Ht) + alpha species,
        #     so the inline loader applies the identical metric conversion via ingest_tree_records!. ---
        base = joinpath(dir, "ld3001_inl")
        db = SQLite.DB(abspath(dbf))
        open(base * ".csv", "w") do io
            println(io, join(F.TREE_CSV_HEADER, ','))
            for r in DBInterface.execute(db, "SELECT * FROM FVS_TreeInit WHERE Stand_ID = 'LD3001'")
                g(k, d) = (v = getproperty(r, k); v === missing ? d : v)
                row = [Int(g(:Plot_ID, 1)), Int(g(:Tree_ID, 0)), Float64(g(:Tree_Count, 1.0)),
                       Int(g(:History, 1)), String(g(:Species, "OT")), Float64(g(:DBH, 0.0)), 0,
                       Float64(g(:Ht, 0.0)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                println(io, join(row, ','))
            end
        end
        SQLite.close(db)
        touch(base * ".tre")                       # TREEDATA resolves base*".tre"; the .csv sibling is read
        inlkey = base * ".key"
        write(inlkey, "STDIDENT\nLD3001\nNUMCYCLE           1\nTREEDATA\nPROCESS\nSTOP\n")

        # Capture the ingested tree state (pre-notre: notre! applies the stand DESIGN expansion, which
        # legitimately differs — DB DESIGN vs inline default — so equivalence is asserted at ingest).
        grab(s) = (n = s.trees.n,
                   sp = copy(s.trees.species[1:s.trees.n]), dbh = copy(s.trees.dbh[1:s.trees.n]),
                   ht = copy(s.trees.height[1:s.trees.n]),  tpa = copy(s.trees.tpa[1:s.trees.n]),
                   plot = copy(s.trees.plot_id[1:s.trees.n]), id = copy(s.trees.tree_id[1:s.trees.n]))

        sdb = first(F.each_stand(dbkey; variant = F.Ontario())); gdb = grab(sdb)
        sin = first(F.each_stand(inlkey; variant = F.Ontario())); gin = grab(sin)

        # 94 raw trees → 91 live (3 dead/non-stockable partitioned out), identically on both paths.
        @test gdb.n == gin.n
        @test gdb.n == 91
        @test gdb.sp   == gin.sp                                             # species indices bit-exact
        @test all(_u32.(gdb.dbh) .== _u32.(gin.dbh))                        # cm→in DBH bit-exact
        @test all(_u32.(gdb.ht)  .== _u32.(gin.ht))                         # m→ft height bit-exact
        @test gdb.id   == gin.id                                            # tree ids
        @test gdb.plot == gin.plot                                          # internal plot index (IPVEC)
        # DB TREE_COUNT is per-hectare; the DB reader pre-scales to per-acre by ACRtoHA. The inline
        # `.tre` PROB is taken raw, so DB-ingested TPA == inline TPA × 0.40468564, bit-exact.
        @test all(_u32.(gdb.tpa) .== _u32.(gin.tpa .* _ONDB_HA2AC))

        # --- Stand attributes the DB reader sets directly (apply_fia_stand!): match the raw DB row. ---
        p = sdb.plot
        @test p.stand_age == Int32(100)                                    # Age
        @test _u32(p.slope)  == _u32(5f0 / 100f0)                          # Slope 5% → 0.05
        @test _u32(p.aspect) == _u32(0f0)                                  # Aspect 0
        @test _u32(p.latitude) == _u32(49f0)                              # Latitude
        # Elevation: DB stores metres (518) with a BLANK ElevFt string; the reader must fall through
        # to ELEVATION and convert m→hundreds-of-ft (518 · 3.280839895 / 100), NOT read the blank
        # ElevFt as 0. (Regression guard for the blank-ELEVFT gate fix.)
        @test _u32(p.elevation) == _u32(518f0 * 3.280839895f0 / 100f0)

        # --- The DB-read path runs end-to-end and emits a cyc0 `.sum` row (oracle DB-path blocked). ---
        sumtext = F.run_keyfile(dbkey; variant = F.Ontario())
        rows = filter(l -> startswith(l, "2004"), split(sumtext, '\n'))
        @test length(rows) == 1                                            # the LD3001 inventory row
        @test occursin("100", rows[1])                                     # AGE 100 rendered
    end
end
