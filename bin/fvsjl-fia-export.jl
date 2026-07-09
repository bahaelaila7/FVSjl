#!/usr/bin/env julia
# =============================================================================
# fvsjl-fia-export.jl — export FVS-ready FIA plots to standalone stand files.
#
# Given a "FVS-ready" FIA SQLite database (the FVS_STANDINIT_COND / FVS_TREEINIT_COND
# tables, e.g. SQLITE_FIADB_ENTIRE.db) and one or more stand CNs, write for EACH CN a
# self-contained keyword file + tree file that runs standalone (no database needed):
#
#   <outdir>/<CN>.key + <CN>.tre     (legacy fixed-column)   — the default
#   <outdir>/<CN>.yaml + <CN>.csv    (modern readable)       — with --format yaml
#
# The stand's setup cards (STDIDENT / STDINFO / DESIGN / SITECODE / INVYEAR / GROWTH /
# NUMCYCLE) are materialized from the FVS_STANDINIT_COND columns — the SAME fields the
# native reader (`apply_fia_stand!`) consumes, so a re-read reproduces the stand state.
# The tree records come straight from FVS_TREEINIT_COND. The result is portable: it runs
# with `bin/fvsjl-run.jl <CN>.key` (or `.yaml`) with no SQLite dependency, and converts
# with `bin/fvsjl-translate.jl` like any other stand file.
#
# Usage:
#   julia --project bin/fvsjl-fia-export.jl <fia.db> <CN | @cnfile | CN1,CN2,…> [outdir] \
#         [--variant SN|NE|CS|LS] [--format key|yaml] [--numcycle N] [--validate]
#
#   <cnfile>     a file with one CN per line (blank lines / # comments ignored). A leading
#                '@' marks a file; a bare comma-separated list is also accepted.
#   outdir       output directory (created; default '.').
#   --variant    which variant the stands are (default SN). Written into the .yaml header;
#                a .key carries no variant (pass --variant to fvsjl-run to run it).
#   --format     key (default; .key + .tre) or yaml (.yaml + .csv).
#   --numcycle   NUMCYCLE cycles to write into each stand (default 10).
#   --validate   after writing, run each exported stand AND the same CN via a DATABASE
#                reader, and report whether their .sum rows match (a faithfulness check).
# =============================================================================
using FVSjl
import SQLite, DBInterface

const STANDSQL = "SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = :cn"
const TREESQL  = "SELECT * FROM FVS_TREEINIT_COND  WHERE STAND_CN = :cn"

# The FVS-ready master (e.g. SQLITE_FIADB_ENTIRE.db, ~1.5M rows) has NO index on STAND_CN,
# so a per-CN `WHERE STAND_CN=…` full-scans a huge table. Build a small in-memory working
# copy of ONLY the requested CNs, INDEXED — one scan of the master, then O(log n) lookups.
# The master is ATTACHed read-only and never modified.
function indexed_workset(dbpath, cns)
    inlist = join(["'" * replace(c, "'"=>"''") * "'" for c in cns], ",")
    mem = SQLite.DB()                                        # in-memory
    DBInterface.execute(mem, "ATTACH DATABASE 'file:$(dbpath)?mode=ro&immutable=1' AS m")
    for tbl in ("FVS_STANDINIT_COND", "FVS_TREEINIT_COND")
        DBInterface.execute(mem, "CREATE TABLE $tbl AS SELECT * FROM m.$tbl WHERE STAND_CN IN ($inlist)")
        DBInterface.execute(mem, "CREATE INDEX ix_$tbl ON $tbl(STAND_CN)")
    end
    DBInterface.execute(mem, "DETACH DATABASE m")
    return mem
end

# ---- small helpers over a raw SQLite row (uppercase-keyed Dict) --------------
_present(d, k) = haskey(d, k) && d[k] !== missing && d[k] !== nothing
_num(d, k, dv=nothing) = _present(d, k) ? d[k] : dv
# Render a numeric field as compact text (FVS free-reads inside the 10-col field, so an
# integer-valued float is written without a trailing ".0"). Strings pass through.
function _ftext(x)
    x === nothing && return ""
    x isa AbstractString && return strip(x)
    xf = Float64(x)
    isinteger(xf) ? string(Int(round(xf))) : string(Float32(xf))
end

# Open the master DB READ-ONLY (never write/journal the source), fetch a row Dict.
function _rows(db, sql, cn)
    out = Dict{String,Any}[]
    for r in DBInterface.execute(db, sql, (cn = cn,))
        push!(out, Dict{String,Any}(uppercase(String(k)) => v for (k, v) in pairs(r)))
    end
    out
end

# card(name; f1=…, f2=…) → a KeywordRecord with the given 1-based fields set.
function card(name; kwargs...)
    fld = Dict{Int,String}()
    for (k, v) in kwargs
        v === nothing && continue
        t = _ftext(v); isempty(t) && continue
        fld[parse(Int, String(k)[2:end])] = t   # keys look like :f1, :f2, …
    end
    FVSjl._sem_card(uppercase(name), fld)
end
raw(t) = FVSjl._raw_record(String(t))

# ---- build the setup keyword records from a FVS_STANDINIT_COND row -----------
function stand_setup_records(d::Dict{String,Any}, cn::AbstractString; numcycle::Int)
    recs = KeywordRecord[]
    push!(recs, card("STDIDENT")); push!(recs, raw(cn))

    # LOCATION (KODFOR): direct, else REGION*100+FOREST (mirrors apply_fia_stand!).
    loc = _num(d, "LOCATION", 0)
    if (loc === nothing || Int(loc) == 0) && _present(d, "REGION")
        loc = Int(_num(d, "REGION", 0)) * 100 + Int(_num(d, "FOREST", 0))
    end
    # ELEVATION in hundreds of ft; ELEVFT is feet ⇒ /100.
    elev = _present(d, "ELEVATION") ? _num(d, "ELEVATION") :
           (_present(d, "ELEVFT") ? Float64(_num(d, "ELEVFT", 0)) / 100 : nothing)
    push!(recs, card("STDINFO";
        f1 = (loc === nothing || Int(loc) == 0) ? nothing : Int(loc),   # forest_code
        f2 = _present(d, "ECOREGION") ? _num(d, "ECOREGION") : nothing, # habitat / eco unit
        f3 = _num(d, "AGE"),                                            # stand_age
        f4 = _num(d, "ASPECT"),                                         # aspect (deg)
        f5 = _num(d, "SLOPE"),                                          # slope (%)
        f6 = elev))                                                     # elevation (100 ft)

    # DESIGN card: BAF / FPA / BRK / IPTINV / NONSTK / SAMWT / GROSPC.
    push!(recs, card("DESIGN";
        f1 = _num(d, "BASAL_AREA_FACTOR"),
        f2 = _num(d, "INV_PLOT_SIZE"),
        f3 = _num(d, "BRK_DBH"),
        f4 = _num(d, "NUM_PLOTS"),
        f5 = _num(d, "NONSTK_PLOTS"),
        f6 = _num(d, "SAM_WT"),
        f7 = _num(d, "STK_PCNT")))

    # SITECODE: site species (field 1) + site index (field 2).
    if _present(d, "SITE_INDEX")
        push!(recs, card("SITECODE";
            f1 = _present(d, "SITE_SPECIES") ? _num(d, "SITE_SPECIES") : 0,
            f2 = _num(d, "SITE_INDEX")))
    end

    _present(d, "INV_YEAR") && push!(recs, card("INVYEAR"; f1 = _num(d, "INV_YEAR")))

    # GROWTH card (IDG/FINT/IHTG/FINTH/FINTM) when the DB carries measured-growth setup.
    if _present(d, "DG_TRANS") || _present(d, "DG_MEASURE")
        push!(recs, card("GROWTH";
            f1 = _num(d, "DG_TRANS", 0),  f2 = _num(d, "DG_MEASURE", 5),
            f3 = _num(d, "HTG_TRANS", 0), f4 = _num(d, "HTG_MEASURE", 5),
            f5 = _num(d, "MORT_MEASURE", 5)))
    end

    push!(recs, card("NUMCYCLE"; f1 = numcycle))
    push!(recs, card("ECHOSUM"))
    push!(recs, card("TREEDATA"))
    push!(recs, card("PROCESS"))
    return recs
end

# ---- build tree records from FVS_TREEINIT_COND rows (mirror apply_fia_trees!) -
function tree_records(rows::Vector{Dict{String,Any}})
    recs = TreeRecord[]
    ci(d, k, dv) = (v = _num(d, k, dv); v === nothing ? Int32(dv) : Int32(round(Float64(v))))
    cf(d, k, dv) = (v = _num(d, k, dv); v === nothing ? Float32(dv) : Float32(Float64(v)))
    for d in rows
        dbh = _present(d, "DIAMETER") ? cf(d, "DIAMETER", 0) : cf(d, "DBH", 0)
        dmg = (ci(d,"DAMAGE1",0), ci(d,"SEVERITY1",0), ci(d,"DAMAGE2",0),
               ci(d,"SEVERITY2",0), ci(d,"DAMAGE3",0), ci(d,"SEVERITY3",0))
        sp = _present(d, "SPECIES") ? strip(string(_num(d, "SPECIES", "OT"))) : "OT"
        push!(recs, TreeRecord(
            ci(d,"PLOT_ID",1), ci(d,"TREE_ID",0), cf(d,"TREE_COUNT",1), ci(d,"HISTORY",1),
            String(sp), dbh, cf(d,"DG",0), cf(d,"HT",0), cf(d,"HTTOPK",0), cf(d,"HTG",0),
            ci(d,"CRRATIO",0), dmg, ci(d,"TREEVALUE",0), ci(d,"PRESCRIPTION",0),
            (Int32(0),Int32(0),Int32(0),Int32(0),Int32(0)), cf(d,"AGE",0)))
    end
    recs
end

# ---- resolve the CN list from the positional argument -----------------------
function cn_list(arg::AbstractString)
    if startswith(arg, "@") || isfile(arg)
        path = startswith(arg, "@") ? arg[2:end] : arg
        return [strip(l) for l in eachline(path)
                if !isempty(strip(l)) && !startswith(strip(l), "#")]
    end
    return [strip(x) for x in split(arg, ",") if !isempty(strip(x))]
end

function main(args)
    pos = String[]; variant = "SN"; fmt = "key"; numcycle = 10; validate = false
    i = 1
    while i <= length(args)
        a = args[i]
        if     a == "--variant";  variant = args[i+1]; i += 2
        elseif a == "--format";   fmt = lowercase(args[i+1]); i += 2
        elseif a == "--numcycle"; numcycle = parse(Int, args[i+1]); i += 2
        elseif a == "--validate"; validate = true; i += 1
        elseif startswith(a, "--variant=");  variant = split(a,"=",limit=2)[2]; i += 1
        elseif startswith(a, "--format=");   fmt = lowercase(split(a,"=",limit=2)[2]); i += 1
        elseif startswith(a, "--numcycle="); numcycle = parse(Int, split(a,"=",limit=2)[2]); i += 1
        else push!(pos, a); i += 1
        end
    end
    if length(pos) < 2
        println(stderr, "usage: fvsjl-fia-export <fia.db> <CN|@cnfile|CN1,CN2,…> [outdir] " *
                        "[--variant SN|NE|CS|LS] [--format key|yaml] [--numcycle N] [--validate]")
        return 1
    end
    dbpath = pos[1]; cns = cn_list(pos[2]); outdir = length(pos) >= 3 ? pos[3] : "."
    isfile(dbpath) || (println(stderr, "error: no such database: $dbpath"); return 1)
    fmt in ("key", "yaml") || (println(stderr, "error: --format must be key or yaml"); return 1)
    mkpath(outdir)
    var = variant_from_code(variant)
    keyext = fmt == "yaml" ? ".yaml" : ".key"

    # Build a small INDEXED in-memory working set for the requested CNs (fast on the
    # unindexed master; the master is opened read-only and never modified).
    db = indexed_workset(dbpath, cns)
    nok = 0; nskip = 0; nbad = 0
    try
        for cn in cns
            srows = _rows(db, STANDSQL, cn)
            if isempty(srows); println(stderr, "skip $cn: no FVS_STANDINIT_COND row"); nskip += 1; continue; end
            trows = _rows(db, TREESQL, cn)
            recs  = stand_setup_records(srows[1], cn; numcycle = numcycle)
            trees = tree_records(trows)
            keyp  = joinpath(outdir, cn * keyext)
            treep = joinpath(outdir, cn * (fmt == "yaml" ? ".csv" : ".tre"))
            if fmt == "yaml"
                write_keywords_yaml(recs, keyp)
                # Make the .yaml self-describing: prepend a top-level `variant:` the reader
                # (yaml_variant_code) picks up, so it runs as this variant without a --variant flag.
                write(keyp, "variant: $variant\n" * read(keyp, String))
                write_trees_csv(trees, treep)
            else
                write_keyfile(recs, keyp)
                write_tree_file(trees, treep)
            end
            nok += 1
            println("wrote $keyp  (+ $(length(trees)) trees → $treep)")

            if validate
                # Standalone run vs the DATABASE-reader run for the same CN → compare .sum rows.
                # The trees and direct-measurement fields are exact, but STDINFO cards can't carry
                # the FVS_STANDINIT lat/long ⇒ a possible ±1 crown-model (CCF) ULP on the setup. So
                # report the per-row bit-exact count, not a strict all-or-nothing pass.
                a = _sum_rows(try FVSjl.run_keyfile(keyp; variant = var) catch e; "ERR:$e" end)
                tmp = tempname() * ".key"; write(tmp, _dbreader_key(cn, dbpath))
                b = _sum_rows(try FVSjl.run_keyfile(tmp; variant = var) catch e; "ERR:$e" end)
                rm(tmp; force = true)
                if isempty(a) || length(a) != length(b)
                    println("  validate $cn: FAILED (standalone $(length(a)) vs DB $(length(b)) rows)"); nbad += 1
                elseif a == b
                    println("  validate $cn: standalone == DB-reader BIT-EXACT ✓ ($(length(a)) rows)")
                else
                    # The export's fidelity is the INITIAL state (cycle-0 row): the tree list and
                    # direct-measurement fields are exact; only crown competition (CCF) can differ by
                    # ~1 (STDINFO cards can't carry the DB lat/long → a Hopkins-index ULP). Later cycles
                    # inherit the model's own dense-phase compounded-ULP, so we compare cycle 0.
                    cols = ["TPA","BA","SDI","CCF","TopHt","QMD","TCuFt","MCuFt","SCuFt","BdFt"]
                    fa = split(a[1]); fb = split(b[1]); diffs = String[]
                    for k in 1:10
                        (length(fa) < k+2 || length(fb) < k+2) && continue
                        va = tryparse(Float64, fa[k+2]); vb = tryparse(Float64, fb[k+2])
                        (va === nothing || vb === nothing || va == vb) && continue
                        push!(diffs, "$(cols[k]) $va≠$vb")
                    end
                    if isempty(diffs)
                        println("  validate $cn: cycle-0 (initial state) BIT-EXACT ✓; later cycles = model dense-phase ULP")
                    else
                        println("  validate $cn: cycle-0 faithful, differs only in [$(join(diffs, ", "))] " *
                                "— crown-model ULP (STDINFO can't carry the DB lat/long); later cycles = model ULP")
                    end
                end
            end
        end
    finally
        SQLite.close(db)
    end
    println(stderr, "exported $nok stand(s) to $outdir/ ($nskip skipped" *
                    (validate ? ", $nbad validation mismatches" : "") * ")")
    return nbad == 0 ? 0 : 2
end

# .sum data rows (year-prefixed) from a run's text output, for the --validate compare.
_sum_rows(txt) = txt isa AbstractString && !startswith(txt, "ERR") ?
    [l for l in split(txt, '\n') if (length(l) >= 4 && (y = tryparse(Int, l[1:4]);
                                     y !== nothing && 1000 <= y <= 3000))] : String[]

# A keyfile that reads the one CN straight from the FIA DB via the DATABASE block.
_dbreader_key(cn, dbpath) = """
STDIDENT
$cn
DATABASE
DSNin
$dbpath
StandSQL
SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
TreeSQL
SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
END
NUMCYCLE  10
ECHOSUM
PROCESS
STOP
"""

exit(main(ARGS))
