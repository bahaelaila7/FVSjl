# =============================================================================
# live_fvs_db.jl — run a .key through the LIVE Fortran FVSsn and read its DBS output.
#
# DISCOVERY (2026): the rebuilt FVSsn (bin/FVSsn_buildDir, assembled via
#   `gfortran -o FVSsn_full $(ls *.o) glibc_shim.o -lpthread -ldl`) is NOT a stripped
# DBS build — it emits the FULL table set (FVS_EconSummary, FVS_BurnReport,
# FVS_PotFire_East, FVS_Carbon, FVS_Fuels, FVS_Mortality, FVS_SnagSum, FVS_Down_Wood_*,
# FVS_Hrv_Carbon, FVS_TreeList, FVS_CutList, FVS_Summary…) when the keyfile's `DataBase`
# block requests them (Summary / TreeliDB / ECONRPTS / CARBRPTS / FuelRpts / BurnRpt …).
# The earlier "stripped binary" belief was a tooling artifact: the `sqlite3` CLI is absent,
# so DBS tables must be read via SQLite.jl (this helper), and the committed golden
# FVSOut.db was generated with a minimal DataBase block (3 tables only).
#
# This makes the gold-standard live-FVS differential feasible for the value-path flags
# (B6 econ, B1 fire flame/scorch, carbon, mortality, …) — no fuller rebuild needed.
# =============================================================================
module LiveFVS
using SQLite

const BUILDDIR = "/workspace/ForestVegetationSimulator/bin/FVSsn_buildDir"
const SHIM     = "/workspace/FVSjl/tmp/glibc_shim.o"
# Persist the live binary under FVSjl/tmp (survives container restarts; /tmp is wiped). ensure_binary()
# rebuilds it from the buildDir object set if it is missing, so a fresh container self-heals on first use.
const BINARY   = "/workspace/FVSjl/tmp/FVSsn_full"

"Build the live FVSsn binary from the buildDir object set (once). Returns the path."
function ensure_binary()
    isfile(BINARY) && return BINARY
    objs = filter(f -> endswith(f, ".o"), readdir(BUILDDIR; join = true))
    run(pipeline(`gfortran -o $BINARY $objs $SHIM -lpthread -ldl`; stdout = devnull, stderr = devnull))
    return BINARY
end

"""
    run_key(keypath) -> dbpath

Run `keypath` through the live FVSsn in an isolated temp dir (FVS writes outputs next
to the key). Returns the path to the SQLite DB it produced (`*.db`), or errors if none.
The key's `DataBase`/`DSNOut` block decides which tables appear.
"""
function run_key(keypath::AbstractString; timeout = 120)
    bin = ensure_binary()
    tmp = mktempdir()
    key = joinpath(tmp, basename(keypath)); cp(keypath, key; force = true)
    tre = replace(keypath, r"\.key$" => ".tre")
    isfile(tre) && cp(tre, joinpath(tmp, basename(tre)); force = true)
    cd(tmp) do
        # FVS exits non-zero (interactive-loop EOF) even on a successful run — ignore the status
        # and judge success by whether a .db was written.
        run(pipeline(ignorestatus(`timeout $timeout $bin`); stdin = IOBuffer(basename(key) * "\n"),
                     stdout = devnull, stderr = devnull); wait = true)
    end
    dbs = filter(f -> endswith(f, ".db"), readdir(tmp; join = true))
    isempty(dbs) && error("live FVSsn produced no .db for $(basename(keypath)) (check the DataBase block)")
    return first(dbs)
end

"List the tables (and row counts) in a DBS file."
function tables(dbpath::AbstractString)
    db = SQLite.DB(dbpath)
    [(name = r.name,
      rows = first(SQLite.DBInterface.execute(db, "SELECT count(*) c FROM \"$(r.name)\"")).c)
     for r in SQLite.DBInterface.execute(db,
        "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")]
end

"Return all rows of `table` from `dbpath` as a Vector of NamedTuples."
function rows(dbpath::AbstractString, table::AbstractString)
    db = SQLite.DB(dbpath)
    [NamedTuple(r) for r in SQLite.DBInterface.execute(db, "SELECT * FROM \"$table\"")]
end

end # module
