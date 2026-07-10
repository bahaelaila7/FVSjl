# sweep_db.jl — DURABLE, cross-session record of the FIA/FVS coverage sweep.
#
# The full-population sweep (run_expand_cycle.sh) compares every real FIA stand's multi-cycle projection
# (FVSjl vs freshly-relinked live FVS) and writes a per-stand differential row to an EPHEMERAL scratchpad ledger
# that vanishes at session end. This module persists that per-stand outcome to a local SQLite database on the
# durable repo volume (NOT the session scratchpad, NOT git — see .gitignore) so the coverage record survives
# sessions AND container restarts. Purpose: a durable WORKLIST of which stands still need a dig / reassessment.
#
# Every stand is bucketed into exactly one `dig_class`:
#   • bit_exact  — FVSjl == live FVS on all 10 .sum cols, every cycle.
#   • ulp_class  — diverges, but the divergence is an ACCEPTED cornered primitive (print/ULP boundary, self-thin
#                  count-straddle, merch/board threshold-crossing, compounded-ULP dense-phase). No action needed.
#   • needs_dig  — diverges in a way the escalation guard could NOT corner (UNCLASSIFIED signature, a MATERIAL
#                  structure move ≥10 abs units & ≥15%, or a threshold-free total-cubic ≥15%). REASSESS / TRACE.
#
# `dig_class` is the SINGLE SOURCE of the "what to look at" question — it mirrors filter_digworthy.jl's
# escalation guard exactly (kept in sync; both key on struct_max_abs so a young/small-base ±1-unit straddle that
# inflates to a big RELATIVE % is ulp_class, not needs_dig — see audit slice 43g / CN 202567027010854).
#
# Usage:
#   julia --project=. test/harness/fia/sweep_db.jl ingest <db> <ledger-or-cycle.csv>   # upsert rows (idempotent)
#   julia --project=. test/harness/fia/sweep_db.jl stats  <db> [variant]               # dig_class breakdown
#   julia --project=. test/harness/fia/sweep_db.jl digs   <db> [variant]               # list needs_dig CNs
# Programmatic (from ledger_fia.jl): open_sweepdb(path) once, then upsert!(db, row) per stand.

import SQLite, DBInterface, Dates

const DEFAULT_SWEEP_DB = "/workspace/FVSjl/data/fia_sweep.db"

# Escalation-guard mirror (filter_digworthy.jl): structure blow-up must be material in ABSOLUTE terms too.
const _STRUCT_ESCALATE_COLS = Set(["TPA","BA","SDI","CCF","QMD"])
const _ESCALATE_REL = 15.0
const _STRUCT_ABS_FLOOR = 10.0
# TCuFt volume net also needs an ABSOLUTE floor: a real volume-equation bug (FORKOD zero-vol) moves 1000s of
# cuft; a 15% on a tiny/degenerate stand (e.g. a 2-tree stand, 62 cuft on 412) is not a bug. 300 cuft clears the
# degenerate case (both-sides-traced CN 218434248010854) while surfacing any real volume-equation error.
const _VOL_ABS_FLOOR = 300.0

# dig_class over the MEASURED facts. Missing *_abs (legacy rows) ⇒ +Inf ⇒ conservative (escalates).
# NOTE (slice 43p): a struct escalation gated on struct_max_rel_pct instead of worst_col OVER-flags — struct_abs
# is TPA-dominated (always huge in dense stands) and struct_max_rel_pct is inflated by TopHt (AVHT40 ULP) and
# small-base columns. A faithful "structure-divergent-but-worst_col=volume" net needs a density-specific metric
# (BA/SDI relative, excluding TPA+TopHt); deferred. The one real case the dig found (CN 209314057, systematic
# growing BA/SDI divergence) is flagged manually. Keeping the conservative worst_col gate here.
function dig_class(bit_exact::Bool, sig::AbstractString, worst_col::AbstractString,
                   max_rel_pct::Real, struct_max_abs::Union{Real,Nothing}, vol_max_abs::Union{Real,Nothing}=nothing,
                   struct_max_rel_pct::Union{Real,Nothing}=nothing)
    bit_exact && return "bit_exact"
    sa = struct_max_abs === nothing ? Inf : float(struct_max_abs)
    va = vol_max_abs === nothing ? Inf : float(vol_max_abs)
    esc = sig == "UNCLASSIFIED" ||
          (sig == "structure_densephase" && worst_col in _STRUCT_ESCALATE_COLS &&
                 max_rel_pct >= _ESCALATE_REL && sa >= _STRUCT_ABS_FLOOR) ||
          (worst_col == "TCuFt" && max_rel_pct >= _ESCALATE_REL && va >= _VOL_ABS_FLOOR)
    return esc ? "needs_dig" : "ulp_class"
end

# Columns: variant,cn,regime keyed; the rest are the measured ledger facts + the derived dig_class worklist key
# (bit_exact | ulp_class | needs_dig) and the ingest timestamp. struct_max_abs is NULL for legacy 15-col rows.
const _SCHEMA = [
    """CREATE TABLE IF NOT EXISTS sweep (
        variant TEXT NOT NULL, cn TEXT NOT NULL, regime TEXT NOT NULL DEFAULT 'none',
        n_cycles INTEGER, bit_exact INTEGER, div_cols TEXT, worst_col TEXT, worst_cycle INTEGER,
        max_rel_pct REAL, max_abs_diff REAL, struct_max_rel_pct REAL, vol_max_rel_pct REAL,
        struct_max_abs REAL, density_bitexact INTEGER, converges INTEGER, signature TEXT,
        dig_class TEXT NOT NULL, swept_at TEXT, vol_max_abs REAL,
        PRIMARY KEY (variant, cn, regime))""",
    # migrate DBs created before vol_max_abs existed (idempotent — the wrapper below swallows "duplicate column")
    "ALTER TABLE sweep ADD COLUMN vol_max_abs REAL",
    "CREATE INDEX IF NOT EXISTS sweep_digclass ON sweep(dig_class)",
    "CREATE INDEX IF NOT EXISTS sweep_variant ON sweep(variant)",
    # per-variant sweep CURSOR (deterministic-order offset into the population) — a self-contained durable
    # snapshot of PROGRESS so a resume after container restart doesn't re-sweep from scratch even if the
    # working cursor file is gone. Mirrors test/harness/fia/expand/<v>.cursor.
    """CREATE TABLE IF NOT EXISTS progress (
        variant TEXT PRIMARY KEY, cursor INTEGER NOT NULL, population INTEGER, updated_at TEXT)""",
]

# Run a statement and DRAIN its cursor — an unconsumed result (e.g. PRAGMA journal_mode returns a row) leaves
# the statement "in progress" and blocks the next transaction/savepoint.
_run(db, sql) = (for _ in DBInterface.execute(db, sql); end)

function open_sweepdb(path::AbstractString=DEFAULT_SWEEP_DB)
    mkpath(dirname(path))
    db = SQLite.DB(path)
    for stmt in _SCHEMA
        try
            _run(db, stmt)
        catch e
            # tolerate the idempotent ALTER (column already present on a fresh or migrated DB)
            startswith(stmt, "ALTER TABLE") || rethrow(e)
        end
    end
    # pragmatic durability/perf: WAL survives an abrupt container stop mid-write better than the default journal.
    _run(db, "PRAGMA journal_mode=WAL")
    db
end

_i(x) = x === nothing ? nothing : Int(x)
_f(x) = x === nothing ? nothing : float(x)

# Upsert one stand. `row` is a NamedTuple with the ledger fields (see COLS order in ledger_fia.jl).
function upsert!(db, row)
    vma = hasproperty(row, :vol_max_abs) ? row.vol_max_abs : nothing
    dc = dig_class(row.bit_exact, row.signature, something(row.worst_col, ""),
                   something(row.max_rel_pct, 0.0), row.struct_max_abs, vma, row.struct_max_rel_pct)
    DBInterface.execute(db, """
        INSERT INTO sweep (variant,cn,regime,n_cycles,bit_exact,div_cols,worst_col,worst_cycle,max_rel_pct,
                           max_abs_diff,struct_max_rel_pct,vol_max_rel_pct,struct_max_abs,density_bitexact,
                           converges,signature,dig_class,swept_at,vol_max_abs)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        ON CONFLICT(variant,cn,regime) DO UPDATE SET
          n_cycles=excluded.n_cycles, bit_exact=excluded.bit_exact, div_cols=excluded.div_cols,
          worst_col=excluded.worst_col, worst_cycle=excluded.worst_cycle, max_rel_pct=excluded.max_rel_pct,
          max_abs_diff=excluded.max_abs_diff, struct_max_rel_pct=excluded.struct_max_rel_pct,
          vol_max_rel_pct=excluded.vol_max_rel_pct, struct_max_abs=excluded.struct_max_abs,
          density_bitexact=excluded.density_bitexact, converges=excluded.converges,
          signature=excluded.signature, dig_class=excluded.dig_class, swept_at=excluded.swept_at,
          vol_max_abs=excluded.vol_max_abs
        """,
        (row.variant, row.cn, row.regime, _i(row.n_cycles), row.bit_exact ? 1 : 0, row.div_cols,
         row.worst_col, _i(row.worst_cycle), _f(row.max_rel_pct), _f(row.max_abs_diff),
         _f(row.struct_max_rel_pct), _f(row.vol_max_rel_pct), _f(row.struct_max_abs),
         row.density_bitexact ? 1 : 0, row.converges ? 1 : 0, row.signature, dc,
         Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"), _f(vma)))
end

_pbool(s) = lowercase(strip(s)) == "true"
_pf(s) = (v = tryparse(Float64, strip(s)); v)
_pi(s) = (v = tryparse(Int, strip(s)); v)

# A ledger row is well-formed only if these hold. Guards against CSV line-concatenation artifacts (a
# timeout-killed / concurrently-appended ledger file can splice two lines ⇒ shifted columns: the variant ends up
# glued to the CN, the signature slot holds a bool/number). Such rows are DROPPED on ingest, never classified.
const _VALID_VARIANTS = Set(["SN","NE","CS","LS"])
const _VALID_SIGS = Set(["bit_exact","print_boundary","volume_persistent","structure_densephase",
                         "threshold_crossing","count_straddle","UNCLASSIFIED"])
_valid_row(variant, cn, sig) =
    variant in _VALID_VARIANTS && sig in _VALID_SIGS && !isempty(cn) && all(isdigit, cn)

# Ingest a ledger/cycle CSV (positional; header may be 15- or 16-col — struct_max_abs optional). Idempotent.
function ingest_csv(dbpath::AbstractString, csvpath::AbstractString)
    db = open_sweepdb(dbpath)
    n = 0
    SQLite.transaction(db) do                                  # one txn for speed; commits/finalizes cleanly
    for (i, l) in enumerate(eachline(csvpath))
        i == 1 && startswith(l, "variant,") && continue        # skip header if present
        f = split(l, ',')
        length(f) < 15 && continue
        # columns: variant,regime,cn,n_cycles,bit_exact,div_cols,worst_col,worst_cycle,max_rel_pct,max_abs_diff,
        #          struct_max_rel_pct,vol_max_rel_pct,density_bitexact,converges,signature[,struct_max_abs]
        variant=String(strip(f[1])); cn=String(strip(f[3])); sig=String(strip(f[15]))
        _valid_row(variant, cn, sig) || continue      # drop mangled/concatenated lines
        row = (variant=variant, regime=String(strip(f[2])), cn=cn,
               n_cycles=_pi(f[4]), bit_exact=_pbool(f[5]), div_cols=String(strip(f[6])),
               worst_col=String(strip(f[7])), worst_cycle=_pi(f[8]), max_rel_pct=_pf(f[9]),
               max_abs_diff=_pf(f[10]), struct_max_rel_pct=_pf(f[11]), vol_max_rel_pct=_pf(f[12]),
               density_bitexact=_pbool(f[13]), converges=_pbool(f[14]), signature=sig,
               struct_max_abs=(length(f) >= 16 ? _pf(f[16]) : nothing),
               vol_max_abs=(length(f) >= 17 ? _pf(f[17]) : nothing))
        upsert!(db, row); n += 1
    end
    end                                                        # SQLite.transaction
    SQLite.close(db)
    n
end

function set_cursor!(dbpath::AbstractString, variant::AbstractString, cursor::Integer, population=nothing)
    db = open_sweepdb(dbpath)
    DBInterface.execute(db, """
        INSERT INTO progress (variant,cursor,population,updated_at) VALUES (?,?,?,?)
        ON CONFLICT(variant) DO UPDATE SET cursor=excluded.cursor, population=excluded.population,
                                           updated_at=excluded.updated_at""",
        (variant, Int(cursor), population === nothing ? nothing : Int(population),
         Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")))
    SQLite.close(db)
end

function get_cursor(dbpath::AbstractString, variant::AbstractString)
    isfile(dbpath) || return nothing
    db = open_sweepdb(dbpath); c = nothing
    for r in DBInterface.execute(db, "SELECT cursor FROM progress WHERE variant=?", (variant,)); c = r.cursor; end
    SQLite.close(db); c
end

# Delete rows that fail _valid_row (CSV-concatenation artifacts from the one-time master-ledger backfill). Prints
# the digit-only CNs found glued inside those rows so they can be re-swept cleanly. Idempotent.
function scrub!(dbpath::AbstractString)
    db = open_sweepdb(dbpath)
    bad = Tuple{String,String,String}[]; recover = Set{String}()
    for r in DBInterface.execute(db, "SELECT variant,cn,signature FROM sweep")
        v=String(r.variant); c=String(r.cn); s=String(r.signature)
        if !_valid_row(v,c,s)
            push!(bad,(v,c,s))
            # salvage any 12+ digit FIA CN embedded in the mangled variant/cn/signature fields for a re-sweep
            for fld in (c,s), m in eachmatch(r"\d{12,}", fld); push!(recover, m.match); end
        end
    end
    SQLite.transaction(db) do
        for (v,c,s) in bad
            DBInterface.execute(db, "DELETE FROM sweep WHERE variant=? AND cn=? AND signature=?", (v,c,s))
        end
    end
    SQLite.close(db)
    (deleted=length(bad), recover=sort(collect(recover)))
end

# CNs manually confirmed (by a both-sides dig) as genuine needs_dig even though the auto-guard scores them
# ulp_class — a structure-divergent stand whose highest-RELATIVE column is a volume-threshold col (see slice 43p).
# One CN per line; committed so the flag survives reclassify. Format: "<cn> # note".
const MANUAL_NEEDSDIG_FILE = "/workspace/FVSjl/docs/fia_manual_needsdig.txt"
# Symmetric to MANUAL_NEEDSDIG: CNs a both-sides dig CORNERED to a named primitive (doctrine bar met) even though
# the conservative auto-guard still trips needs_dig on them. reclassify forces these to ulp_class. One CN/line.
const CORNERED_STANDS_FILE = "/workspace/FVSjl/docs/fia_cornered_stands.txt"
function _read_cn_set(path)
    s = Set{String}()
    isfile(path) || return s
    for l in eachline(path)
        t = strip(split(l, '#')[1]); isempty(t) || push!(s, String(t))
    end
    s
end
_manual_needsdig() = _read_cn_set(MANUAL_NEEDSDIG_FILE)
_cornered_stands() = _read_cn_set(CORNERED_STANDS_FILE)

# Recompute dig_class for every row from the STORED facts (no FVS re-run) — apply a refined guard to the whole DB.
function reclassify!(dbpath::AbstractString)
    db = open_sweepdb(dbpath)
    manual = _manual_needsdig(); cornered = _cornered_stands()
    rows = Tuple[]
    for r in DBInterface.execute(db, """SELECT variant,cn,regime,bit_exact,signature,worst_col,max_rel_pct,
                                        struct_max_abs,vol_max_abs,struct_max_rel_pct,dig_class FROM sweep""")
        mn(x) = (x === missing || x === nothing) ? nothing : x
        dc = dig_class(r.bit_exact == 1, String(something(r.signature,"")), String(something(r.worst_col,"")),
                       something(mn(r.max_rel_pct), 0.0), mn(r.struct_max_abs), mn(r.vol_max_abs), mn(r.struct_max_rel_pct))
        String(r.cn) in manual && (dc = "needs_dig")   # force manually-confirmed genuine finds
        # cornered-to-named-primitive (dig complete) overrides needs_dig → ulp_class; keep bit_exact as-is
        String(r.cn) in cornered && dc == "needs_dig" && (dc = "ulp_class")
        dc == r.dig_class || push!(rows, (dc, r.variant, r.cn, r.regime))
    end
    SQLite.transaction(db) do
        for (dc,v,cn,rg) in rows
            DBInterface.execute(db, "UPDATE sweep SET dig_class=? WHERE variant=? AND cn=? AND regime=?", (dc,v,cn,rg))
        end
    end
    SQLite.close(db)
    length(rows)
end

function _stats(dbpath, variant)
    db = open_sweepdb(dbpath)
    where = variant === nothing ? "" : "WHERE variant='$(variant)'"
    println("dig_class breakdown ", variant === nothing ? "(all variants)" : "($variant)", ":")
    for r in DBInterface.execute(db, "SELECT dig_class, COUNT(*) c FROM sweep $where GROUP BY dig_class ORDER BY c DESC")
        println("  ", rpad(r.dig_class, 12), r.c)
    end
    for r in DBInterface.execute(db, "SELECT COUNT(*) c, COUNT(DISTINCT cn) d FROM sweep $where")
        println("  total rows=", r.c, "  distinct CNs=", r.d)
    end
    SQLite.close(db)
end

function _digs(dbpath, variant)
    db = open_sweepdb(dbpath)
    where = variant === nothing ? "WHERE dig_class='needs_dig'" : "WHERE dig_class='needs_dig' AND variant='$(variant)'"
    println("variant,cn,regime,signature,worst_col,worst_cycle,max_rel_pct,struct_max_abs,swept_at")
    for r in DBInterface.execute(db, "SELECT variant,cn,regime,signature,worst_col,worst_cycle,max_rel_pct,struct_max_abs,swept_at FROM sweep $where ORDER BY variant,cn")
        println(join([r.variant,r.cn,r.regime,r.signature,r.worst_col,r.worst_cycle,r.max_rel_pct,r.struct_max_abs,r.swept_at], ","))
    end
    SQLite.close(db)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 1 || error("usage: sweep_db.jl {ingest <db> <csv> | stats <db> [variant] | digs <db> [variant]}")
    cmd = ARGS[1]
    if cmd == "ingest"
        length(ARGS) >= 3 || error("ingest needs <db> <csv>")
        n = ingest_csv(ARGS[2], ARGS[3]); println("ingested/updated $n stands into $(ARGS[2])")
    elseif cmd == "stats"
        _stats(ARGS[2], length(ARGS) >= 3 ? ARGS[3] : nothing)
    elseif cmd == "digs"
        _digs(ARGS[2], length(ARGS) >= 3 ? ARGS[3] : nothing)
    elseif cmd == "setcursor"      # setcursor <db> <variant> <offset> [population]
        set_cursor!(ARGS[2], ARGS[3], parse(Int, ARGS[4]), length(ARGS) >= 5 ? parse(Int, ARGS[5]) : nothing)
    elseif cmd == "getcursor"      # getcursor <db> <variant>  → prints the offset (empty if none)
        c = get_cursor(ARGS[2], ARGS[3]); c === nothing || println(c)
    elseif cmd == "cns"            # cns <db> <variant>  → print every already-recorded CN (one per line)
        db = open_sweepdb(ARGS[2])
        for r in DBInterface.execute(db, "SELECT cn FROM sweep WHERE variant=?", (ARGS[3],)); println(r.cn); end
        SQLite.close(db)
    elseif cmd == "reclassify"     # reclassify <db>  → recompute dig_class for all rows from stored facts
        n = reclassify!(ARGS[2]); println("reclassified $n rows")
    elseif cmd == "scrub"          # scrub <db>  → delete malformed rows; print recoverable CNs (one per line)
        res = scrub!(ARGS[2])
        print(stderr, "scrubbed $(res.deleted) malformed rows; $(length(res.recover)) recoverable CNs\n")
        for cn in res.recover; println(cn); end
    else
        error("unknown command $cmd")
    end
end
