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

# dig_class over the MEASURED facts. struct_abs missing (legacy 15-col rows) ⇒ +Inf ⇒ conservative (escalates).
function dig_class(bit_exact::Bool, sig::AbstractString, worst_col::AbstractString,
                   max_rel_pct::Real, struct_max_abs::Union{Real,Nothing})
    bit_exact && return "bit_exact"
    sa = struct_max_abs === nothing ? Inf : float(struct_max_abs)
    esc = sig == "UNCLASSIFIED" ||
          (sig == "structure_densephase" && worst_col in _STRUCT_ESCALATE_COLS &&
                 max_rel_pct >= _ESCALATE_REL && sa >= _STRUCT_ABS_FLOOR) ||
          (worst_col == "TCuFt" && max_rel_pct >= _ESCALATE_REL)
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
        dig_class TEXT NOT NULL, swept_at TEXT,
        PRIMARY KEY (variant, cn, regime))""",
    "CREATE INDEX IF NOT EXISTS sweep_digclass ON sweep(dig_class)",
    "CREATE INDEX IF NOT EXISTS sweep_variant ON sweep(variant)",
]

# Run a statement and DRAIN its cursor — an unconsumed result (e.g. PRAGMA journal_mode returns a row) leaves
# the statement "in progress" and blocks the next transaction/savepoint.
_run(db, sql) = (for _ in DBInterface.execute(db, sql); end)

function open_sweepdb(path::AbstractString=DEFAULT_SWEEP_DB)
    mkpath(dirname(path))
    db = SQLite.DB(path)
    for stmt in _SCHEMA
        _run(db, stmt)
    end
    # pragmatic durability/perf: WAL survives an abrupt container stop mid-write better than the default journal.
    _run(db, "PRAGMA journal_mode=WAL")
    db
end

_i(x) = x === nothing ? nothing : Int(x)
_f(x) = x === nothing ? nothing : float(x)

# Upsert one stand. `row` is a NamedTuple with the ledger fields (see COLS order in ledger_fia.jl).
function upsert!(db, row)
    dc = dig_class(row.bit_exact, row.signature, something(row.worst_col, ""),
                   something(row.max_rel_pct, 0.0), row.struct_max_abs)
    DBInterface.execute(db, """
        INSERT INTO sweep (variant,cn,regime,n_cycles,bit_exact,div_cols,worst_col,worst_cycle,max_rel_pct,
                           max_abs_diff,struct_max_rel_pct,vol_max_rel_pct,struct_max_abs,density_bitexact,
                           converges,signature,dig_class,swept_at)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        ON CONFLICT(variant,cn,regime) DO UPDATE SET
          n_cycles=excluded.n_cycles, bit_exact=excluded.bit_exact, div_cols=excluded.div_cols,
          worst_col=excluded.worst_col, worst_cycle=excluded.worst_cycle, max_rel_pct=excluded.max_rel_pct,
          max_abs_diff=excluded.max_abs_diff, struct_max_rel_pct=excluded.struct_max_rel_pct,
          vol_max_rel_pct=excluded.vol_max_rel_pct, struct_max_abs=excluded.struct_max_abs,
          density_bitexact=excluded.density_bitexact, converges=excluded.converges,
          signature=excluded.signature, dig_class=excluded.dig_class, swept_at=excluded.swept_at
        """,
        (row.variant, row.cn, row.regime, _i(row.n_cycles), row.bit_exact ? 1 : 0, row.div_cols,
         row.worst_col, _i(row.worst_cycle), _f(row.max_rel_pct), _f(row.max_abs_diff),
         _f(row.struct_max_rel_pct), _f(row.vol_max_rel_pct), _f(row.struct_max_abs),
         row.density_bitexact ? 1 : 0, row.converges ? 1 : 0, row.signature, dc,
         Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")))
end

_pbool(s) = lowercase(strip(s)) == "true"
_pf(s) = (v = tryparse(Float64, strip(s)); v)
_pi(s) = (v = tryparse(Int, strip(s)); v)

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
        row = (variant=String(strip(f[1])), regime=String(strip(f[2])), cn=String(strip(f[3])),
               n_cycles=_pi(f[4]), bit_exact=_pbool(f[5]), div_cols=String(strip(f[6])),
               worst_col=String(strip(f[7])), worst_cycle=_pi(f[8]), max_rel_pct=_pf(f[9]),
               max_abs_diff=_pf(f[10]), struct_max_rel_pct=_pf(f[11]), vol_max_rel_pct=_pf(f[12]),
               density_bitexact=_pbool(f[13]), converges=_pbool(f[14]), signature=String(strip(f[15])),
               struct_max_abs=(length(f) >= 16 ? _pf(f[16]) : nothing))
        upsert!(db, row); n += 1
    end
    end                                                        # SQLite.transaction
    SQLite.close(db)
    n
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
    else
        error("unknown command $cmd")
    end
end
