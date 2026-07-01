# =============================================================================
# event_monitor.jl — FVS event monitor (IF/THEN/ENDIF conditional scheduling)
#
# Ported (semantics) from: base/evmon.f + base/algmon.f (ALGCMP/ALGEVL).
#
# An `IF <expr> THEN <activities> ENDIF` block schedules its activities only in the
# cycles where the algebraic condition <expr> evaluates true. The Fortran compiles
# <expr> to a postfix opcode stream (ALGCMP) and evaluates it each cycle (ALGEVL);
# we parse <expr> once into an expression AST (`EvNode`) and evaluate it per cycle
# (`eval_event`) — same results, idiomatic, and trim-friendly (no closures).
#
# Grammar (ALGKEY token set): arithmetic + - * /, comparison GT GE LT LE EQ NE,
# logical AND OR NOT (truth = nonzero), functions FRAC INT MOD EXP SQRT ALOG ABS,
# variables resolved by `_event_var` (CYCLE/YEAR + stand vars, extensible).
# =============================================================================

# --- AST node types (EvNode declared in core/state.jl) ----------------------
struct EvNum <: EvNode; v::Float32; end
struct EvVar <: EvNode; name::String; end
struct EvUn  <: EvNode; op::Symbol; a::EvNode; end
struct EvBin <: EvNode; op::Symbol; a::EvNode; b::EvNode; end
struct EvFun <: EvNode; name::Symbol; a::EvNode; b::EvNode; end   # b unused for 1-arg
struct EvTime <: EvNode; args::Vector{EvNode}; end                # TIME(v0,y1,v1,…) variadic year-step fn

"Per-cycle context the condition reads (extend as variables are needed)."
struct EventCtx
    cycle::Int        # FVS CYCLE: 1-based cycle index
    year::Int         # calendar year of this cycle
    state             # StandState, for stand variables (BBA/SDI/TPA/…)
end

# --- evaluator (single recursive function, no closures) ---------------------
function eval_event(n::EvNode, ctx::EventCtx)::Float32
    if n isa EvNum
        return n.v
    elseif n isa EvVar
        return _event_var(n.name, ctx)
    elseif n isa EvUn
        x = eval_event(n.a, ctx)
        return n.op === :neg ? -x : (x == 0f0 ? 1f0 : 0f0)   # :neg | :not
    elseif n isa EvFun
        x = eval_event(n.a, ctx)
        n.name === :FRAC && return x - trunc(x)
        n.name === :INT  && return trunc(x)
        n.name === :ABS  && return abs(x)
        n.name === :EXP  && return exp(x)
        n.name === :SQRT && return sqrt(x)
        n.name === :ALOG && return log(x)
        n.name === :MOD  && (y = eval_event(n.b, ctx); return x - trunc(x / y) * y)
        return x   # RANN etc.: not yet exercised
    elseif n isa EvTime
        # algevl.f:303 — TIME(v0,y1,v1,y2,v2,…): result v0 while year<y1, else the value vk for the
        # largest critical year yk ≤ year (years assumed increasing; stop at the first yk>year).
        # ≤2 args ⇒ v0 (IF(J.LE.2)). NOT the current year — that is the YEAR variable (evtstv.f:101).
        a = n.args
        result = eval_event(a[1], ctx)
        k = 2
        while k + 1 <= length(a)
            ctx.year >= trunc(Int, eval_event(a[k], ctx)) || break   # IYRCUR.GE.IFIX(year)
            result = eval_event(a[k + 1], ctx)
            k += 2
        end
        return result
    else  # EvBin
        b = n::EvBin
        a = eval_event(b.a, ctx); c = eval_event(b.b, ctx)
        op = b.op
        op === :add && return a + c
        op === :sub && return a - c
        op === :mul && return a * c
        # Divide by zero: FVS flags the result UNDEFINED (algevl.f:332-333 LREG=.TRUE.), not ±Inf — an
        # undefined operand makes the enclosing IF condition false (action skipped). jl has no undefined-flag
        # stack, so return NaN: a NaN comparison is false, approximating FVS's "condition does not fire" for the
        # common `IF a/0 > k` case. (Full LREG propagation — e.g. NOT(undefined) — is not modeled; niche.)
        op === :div && return c == 0f0 ? NaN32 : a / c
        op === :pow && return a ^ c                          # ** (algevl.f:339 XREG**XREG)
        op === :and && return (a != 0f0 && c != 0f0) ? 1f0 : 0f0
        op === :or  && return (a != 0f0 || c != 0f0) ? 1f0 : 0f0
        r = op === :gt ? a >  c : op === :ge ? a >= c :
            op === :lt ? a <  c : op === :le ? a <= c :
            op === :eq ? a == c : a != c               # :ne
        return r ? 1f0 : 0f0
    end
end

"""
    snapshot_compute!(s, year, cycle) -> Vector{Tuple{String,Float32}}

Evaluate the COMPUTE definitions active at `year` (date ≤ year) against the current
(start-of-cycle) stand, updating `control.compute_vars`, and return the active variables in
declaration order — the per-cycle row for the DBS FVS_Compute table (dbscmpu.f). Mirrors the
COMPUTE evaluation in `cuts!` (same start-of-cycle values), so the IF conditions and the table
agree. `cycle` is the 0-based FVSjl cycle (the event monitor's CYCLE is `cycle + 1`).
"""
function snapshot_compute!(s::StandState, year::Integer, cycle::Integer)
    defs = s.control.compute_defs
    isempty(defs) && return Tuple{String,Float32}[]
    ctx = EventCtx(Int(cycle) + 1, Int(year), s)
    fvscyc = Int(cycle) + 1
    snap = Tuple{String,Float32}[]
    for (cd, nm, ast) in defs
        # A COMPUTE block is a scheduled activity (EVUSRV → OPNEW act 33) that fires ONLY at its
        # date (IDT default 1 = cycle 1; IDT=0 = all cycles), NOT every cycle. So a default
        # `COMPUTE  MYCYC = CYCLE` evaluates once (MYCYC frozen at 1), it does NOT track the cycle.
        # The variable then persists in compute_vars for later IF conditions. (Debug-stamp of live
        # evmon proved MYCYC stays 1 ⇒ FRAC(MYCYC/3)=0.333 ⇒ its THEN thin never fires.)
        _compute_due(Int(cd), s, Int(year), fvscyc) || continue
        v = eval_event(ast, ctx)
        s.control.compute_vars[nm] = v
        push!(snap, (nm, v))
    end
    return snap
end

"""
    _compute_due(cd, s, yr, fvscyc) -> Bool

Whether a COMPUTE def scheduled for date `cd` fires this cycle (mirrors OPNEW/OPCYCL): `cd == 0`
= all cycles; `0 < cd < 1000` = a 1-based CYCLE number (fire when the current cycle == cd); else a
calendar year (fire in the cycle whose [start,end) range contains it).
"""
function _compute_due(cd::Integer, s::StandState, yr::Integer, fvscyc::Integer)::Bool
    cd == 0 && return true
    (0 < cd < 1000) && return Int(cd) == Int(fvscyc)
    cyc0 = Int(s.control.cycle)
    cs = cycle_year_at(s.control, cyc0); ce = cycle_year_at(s.control, cyc0 + 1)
    ce <= cs && (ce = cs + 1)
    return cs <= Int(cd) < ce
end

"Resolve an event-monitor variable to its current value. Extend as scenarios need."
function _event_var(name::AbstractString, ctx::EventCtx)::Float32
    if ctx.state !== nothing                       # COMPUTE user variables (checked first)
        cv = ctx.state.control.compute_vars
        (!isempty(cv) && haskey(cv, name)) && return cv[name]
    end
    name == "CYCLE" ? Float32(ctx.cycle) :
    name == "YEAR"  ? Float32(ctx.year) :
    (name == "NO" || name == "ALL") ? 0f0 :    # constant 0.0 (evtstv.f:82,281, code 112)
    name == "YES" ? 1f0 :                       # constant 1.0 (evtstv.f:81, code 111)
    name == "BBA"  ? stand_ba(ctx.state) / ctx.state.plot.gross_space :
    name == "BSDI" ? _event_bsdi(ctx.state) :                       # SDIBC (Reineke SDI, evtstv.f:116)
    # SSTAGE structural-stage event vars (evtstv.f:203-229; need STRCLASS on). The conditions
    # evaluate pre-thin, so the before/after-thin pairs read the same current stand here.
    (name == "BSCLASS" || name == "ASCLASS") ? Float32(structure_class(ctx.state).class) :
    (name == "BSTRDBH" || name == "ASTRDBH") ? Float32(structure_class(ctx.state).strdbh) :
    (name == "BCANCOV" || name == "ACANCOV") ? Float32(structure_class(ctx.state).cover) :
    name == "TPA"  ? stand_tpa(ctx.state) / ctx.state.plot.gross_space :
    # AGE = initial stand age + ELAPSED years (evtstv.f:260 TSTV1(2)=IAGE+IY(ICYC)−IY(1)); `stand_age` is
    # the fixed inventory age, so add (current year − inventory year), exactly as the .sum age does
    # (summary.jl). Was the bare `stand_age` — it omitted the elapsed term (a GAP; untested keyword path).
    name == "AGE"  ? Float32(Int(ctx.state.plot.stand_age) + (ctx.year - Int(ctx.state.control.cycle_year[1]))) :
    error("event-monitor variable not yet ported: $name")
end

# BSDI = SDIBC (evtstv.f:116) = the before-cut REINEKE stand density index — the
# `SDIC = SPROB·A + B·SDSQ` Taylor form (sdical.f:281-327, via SDICLS). This is a DIFFERENT
# SDI from the reported `.sum` column (SDIBC2 = Zeide summation) AND the mortality SDImax —
# FVS keeps all three. The event monitor reports SDIBC RAW (NOT divided by GROSPC, unlike
# BBA/TPA — evtstv.f:285 has no /GROSPC), so `t.tpa` is used directly (already gross-scaled).
# (Earlier this returned the BA — a copy-paste bug from BBA.)
function _event_bsdi(s::StandState)::Float32
    t = s.trees
    sprob = 0f0; sdsq = 0f0
    @inbounds for i in 1:t.n
        sprob += t.tpa[i]; sdsq += t.dbh[i]^2 * t.tpa[i]
    end
    sprob <= 0f0 && return 0f0
    mdsq = sdsq / sprob
    a = 10f0^(-1.605f0) * (1f0 - 1.605f0 / 2f0) * mdsq^(1.605f0 / 2f0)
    b = 10f0^(-1.605f0) * (1.605f0 / 2f0) * mdsq^(1.605f0 / 2f0 - 1f0)
    return sprob * a + b * sdsq
end

# --- tokenizer --------------------------------------------------------------
function _ev_tokens(s::AbstractString)
    toks = String[]; i = 1; n = lastindex(s)
    while i <= n
        c = s[i]
        if isspace(c)
            i = nextind(s, i)
        elseif c == '*' && i < n && s[nextind(s, i)] == '*'   # ** exponentiation (algcmp.f:103, precedence 8)
            push!(toks, "**"); i = nextind(s, nextind(s, i))
        elseif c in ('(', ')', ',', '+', '-', '*', '/')
            push!(toks, string(c)); i = nextind(s, i)
        elseif isdigit(c) || c == '.'
            j = i
            while j <= n && (isdigit(s[j]) || s[j] == '.' || s[j] == 'E' || s[j] == 'e'); j = nextind(s, j); end
            push!(toks, s[i:prevind(s, j)]); i = j
        else
            j = i
            while j <= n && (isletter(s[j]) || isdigit(s[j])); j = nextind(s, j); end
            j == i && (j = nextind(s, j))
            push!(toks, uppercase(s[i:prevind(s, j)])); i = j
        end
    end
    return toks
end

const _EV_CMP = Dict("GT"=>:gt, "GE"=>:ge, "LT"=>:lt, "LE"=>:le, "EQ"=>:eq, "NE"=>:ne)
const _EV_FUN = Dict("FRAC"=>:FRAC, "INT"=>:INT, "MOD"=>:MOD, "EXP"=>:EXP,
                     "SQRT"=>:SQRT, "ALOG"=>:ALOG, "ABS"=>:ABS, "RANN"=>:RANN)

# --- recursive-descent parser → AST -----------------------------------------
mutable struct _EvP; toks::Vector{String}; pos::Int; end
_peek(p) = p.pos <= length(p.toks) ? p.toks[p.pos] : ""
_next!(p) = (t = _peek(p); p.pos += 1; t)

function _ev_or(p)
    a = _ev_and(p)
    while _peek(p) == "OR"; _next!(p); a = EvBin(:or, a, _ev_and(p)); end
    return a
end
function _ev_and(p)
    a = _ev_not(p)
    while _peek(p) == "AND"; _next!(p); a = EvBin(:and, a, _ev_not(p)); end
    return a
end
function _ev_not(p)
    # Only NOT is the negation operator (algkey.f CTAB3). `NO`/`ALL` are NOT operators — they are the
    # constant 0.0 (test-var code 112, evtstv.f:82,281), resolved as a variable below.
    _peek(p) == "NOT" && (_next!(p); return EvUn(:not, _ev_not(p)))
    return _ev_cmp(p)
end
function _ev_cmp(p)
    a = _ev_add(p); op = _peek(p)
    haskey(_EV_CMP, op) && (_next!(p); return EvBin(_EV_CMP[op], a, _ev_add(p)))
    return a
end
function _ev_add(p)
    a = _ev_mul(p)
    while _peek(p) in ("+", "-"); op = _next!(p); a = EvBin(op == "+" ? :add : :sub, a, _ev_mul(p)); end
    return a
end
function _ev_mul(p)
    a = _ev_unary(p)
    while _peek(p) in ("*", "/"); op = _next!(p); a = EvBin(op == "*" ? :mul : :div, a, _ev_unary(p)); end
    return a
end
function _ev_unary(p)
    _peek(p) == "-" && (_next!(p); return EvUn(:neg, _ev_unary(p)))
    _peek(p) == "+" && (_next!(p); return _ev_unary(p))
    return _ev_pow(p)
end
# Exponentiation `**` (algcmp.f:103 precedence 8 — binds TIGHTER than unary minus `7`, so `-a**b` = `-(a**b)`;
# RIGHT-associative like Fortran, and the exponent may carry a unary sign, e.g. `a**-2`). algevl.f:339.
function _ev_pow(p)
    a = _ev_atom(p)
    _peek(p) == "**" && (_next!(p); return EvBin(:pow, a, _ev_unary(p)))
    return a
end
function _ev_atom(p)
    t = _next!(p)
    if t == "("
        a = _ev_or(p); _peek(p) == ")" && _next!(p); return a
    elseif t == "TIME"                       # variadic year-step fn (algevl.f:303) — keep ALL args
        _peek(p) == "(" && _next!(p)
        args = EvNode[_ev_or(p)]
        while _peek(p) == ","; _next!(p); push!(args, _ev_or(p)); end
        _peek(p) == ")" && _next!(p)
        return EvTime(args)
    elseif haskey(_EV_FUN, t)
        _peek(p) == "(" && _next!(p)
        a = _ev_or(p); b = EvNum(0f0)
        while _peek(p) == ","; _next!(p); b = _ev_or(p); end
        _peek(p) == ")" && _next!(p)
        return EvFun(_EV_FUN[t], a, b)
    else
        num = tryparse(Float32, t)
        return num !== nothing ? EvNum(num) : EvVar(t)
    end
end

"Parse an IF condition expression into an `EvNode` AST."
parse_event_condition(expr::AbstractString) = _ev_or(_EvP(_ev_tokens(expr), 1))
