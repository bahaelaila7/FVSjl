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
        n.name === :TIME && return Float32(ctx.year)
        n.name === :MOD  && (y = eval_event(n.b, ctx); return x - trunc(x / y) * y)
        return x   # RANN etc.: not yet exercised
    else  # EvBin
        b = n::EvBin
        a = eval_event(b.a, ctx); c = eval_event(b.b, ctx)
        op = b.op
        op === :add && return a + c
        op === :sub && return a - c
        op === :mul && return a * c
        op === :div && return a / c
        op === :and && return (a != 0f0 && c != 0f0) ? 1f0 : 0f0
        op === :or  && return (a != 0f0 || c != 0f0) ? 1f0 : 0f0
        r = op === :gt ? a >  c : op === :ge ? a >= c :
            op === :lt ? a <  c : op === :le ? a <= c :
            op === :eq ? a == c : a != c               # :ne
        return r ? 1f0 : 0f0
    end
end

"Resolve an event-monitor variable to its current value. Extend as scenarios need."
function _event_var(name::AbstractString, ctx::EventCtx)::Float32
    if ctx.state !== nothing                       # COMPUTE user variables (checked first)
        cv = ctx.state.control.compute_vars
        (!isempty(cv) && haskey(cv, name)) && return cv[name]
    end
    name == "CYCLE" ? Float32(ctx.cycle) :
    name == "YEAR"  ? Float32(ctx.year) :
    name == "BBA" || name == "BSDI" ? stand_ba(ctx.state) / ctx.state.plot.gross_space :
    name == "TPA"  ? stand_tpa(ctx.state) / ctx.state.plot.gross_space :
    name == "AGE"  ? Float32(ctx.state.plot.stand_age) :
    error("event-monitor variable not yet ported: $name")
end

# --- tokenizer --------------------------------------------------------------
function _ev_tokens(s::AbstractString)
    toks = String[]; i = 1; n = lastindex(s)
    while i <= n
        c = s[i]
        if isspace(c)
            i = nextind(s, i)
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
                     "SQRT"=>:SQRT, "ALOG"=>:ALOG, "ABS"=>:ABS, "TIME"=>:TIME, "RANN"=>:RANN)

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
    (_peek(p) == "NOT" || _peek(p) == "NO") && (_next!(p); return EvUn(:not, _ev_not(p)))
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
    return _ev_atom(p)
end
function _ev_atom(p)
    t = _next!(p)
    if t == "("
        a = _ev_or(p); _peek(p) == ")" && _next!(p); return a
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
