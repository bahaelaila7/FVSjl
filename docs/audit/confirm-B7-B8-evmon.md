# Confirm B7/B8 — Event Monitor: TIME and NO

Two independent sub-flags in `src/engine/event_monitor.jl`. Both **CONFIRMED** as
semantic bandaids. Both currently **unreachable** by the keyword-coverage suite
(no test `.key` uses `TIME(` or bare `NO`/`ALL` in an IF/COMPUTE expression).

---

## B7 — `TIME` mapped to current year vs. the 103NN year-indexed step function

### Julia logic
`src/engine/event_monitor.jl`
- Tokenizer/funtab register `TIME` as a function: `:48` (eval) + `:155` (`_EV_FUN["TIME"=>:TIME]`).
- Evaluation (`:48`):
  ```julia
  n.name === :TIME && return Float32(ctx.year)
  ```
  i.e. TIME() unconditionally returns the calendar year of the current cycle and
  **ignores all of its arguments**.
- Parser (`_ev_atom`, `:200-205`): an `EvFun` keeps only `a` (first arg) and `b`
  (last comma arg — each `,` overwrites `b`, `:203`). So a multi-arg `TIME(v0,y1,v1,y2,v2,…)`
  cannot even be *represented*, let alone evaluated. The arguments are structurally discarded.

Note: the same value (`ctx.year`) is what the genuine `YEAR` variable returns
(`:98`, and the `YEAR` opcode 101 → `TSTV1(1)=IY(ICYC)` in evtstv.f:259). So Julia's
`TIME` is just an alias of `YEAR`.

### FVS logic
- `algkey.f:345/350`: `TIME` is in CTAB4 → `IOPT4 = 10300`. The `103NN` family
  (comment `algkey.f:289`) is a multi-argument function, dispatched in ALGEVL at
  label `303`.
- `algevl.f:429-460` (label 303) — `TIME(v0, y1, v1, y2, v2, …)` is a **year-indexed
  step interpolation**:
  - Start `XREG = arg1` (the default value, `:432`).
  - Walk the `(year, value)` pairs (`DO 310 … ,2`): when the current calendar year
    `IYRCUR >= IFIX(year)` (`:451`), advance the result to the value following that
    year (`XREG(IXSTK)=XREG(NDC+1)`, `:452`).
  - Result = the value whose year-threshold the current year has reached
    (a left-continuous step function of `IYRCUR`). If only one arg (`J<=2`), it
    returns that single default (`:430`).
  - `IYRCUR` is the **selector**, never the return value.

### Semantic diff
FVS `TIME` returns a user-specified *value schedule* keyed by year thresholds
(e.g. `TIME(0, 2010, 1, 2030, 2)` → 0 before 2010, 1 in 2010–2029, 2 from 2030).
Julia returns the year integer itself and throws every argument away. For any real
use the magnitudes differ by orders (years ~2000 vs scheduled values ~0–N), and the
intended time-varying behavior is entirely absent. This is a stub, not a port.

### Faithful fix
Two changes (parser + eval), because the current `EvFun{a,b}` cannot hold the
variadic arg list:
1. Make `EvFun` carry all arguments (e.g. `args::Vector{EvNode}`); in `_ev_atom`
   collect every comma-separated sub-expression instead of overwriting `b`.
2. Evaluate TIME per algevl 303:
   ```julia
   # args = [v0, y1, v1, y2, v2, ...]
   res = eval_event(args[1], ctx)
   i = 2
   while i + 1 <= length(args)
       eval_event(args[i], ctx) <= Float32(ctx.year) && (res = eval_event(args[i+1], ctx))
       i += 2
   end
   return res
   ```
   (`year >= threshold` ⇔ `threshold <= ctx.year`; matches `IYRCUR.GE.IFIX(XREG)`.)
   Keep `DECADE` (10200, label 302) in mind as the sibling if later needed.

### Upstream rank — MID
Consumed by `eval_event` → IF-condition truth and COMPUTE values
(`snapshot_compute!`, and the `cuts!`/scheduling path). When a keyword uses TIME to
gate an activity (thin/fertilize) the wrong value mis-fires the activity, altering
the stand trajectory (upstream). When TIME feeds only a reported COMPUTE var it is
LEAF. Net: MID.

### Reachability
Not exercised. No suite `.key` contains `TIME(`. A faithful fix is silent on the
current suite; it would only move results once a scenario using TIME is added.

---

## B8 — `NO` parsed as `NOT` vs. the constant 0.0 variable (opcode 112)

### Julia logic
`src/engine/event_monitor.jl:172-175` (`_ev_not`):
```julia
(_peek(p) == "NOT" || _peek(p) == "NO") && (_next!(p); return EvUn(:not, _ev_not(p)))
```
`NO` is treated as a **prefix logical-negation operator** identical to `NOT`
(`EvUn(:not, …)`, evaluated at `:39` as `x==0 ? 1 : 0`). It consumes and negates the
following expression.

### FVS logic
- `algkey.f:331-332`: `NO` is in CTAB2 → `IOPT2 = 112`. Code 112 is a **variable**
  ("NO" and "ALL" both equate to the constant 0.0 — `algkey.f:107-109`).
- `evtstv.f:281`: `TSTV1(12) = 0.0`; opcode `10N → TSTV1(N)`, so 112 → `TSTV1(12) = 0.0`.
  `NO` (and `ALL`) is therefore a leaf **constant 0.0**, not an operator.
- `NOT` is the operator: `algkey.f:335/339`, CTAB3 → `IOPT3 = 001` (opcode 1 = NOT,
  `algkey.f:70`). It is a distinct token.

### Semantic diff
FVS: `NO` is an atom equal to `0.0` (and `ALL` likewise; `YES` = 1.0 via
`TSTV1(11)=1.0`, opcode 111). Julia: `NO` is a unary negation that consumes the next
operand. E.g. `IF FLAG EQ NO` → FVS compares FLAG to 0.0; Julia would mis-parse `NO`
as expecting an operand to negate. And `IF NO BBA GT 100` → Julia yields
`NOT(BBA GT 100)`; FVS yields the atom `0.0` followed by a syntax/eval that never
negates. Different semantics in every case.

### Faithful fix
- Remove the `|| _peek(p) == "NO"` clause from `_ev_not` (`:173`); leave only `NOT`.
- Treat `NO` and `ALL` as constant `0.0`, and `YES` as constant `1.0`, in
  `_ev_atom` (or pre-seed them as numeric tokens):
  ```julia
  t == "NO"  || t == "ALL" ? EvNum(0f0) :
  t == "YES"               ? EvNum(1f0) : ...
  ```
  (Currently `YES`/`ALL` also fall through to `EvVar` and would `error` in
  `_event_var` — same family of miss; fixing all three matches algkey opcodes
  111/112.)

### Upstream rank — MID
Same consumer chain as B7: `eval_event` → IF truth / COMPUTE. A condition using `NO`
or `ALL` as a constant would gate activities differently (upstream when it drives a
cut/treatment, LEAF when only reported). Net: MID.

### Reachability
Not exercised. No suite `.key` uses `NO`/`ALL` as an event-monitor expression token.
Silent path; fix won't move the current suite.

---

## Masked-bug watch
Neither fix can regress a current test (both paths are dead in the suite). The B7
parser change (variadic `EvFun`) touches the shared function-arg path used by `MOD`
and any future multi-arg function; it must preserve the existing 1-/2-arg behavior
(`FRAC/INT/ABS/EXP/SQRT/ALOG`, and `MOD(x,y)`), so verify those keep `args[1]`
(and `args[2]` for MOD) semantics after the refactor. No masked downstream bug is
implied — these are pure stubs, not compensations.
