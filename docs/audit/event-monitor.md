# Bandaid audit — Event monitor (`src/engine/event_monitor.jl`)

Audited against FVS `base/algevl.f`, `base/algcmp.f`, `base/algkey.f`, `base/evmon.f`,
and `bin/FVSsn_buildDir/{evtstv.f, sdical.f, grinit.f}`.

The core evaluator and the recursive-descent parser are mostly a faithful re-expression of
ALGCMP/ALGEVL. The arithmetic op set (`add/sub/mul/div`), the comparison set
(`gt/ge/lt/le/eq/ne` = opcodes 6/7/8/9/4/5), `FRAC = x-trunc(x)` (algevl.f:269),
`INT = trunc` (algevl.f:272), `MOD = x-trunc(x/y)*y` (= Fortran `AMOD`, algevl.f:409),
`ALOG = natural log` (algevl.f:263), the operator-precedence ladder
`OR<AND<NOT<CMP<ADD<MUL<unary‑minus` (matches `IPRMAP`/comments algcmp.f:86–112), and the
unary `+`/`-` handling (algcmp.f:188–218) all check out. `BBA=BA/GROSPC` (evtstv.f:268),
`TPA=TPROB/GROSPC` (evtstv.f:263), `CYCLE`/`YEAR`, and the start-of-cycle COMPUTE snapshot
are all faithful. Six items below are flagged.

---

## 1. BANDAID — `TIME` is mapped to "current year"

- **jl:** `eval_event`, line 48: `n.name === :TIME && return Float32(ctx.year)`; `TIME` is also
  registered as a 1-arg function in `_EV_FUN` (line 155).
- **Claim:** that `TIME(...)` yields the calendar year. No source citation.
- **FVS checked:** `algkey.f:345/350` — `TIME` is opcode `10300`, i.e. a *multi-argument*
  function (the `103NN` family). `algevl.f:429–460` implements it as a **year-indexed step
  interpolation**: `TIME(v0, y1, v1, y2, v2, …)` returns `v0` while `IYRCUR<y1`, `v1` once
  `IYRCUR>=y1`, etc. With `J<=2` args (algevl.f:430) it returns the *first argument*, never the
  year. FVS has no scalar "TIME" variable that returns the current year.
- **Severity:** BANDAID. The mapping contradicts the FVS source. `TIME(...)` will silently
  return the calendar year regardless of its arguments. (The "current year" intent is what the
  `YEAR` variable already provides, evtstv.f:259.) Likely unexercised in current SN scenarios,
  but the behavior is wrong if any keyfile uses `TIME`.

## 2. BANDAID — `"NO"` parsed as logical NOT

- **jl:** `_ev_not`, line 173: `(_peek(p) == "NOT" || _peek(p) == "NO") && (... EvUn(:not, …))`.
- **Claim:** that `NO` is a synonym for the logical negation operator.
- **FVS checked:** `algkey.f:331–332` — `CTAB2` lists `'NO'` with `IOPT2 = 112`, i.e. `NO` is
  **test-variable code 112 = the constant 0.0** (evtstv.f:81–82, "112 NO A CONSTANT EQUAL TO 1/0";
  `TSTV1(12)=0.0` at evtstv.f:281). The negation operator is the distinct token `'NOT'`
  (`CTAB3`, opcode 1, algkey.f:335/339). `NO` is the boolean complement *value* of `YES` (=1.0),
  not an operator.
- **Severity:** BANDAID. Treating `NO` as NOT directly contradicts algkey.f. An expression using
  `NO` as the constant 0 would be mis-parsed into a unary negation of the following sub-expression.

## 3. GAP — exponentiation `**` is unsupported and silently mis-parsed

- **jl:** tokenizer `_ev_tokens` (line 137) only treats `+ - * / ( ) ,` as operator chars; there
  is no `**` (opcode 15) or `^` handling, and the parser has no power level.
- **FVS checked:** `algcmp.f:234–250` explicitly scans `*` followed by `*` and emits opcode **15
  (RAISE TO A POWER)** with precedence **8** (algcmp.f:103, higher than unary minus); `algevl.f:338–340`
  evaluates `XREG**XREG`. Exponentiation is a first-class operator in both IF conditions and
  COMPUTE expressions.
- **Effect:** `a**b` tokenizes to `a`,`*`,`*`,`b`; `_ev_mul` consumes the first `*`, `_ev_unary`/`_ev_atom`
  then consumes the second `*` as a token and (failing numeric parse) builds `EvVar("*")` — i.e.
  silent garbage / a spurious "variable not ported" error, not the intended power.
- **Severity:** GAP. The default tested SN conditions evidently avoid `**`, but any COMPUTE/IF
  using exponentiation is silently mishandled.

## 4. GAP — `AGE` omits the elapsed-years term (wrong for cycle > 1)

- **jl:** `_event_var`, line 107: `name == "AGE" ? Float32(ctx.state.plot.stand_age)`. `plot.stand_age`
  is the *input* stand age `IAGE` (core/state.jl:345), never advanced.
- **FVS checked:** `evtstv.f:260` — `TSTV1(2) = IAGE + IY(MAX(1,ICYC)) - IY(1)`, i.e. input age
  **plus elapsed calendar years** since the first cycle. The FVSjl codebase already knows this
  formula: `io/summary.jl:224` computes `age = stand_age + (yr - cycle_year[1])`.
- **Effect:** correct only at cycle 1 (elapsed = 0); for every later cycle the event-monitor `AGE`
  stays frozen at the input age instead of growing.
- **Severity:** GAP. Faithful for the first cycle, wrong thereafter; contradicts evtstv.f:260.

## 5. GAP — division/MOD by zero yields Inf/NaN instead of FVS "undefined"

- **jl:** `eval_event` `:div` (line 58) returns `a / c` (→ `Inf` when `c==0`); `:MOD` (line 49)
  returns `x - trunc(x/y)*y` (→ `NaN` when `y==0`). There is no "defined/undefined" status
  tracking at all.
- **FVS checked:** `algevl.f:332–336` (divide) and `algevl.f:406–410` (AMOD) set the operand's
  LREG status to "undefined" on a zero divisor, propagating an undefined result and returning
  `IRC=1`. `evmon.f:109` fires an event only when `LREG(1) .AND. IRC.EQ.0` — so an undefined
  condition is treated as **not occurring**.
- **Effect:** a condition like `X/0 GT 5` (or an undefined COMPUTE var) evaluates to `Inf>5 = true`
  in FVSjl and would *fire*, whereas FVS suppresses it. Also stores `Inf/NaN` into
  `compute_vars` where FVS would mark it unset.
- **Severity:** GAP (edge case; whole undefined-propagation subsystem is absent).

## 6. GAP — `_event_bsdi` ignores the DBHSTAGE / dead-tree exclusions of the SDIBC sum

- **jl:** `_event_bsdi` (lines 117–128) sums the Stage-SDI Taylor form over **all** tree records.
  The header comment (lines 111–115) correctly ties `BSDI` to `SDIBC` and to the raw (no `/GROSPC`)
  report at evtstv.f:285.
- **FVS checked:** the formula matches `sdical.f:281–283` exactly (`A`, `B`, `SDIC=SPROB*A+B*SDSQ`)
  — faithful for the default stand. **But** the returned `SDIBC` is actually the re-summed value at
  `sdical.f:292–327`, which (a) excludes trees with `DBH < DBHSTAGE` (sdical.f:269, 327) and
  (b) excludes recent/older dead-at-inventory records `I >= IREC2` (sdical.f:298). The A/B
  coefficients are likewise built only from trees passing `DBHSTAGE` (sdical.f:269).
- **Effect:** `DBHSTAGE` defaults to **0** (grinit.f:263) and is only nonzero under the `SDIMINDBH`
  keyword (initre.f:5959–5960), and dead-at-inventory records are normally absent — so the default
  path is bit-faithful. Divergence arises only when `SDIMINDBH` is set or dead-inventory trees
  exist.
- **Severity:** GAP (low; default SN path is faithful).

---

### Notes (not flagged)
- Multi-arg functions other than `MOD` (`MIN/MAX/BOUND/LININT/DECADE/NORMAL`) and the extra unary
  functions (`ALOG10`, `SIN/COS/TAN`, `ARCSIN/ARCCOS/ARCTAN`) are absent from `_EV_FUN`; they
  resolve to `EvVar` and raise an explicit "not yet ported" error rather than silently miscomputing,
  so they are coverage gaps, not bandaids.
- The tokenizer splits `%`-suffixed variable names (`DTPA%`, `DBA%`, …) into `NAME` + `%`; those
  variables are not in `_event_var` and would error anyway, so no silent miscompute.
