# test_event_monitor.jl — FVS event monitor (IF/THEN/ENDIF). The algebraic-condition
# evaluator (parse_event_condition + eval_event over the ALGKEY grammar) and the
# per-cycle firing of an IF block's activities (cuts!). snt01 stand 2 thins with a
# `THINDBH` block gated by `IF (FRAC(CYCLE/3.0) EQ 0.0) THEN … ENDIF`.

using Test, FVSjl

# minimal context: only CYCLE matters for these conditions
_ctx(cycle) = FVSjl.EventCtx(cycle, 1990 + (cycle - 1) * 5, nothing)
_evalcond(expr, cycle) = FVSjl.eval_event(FVSjl.parse_event_condition(expr), _ctx(cycle)) != 0f0

@testset "event monitor — algebraic condition evaluator" begin
    # snt01 stand 2's condition: fire every 3rd cycle (FVS CYCLE is 1-based)
    fires = [c for c in 1:10 if _evalcond("(FRAC(CYCLE/3.0) EQ 0.0)", c)]
    @test fires == [3, 6, 9]

    # operators / functions / precedence
    @test _evalcond("CYCLE GE 3", 3) && !_evalcond("CYCLE GE 3", 2)
    @test _evalcond("CYCLE GT 2 AND CYCLE LT 5", 3)
    @test !_evalcond("CYCLE GT 2 AND CYCLE LT 5", 6)
    @test _evalcond("CYCLE EQ 2 OR CYCLE EQ 4", 4)
    @test _evalcond("NOT CYCLE EQ 3", 2) && !_evalcond("NOT CYCLE EQ 3", 3)
    @test _evalcond("MOD(CYCLE,2) EQ 0", 4) && !_evalcond("MOD(CYCLE,2) EQ 0", 3)
    @test _evalcond("INT(CYCLE/2) EQ 2", 4)            # 4/2=2
    @test _evalcond("CYCLE * 2 GT 5 + 1", 4)           # 8 > 6
end

const _EM_KEY = "/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key"

@testset "event monitor — snt01 stand 2 IF/THEN fires the THINDBH block" begin
    if !isfile(_EM_KEY)
        @test_skip "snt01.key not available"
    else
        stands = each_stand(_EM_KEY)
        s2 = stands[2]
        @test !isempty(s2.control.conditionals)        # the IF block was parsed
        c = s2.control.conditionals[1]
        @test !isempty(c.acts)                          # THINDBH activities captured
        # the condition fires at cycles 3/6/9 (years 2000/2015/2030)
        @test FVSjl.eval_event(c.cond, _ctx(3)) != 0f0
        @test FVSjl.eval_event(c.cond, _ctx(4)) == 0f0
        # project: stand 2 must thin (TPA at 2030 well below the unthinned control)
        notre!(s2); FVSjl.setup_growth!(s2); FVSjl.compute_volumes!(s2)
        io = IOBuffer(); FVSjl.write_sum_file(io, s2)
        tpa2030 = 0.0
        for ln in eachline(IOBuffer(String(take!(io))))
            t = split(strip(ln)); length(t) >= 3 && t[1] == "2030" && (tpa2030 = parse(Float64, t[3]))
        end
        # the IF/THEN THINDBH thins at cycles 3/6/9 → managed density 257 @2030 (vs the
        # self-thinned unthinned control's 171), validated bit-exact for the 1st two thins
        @test isapprox(tpa2030, 257.0; atol = 3)
    end
end

@testset "event monitor — BSDI is the Reineke SDIBC, not the BA" begin
    # Regression for the copy-paste bug where BSDI returned the basal area. BSDI = SDIBC
    # (the raw Reineke SDIC = SPROB·A + B·SDSQ, sdical.f:281-327) — a DIFFERENT SDI from the
    # `.sum` Zeide column and the mortality SDImax. Reported RAW (no /GROSPC, unlike BBA/TPA).
    key = "/workspace/FVSjl/test/harness/scenarios/fire_early.key"
    if !isfile(key)
        @test_skip "fire_early.key not available"
    else
        s, _ = FVSjl.initialize(key)
        notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s); FVSjl.compute_density!(s)
        bsdi = FVSjl._event_bsdi(s)
        ctx = FVSjl.EventCtx(1, 1990, s)
        @test FVSjl._event_var("BSDI", ctx) ≈ bsdi          # the variable dispatches to BSDI
        @test FVSjl._event_var("BSDI", ctx) != FVSjl._event_var("BBA", ctx)   # not the BA
        @test isapprox(bsdi, 202.9; atol = 0.2)             # bit-exact vs live Fortran COMPUTE BSDI
    end
end

@testset "event monitor — NO/YES constants + TIME step function (algkey/algevl/evtstv)" begin
    _val(expr, cycle) = FVSjl.eval_event(FVSjl.parse_event_condition(expr), _ctx(cycle))
    # B8: NO/ALL = constant 0.0 (evtstv.f:82,281), YES = 1.0 (evtstv.f:81) — NOT a logical-NOT operator.
    @test _val("NO", 1) == 0f0
    @test _val("ALL", 1) == 0f0
    @test _val("YES", 1) == 1f0
    @test _evalcond("YES", 1) && !_evalcond("NO", 1)         # used as truth values
    @test _evalcond("NOT CYCLE EQ 3", 2)                     # NOT is still the negation operator
    # B7: TIME(v0,y1,v1,…) year-indexed step fn (algevl.f:303), NOT the current year (that is YEAR).
    @test _val("TIME(5)", 1) == 5f0                          # ≤2 args ⇒ v0
    # _ctx years: cycle1=1990, cycle3=2000, cycle5=2010, cycle7=2020
    step = "TIME(1, 2000, 2, 2010, 3)"
    @test _val(step, 1) == 1f0                               # 1990 < 2000 ⇒ v0
    @test _val(step, 3) == 2f0                               # 2000 ≤ yr < 2010 ⇒ v1
    @test _val(step, 5) == 3f0                               # yr ≥ 2010 ⇒ v2
    @test _val(step, 7) == 3f0                               # stays at last
    @test _val("YEAR", 3) == 2000f0                          # current year is YEAR, not TIME
end

@testset "event monitor — ** exponentiation (algcmp.f:103, algevl.f:339)" begin
    _val(expr) = FVSjl.eval_event(FVSjl.parse_event_condition(expr), _ctx(1))
    # FVS supports `**` (precedence 8 — higher than unary minus 7, higher than * / 6; RIGHT-associative like
    # Fortran). jl previously lacked it (tokenized `**` as two `*` ⇒ wrong). The prior "FVS has no **" audit
    # verdict was a misread of algcmp.f.
    @test _val("2 ** 3")      == 8f0
    @test _val("3 ** 2")      == 9f0
    @test _val("-2 ** 2")     == -4f0       # ** binds tighter than unary minus ⇒ −(2²)
    @test _val("2 ** -1")     == 0.5f0      # exponent may carry a unary sign
    @test _val("2 ** 3 ** 2") == 512f0      # right-associative ⇒ 2^(3^2) = 2^9
    @test _val("10 ** 2 + 1") == 101f0      # ** > +
    @test _val("2 * 3 ** 2")  == 18f0       # ** > *
end
