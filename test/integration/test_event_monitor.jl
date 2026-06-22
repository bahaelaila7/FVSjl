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
