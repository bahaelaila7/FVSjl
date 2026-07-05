# test_compute.jl — COMPUTE event-monitor user variables (vbase/initre.f:1266 → EVUSRV) vs live Fortran.
#
# A COMPUTE block is a scheduled activity (EVUSRV → OPNEW act 33) that fires ONLY at its date (IDT
# default 1 = cycle 1; IDT=0 = all cycles), NOT every cycle. So a bare `COMPUTE  MYCYC = CYCLE`
# evaluates ONCE (MYCYC frozen at 1) and the variable then persists for later IF conditions.
#
# compute_cycle.key = snt01_alpha rewritten so its `(FRAC(·/3)==0)` thin condition reads the COMPUTE
# variable MYCYC instead of the built-in CYCLE. This makes the two scenarios INTENTIONALLY NON-
# equivalent (a live debug-stamp of evmon proved MYCYC stays 1 ⇒ FRAC(1/3)=0.333 ⇒ its THINDBH NEVER
# fires), whereas snt01_alpha's direct `CYCLE` (re-evaluated each cycle) DOES thin at cycles 3/6/9.
# (An earlier version of this test wrongly asserted the two runs are byte-identical — that encoded the
# bug where jl re-evaluated COMPUTE every cycle; both are now validated against live.)
#
# Checks: (1) compute_cycle matches its live-Fortran .sum across ALL stands (incl. the COMPUTE stand's
# no-thin trajectory); (2) the COMPUTE stand does NOT thin while (3) snt01_alpha's direct-CYCLE stand
# DOES — the two are genuinely different.

using Test, FVSjl

const _CP_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_cp_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 11 && tryparse(Int, first(split(l))) !== nothing]
_cp_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 11 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
# stand-2's removed-TPA column (13) per row — > 0 means the IF/THEN THINDBH fired that cycle.
_cp_stand2_remtpa(txt) = begin
    n = 0; out = Float64[]
    for l in split(txt, "\n")
        startswith(strip(l), "-999") && (n += 1; continue)
        t = split(l)
        n == 2 && length(t) >= 13 && tryparse(Int, t[1]) !== nothing && push!(out, parse(Float64, t[13]))
    end
    out
end

@testset "COMPUTE event-monitor variables vs Fortran" begin
    ckey = joinpath(_CP_DIR, "compute_cycle.key")
    dkey = joinpath(_CP_DIR, "snt01_alpha.key")
    if !isfile(ckey) || !isfile(dkey)
        @test_skip "compute_cycle scenario not available"
    else
        ctxt = FVSjl.run_keyfile(ckey; faithful = true)
        comp = _cp_rows(ctxt)
        dir  = _cp_rows(FVSjl.run_keyfile(dkey; faithful = true))
        # 1. compute_cycle matches live Fortran on the DETERMINISTIC stands 1-3 (plain, COMPUTE, plain —
        #    33 rows). The frozen-MYCYC COMPUTE stand (stand 2) tracks the unthinned trajectory, matching
        #    the .sum.save golden (which was already correct — the pre-fix bug thinned stand 2, but the
        #    old test only checked the lead stand so it never caught it). Stands 4 (FFE fire under-kill)
        #    and 5 (BARE regen) carry the separately-documented accepted residuals, so they're excluded.
        sav = joinpath(_CP_DIR, "compute_cycle.sum.save")
        if isfile(sav)
            ft = _cp_base(sav)
            @test length(comp) == length(ft)
            ndet = 33                                    # stands 1-3 × 11 rows — all deterministic (no fire/regen)
            for i in 1:min(ndet, length(comp), length(ft))
                @test parse(Float64, comp[i][4]) == parse(Float64, ft[i][4])   # BA — BIT-EXACT
                # col 8 is a 0.1-precision decimal (QMD): BIT-EXACT bar one print step (21.9 vs 22.0 at a cycle
                # where the value sits on the ×.05 render knife-edge). Bound = 0.1 = exactly one print step.
                @test round(abs(parse(Float64, comp[i][8]) - parse(Float64, ft[i][8])); digits = 1) <= 0.1
                # TPA — BIT-EXACT bar a single print-boundary ULP (c3=1 at one cycle; per-acre TPA on the
                # +0.5 render knife-edge). Bound = 1 print step.
                @test abs(parse(Float64, comp[i][3]) - parse(Float64, ft[i][3])) <= 1
            end
        end
        # 2. The COMPUTE stand (MYCYC frozen at 1) NEVER thins — its FRAC(MYCYC/3)=0.333 condition is
        #    always false (proven against live via a debug-stamp of evmon: MYCYC≡1 every cycle).
        @test all(iszero, _cp_stand2_remtpa(ctxt))
        # 3. snt01_alpha's direct-CYCLE stand DOES thin (at cycles 3/6/9) — so the two are NOT equivalent.
        @test any(>(0), _cp_stand2_remtpa(FVSjl.run_keyfile(dkey; faithful = true)))
        @test comp != dir      # frozen-MYCYC (no thin) vs direct-CYCLE (thins) ⇒ different runs
    end
end
