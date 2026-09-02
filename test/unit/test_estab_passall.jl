# =============================================================================
# test_estab_passall.jl — PASSALL (esin.f opt 18 → CONFID → estab.f PASMAX) ESTAB-packet keyword.
#
# PASSALL sets CONFID, used as PASMAX = "the MAXIMUM number of EXCESS regeneration trees passed per plot per
# species" (esin.f:551-558 echo; estab.f:171 PASMAX=CONFID). When a plot's per-species EXCESS count is passed to
# the tree list it is broken into IBRKUP=INT(EXCESS/5+1) records each carrying XCSMAX=EXCESS/BRKUP; PASMAX caps
# that to min(EXCESS/BRKUP, PASMAX/BRKUP) (estab.f:1318-1321). The record PROB scales LINEARLY with XCSMAX, so
# the cap is a purely DETERMINISTIC post-draw scaling — it consumes NO RNG draw. Default CONFID=5 (esinit.f:50);
# esin.f:557 raises CONFID<1 to 1. The AUTOES ingrowth path (estab.f:249, NTALLY==99) OVERRIDES PASMAX=15.
#
# STATUS (measured 2026-09-02 vs live FVSie_g16, single-.o estab.f dump swap):
#   • PASSALL IS model-affecting on the oracle. On a dated-ESTAB bare stand (DESIGN 11, NOTREES, ESTAB 1992)
#     PASSALL 1 → 2002 TPA 477 vs default(5)/PASSALL 100 → 555 (default==100 there: max EXCESS≤5, so only
#     PASMAX 1 caps). At PLANT densities the cap fires at EXCESS≥6 (default PASMAX=5).
#   • The cap KERNEL (es_pasmax_xcsmax) is a faithful transcription of estab.f:1288/1318-1321 and is BIT-EXACT
#     vs the instrumented oracle across EXCESS 1..6 × PASMAX {1,5} (numbers asserted below).
#   • WIRED 2026-09-02 (estab-ie-discrete): the cap is now applied in ie_autoes_tally, gated to InlandEmpire.
#     Per plot the NOTE best/excess split (estab.f:1079-1145) uses the MEASURED bit-exact rule NBEST =
#     min(ITP,max(4,NUMSPE)) (=54/54 plots on the FVSie dump); the (NBEST−NUMSPE) tallest excess are promoted to
#     best (generation-order proxy — EXACT for single-species-excess plots, cornered for the ~6% multi-species
#     ones), and each residual per-species EXCESS is capped via es_pasmax_xcsmax. The uncapped path (PASMAX ≥
#     EXCESS) is left byte-identical, so the 339/11 gate + every non-IE variant are untouched. VALIDATED on the
#     bare dated-ESTAB stand: default(5)==PASSALL 100 (max EXCESS≤5 ⇒ no cap), PASSALL 1 moves 2002 TPA 572→496
#     (oracle 555→477: cap delta −76 vs −78, i.e. −13.3% vs −14.1%). The absolute residual (jl 572 vs oracle 555
#     at default) is the PRE-EXISTING IE establishment placeholder-height/growth approximation, NOT the cap.
# =============================================================================

using Test
using FVSjl
using FVSjl: KeywordReader, read_keyword!, kw_estab!, StandState, InlandEmpire, init_blockdata!,
             es_pasmax_xcsmax

kwline(kw, fields...) = rpad(kw, 10) * join(rpad.(string.(fields), 10)) * "\n"

@testset "PASSALL ESTAB-packet keyword (esin.f opt 18 → CONFID → PASMAX)" begin
    @testset "parse + clamp (default 5.0; CONFID<1 → 1; blank → 1)" begin
        # default (no PASSALL): esinit.f CONFID=5.0
        s = StandState(InlandEmpire()); init_blockdata!(s, s.variant)
        @test s.estab.pasmax == 5.0f0

        # PASSALL 100 → PASMAX=100
        s = StandState(InlandEmpire()); init_blockdata!(s, s.variant)
        kr = KeywordReader(IOBuffer("ESTAB\n" * kwline("PASSALL", 100) * "END\n"))
        kw_estab!(s, read_keyword!(kr), kr)
        @test s.estab.pasmax == 100.0f0

        # PASSALL 1 → PASMAX=1
        s = StandState(InlandEmpire()); init_blockdata!(s, s.variant)
        kr = KeywordReader(IOBuffer("ESTAB\n" * kwline("PASSALL", 1) * "END\n"))
        kw_estab!(s, read_keyword!(kr), kr)
        @test s.estab.pasmax == 1.0f0

        # PASSALL 0.5 clamps to the 1.0 floor (esin.f:557)
        s = StandState(InlandEmpire()); init_blockdata!(s, s.variant)
        kr = KeywordReader(IOBuffer("ESTAB\n" * kwline("PASSALL", 0.5) * "END\n"))
        kw_estab!(s, read_keyword!(kr), kr)
        @test s.estab.pasmax == 1.0f0

        # blank PASSALL ⇒ ARRAY(1)=0 ⇒ CONFID=0<1 ⇒ 1.0 (esin.f has no LNOTBK guard on field 1)
        s = StandState(InlandEmpire()); init_blockdata!(s, s.variant)
        kr = KeywordReader(IOBuffer("ESTAB\nPASSALL\nEND\n"))
        kw_estab!(s, read_keyword!(kr), kr)
        @test s.estab.pasmax == 1.0f0
    end

    @testset "es_pasmax_xcsmax cap kernel is BIT-EXACT vs live FVSie (estab.f:1288/1318-1321)" begin
        # (excess, pasmax) => (xcsmax, ibrkup) — every value read off the instrumented oracle dump.
        # PASMAX=5 (default): the cap fires only at EXCESS=6 (uncapped 3.0 → 2.5).
        for (exc, xcs, brk) in ((1, 1.0f0, 1), (2, 2.0f0, 1), (3, 3.0f0, 1),
                                (4, 4.0f0, 1), (5, 2.5f0, 2), (6, 2.5f0, 2))
            x, b = es_pasmax_xcsmax(exc, 5f0)
            @test x == xcs && b == brk
        end
        # PASMAX=1: the cap fires from EXCESS=2 up (BRK=1 → cap to 1.0; BRK=2 → cap to 0.5).
        for (exc, xcs, brk) in ((1, 1.0f0, 1), (2, 1.0f0, 1), (3, 1.0f0, 1),
                                (4, 1.0f0, 1), (5, 0.5f0, 2), (6, 0.5f0, 2))
            x, b = es_pasmax_xcsmax(exc, 1f0)
            @test x == xcs && b == brk
        end
        # PASMAX=15 (the AUTOES ingrowth override): inert until EXCESS>15 (all BRK small here).
        @test es_pasmax_xcsmax(6, 15f0)  == (3.0f0, 2)   # IBRKUP=INT(6/5+1)=2; 6/2=3.0 ≤ 15/2=7.5 ⇒ uncapped
        @test es_pasmax_xcsmax(20, 15f0) == (3.0f0, 5)   # IBRKUP=INT(20/5+1)=5; 20/5=4.0 > 15/5=3.0 ⇒ cap to 3.0
    end

    @testset "PASSALL is WIRED into the IE tally seam (caps EXCESS; default/PASSALL-100 stay uncapped)" begin
        # A dated-ESTAB bare stand where the ORACLE .sum DIFFERS for PASSALL 1 (477 TPA) vs default (555). The cap
        # is now wired (ie_autoes_tally, IE-gated): PASSALL 1 MOVES the 2002 TPA down (jl 572→496, oracle 555→477),
        # while default(PASMAX=5) and PASSALL 100 stay identical (this stand's max per-plot EXCESS is ≤5 ⇒ neither
        # caps). This test FIRES the ported branch — deleting the cap wiring makes `sum_def == sum_p1` hold and the
        # inequality assertion fail.
        tpa2002(t) = begin
            v = nothing
            for l in split(t, "\n")
                f = split(l); (length(f) >= 3 && f[1] == "2002") && (v = parse(Int, f[3]))
            end
            v
        end
        mktempdir() do dir
            # column-sensitive keyfile — NO leading indentation (STDINFO/DESIGN fields are fixed-column)
            base = "STDIDENT\n" *
                   "PASSALLT\n" *
                   "DESIGN                                        11.0       1.0\n" *
                   "NOTREES\n" *
                   "STDINFO        101.0     520.0       0.0     315.0      30.0      45.0\n" *
                   "INVYEAR         1992\n" *
                   "NUMCYCLE           3\n" *
                   "TREELIST           0\n" *
                   "ESTAB           1992\n" *
                   "%PASSALL%END\n" *
                   "ECHOSUM\n" *
                   "PROCESS\n" *
                   "STOP\n"
            function runsum(passall)
                kf = joinpath(dir, "k.key")
                write(kf, replace(base, "%PASSALL%" => passall))
                FVSjl.run_keyfile(kf; variant = InlandEmpire(), output = :sum)
            end
            sum_def = runsum("")
            sum_p1  = runsum("PASSALL            1.0\n")
            sum_p100 = runsum("PASSALL          100.0\n")
            # strip the timestamp header line (col-varying date) before comparing the projection rows
            body(t) = join([l for l in split(t, "\n") if !occursin("PASSALLT", l)], "\n")
            # PASSALL 100 == default: this stand's max per-plot EXCESS ≤ 5 ⇒ min(EXCESS,PASMAX)=EXCESS for both.
            @test body(sum_def) == body(sum_p100)
            # PASSALL 1 caps the excess ⇒ the projection DIFFERS and 2002 TPA drops materially (jl 572→496).
            @test body(sum_def) != body(sum_p1)
            @test tpa2002(sum_p1) < tpa2002(sum_def)
            @test tpa2002(sum_def) - tpa2002(sum_p1) > 30      # the cap is a large, model-affecting reduction
        end
    end
end
