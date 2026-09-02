# STOCKADJ field-position fix (esin.f opt 13 → activity 440) — validated against the live FVSie echo.
# The oracle (FVSie_clean --keywordfile) prints the parsed fields, which is the ground truth for the layout:
#   STOCKADJ 2000. 0.5  →  "STOCKADJ   DATE/CYCLE= 2000; MULTIPLIER=  0.50"
#   STOCKADJ 0.5        →  "STOCKADJ   DATE/CYCLE=    0; MULTIPLIER=  0.00"   (0.5 read as the DATE; field-2 blank ⇒ 0)
# esin.f opt 13: IDT = IFIX(ARRAY(1)) (field 1 = date/cycle); OPNEW(KODE,IDT,440,1,ARRAY(2)) schedules the
# field-2 multiplier at activity 440 for that date. jl PREVIOUSLY read the multiplier from field 1 and dropped
# the date — the confirmed field-position bug this test locks down (the earlier fixture masked it by putting the
# value in field 1, where jl happened to read it but the oracle read it as the date and cancelled regen).
using Test
using FVSjl
using FVSjl: KeywordReader, read_keyword!, kw_estab!, StandState, InlandEmpire, init_blockdata!

kwline(kw, fields...) = rpad(kw, 10) * join(rpad.(string.(fields), 10)) * "\n"
sched440(s) = filter(x -> x.icflag == Int32(440), s.control.schedule)

@testset "STOCKADJ field position (esin.f opt 13 → activity 440) vs live FVSie echo" begin
    @testset "dated card: STOCKADJ 2000. 0.5 → date 2000, mult 0.5" begin
        s = StandState(InlandEmpire()); init_blockdata!(s, s.variant)
        kr = KeywordReader(IOBuffer("ESTAB\n" * kwline("STOCKADJ", "2000.", "0.5") * "END\n"))
        kw_estab!(s, read_keyword!(kr), kr)
        a = only(sched440(s))
        @test a.year == Int32(2000)          # field 1 = date/cycle (was: read as the multiplier)
        @test a.params[1] == 0.5f0           # field 2 = multiplier
    end

    @testset "single-field card: STOCKADJ 0.5 → date IFIX(0.5)=0 (→1), mult 0.0 (blank field 2)" begin
        s = StandState(InlandEmpire()); init_blockdata!(s, s.variant)
        kr = KeywordReader(IOBuffer("ESTAB\n" * kwline("STOCKADJ", "0.5") * "END\n"))
        kw_estab!(s, read_keyword!(kr), kr)
        a = only(sched440(s))
        @test a.year == Int32(1)             # IFIX(0.5)=0, clamped to the ScheduledActivity year-≥1 floor
        @test a.params[1] == 0.0f0           # blank field 2 ⇒ 0.0 (oracle MULTIPLIER=0.00) — NOT the old 0.5
    end

    @testset "no STOCKADJ ⇒ no activity-440 scheduled (inert)" begin
        s = StandState(InlandEmpire()); init_blockdata!(s, s.variant)
        kr = KeywordReader(IOBuffer("ESTAB\n" * kwline("MINPLOTS", 100) * "END\n"))
        kw_estab!(s, read_keyword!(kr), kr)
        @test isempty(sched440(s))
    end
end
