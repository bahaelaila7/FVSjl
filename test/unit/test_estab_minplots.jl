# MINPLOTS (esin.f opt 20 → MINREP) ESTAB-packet keyword (2026-08-23). MINREP is the minimum plot-replication
# target; it drives the estab.f:199-207 IDUP loop (IDUP = smallest I with NPTIDS·I ≥ MINREP = ceil(MINREP/NPTIDS))
# so DUPNPT = NPTIDS·IDUP sets the number of establishment plots looped — hence the ESRANN draw-stream length
# (the one-time WK6 site-prep fill DO 183 I=1,IDUP*NPTIDS + the per-plot bodies + the last plot's ESAVE that
# becomes the continuing stream). WIRED bit-exact once AUTOES #143 aligned the draw stream.
#
# Bit-exact evidence (live instrumented oracle FVSie_estab2, under-stocked IE fixture 753189105290487, ihab 9,
# NPTIDS 50, per-plot body 16+3·23+2·21 = 127): the per-tally ESDRAW seed chain shifts with MINPLOTS exactly as
# the DUPNPT-scaled ESAVE chain predicts —
#   MINPLOTS default (50→idup 1, dupnpt 50):  43303 / 61997 / 49053
#   MINPLOTS 100     (idup 2, dupnpt 100):     43303 / 31334 / 6152
#   MINPLOTS 150     (idup 3, dupnpt 150):     43303 / 53059 / 15474
# jl reproduces all three byte-for-byte (scratchpad/estab/task1). The critical fix was passing wk6 = DUPNPT (not
# the default 50) so the one-time WK6 site-prep prefix scales with IDUP; with the stale wk6=50 the idup>1 chain
# desynced (jl 4681/17715 vs oracle 31334/6152).

using Test
using FVSjl
using FVSjl: KeywordReader, read_keyword!, kw_estab!, StandState, InlandEmpire, init_blockdata!,
             ie_autoes_plot_seeds, ie_esrann!, IEEstabRNG

kwline(kw, fields...) = rpad(kw, 10) * join(rpad.(string.(fields), 10)) * "\n"

@testset "MINPLOTS ESTAB-packet keyword (esin.f opt 20 → MINREP)" begin
    @testset "parse + clamp (default 50; MINREP<20 → 20)" begin
        # default (no MINPLOTS): esinit.f MINREP=50
        s = StandState(InlandEmpire()); init_blockdata!(s, s.variant)
        @test s.estab.minrep == Int32(50)

        # MINPLOTS 100 sets MINREP=100
        s = StandState(InlandEmpire()); init_blockdata!(s, s.variant)
        kr = KeywordReader(IOBuffer("ESTAB\n" * kwline("MINPLOTS", 100) * "END\n"))
        kw_estab!(s, read_keyword!(kr), kr)
        @test s.estab.minrep == Int32(100)

        # MINPLOTS 5 clamps to the 20-plot floor (esin.f:588)
        s = StandState(InlandEmpire()); init_blockdata!(s, s.variant)
        kr = KeywordReader(IOBuffer("ESTAB\n" * kwline("MINPLOTS", 5) * "END\n"))
        kw_estab!(s, read_keyword!(kr), kr)
        @test s.estab.minrep == Int32(20)
    end

    @testset "DUPNPT-scaled ESAVE chain is bit-exact vs the live oracle" begin
        # under-stocked IE fixture: NPTIDS=50, per-plot body=127. seed0 (tally-1) = 43303 for every MINPLOTS.
        body = 127
        # (MINREP, idup, dupnpt, oracle es_stream endpoint, oracle NEXT-tally seed0)
        cases = ((50, 1, 50, 79215, 61997),
                 (100, 2, 100, 40037, 31334),
                 (150, 3, 150, 67795, 53059))
        for (minrep, idup, dupnpt, es_expect, seed_next_expect) in cases
            @test dupnpt == 50 * idup                       # DUPNPT = NPTIDS·ceil(MINREP/NPTIDS)
            # es_stream after the tally = the last plot's ESAVE; the WK6 site-prep prefix = DUPNPT (not 50).
            es = ie_autoes_plot_seeds(43303, dupnpt + 1; wk6 = dupnpt, body = body)[end]
            @test es == es_expect
            # the NEXT disturbance/ingrowth tally draws seed0 = round(ESRANN(es_stream)·100000+0.5)
            dr = ie_esrann!(IEEstabRNG(Float64(es)))
            seed_next = floor(Int, dr * 100000f0 + 0.5f0)
            @test seed_next == seed_next_expect
        end
    end
end
