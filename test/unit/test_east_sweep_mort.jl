# test_east_sweep_mort.jl — regression fixtures from the eastern FIA full-population re-sweep (2026-10).
# Each fixture = a real FIA stand (test/fixtures/east_sweep/<var>_<CN>.db, extracted from the FIA master) + the
# golden .sum captured from the LIVE stock oracle (FVS{sn,ne,cs,ls}_new, "NONE" regime, NUMCYCLE 5). Every
# .sum cycle row must match bit-for-bit (TPA..QMD + volumes + mortality + accretion columns).
#
# LS 18406903010661: morts.f IPASS QMD-convergence loop. Cycle 2 enters with the density rate OFF (RN=0, TN10=T);
# the first (background) pass still computes D10N, finds |D10-D10N|>0.1 and re-enters with D10=D10N, where the
# density rate turns ON; the loop then limit-cycles bg/density for all 10 passes and keeps the pass-10 state. jl
# used to `break` after any background pass (one pass), keeping the background kill (-13 TPA at 1998, −10% of
# every small tree's survival).
# SN 258567874010854: sn/crown.f LSTART dub reads RELDEN (CCF) from CRATET's DENSE that still holds the 10 cycle-0 DEAD
# records (live 94.8 vs live-FVS 554.2 ⇒ SCALE 0.84 vs 0.30); jl dubbed every seedling crown 7-10 points high.
# SN 1893364809290487: the same IPASS loop as the LS stand above (1011 TPA at 2030).
using Test
using FVSjl

const ES_DIR = joinpath(@__DIR__, "..", "fixtures", "east_sweep")
const ES_VAR = Dict("sn"=>"SN", "ne"=>"NE", "cs"=>"CS", "ls"=>"LS")

_es_rows(text) = [split(strip(ln)) for ln in split(text, '\n')
                  if (f = split(strip(ln)); length(f) >= 20 && (y = tryparse(Int, f[1])) !== nothing && y >= 1000)]

@testset "eastern FIA sweep regression fixtures vs live .sum" begin
    for sumf in sort(filter(f -> endswith(f, ".live.sum"), readdir(ES_DIR)))
        base = replace(sumf, ".live.sum" => "")
        vcode, cn = split(base, '_'; limit = 2)
        live = _es_rows(read(joinpath(ES_DIR, sumf), String))
        ncyc = length(live) - 1
        dir = mktempdir()
        key = joinpath(dir, "s.key")
        write(key, """
STDIDENT
$cn
DATABASE
DSNin
$(abspath(joinpath(ES_DIR, base * ".db")))
StandSQL
SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
TreeSQL
SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
END
NUMCYCLE         $(ncyc).0
ECHOSUM
PROCESS
STOP
""")
        jl = _es_rows(FVSjl.run_keyfile(key; variant = FVSjl.variant_from_code(ES_VAR[vcode])))
        @testset "$base" begin
            @test length(jl) == length(live)
            for (a, b) in zip(live, jl)
                @test a[1:21] == b[1:21]
            end
        end
    end
end
