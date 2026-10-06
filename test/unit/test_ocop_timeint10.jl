# test_ocop_timeint10.jl — OC/OP at TIMEINT 10 (FINT = 10 != the 5-yr YR): every .sum row + every DBS table vs the LIVE oracle.
# Goldens test/fixtures/tiered_ti10/{oc,op}/<stand>_<regime>_ti10.* are the live runs of the private crash-fixed relinks
#   /workspace/.ocwork/tiered/FVSoc_tiered  and  /workspace/.opwork/tiered/FVSop_tiered   (ORACLE_SOURCE_AUDIT §8)
# produced by `test/harness/ocop_ti10/live_ti10.sh V STAND REGIME` (the committed tiered fixture key with `NUMCYCLE 4.0 /
# TIMEINT 0.0 10.0 / TREELIST 0 / TREELIDB`, `FVS<v>_tiered --keywordfile=s.key`) and turned into goldens by
# `julia --project test/harness/ocop_ti10/build_ti10_fixture.jl <OUT> V:STAND:REGIME ...`.
# 645182389126144 / 645142535126144 carry FIA DF dwarf-mistletoe ratings.
module OcOpTimeint10Test
using Test
using FVSjl
include(joinpath(@__DIR__, "..", "harness", "tiered", "tiered_runner.jl"))

FX(v) = joinpath(@__DIR__, "..", "fixtures", "tiered_ti10", lowercase(v))   # stands.db + <cn>_<regime>_ti10.{key,live.sum,tables,<Table>.csv}
const CASES = [
    ("OC", "24061836010900", "none"),
    ("OC", "41133570010497", "none"),
    ("OC", "645142535126144", "none"),
    ("OC", "24061836010900", "thinbba"),
    ("OC", "645182389126144", "thinbba")
]

@testset "OC/OP TIMEINT 10 vs live: $v $cn $r" for (v, cn, r) in CASES
    d = mktempdir()
    txt, db, crashed, err = run_case(v, cn, r * "_ti10"; dir = d, fx = FX(v))
    @test !crashed
    ms = compare_case(v, cn, r * "_ti10", txt, db; fx = FX(v))
    isempty(ms) || foreach(m -> println("  ", key_str(m), " gold=", m.gold, " got=", m.got), first(ms, 8))
    @test isempty(ms)
end
end
