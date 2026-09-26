# test_em_regent_dk.jl — EM REGENT (em/regent.f) on every sub-model, vs FVSem_g16.
#
# Fixture: the shipped PN tree list plus seedlings/saplings of every EM small-tree sub-model — EMVAR (AF/WL/PP/LP/WB),
# NIVAR (LL, including 3-10" trees in the height-blend range), TTVAR (LM), CRVAR (CW/GA/OH/BA), UTVAR (AS/PB aspen and
# RM juniper) — run with NOTRIPLE and DGSTDEV 0. Golden = the FVSem_g16 .sum rows (live, 2026-09-26). Cycle 1 is
# bit-exact on every printed field: the Fortran-shaped REGENT (one subcycle loop + one DO-30 assembly), em/bratio.f's
# TEMD floor, the cycle-1 mortality WK1 dub and the LL crown OBA/RDM1/OLDPCT all have to be right for that row.
using Test
using FVSjl

@testset "EM REGENT all sub-models (swap fixture) vs FVSem_g16" begin
    fx = joinpath(@__DIR__, "..", "fixtures", "easternmontana", "regent_swap")
    rows = cd(fx) do
        txt = FVSjl.run_keyfile("regent_swap.key"; variant = FVSjl.EasternMontana(), output = :sum)
        [split(strip(l)) for l in split(txt, '\n') if occursin(r"^\d{4}\s", strip(l))]
    end
    gold = ["1990  60  2856 127  382 201  63  2.9  2995  1121     0  5802     0     0     0     0     0 127  382 201  63  2.9      10   78    24    18.7 266 22",
            "2000  70  2572 156  442 239  68  3.3  3540  1553     0  7754     0     0     0     0     0 156  442 239  68  3.3      10   62    25    22.2 266 22"]
    @test length(rows) >= 2
    for (k, g) in enumerate(gold)
        @test rows[k] == split(g)
    end
end
