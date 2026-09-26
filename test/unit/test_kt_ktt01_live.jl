# KT ktt01 (S248112) vs the live FVSkt_clean oracle: every .sum data row 1990-2050 must equal live.
# Locks the KT chain measured bit-exact on 2026-09-26: kt/cratet.f √D height dub + SMHTRG, htgf.f tripled-copy HTG,
# REGENT subcycle density feedback / IND1 walk / DNEW≥3 DUBSCR re-dub, crown.f OLDPCT/OBA/RDM1, cftopk.f KTFCTR.
using Test, FVSjl

@testset "KT ktt01 .sum == live FVSkt" begin
    dir = joinpath(@__DIR__, "..", "fixtures", "kootenai", "ktt01live")
    gold = [l for l in readlines(joinpath(dir, "ktt01.live.sum")) if !startswith(l, "-999")]
    out = mktempdir()
    cp(joinpath(dir, "ktt01.key"), joinpath(out, "ktt01.key")); cp(joinpath(dir, "ktt01.tre"), joinpath(out, "ktt01.tre"))
    jl = cd(out) do
        [l for l in split(FVSjl.run_keyfile("ktt01.key"; variant = FVSjl.Kootenai(), output = :sum), '\n')
         if !isempty(strip(l)) && !startswith(l, "-999")]
    end
    @test length(jl) == length(gold)
    for (g, j) in zip(gold, jl)
        @test rstrip(j) == rstrip(g)
    end
end
