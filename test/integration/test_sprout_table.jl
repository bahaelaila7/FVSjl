# test_sprout_table.jl — SPROUT keyword per-species / stump-DBH-range multiplier table (esuckr.f activity 450).
#
# The SPROUT keyword (esin.f opt 26) carries, per species, a sprout-COUNT multiplier (SMULT), a HEIGHT
# multiplier (HMULT) and a stump-DBH window [DLO, DHI). esuckr.f:197-205 looks these up by the PARENT
# stump's species and DBH (DO 450: default 1/1, each matching activity with DSTMP ∈ [DLO,DHI) overwrites).
# FVSjl previously applied a single global SMULT/HMULT to ALL species and ignored the window — so e.g.
# `SPROUT 22 3` tripled every species instead of just 22. Now ported to a per-species override table.
#
# Validated BIT-EXACT vs live FVSsn three ways (all on the sprout.key thin-then-sprout stand):
#   • smult=3 for species 22         → 2005 TPA 491→729  (sprout_smult.sum.save)
#   • smult=3 for species 22, [8,99) → 2005 TPA = 509, strictly between 491 (all 1×) and 729 (all 3×),
#                                       proving only the ≥8" stumps took the 3× override (sprout_win3.sum.save)
#   • the default form (sprout.key, SMULT=HMULT=1) stays bit-exact (covered by test_sprout_regen.jl).

using Test, FVSjl

const _SPT_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_spt_rows(txt) = [split(l) for l in split(txt, "\n")
                  if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                                y !== nothing && 1980 < y < 2110)]
_spt_base(path) = [split(l) for l in eachline(path)
                   if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                                 y !== nothing && 1980 < y < 2110)]

@testset "SPROUT per-species/DBH-range multiplier table" begin
    # 1. PARSE — the SPROUT keyword (read inside the ESTAB block) populates the per-species override table
    #    with the right fields. sprout_win3.key has `SPROUT 2000. 22.0 3.0 1.0 8.0 99.0`.
    s = first(FVSjl.each_stand(joinpath(_SPT_DIR, "sprout_win3.key"); variant = FVSjl.Southern(), faithful = true))
    @test s.control.lsprut
    @test length(s.control.sprout_overrides) == 1
    code, sm, hm, dlo, dhi = s.control.sprout_overrides[1]
    @test code == 22f0 && sm == 3f0 && hm == 1f0 && dlo == 8f0 && dhi == 99f0

    # 2. FORTRAN — per-species multiplier (smult=3 for species 22) bit-exact vs live FVS.
    for (key, save) in (("sprout_smult.key", "sprout_smult.sum.save"),   # smult=3, no window
                        ("sprout_win3.key",  "sprout_win3.sum.save"))    # smult=3, window [8,99)
        base = _spt_base(joinpath(_SPT_DIR, save))
        got  = _spt_rows(FVSjl.run_keyfile(joinpath(_SPT_DIR, key)))
        @test length(got) == length(base) && !isempty(base)
        for (g, b) in zip(got, base)
            @test g[1] == b[1]            # year
            @test g[3] == b[3]            # TPA — bit-exact (the multiplied/windowed sprout count)
            @test g[5] == b[5]            # BA
            @test g[8] == b[8]            # QMD
            @test g[9] == b[9]           # TopHt — BIT-EXACT (measured Δ0 both scenarios; was padded ≤2)
        end
        # TCuFt (col 10): sprout_smult is BIT-EXACT every row; sprout_win3 straddles the print/tree-sum ±1
        # boundary at 2020 (jl 1908 / ft 1907) ⇒ exposed @test_broken (doctrine #9), not a passing ±1.
        if key == "sprout_win3.key"
            @test_broken all(parse(Float32, g[10]) == parse(Float32, b[10]) for (g, b) in zip(got, base))  # tree-sum order
        else
            @test all(parse(Float32, g[10]) == parse(Float32, b[10]) for (g, b) in zip(got, base))  # BIT-EXACT
        end
    end

    # 3. the window must actually be honored: smult=3 windowed (509) lies strictly between the
    #    all-1× baseline (491) and the all-3× result (729) at 2005 — only the ≥8" stumps got 3×.
    row2005(key) = begin
        r = first(filter(x -> x[1] == "2005", _spt_rows(FVSjl.run_keyfile(joinpath(_SPT_DIR, key)))))
        parse(Int, r[3])
    end
    all1 = row2005("sprout.key"); all3 = row2005("sprout_smult.key"); win = row2005("sprout_win3.key")
    @test all1 < win < all3
end
