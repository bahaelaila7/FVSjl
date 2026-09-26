# test_cr_cyc1_heights.jl — CR cycle-1 per-tree height growth vs FVScr_clean (2026-09-26), on the REGCAL fixture
# (test/fixtures/centralrockies/regcal, IMODTY 5 lodgepole model, DGSTDEV 0, NOTRIPLE).
#
# Two faithful fixes are pinned here:
#  * cratet.f:535-552 FINDAG runs AFTER :522 CROWN and reads COMMON BA/RELDEN as the :175 backdating DENSE left
#    them (live + cycle-0 dead) — fndag.f IMODTY 5 uses CCFTEM=RELDEN−125. A live-only recompute aged the DF
#    D=10.4 H=55 at 55 instead of live's 60 ⇒ HTGF AP 55 ⇒ HTG 7.764 vs live 7.463.
#  * htgf.f:246-256: the IMODTY 5 ABIRTH≤31 H30 replacement also resets HHE2=H30, HHE1=HT, which the AGERNG>40
#    uneven-aged blend reads (HGE); jl only replaced HTG ⇒ the five 1.6–2.1" saplings in the XWT window grew
#    0.02–0.18 ft short (AF D=1.2: large-tree HTG 3.915 vs live 4.1776 ⇒ blended 2.544 vs 2.666).
# Golden: live FVS_TreeList (TREELIST 0 + TREELIDB) and .sum. Volumes are NOT pinned (FW2 cubic precision residual).
using Test
using FVSjl

const _CRC1_KEY = joinpath(@__DIR__, "..", "fixtures", "centralrockies", "regcal", "regcal.key")

@testset "CR cycle-1 heights (FINDAG ABIRTH on CRATET-DENSE RELDEN; H30 resets HHE1/HHE2)" begin
    s, _ = FVSjl.initialize(_CRC1_KEY; variant = FVSjl.CentralRockies())
    FVSjl.notre!(s); FVSjl.setup_growth!(s)
    t = s.trees
    # live FINDAG EFFAGE (DEBUG CRATET "IN FINDAG I,EFFAGE,…"): the DF D=10.4 H=55 → 60, DF D=9.4 H=60 → 70,
    # LP D=11.5 H=73 → 130, WL(10) D=7.9 H=75 → 200 (AGEMAX), ES D=5.3 H=27 → 35.
    for (sp, d, h, age) in ((3, 10.4f0, 55f0, 60f0), (3, 9.4f0, 60f0, 70f0), (11, 11.5f0, 73f0, 130f0),
                            (10, 7.9f0, 75f0, 200f0), (5, 5.3f0, 27f0, 35f0))
        i = findfirst(k -> Int(t.species[k]) == sp && t.dbh[k] == d && t.height[k] == h, 1:t.n)
        @test i !== nothing
        i === nothing || @test t.birth_age[i] == age
    end

    # .sum growth columns (TPA BA SDI CCF TopHt QMD) — every row equal to live; the volume columns carry the FW2
    # total-cubic residual and are excluded.
    live = [(1990, 7409, 138, 400, 190, 63, 1.8), (2000, 4392, 169, 449, 223, 69, 2.7),
            (2010, 2821, 186, 454, 233, 75, 3.5), (2020, 1989, 202, 461, 244, 79, 4.3)]
    out = FVSjl.run_keyfile(_CRC1_KEY; variant = FVSjl.CentralRockies(), output = :sum)
    rows = Dict{Int,Vector{String}}()
    for l in split(out, '\n')
        f = split(l)
        (length(f) > 8 && !startswith(l, "-999")) || continue
        y = tryparse(Int, f[1]); y === nothing && continue
        rows[y] = f
    end
    for (y, tpa, ba, sdi, ccf, toph, qmd) in live
        f = rows[y]
        @test parse(Int, f[3]) == tpa
        @test parse(Int, f[4]) == ba
        @test parse(Int, f[5]) == sdi
        @test parse(Int, f[6]) == ccf
        @test parse(Int, f[7]) == toph
        @test parse(Float64, f[8]) == qmd
    end
end
