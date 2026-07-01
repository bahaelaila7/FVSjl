# =============================================================================
# test_r8_intl_board.jl — R8 International ¼" board feet for specific National Forests (D11) vs live FVSsn.
#
# The R8-CLK volume path reports INTERNATIONAL ¼" board feet (FVS volinit2.f:269-272 `VOL(2)=VOL(10)`)
# for GW/JF (IFORST 8), Ouachita (9), Ozark-St Francis (10), and Francis Marion & Sumter (12) except the
# Andrew Pickens district (IDIST 2); every other R8 forest keeps Scribner. jl previously always used
# Scribner (`_r8_scribner_bf`), so a stand on one of those forests under-reported board feet.
#
# Diagnosis (9-layer live debug-stamp trace): the SM tree (DBH 12.7, eq 831CLKE318) on forest 808 has live
# board 85 bd but jl's Scribner = 69; the 85 = International (`_r9_intl_log`: saw logs 30+30+15+10 = 85),
# selected by the IFORST=8 branch. Fix: `_r8_intlqtr_bf` + the IFORST/IDIST gate in compute_volumes!.
#
# Golden = live FVSsn: s07_forest_808 (IFORST 8) and s22_forest_809 (IFORST 9) Bdft@1990 = 351 (bit-exact
# across all cycles); the homogeneous all_GA/PC/BY (IFORST 1) stay Scribner (174/861/1362), unchanged.
# =============================================================================

using Test
using FVSjl

@testset "R8 International ¼\" board (D11) vs live FVSsn" begin
    sdir = joinpath(@__DIR__, "..", "harness", "scenarios")
    _bd(scn) = begin
        txt = FVSjl.run_keyfile(joinpath(sdir, scn * ".key"); variant = FVSjl.Southern(), output = :sum)
        row = nothing
        for ln in split(txt, '\n')
            t = split(strip(ln)); length(t) >= 12 && t[1] == "1990" && (row = t; break)
        end
        row === nothing ? nothing : parse(Int, row[12])
    end
    # International forests (IFORST 8 / 9) — the fix: 351, was 285 (Scribner) before.
    for scn in ("s07_forest_808", "s22_forest_809")
        if isfile(joinpath(sdir, scn * ".key"))
            @test _bd(scn) == 351
        else
            @test_skip "$scn not available"
        end
    end
    # Scribner forests (IFORST 1) — unchanged, must stay bit-exact.
    for (scn, bd) in (("all_GA", 174), ("all_PC", 861), ("all_BY", 1362))
        if isfile(joinpath(sdir, scn * ".key"))
            @test _bd(scn) == bd
        else
            @test_skip "$scn not available"
        end
    end
end
