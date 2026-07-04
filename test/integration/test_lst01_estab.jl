# =============================================================================
# test_lst01_estab.jl — Lake States (LS) BARE-GROUND establishment (PLANT) stand
#
# The lst01 stand-5 scenario: bare ground, plant 400 jack pine (sp1) + 400 sp3
# in 1992, grow 10 cycles. Validates the LS establishment/regent seedling growth
# BIT-EXACT vs live FVSls. The fix that made this bit-exact: the planted-seedling
# default-height RAN acceptance window (estab.f:483/490) is [-2.5,2.5] for NE/CS/LS
# (only SN uses [0,1.5]); jl had wrongly given CS+LS the SN window, which rejects the
# low RAN tail and biased the seedling heights HIGH ("BARE-PLANT over-sizing"), which
# in turn shifted the merch report AND the dense-stand self-thinning.
#
# Live FVSls (estab.key): 2022 TPA 708 / BA 114 / SDI 266 / MCuFt 1325; 2042 TPA 503 /
# BA 192 / SDI 373 / MCuFt 4491; 2092 TPA 169 / BA 209.
# =============================================================================

using Test
using FVSjl

@testset "LS BARE-plant establishment (vs live FVSls)" begin
    key = joinpath(@__DIR__, "..", "fixtures", "ls", "estab.key")
    if !isfile(key)
        @info "estab.key fixture not present; skipping LS establishment test"
    else
        txt = FVSjl.run_keyfile(key; variant = LakeStates(), output = :sum)
        rows = Dict{Int,Vector{Int}}()
        for l in split(txt, "\n")
            m = match(r"^(20\d\d)\s+\d+\s+(\d+)\s+(\d+)\s+(\d+)\s+\d+\s+\d+\s+[\d.]+\s+\d+\s+(\d+)", l)
            m === nothing && continue
            rows[parse(Int, m[1])] = [parse(Int, m[2]), parse(Int, m[3]), parse(Int, m[4]), parse(Int, m[5])]
        end
        # planting: 800 TPA in 2002 (bit-exact)
        @test rows[2002][1] == 800
        # BIT-EXACT vs live across the whole trajectory (TPA / BA / SDI / MCuFt): the RAN-window fix.
        @test rows[2022] == [708, 114, 266, 1325]
        @test rows[2032][1] == 597 && rows[2032][4] == 2990
        @test rows[2042] == [503, 192, 373, 4491]   # dense self-thinning + merch, all bit-exact
        @test rows[2092][1] == 169 && rows[2092][2] == 209
    end
end
