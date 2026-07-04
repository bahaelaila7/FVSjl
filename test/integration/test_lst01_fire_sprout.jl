# =============================================================================
# test_lst01_fire_sprout.jl — LS post-fire STUMP SPROUTING (vs live FVSls)
#
# The default (no-NOAUTOES) LS establishment behaviour is dominated NOT by the
# STEMS stocking model (that path is LAUTAL-gated and LAUTAL is never set .TRUE.
# in the LS build) but by STUMP SPROUTING of FIRE-KILLED trees:
#   fmkill.f:80 (ICALL=1, GRADD, after the fire) calls
#       ESTUMP(ISP,DBH,FIRKIL,ITRE,ISHAG)
#   for every record with FIRKIL>0.00001 — feeding fire kills into the SAME
#   sprout pool as cutting (cuts.f:1713). ESUCKR then sprouts them, gated by
#   LSPRUT (default .TRUE.; NOAUTOES/NOSPROUT turn it off). ISHAG = IY(ICYC+1) −
#   BURNYR; the SIMFIRE is booked at BURNYR = cycle-start, so ISHAG = cyclen.
#
# This is the mechanism behind the "regen after a fire without NOAUTOES" the
# sweep surfaced (live 456 vs a jl-without-this-fix 177 at 2010). jl already had
# the full LS sprout machinery (esuckr!/nsprec_ls/essprt_ls/sprtht_ls, aspen
# sp41) and lsprut default ON; it just never fed FIRE kills into cut_log (only
# harvest cuts). The fix (fmburn.jl) pushes fire-killed sprouting records into
# cut_log, and esuckr! drains the pool once (empty! at its end, esuckr.f:380
# ITRNRM=0) so a cut-free fire cycle does not re-sprout every subsequent cycle.
#
# Live FVSls (ffe_sprout.key, SIMFIRE 2005, NO NOAUTOES): the 2005 fire kills the
# overstory and the fire-killed hardwoods sprout ⇒ 2010 TPA 456 / BA 89 / SDI 164
# (vs 177/88/158 WITH NOAUTOES); 2020 TPA 446. The fire cycle and the next are
# BIT-EXACT; later cycles track live within the documented tripling/DGSCOR tail.
# =============================================================================

using Test
using FVSjl

@testset "LS post-fire stump sprouting (vs live FVSls)" begin
    key = joinpath(@__DIR__, "..", "fixtures", "ls", "ffe_sprout.key")
    if !isfile(key)
        @info "ffe_sprout.key fixture not present; skipping LS fire-sprout test"
    else
        txt = FVSjl.run_keyfile(key; variant = LakeStates(), output = :sum)
        rows = Dict{Int,NTuple{3,Int}}()
        for l in split(txt, "\n")
            m = match(r"^(20\d\d)\s+\d+\s+(\d+)\s+(\d+)\s+(\d+)", l)
            m === nothing && continue
            rows[parse(Int, m[1])] = (parse(Int, m[2]), parse(Int, m[3]), parse(Int, m[4]))
        end
        # Pre-fire cycles are unaffected (identical with/without NOAUTOES).
        @test rows[2000][1] == 524
        # THE fire cycle: post-fire sprouts make TPA/BA/SDI BIT-EXACT vs live FVSls
        # (without the fire-kill→sprout feed this was 177/88/158 — survivors only).
        @test rows[2010] == (456, 89, 164)
        # Next cycle TPA is bit-exact too (sprout cohort grows in; live 446).
        @test rows[2020][1] == 446
    end
end
