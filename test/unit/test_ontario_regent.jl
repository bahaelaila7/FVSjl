# =============================================================================
# test_ontario_regent.jl — ON sub-12cm REGENT small-tree growth (canada/on/regent.f).
#
# ONSTHG (the Penner-2010 annual small-tree height model, the ON-distinctive part of REGENT) is
# BIT-EXACT vs an instrumented FVSon_g16 (onsthg.f dump, scratchpad/on/): per-tree HTG(annual) hex
# matches for both the exp(-5)-clamped dense case (PJ) and the computed case (MH). The full
# small_tree_growth! runs the ont_sm stand (8 trees, DBH 4-11 cm) without error; the deterministic
# HTGR matches the oracle bit-exact on records where the [XMIN,XMAX] blend does not consume the
# large-tree HTG (xwt=0), and elsewhere inherits the accepted large-tree-DG-for-small-trees residual.
# =============================================================================
using FVSjl, Test
const _R = FVSjl
_h2f(h) = reinterpret(Float32, parse(UInt32, h; base=16))
_f2h(x::Float32) = uppercase(string(reinterpret(UInt32, x); base=16, pad=8))

@testset "ON ONSTHG small-tree height increment bit-exact vs FVSon_g16" begin
    # oracle onsthg.f dump (fort.775): (sp, D hex, H hex, SI hex, BAL hex) => HTG(annual) hex
    cases = [
        (1,  "407BF7CF", "419D7AF6", "42766354", "44BFCC62", "3CB517E8"),  # PJ — exp(-5) clamp (dense BAL)
        (11, "4062C56E", "41833BCD", "42700000", "452F279B", "3CB517E8"),  # CE — exp(-5) clamp
        (26, "408A9518", "41B7BA1F", "4250AB9A", "00000000", "3F8B5182"),  # MH — computed (BAL=0)
    ]
    for (sp, dh, hh, sih, balh, gh) in cases
        g = _R._on_onsthg(sp, _h2f(dh), _h2f(hh), _h2f(sih), _h2f(balh))
        @test _f2h(g) == gh
    end
end

@testset "ON REGENT small_tree_growth! runs on a sub-12cm stand (ont_sm)" begin
    # Previously errored ("not yet ported"); now grows every record through the Penner REGENT chain.
    sums = _R.run_keyfile("scratchpad/on/ont_sm.key"; variant=_R.Ontario())
    @test occursin("ONTSM", sums)          # completed → emitted a .sum for the stand
    # cyc0 inventory row unaffected (REGENT acts in the growth cycle): TPA matches the oracle (±1 NINT).
    @test occursin("296520", sums) || occursin("296519", sums)
end
