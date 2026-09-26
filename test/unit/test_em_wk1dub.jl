# test_em_wk1dub.jl — EM cycle-1 mortality WK1 for an UNMEASURED tree of a CALIBRATED added species, vs FVSem_g16.
#
# em/dgdriv.f DO 220 dubs WK1 = DG for records without a measured increment from the SECOND calibration DGF call
# (:770 CALL DGF(WK3)), which runs after the correction terms are final — em/dgf.f adds COR(ISPC) into CONSPP/DDS. jl
# captured the first call's WK2 (COR still 0), so a calibrated added species' unmeasured tree got the wrong Hamilton
# vigor G. Fixture: the regcal stand with 5 LM (sp4) carrying a measured DG at 5-9" (LM COR scale 2.80) plus one
# unmeasured 1.2" LM. Live MORTS DEBUG: that tree's G 0.1547, 6.462 TPA killed (first-call capture gave 0.0857 /
# 13.12). Golden = the FVSem_g16 .sum rows (2026-09-26).
using Test
using FVSjl

@testset "EM DO-220 WK1 dub from the post-COR DGF (calibrated added species) vs FVSem_g16" begin
    fx = joinpath(@__DIR__, "..", "fixtures", "easternmontana", "wk1dub")
    rows = cd(fx) do
        txt = FVSjl.run_keyfile("wk1dub.key"; variant = FVSjl.EasternMontana(), output = :sum)
        [split(strip(l)) for l in split(txt, '\n') if occursin(r"^\d{4}\s", strip(l))]
    end
    gold = ["1990  60  6250 180  590 338  63  2.3  5921  1511     0  7672     0     0     0     0     0 180  590 338  63  2.3      10   98    82    25.2 366 22",
            "2000  70  5231 204  631 380  67  2.7  6082  1644     0  8013     0     0     0     0     0 204  631 380  67  2.7      10  118    93    23.5 366 21",
            "2010  80  4320 220  644 386  55  3.1  6331  1665     0  7896     0     0     0     0     0 220  644 386  55  3.1       0    0     0    20.8 366 21"]
    @test length(rows) >= 3
    for (k, g) in enumerate(gold)
        @test rows[k] == split(g)
    end
end
