# test_svs_snag.jl — SVS mortality→snag data path, sub-chunk 1: standing-dead snag AGING + DISPLAY.
#
# Validates the SVS snag aging (SVSNAGE) + emission (SVOUT IOBJTP=2 branch) bit-exact vs the LIVE
# FVSkt snag record, MEASURED from an instrumented svout.f dump (WRITE of the SVDEAD arrays right
# after the SVSNAGE call) on stand S248112, kt2c (2-cycle) at the year-2010 End-of-projection
# picture. That picture's one snag is the tripled-upper WL that died in cycle 2:
#
#   FROZEN AT DEATH (SVSNAD): sp=WL ODIA=10.0664234 OLEN=86.4655304 SNGDIA=10.0664234
#                             SNGLEN=86.4655304 CRNDIA=13.0649385 CRNRTO=27.0 IYRCOD=2009 ISTATUS=2
#   AGED TO 2010 (SVSNAGE):   SNDI=10.0664234 SNHT=84.6912613 SNCRDI=11.7584448 SNCRTO=26.9999981
#   DISPLAYED (SVOUT):        WL tree#=1 class=98 dbh=10.1 ht=85. crad=5.9 cr=0.27 at (201.15,181.52)
#
# The oracle snag line is kt2c_003.svs line 10. This sub-chunk ports SVSSnag (SVDEAD), the KT
# snag params (fmvinit: HTR1=0.0228, HTX(WL)=0.9), FMSNGHT/FMSNGDK (KT CASE DEFAULT), svsnage_standing!
# and svs_write_snag!. The mortality object-SELECTION (svrmov) + snag CREATION (svsnad) + engine
# SVMORT seam are the NEXT sub-chunk — so this validates against the measured frozen record, not a
# full run. (The displayed HT/CR are the same pre-existing tripled-record growth residual as _002.)

using Test, FVSjl

@testset "SVS snag — standing-dead aging + display bit-exact vs live FVSkt (measured)" begin
    # A minimal KT stand just for its species table (code2[2] == "WL").
    s = FVSjl.each_stand(joinpath(@__DIR__, "..", "fixtures", "svs", "kt0.key");
                         variant = FVSjl.Kootenai())[1]
    FVSjl.notre!(s); FVSjl.setup_growth!(s)
    @test rstrip(String(s.species.code2[2])) == "WL"

    # Reconstruct the frozen snag record measured from the live FVSkt svout dump.
    sn = FVSjl.SVSSnag(2,                 # sp = WL
                       10.0664234f0,      # odia
                       86.4655304f0,      # olen
                       10.0664234f0,      # sngdia
                       86.4655304f0,      # snglen
                       13.0649385f0,      # crndia
                       27.0f0,            # crnrto
                       2009,              # iyrcod
                       2,                 # istatus (red hard snag)
                       -1.0f0,            # falldir (standing)
                       0)                 # oidtre (unused here)

    # SVSNAGE aging to 2010 (ILYEAR=2000, XMOD=1) — bit-exact vs the measured live SVSNAGE outputs.
    snw = deepcopy(sn)
    (sndi, snht, sncrdi, sncrto) = FVSjl.svsnage_standing!(snw, 2010, 2000, 1.0f0)
    @test sndi   ≈ 10.0664234f0 atol=1f-4
    @test snht   ≈ 84.6912613f0 atol=1f-3      # FMSNGHT: 86.4655·(1−0.0228·0.9)
    @test sncrdi ≈ 11.7584448f0 atol=1f-3      # CRNDIA·0.90
    @test sncrto ≈ 26.9999981f0 atol=1f-3      # (OLEN/SNGLEN)·(CRNRTO·.01−1)+1, ×100
    @test snw.istatus == 2                     # ITIDIF=1 (<2) ⇒ stays red hard snag (class 98)
    @test FVSjl.svs_snag_class(snw.istatus) == 98

    # SVOUT emission — the whole snag object record must be BYTE-IDENTICAL to the oracle line.
    io = IOBuffer()
    ok = FVSjl.svs_write_snag!(io, sn, s, 1, 201.146072f0, 181.524506f0; iyear = 2010, ilyear = 2000)
    @test ok
    line = rstrip(String(take!(io)), '\n')
    oracle = "WL                 1 98 0 1  10.1   85. 0   0 0   5.9 0.27   5.9 0.27   5.9 0.27   5.9 0.27 1 0  201.15  181.52 0"
    @test line == oracle
end
