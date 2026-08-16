# Western Root Disease (WRD) — Chunk −1 (infrastructure + keyword reader + inert seam).
#
# Validates the RD keyword reader / default-init against the LIVE FVSkt oracle's RDIN
# option echo, and that the engine seam is INERT (a stand with an RDIN block produces
# byte-identical summary data rows to the same stand without one — the Chunk-0 mortality
# body is not yet ported). Oracle reference: FVSkt_clean on the same S248112 stand with
#   RDin / RRType 3 / RRInit 0 10 10 20 0.1 10 3 / SArea 100 / BBClear / End
# echoes (rd.out): disease=ARMILLARIA, #centers=10, infected TPA=10, uninfected TPA=20,
# proportion roots infected=0.10, disease area=10 ac, stand area=100 ac.

using Test
using FVSjl

const _KT_TRE = """
   1      248112       0101   011LP 11510   0734   00111     0  0
   2      248112       0101   031DF 001     0026   00222     0  0
   3      248112       0102   011WH 06523   0308   00111     0  0
   4      248112       0102   011WL 07906   0753   00111     0  0
   5      248112       0102   018WL 346            10322     0  0
   6      248112       0103   011WL 08007   0633   96222     0 56
   7      248112       0103   011GF 06220   0385   00111     0  0
   8      248112       0103   011WL 084       54   00111     0  0
   9      248112       0103   011LP 09511   0603   00111     0  0
  10      248112       0104   011DF 040     0203   00111    50  0
  11      248112       0104   011WL 08212   0655   50111     0  0
  12      248112       0105   011DF 012     0116   00222    42  0
  13      248112       0105   011DF 019     0135   00222    47  0
  14      248112       0105   016LP 072            11322     0  0
  15      248112       0105   031GF 001     0037   00222     0  0
  16      248112       0105   011GF 05309   0277   00111     0  0
  17      248112       0106   011DF 10010   0654   00111     0  0
  18      248112       0106   011GF 06112   0388   00111     0  0
  19      248112       0106   011DF 12716   0674   00111     0  0
  20      248112       0107                          800
  21      248112       0108   011LP 09605   0603   00222     0  0
  22      248112       0108   011DF 10409   0555   97222     0 49
  23      248112       0108   011LP 085       03   00111     0  0
  24      248112       0109   011GF 10910   0657   00111     0  0
  25      248112       0109   011DF 09418   0604   00111     0  0
  26      248112       0110   011ES 03206   0175   00222    32  0
  27      248112       0110   011ES 001     0027   00222     0  0
  28      248112       0110   011ES 05810   0287   00111     0  0
  29      248112       0110   011ES 05010   0253   00111    37  0
  30      248112       0111   011GF 06614   0307   00111     0  0
"""

_kt_head(title) = """
SCREEN
NOAUTOES
NOTRIPLE
STATS
STDIDENT
S248112  $title
DESIGN                                        11.0       1.0
STDINFO     11406001     570.0      60.0     315.0      30.0      34.0
INVYEAR       1990.0
NUMCYCLE        10.0
TREEFMT
(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,
T52,I2,T66,5I1,T54,7I1,T75,F3.0)
TREEDATA
"""

const _RDIN_BLOCK = """
RDIN
RRTYPE             3
RRINIT             0        10        10        20       0.1        10         3
SAREA            100
BBCLEAR
END
"""

_datarows(sumtext) = filter(l -> !startswith(l, "-999"), split(strip(sumtext), '\n'))

@testset "Western Root Disease (WRD) — Chunk −1" begin

    @testset "rdinit.f defaults + rdatv gate" begin
        rd = FVSjl.RootDiseaseState()
        @test rd.iroot == 0
        @test rd.rrman == false && rd.rrtinv == false
        @test rd.minrr == 1 && rd.maxrr == 2      # rdinit annosus default (overwritten by RRTYPE)
        @test rd.sarea == 100.0f0
        @test isapprox(rd.dimen, 2087.0f0; atol = 0.5f0)
        @test all(==(Int32(20)), rd.ncents)       # rdinit NCENTS=20
        @test all(==(0.1f0), rd.rrincs)           # rdinit RRINCS=0.1
        @test rd.xxinf == (0.0f0, 3.9f0, 35.4f0, 0.0f0, 0.0f0)
        @test rd.yyinf == (0.0f0, 5.0f0, 40.0f0, 0.0f0, 0.0f0)
        @test rd.nninf == 3
        @test rd.xminlf == Float32[1.0, 1.0, 1.0, 20.0]
        @test FVSjl.rd_active(rd) == false        # rdatv: L = RRTINV .OR. RRMAN
        @test FVSjl.rd_active(nothing) == false
    end

    # end-to-end through the real keyword dispatch + engine
    dir = mktempdir()
    write(joinpath(dir, "shared.tre"), _KT_TRE)
    ctrl_key = joinpath(dir, "ctrl.key")
    rd_key   = joinpath(dir, "rd.key")
    write(ctrl_key, _kt_head("RD CONTROL") * "ECHOSUM\nPROCESS\nSTOP\n")
    write(rd_key,   _kt_head("RD ACTIVE ") * _RDIN_BLOCK * "ECHOSUM\nPROCESS\nSTOP\n")
    # run_keyfile locates <base>.tre alongside the keyfile
    cp(joinpath(dir, "shared.tre"), joinpath(dir, "ctrl.tre"))
    cp(joinpath(dir, "shared.tre"), joinpath(dir, "rd.tre"))
    v = FVSjl.Kootenai()

    @testset "RDIN reader vs live FVSkt RDIN echo" begin
        rd = nothing
        for s in FVSjl.each_stand(rd_key; variant = v)
            rd = s.root_disease
            break
        end
        @test rd !== nothing                       # RDIN block was dispatched + parsed
        @test FVSjl.rd_active(rd) == true          # rdatv gate open (RRMAN set by RRINIT)
        @test rd.iroot == 1
        @test rd.rrman == true && rd.rrtinv == false
        @test rd.lrtype == true
        @test rd.minrr == 3 && rd.maxrr == 3       # RRTYPE 3
        @test FVSjl.RD_DISEASE_NAMES[Int(rd.maxrr)] == "ARMILLARIA"
        @test rd.ncents[3] == 10                   # NUMBER OF CENTERS = 10
        @test rd.prkill[3] == 10.0f0               # INFECTED TREES/ACRE = 10.00
        @test rd.prun[3]   == 20.0f0               # UNINFECTED TREES/ACRE = 20.00
        @test rd.rrincs[3] == 0.1f0                # PROPORTION ROOTS INFECTED = 0.10
        @test rd.parea[3]  == 10.0f0               # DISEASE AREA = 10.00 ACRES
        @test rd.lparea[3] == true
        @test rd.sarea == 100.0f0                  # STAND AREA = 100.00 ACRES
        @test isapprox(rd.dimen, 2087.0f0; atol = 0.5f0)
        @test rd.bbclear == true && rd.lbbon == false
    end

    @testset "engine seam is inert (Chunk-0 mortality pending)" begin
        out_ctrl = FVSjl.run_keyfile(ctrl_key; variant = v, output = :sum)
        out_rd   = FVSjl.run_keyfile(rd_key;   variant = v, output = :sum)
        # header (-999) carries a wall-clock timestamp; compare the data rows only
        @test _datarows(out_ctrl) == _datarows(out_rd)
    end
end
