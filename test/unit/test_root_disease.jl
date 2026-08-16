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

    # -------------------------------------------------------------------------
    # Chunk 0a — RD stochastic primitives (rd/rdrani.f + rd/rdrann.f + rd/rdranp.f).
    # Golden values captured from the LIVE relinked FVSkt oracle (a single-.o
    # instrumentation swap dumping the RD RNG stream and every RDRANP call; the
    # instrumented .sum was byte-identical to FVSkt_clean). Full runs were checked
    # bit-exact: 2000 RNG draws (0 mismatch) and 297 RDRANP calls (0 mismatch).
    # -------------------------------------------------------------------------
    @testset "RD RNG (rdrani/rdrann) — bit-exact vs live FVSkt" begin
        rd = FVSjl.RootDiseaseState()
        FVSjl.rd_rani!(rd, 889347.0)              # rd/rdinit.f DSEED default
        @test rd.rd_ss == 889347.0                # INT(889347) already odd
        @test rd.rd_s0 == 889347.0
        # (integer state S1, Float32 return) for the first 8 draws, seed 889347
        golden = [
            (2062353147.0, 0.9603580236434937),
            (1583279049.0, 0.7372717857360840),
            ( 701106566.0, 0.3264781832695007),
            ( 255283673.0, 0.1188757270574570),
            (2027849052.0, 0.9442908167839050),
            (1493539074.0, 0.6954833269119263),
            (2122350582.0, 0.9882965087890625),
            ( 642855004.0, 0.2993526756763458),
        ]
        for (s1, val) in golden
            v_jl = FVSjl.rd_rann!(rd)
            @test rd.rd_s1 == s1                  # exact integer generator state
            @test v_jl === Float32(val)           # exact Float32 return
        end
        # default-constructed state is already seeded from DSEED
        rd2 = FVSjl.RootDiseaseState()
        @test rd2.dseed == 889347.0 && rd2.rd_s0 == 889347.0
    end

    @testset "RD binomial proportion (rdranp) — bit-exact vs live FVSkt" begin
        rd = FVSjl.RootDiseaseState()
        # (prop, rannum, expected) triples straight from the oracle dump.
        triples = [
            (0.1,  0.2924090325832367, 0.07999999821186066),
            (0.1,  0.8925243616104126, 0.15999999642372131),
            (0.1,  0.04066308960318565, 0.03999999910593033),
            (0.6,  0.2924090325832367, 0.6999999880790710),   # LREV (prop>0.5)
            (0.6,  0.6669166088104248, 0.5000000000000000),
            (0.6,  0.8678408265113831, 0.3999999761581421),
            (0.05, 0.2255645245313644, 0.02999999932944775),
            (0.05, 0.06300469487905502, 0.01999999955296516),
            (0.05, 0.9199229478836060, 0.07999999821186066),
        ]
        for (prop, rannum, exp) in triples
            got = FVSjl.rd_ranp_from(rd, Float32(prop), Float32(rannum))
            @test got === Float32(exp)
        end
    end

    # -------------------------------------------------------------------------
    # Chunk 0b-1 — disease-center placement (rd/rdcloc.f + rd/rdarea.f). Golden
    # PCENTS/PAREA captured from the LIVE relinked FVSkt oracle (single-.o swap of
    # rdcloc/rdarea dumping PCENTS after placement and PAREA after the grid-area
    # match; the instrumented .sum was byte-identical to FVSkt_clean). Dominant-
    # signal keyfile: RDin/RRType 2 (S-annosus)/RRInit 0 20 300 300 1.0 100 2/SArea
    # 100, seed 889347. All 20 centers + PAREA reproduced bit-exact.
    # -------------------------------------------------------------------------
    @testset "RD center placement (rdcloc/rdarea) — bit-exact vs live FVSkt" begin
        rd = FVSjl.RootDiseaseState()
        rd.minrr = Int32(2); rd.maxrr = Int32(2)     # S-annosus only
        rd.ncents[2] = Int32(20)
        rd.parea[2]  = 100.0f0
        rd.ipcflg[2] = Int32(0)                      # random centers
        rd.lonect[2] = Int32(0)
        rd.sarea = 100.0f0
        FVSjl.rd_rani!(rd, 889347.0)                 # placement is the first RD-stream consumer
        FVSjl.rd_place_centers!(rd, true)            # LSTART init placement

        # (x, y) of all 20 centers, straight from the oracle PCENTS dump.
        xy = [
            (2004.267211914062, 1538.686157226562),
            ( 681.3599853515625, 248.0936431884766),
            (1122.115234375000, 1258.004516601562),
            (1442.703369140625,  748.8413085937500),
            ( 187.6241149902344, 2028.444458007812),
            (2049.195068359375, 1148.000610351562),
            ( 727.5288085937500, 1931.427978515625),
            ( 798.1463623046875, 1296.653808593750),
            ( 262.4423828125000, 1037.846435546875),
            ( 470.4659729003906, 1565.640991210938),
            (1903.835449218750, 1966.276367187500),
            ( 626.0647583007812, 1703.881835937500),
            (1977.244140625000,  240.7441101074219),
            (1673.016357421875,  235.7854003906250),
            (1702.945556640625,  288.0748291015625),
            (2036.563476562500, 1721.330932617188),
            ( 536.5933227539062,  596.2082519531250),
            (2071.151855468750,  778.6830444335938),
            (1288.900634765625, 1580.191650390625),
            (1161.100952148438, 1173.102539062500),
        ]
        for (i, (x, y)) in enumerate(xy)
            @test rd.pcents[2, i, 1] === Float32(x)
            @test rd.pcents[2, i, 2] === Float32(y)
        end
        # rdarea grew every radius to hit the target area; final PAREA/OOAREA match.
        @test rd.pcents[2, 1, 3] === 525.50903f0     # grown radius (init 263.289)
        @test rd.parea[2]  === 99.022224f0           # SAREA·IN/75²
        @test rd.ooarea[2] === 99.022224f0
    end

    # -------------------------------------------------------------------------
    # Chunk 0b-2 — initial infection state (rd/rdinit.f host tables + rd/rdslp.f +
    # rd/rdiprp.f + the rd/rdsetp.f tree loop + rd/rdinoc.f). Golden PROPN/PROPI/
    # PROBI/PROBIU/FPROB captured from the LIVE relinked FVSkt oracle via a single-.o
    # swap of rdsetp.f dumping every per-record value at LSTART (RRType 3 Armillaria,
    # RRInit 0 10 10 20 0.1 10 3, SArea 100, seed 889347; the instrumented .sum was
    # byte-identical to FVSkt_clean). All 27 records × 5 fields verified 0-ULP
    # bit-exact (parsing the oracle's E18.10 dump to Float32). Not yet .sum-visible
    # (mortality is 0b-3), so the g16 dump IS the oracle here.
    # -------------------------------------------------------------------------
    @testset "RD initial infection (rdsetp/rdiprp) — bit-exact vs live FVSkt" begin
        s = nothing
        for st in FVSjl.each_stand(rd_key; variant = v)
            s = st; break
        end
        FVSjl.notre!(s)
        FVSjl.setup_growth!(s)              # runs root_disease_setup! → rd_setp!
        rd = s.root_disease; t = s.trees

        # RDSETP header scalars (oracle HDR dump).
        @test rd.yincpt === 2.2177f0                     # 1 - SDNORM·SDISLP (369·-0.0033)
        @test isapprox(rd.dimen, 2087.0f0; atol = 0.5f0)
        @test rd.rrgen1[3] === 1.0f0                     # RRGEN(3,1) after normalization
        @test rd.prkill[3] === 0.33333334f0              # 10/(10+20)
        @test rd.prun[3]   === 0.6666667f0               # 20/(10+20)
        @test rd.parea[3]  === 10.026667f0               # grid-matched disease area (rdarea)

        # 27 projectable records (ITRN = IREC1, 0 dead in FVS partition).
        @test t.n == 27

        # Per-species PROPN (rdiprp) — same value for every record of a species.
        propn_sp = Dict(2 => 0.18032679f0, 3 => 0.27565289f0, 4 => 0.32363024f0,
                        5 => 0.13220361f0, 7 => 0.75632572f0, 8 => 0.42801791f0)
        for i in 1:t.n
            @test rd.propn[i] === propn_sp[Int(t.species[i])]
        end

        # Full per-record golden vectors (oracle rdsetp dump), asserted bit-exact.
        GOLD_PROPI = Float32[0.07999999821186066, 0.10000000149011612, 0.07999999821186066, 0.11999999731779099, 0.14000000059604645, 0.14000000059604645, 0.07999999821186066, 0.10000000149011612, 0.11999999731779099, 0.03999999910593033, 0.07999999821186066, 0.03999999910593033, 0.14000000059604645, 0.07999999821186066, 0.1599999964237213, 0.11999999731779099, 0.20000000298023224, 0.10000000149011612, 0.03999999910593033, 0.20000000298023224, 0.14000000059604645, 0.0139999995008111, 0.05999999865889549, 0.07999999821186066, 0.07999999821186066, 0.07999999821186066, 0.03999999910593033]
        GOLD_PROBI = Float32[42.053104400634766, 248.74668884277344, 23.009193420410156, 21.24666976928711, 20.71882438659668, 61.90854263305664, 18.792587280273438, 61.62351989746094, 82.91556549072266, 19.720476150512695, 82.91556549072266, 82.91556549072266, 292.041015625, 84.71927642822266, 20.269702911376953, 63.95496368408203, 12.567241668701172, 60.34638214111328, 18.7404842376709, 76.9760971069336, 20.03000259399414, 22.939910888671875, 128.7465057373047, 128.7465057373047, 93.56011962890625, 125.8945083618164, 54.63187789916992]
        GOLD_PROBIU = Float32[13.548739433288574, 653.6443481445312, 151.034423828125, 96.57648468017578, 94.17716217041016, 129.38552856445312, 85.42147064208984, 19.853967666625977, 217.88145446777344, 89.63917541503906, 217.88145446777344, 217.88145446777344, 610.3500366210938, 177.0587158203125, 53.263729095458984, 133.6624298095703, 33.02357864379883, 19.44249725341797, 49.24531936645508, 24.8002872467041, 41.86162185668945, 60.280372619628906, 172.05050659179688, 172.05050659179688, 125.0291519165039, 168.23924255371094, 114.17767333984375]
        GOLD_FPROB = Float32[5.545452117919922, 90.00000762939453, 17.358247756958008, 11.751097679138184, 11.459156036376953, 19.07872200012207, 10.393794059753418, 8.126160621643066, 30.000001907348633, 10.906990051269531, 30.000001907348633, 30.000001907348633, 90.00000762939453, 26.108436584472656, 7.333859920501709, 19.709379196166992, 4.54700231552124, 7.957746982574463, 6.780566692352295, 10.15067195892334, 6.172763824462891, 8.299978256225586, 30.000001907348633, 30.000001907348633, 21.801010131835938, 29.335439682006836, 16.836227416992188]
        for i in 1:t.n
            @test rd.propi[i]  === GOLD_PROPI[i]
            @test rd.probi[i]  === GOLD_PROBI[i]
            @test rd.probiu[i] === GOLD_PROBIU[i]
            @test rd.fprob[i]  === GOLD_FPROB[i]
            @test rd.probl[i]  === GOLD_FPROB[i]         # PROBL(I) = PROB(I) = FPROB(I)
        end

        # rd_slp (rd/rdslp.f) knot interpolation — spot values on the YTKILL curve.
        @test FVSjl.rd_slp(0.0f0,  rd.xxinf, rd.yyinf, 3) === 0.0f0    # below XX(1)
        @test FVSjl.rd_slp(3.9f0,  rd.xxinf, rd.yyinf, 3) === 5.0f0    # exactly XX(2)
        @test FVSjl.rd_slp(50.0f0, rd.xxinf, rd.yyinf, 3) === 40.0f0   # above XX(3) → YY(3)
    end
end
