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

    # -------------------------------------------------------------------------
    # Chunk 0b-3 — the per-cycle mortality kernel (rd/rdmort.f + rd/rdsum.f).
    # Golden entry/exit state captured from the LIVE relinked FVSkt oracle: a
    # single-.o instrumentation swap of rdmort.f dumped every record's PROBI/PROPI
    # entry state and RRKILL/RDKILL/PROPI exit state for all 10 cycles of the
    # turnkey scenario (RRType 3 Armillaria, RRInit 0 10 10 20 0.1 10 3, SArea 100;
    # the instrumented .sum stayed byte-identical to FVSkt_clean). Driving
    # rd_mort_kernel! with the oracle's ENTRY state reproduced every EXIT value
    # bit-exact across the full run — 270 records, 254 nonzero-kill, 0 mismatches.
    # This testset asserts the complete CYCLE-1 block (ISTEP=2 ⇒ slots (1,1) initial
    # infection from rdsetp, (2,1) RDINSD inside-patch aging, (2,2) RDINF new
    # infection). The upstream chain that BUILDS the (2,*) entry slots and the RDEND
    # step that makes RRKILL .sum-visible are later sub-chunks — see the port report.
    # -------------------------------------------------------------------------
    @testset "RD mortality kernel (rdmort/rdsum) — bit-exact vs live FVSkt" begin
        G1_ISP = Int[7, 3, 5, 2, 2, 4, 2, 7, 3, 2, 3, 3, 4, 4, 3, 4, 3, 7, 3, 7, 4, 3, 8, 8, 8, 8, 4]
        G1_DBH = Float32[11.5, 0.107423656, 6.5, 7.9, 8.0, 6.2, 8.4, 9.5, 4.0, 8.2, 1.2, 1.9, 0.1, 5.3, 10.0, 6.1, 12.7, 9.6, 10.4, 8.5, 10.9, 9.4, 3.2, 0.1, 5.8, 5.0, 6.6]
        G1_PROBI_11 = Float32[42.053104, 248.74669, 23.009193, 21.24667, 20.718824, 61.908543, 18.792587, 61.62352, 82.915565, 19.720476, 82.915565, 82.915565, 292.04102, 84.71928, 20.269703, 63.954964, 12.567242, 60.346382, 18.740484, 76.9761, 20.030003, 22.93991, 128.7465, 128.7465, 93.56012, 125.89451, 54.631878]
        G1_PROPI_11 = Float32[0.08, 0.1, 0.08, 0.12, 0.14, 0.14, 0.08, 0.1, 0.12, 0.04, 0.08, 0.04, 0.14, 0.08, 0.16, 0.12, 0.2, 0.1, 0.04, 0.2, 0.14, 0.0139999995, 0.06, 0.08, 0.08, 0.08, 0.04]
        G1_PROBI_21 = Float32[7.3682227, 17.767488, 25.147968, 11.330256, 12.265579, 81.90156, 11.054939, 11.137711, 66.64456, 10.967106, 24.348478, 29.870813, 22.45291, 89.70798, 42.717327, 86.674324, 31.521011, 8.642616, 40.65856, 10.926502, 39.378014, 50.524006, 33.259735, 6.001759, 66.046585, 73.94884, 79.15031]
        G1_PROPI_21 = Float32[-0.6544206, -0.45166653, -0.60440654, -0.6084938, -0.6205581, -0.6337918, -0.61797094, -0.5690732, -0.62119454, -0.64291286, -0.6099856, -0.58321816, -0.5647853, -0.63445, -0.65258735, -0.6492604, -0.62875396, -0.62178385, -0.6483234, -0.65945894, -0.60518396, -0.64650244, -0.5804697, -0.4810371, -0.62053853, -0.6275188, -0.6224158]
        G1_PROBI_22 = Float32[1.8139817, 73.60012, 2.8390431, 0.9609802, 0.93710595, 18.722614, 0.84998286, 2.6581614, 24.533373, 0.891951, 24.533373, 24.533373, 88.320145, 25.62112, 5.9974766, 19.3415, 3.718443, 2.6030712, 5.545005, 3.3204024, 6.0575476, 6.7875476, 24.533371, 24.533371, 17.828411, 23.989906, 16.521976]
        G1_PROPI_22 = Float32[0.05, 0.04, 0.06, 0.03, 0.02, 0.02, 0.08, 0.04, 0.05, 0.03, 0.03, 0.03, 0.04, 0.06, 0.06, 0.04, 0.03, 0.05, 0.04, 0.02, 0.04, 0.04, 0.03, 0.02, 0.06, 0.03, 0.04]
        G1_RRKILL = Float32[0.0, 340.11432, 25.848236, 0.0, 0.0, 80.63116, 0.0, 0.0, 174.09349, 0.0, 131.79742, 137.31975, 402.8141, 200.04837, 20.269703, 83.29646, 0.0, 0.0, 0.0, 76.9761, 26.087551, 0.0, 186.53961, 159.28163, 111.388535, 149.88441, 71.153854]
        G1_APROPI_11 = Float32[0.20396695, 58.187775, 1.2067606, 0.64941174, 0.6632558, 1.5517648, 0.58000004, 0.24851489, 1.6852175, 0.55136365, 5.2799997, 3.3242106, 83.340004, 1.7071187, 0.8392453, 1.552836, 0.74135345, 0.24705881, 0.6945455, 0.36483517, 0.9747827, 0.734, 1.7219319, 53.26182, 1.0388069, 1.1757792, 1.3733336]
        G1_APROPI_21 = Float32[-0.5304537, 57.6361, 0.5223541, -0.07908207, -0.09730226, 0.777973, -0.117970906, -0.42055836, 0.94402283, -0.13154915, 4.5900145, 2.7009926, 82.63522, 0.9926687, 0.026658043, 0.78357553, -0.08740053, -0.474725, 0.006222132, -0.49462375, 0.22959876, 0.07349758, 1.0814621, 52.700783, 0.3382683, 0.46826038, 0.7109177]
        G1_APROPI_22 = Float32[0.17396696, 58.127773, 1.1867608, 0.55941176, 0.5432558, 1.4317648, 0.58000004, 0.18851486, 1.6152174, 0.54136366, 5.23, 3.3142107, 83.240005, 1.6871188, 0.73924536, 1.4728359, 0.57135344, 0.19705881, 0.6945455, 0.18483518, 0.8747828, 0.76000005, 1.691932, 53.201817, 1.0188068, 1.1257793, 1.3733336]
        PAREA3 = 11.662223f0

        n = 27; istep = 2
        rd = FVSjl.RootDiseaseState()
        rd.minrr = Int32(3); rd.maxrr = Int32(3); rd.irhab = Int32(1)
        rd.parea[3] = PAREA3
        probi = zeros(Float32, n, istep, 2)
        propi = zeros(Float32, n, istep, 2)
        for i in 1:n
            probi[i,1,1] = G1_PROBI_11[i]; propi[i,1,1] = G1_PROPI_11[i]
            probi[i,2,1] = G1_PROBI_21[i]; propi[i,2,1] = G1_PROPI_21[i]
            probi[i,2,2] = G1_PROBI_22[i]; propi[i,2,2] = G1_PROPI_22[i]
        end
        rrkill = zeros(Float32, n); rdkill = zeros(Float32, n)
        FVSjl.rd_mort_kernel!(rd, probi, propi, rrkill, rdkill,
                              G1_DBH, G1_ISP, istep, 10.0f0)

        for i in 1:n
            @test rrkill[i] === G1_RRKILL[i]          # killed infected TPA (RRKILL)
            @test rdkill[i] === G1_RRKILL[i]          # RDKILL == RRKILL within a cycle
            @test propi[i,1,1] === G1_APROPI_11[i]    # aged proportion-of-roots-infected
            @test propi[i,2,1] === G1_APROPI_21[i]
            @test propi[i,2,2] === G1_APROPI_22[i]
            # a killed slot has its PROBI zeroed; a survivor keeps it unchanged
            @test probi[i,1,1] === (G1_APROPI_11[i] < FVSjl.RD_PKILLS[Int(rd.irtspc[G1_ISP[i]]), 3] ? G1_PROBI_11[i] : 0.0f0)
        end

        # rd/rdsum.f — PROBIT is the per-record sum over all (IT,IP) slots.
        probit = zeros(Float32, n)
        FVSjl.rd_sum!(probit, probi, istep)
        for i in 1:n
            @test probit[i] === probi[i,1,1] + probi[i,1,2] + probi[i,2,1] + probi[i,2,2]
        end

        # TAREA ≤ 0 (no disease area) ⇒ kernel is a no-op that only zeroes RDKILL.
        rd0 = FVSjl.RootDiseaseState(); rd0.minrr = Int32(3); rd0.maxrr = Int32(3)
        rd0.parea .= 0.0f0
        p2 = copy(probi); pr2 = copy(propi); rk = fill(9.0f0, n); rd2 = fill(9.0f0, n)
        FVSjl.rd_mort_kernel!(rd0, p2, pr2, rk, rd2, G1_DBH, G1_ISP, istep, 10.0f0)
        @test p2 == probi && pr2 == propi        # untouched
        @test all(==(0.0f0), rd2)                # RDKILL zeroed
        @test all(==(9.0f0), rk)                 # RRKILL accumulator left alone
    end

    @testset "RD mortality→WK2 reconciliation (rdend) — bit-exact vs live FVSkt" begin
        # (WK2_before, PROB, RRKILL, PROBIT, WK2_after) captured from a single-.o
        # rd/rdend.f instrumentation swap on the live FVSkt oracle (cycle 1, all 27
        # host records; the instrumented .sum stayed byte-identical to FVSkt_clean).
        # OAKL=0 in the turnkey scenario ⇒ TDIUN=TDIEOU=0 ⇒ WMESS=TDIEN/SAREA.
        RDEND_G1 = [
(0.3851258158684f0,0.5545452117920f1,0.0000000000000f0,0.5123530578613f2, 0.3851258158684f0),
(0.2900475120544f2,0.9000000762939f2,0.3401143188477f3,0.0000000000000f0, 0.2900475120544f2),
(0.5265119075775f0,0.1735824775696f2,0.2584823608398f2,0.2514796829224f2, 0.5265119075775f0),
(0.8497061133385f0,0.1175109767914f2,0.0000000000000f0,0.3353790664673f2, 0.8497061133385f0),
(0.4544001817703f0,0.1145915603638f2,0.0000000000000f0,0.3392151260376f2, 0.4544001817703f0),
(0.9149742722511f0,0.1907872200012f2,0.8063115692139f2,0.8190155792236f2, 0.9149742722511f0),
(0.7114505767822f0,0.1039379405975f2,0.0000000000000f0,0.3069750976562f2, 0.7114505767822f0),
(0.5937915444374f0,0.8126160621643f1,0.0000000000000f0,0.7541939544678f2, 0.5937915444374f0),
(0.4530912876129f1,0.3000000190735f2,0.1740934906006f3,0.0000000000000f0, 0.4530912876129f1),
(0.7029584050179f0,0.1090699005127f2,0.0000000000000f0,0.3157953262329f2, 0.7029584050179f0),
(0.7142680168152f1,0.3000000190735f2,0.1317974243164f3,0.0000000000000f0, 0.7142680168152f1),
(0.5775949954987f1,0.3000000190735f2,0.1373197479248f3,0.0000000000000f0, 0.5775949954987f1),
(0.3088256645203f2,0.9000000762939f2,0.4028140869141f3,0.0000000000000f0, 0.3088256645203f2),
(0.7539042830467f0,0.2610843658447f2,0.2000483703613f3,0.0000000000000f0, 0.2000483751297f1),
(0.3027817308903f0,0.7333859920502f1,0.2026970291138f2,0.4871480560303f2, 0.3027817308903f0),
(0.7808734774590f0,0.1970937919617f2,0.8329646301270f2,0.8667432403564f2, 0.8329646587372f0),
(0.2041696161032f0,0.4547002315521f1,0.0000000000000f0,0.4780669403076f2, 0.2041696161032f0),
(0.8079369068146f0,0.7957746982574f1,0.0000000000000f0,0.7159207153320f2, 0.8079369068146f0),
(0.2555747926235f0,0.6780566692352f1,0.0000000000000f0,0.6494405364990f2, 0.2555747926235f0),
(0.4972422122955f0,0.1015067195892f2,0.7697609710693f2,0.1424690437317f2, 0.7697609663010f0),
(0.8882896602154f-1,0.6172763824463f1,0.2608755111694f2,0.3937801361084f2, 0.2608755230904f0),
(0.1456833481789f0,0.8299978256226f1,0.0000000000000f0,0.8025145721436f2, 0.1456833481789f0),
(0.1683012962341f1,0.3000000190735f2,0.1865396118164f3,0.0000000000000f0, 0.1865396142006f1),
(0.6447217464447f1,0.3000000190735f2,0.1592816314697f3,0.0000000000000f0, 0.6447217464447f1),
(0.5072642564774f0,0.2180101013184f2,0.1113885345459f3,0.6604658508301f2, 0.1113885402679f1),
(0.2338130474091f1,0.2933543968201f2,0.1498844146729f3,0.7394883728027f2, 0.2338130474091f1),
(0.3865458667278f0,0.1683622741699f2,0.7115385437012f2,0.7915030670166f2, 0.7115385532379f0)
]
        for (wk2b, prob, rrkill, probit, wk2a) in RDEND_G1
            @test FVSjl.rd_end_newwk2(wk2b, prob, rrkill, probit, 0.0f0, 0.0f0, 100.0f0) === wk2a
        end

        # Full rd_end_kernel! (the array-level RDEND body): drive it with the golden
        # entry state (oakl/bbkill all zero) and confirm every WK2 matches bit-exact.
        n = length(RDEND_G1)
        rd = FVSjl.RootDiseaseState()
        rd.minrr = Int32(3); rd.maxrr = Int32(3); rd.sarea = 100.0f0
        rd.parea[3] = 11.66222f0                       # cycle-1 PAREA (any >0 works; OAKL=0)
        wk2  = Float32[r[1] for r in RDEND_G1]
        prob = Float32[r[2] for r in RDEND_G1]
        rrk  = Float32[r[3] for r in RDEND_G1]
        pbit = Float32[r[4] for r in RDEND_G1]
        exp  = Float32[r[5] for r in RDEND_G1]
        rdk  = zeros(Float32, n); pbiu = fill(1.0f0, n); fpr = fill(1.0f0, n)
        isp  = fill(Int32(1), n)                       # host sp 1 → base-RD sp 1
        rootl = fill(3.0f0, n); wk22 = zeros(Float32, n); rroott = zeros(Float32, n)
        probi = zeros(Float32, n, 2, 2); propi = zeros(Float32, n, 2, 2)
        dprob = zeros(Float32, n, 4, 3)
        oakl = zeros(Float32, 3, n); bbkill = zeros(Float32, 3, n)
        FVSjl.rd_end_kernel!(rd, wk2, prob, rrk, rdk, pbit, pbiu, fpr, probi, propi,
                             isp, rootl, wk22, rroott, dprob, oakl, bbkill)
        for i in 1:n
            @test wk2[i] === exp[i]
        end
    end

    @testset "RD growth-loss (rdgrow) — bit-exact vs live FVSkt" begin
        # (DG_b, HTG_b, DGTOT, HTTOT, OUTNUM, PROBIU, PROBIT, PAREA, DG_a, HTG_a)
        # from a single-.o rd/rdgrow.f instrumentation swap (cycle 1, 27 records).
        RDGROW_G1 = [
(0.6511259078979f0,0.5735917091370f1,0.6528559923172f0,0.6528559923172f0,0.4898734436035f3,0.1250329780579f2,0.4767707061768f2,0.1166222286224f2, 0.6315338611603f0,0.5563326358795f1),
(0.0000000000000f0,0.2423654317856f1,0.0000000000000f0,0.0000000000000f0,0.7950409179688f4,0.7094769287109f3,0.0000000000000f0,0.1166222286224f2, 0.9999999747379f-4,0.2423654317856f1),
(0.1177871227264f1,0.6388757705688f1,0.0000000000000f0,0.0000000000000f0,0.1533390869141f4,0.1514378356934f3,0.2514796829224f2,0.1166222286224f2, 0.1160548686981f1,0.6294800758362f1),
(0.7816820144653f0,0.7996763706207f1,0.3378342986107f0,0.3378342986107f0,0.1038066894531f4,0.9602055358887f2,0.3111282539368f2,0.1166222286224f2, 0.7678611874580f0,0.7855373859406f1),
(0.1189398288727f1,0.1017678070068f2,0.3615870773792f0,0.3615870773792f0,0.1012277526855f4,0.9576243591309f2,0.3257638931274f2,0.1166222286224f2, 0.1167711615562f1,0.9991224288940f1),
(0.9328789710999f0,0.7552839756012f1,0.0000000000000f0,0.0000000000000f0,0.1685373779297f4,0.5996571350098f2,0.8190155792236f2,0.1166222286224f2, 0.8910650014877f0,0.7214302062988f1),
(0.7851281166077f0,0.5675306797028f1,0.3601249456406f0,0.3601249456406f0,0.9181656494141f3,0.8432041168213f2,0.2859627914429f2,0.1166222286224f2, 0.7711948752403f0,0.5574590206146f1),
(0.7448329925537f0,0.6368978023529f1,0.5805994868279f0,0.5805994868279f0,0.7178477783203f3,0.1793505477905f2,0.6990837860107f2,0.1166222286224f2, 0.7177280783653f0,0.6137206554413f1),
(0.5277636051178f0,0.7001430511475f1,0.0000000000000f0,0.0000000000000f0,0.2650136474609f4,0.1757702636719f3,0.0000000000000f0,0.1166222286224f2, 0.5277636051178f0,0.7001430988312f1),
(0.8311281204224f0,0.9308626174927f1,0.3472852110863f0,0.3472852110863f0,0.9635003662109f3,0.8945645904541f2,0.2954422378540f2,0.1166222286224f2, 0.8163222074509f0,0.9142800331116f1),
(0.5064020156860f0,0.4937581062317f1,0.0000000000000f0,0.0000000000000f0,0.2650136474609f4,0.2180663452148f3,0.0000000000000f0,0.1166222286224f2, 0.5064020156860f0,0.4937581062317f1),
(0.5526168346405f0,0.5379607677460f1,0.0000000000000f0,0.0000000000000f0,0.2650136474609f4,0.2125440063477f3,0.0000000000000f0,0.1166222286224f2, 0.5526168346405f0,0.5379607677460f1),
(0.3026429414749f0,0.2824449539185f1,0.0000000000000f0,0.0000000000000f0,0.7950410156250f4,0.6467772216797f3,0.0000000000000f0,0.1166222286224f2, 0.3026429414749f0,0.2824449539185f1),
(0.1438930988312f1,0.8416152000427f1,0.0000000000000f0,0.0000000000000f0,0.2306364013672f4,0.1044314880371f3,0.0000000000000f0,0.1166222286224f2, 0.1438930988312f1,0.8416152000427f1),
(0.7862329483032f0,0.8269799232483f1,0.8301337361336f0,0.8301337361336f0,0.6478575439453f3,0.1654387855530f2,0.4871480560303f2,0.1166222286224f2, 0.7771095037460f0,0.8173836708069f1),
(0.1098357677460f1,0.8128564834595f1,0.0000000000000f0,0.0000000000000f0,0.1741084594727f4,0.5988243865967f2,0.8667432403564f2,0.1166222286224f2, 0.1047924637794f1,0.7755327701569f1),
(0.5648336410522f0,0.6642565250397f1,0.6593430638313f0,0.6593430638313f0,0.4016724853516f3,0.4986576080322f1,0.4566007614136f2,0.1166222286224f2, 0.5454100370407f0,0.6414139747620f1),
(0.5176057815552f0,0.5613164424896f1,0.5691684484482f0,0.5691684484482f0,0.7029704589844f3,0.1905853271484f2,0.6432344055176f2,0.1166222286224f2, 0.4993643760681f0,0.5415345668793f1),
(0.8411388397217f0,0.8361371994019f1,0.6182643771172f0,0.6182643771172f0,0.5989808349609f3,0.1359910488129f2,0.6249616622925f2,0.1166222286224f2, 0.8114132285118f0,0.8065882682800f1),
(0.1098293304443f1,0.7464504718781f1,0.9138440489769f0,0.9138440489769f0,0.8966887817383f3,0.2715539550781f2,0.1424690437317f2,0.1166222286224f2, 0.1096856236458f1,0.7454737663269f1),
(0.2426842689514f1,0.1242214488983f2,0.5408024787903f0,0.5408024787903f0,0.5452887573242f3,0.6521972656250f1,0.3937801361084f2,0.1166222286224f2, 0.2352614641190f1,0.1204219722748f2),
(0.1816849708557f1,0.1063282775879f2,0.5370271801949f0,0.5370271801949f0,0.7332023925781f3,0.1625353240967f2,0.7884287261963f2,0.1166222286224f2, 0.1736783266068f1,0.1016425228119f2),
(0.9326717853546f0,0.6665924072266f1,0.0000000000000f0,0.0000000000000f0,0.2650136230469f4,0.1633241424561f3,0.0000000000000f0,0.1166222286224f2, 0.9326717853546f0,0.6665924549103f1),
(0.2000000178814f0,0.2868623256683f1,0.0000000000000f0,0.0000000000000f0,0.2650136230469f4,0.1905821228027f3,0.0000000000000f0,0.1166222286224f2, 0.2000000029802f0,0.2868623256683f1),
(0.1354193687439f1,0.8921731948853f1,0.3234633803368f0,0.3234633803368f0,0.1925855102539f4,0.7681097412109f2,0.6604658508301f2,0.1166222286224f2, 0.1324943780899f1,0.8729028701782f1),
(0.5371613502502f0,0.7228552818298f1,0.6347924470901f-1,0.6347924470901f-1,0.2591430175781f4,0.1182803115845f3,0.7394883728027f2,0.1166222286224f2, 0.5237973928452f0,0.7048713684082f1),
(0.1633671283722f1,0.9058808326721f1,0.0000000000000f0,0.0000000000000f0,0.1487276611328f4,0.4604201507568f2,0.7915030670166f2,0.1166222286224f2, 0.1553480148315f1,0.8614144325256f1)
]
        for (dgb,htgb,dgtot,httot,outnum,probiu,probit,parea, dga,htga) in RDGROW_G1
            diff = 100.0f0 - parea
            @test FVSjl.rd_grow_newg(dgb, dgtot, outnum, probiu, probit, diff) === dga
            @test FVSjl.rd_grow_newg(htgb, httot, outnum, probiu, probit, diff) === htga
        end

        # rd_grow_gtot: a single infection slot with proportion PROPI contributes
        # RDSLP(max(PROPI,0),XDBH,YDBH,2)·PROBI/(PROBIT+1e-6); a NEGATIVE PROPI
        # (infection not yet at the tree center) contributes nothing (PRADI clamps to 0
        # ⇒ RDSLP=1, i.e. no growth loss factor <1 — full growth for that slot).
        probi = zeros(Float32, 1, 2, 2); propi = zeros(Float32, 1, 2, 2)
        probi[1,1,1] = 10.0f0; propi[1,1,1] = 0.25f0     # RDSLP(0.25) = 0.5
        gt = FVSjl.rd_grow_gtot(probi, propi, 1, 2, 10.0f0, FVSjl.RD_XDBH, FVSjl.RD_YDBH, 1.0f0)
        @test gt === Float32(0.5f0 * 10.0f0 / (10.0f0 + 1.0f-6))
        propi[1,1,1] = -0.3f0                             # negative ⇒ PRADI=0 ⇒ RDSLP=1
        gt2 = FVSjl.rd_grow_gtot(probi, propi, 1, 2, 10.0f0, FVSjl.RD_XDBH, FVSjl.RD_YDBH, 1.0f0)
        @test gt2 === Float32(1.0f0 * 10.0f0 / (10.0f0 + 1.0f-6))
    end

    @testset "RD cycle advance (rdmn2)" begin
        rd = FVSjl.RootDiseaseState()
        rd.minrr = Int32(3); rd.maxrr = Int32(3); rd.parea[3] = 10.0f0
        rd.iroot = Int32(1); rd.istep = Int32(1)
        n = 3
        probi = zeros(Float32, n, 2, 2)
        probi[1,1,1] = 5.0f0; probi[2,2,1] = 2.0f0; probi[2,2,2] = 1.0f0
        rrkill = fill(7.0f0, n); probit = zeros(Float32, n)
        FVSjl.rd_mn2_advance!(rd, rrkill, probit, probi, 10.0f0)
        @test rd.istep == Int32(2)                       # ISTEP advanced
        @test all(==(0.0f0), rrkill)                     # RRKILL zeroed
        @test probit[1] === 5.0f0                         # PROBIT re-summed over slots
        @test probit[2] === 3.0f0
        # inert when RD off
        rd0 = FVSjl.RootDiseaseState(); rd0.iroot = Int32(0)
        s0 = Int(rd0.istep)
        FVSjl.rd_mn2_advance!(rd0, Float32[], Float32[], zeros(Float32,0,2,2), 10.0f0)
        @test Int(rd0.istep) == s0
    end

    @testset "RD live-tree root radius (rdroot) — bit-exact vs live FVSkt" begin
        # Upstream spread-chain routine rd/rdroot.f (RDTREG DO-1001), the deterministic
        # ROOTL(record) that feeds RDSPRD/RDJUMP/RDINSD. Golden per-tree inputs+ROOTL
        # dumped at E22.13 from the g16-instrumented FVSkt (rd/rdtreg.f, turnkey rd.key,
        # stand S248112, RRType 3 Armillaria). Each tuple:
        #   (DBH, HT, PROOT, RSLOP, expected ROOTL); scalars are (SDISLP,YINCPT,OLDTPA,
        #   GROSPC,ORMSQD,BA) for that cycle. Full Float32 identity (===) incl. ^1.605.
        sc1 = (-0.329999998211861f-02, 2.21770000457764f0, 589.652709960938f0,
               1.10000002384186f0, 5.14496803283691f0, 85.1312713623047f0)
        cyc1 = [
            (11.5f0,              73.0f0, 14.2600002288818f0, 1.0f0, 20.4987487792969f0),
            (0.107423655688763f0,  2.0f0, 14.2600002288818f0, 1.0f0,  1.22190737724304f0),  # sub-3.5" HT/BA allometry
            (6.5f0,               30.0f0, 14.5f0,             1.0f0, 11.78125f0),
            (7.9f0,               75.0f0, 14.2600002288818f0, 1.0f0, 14.0817499160767f0),
            (8.0f0,               63.0f0, 14.2600002288818f0, 1.0f0, 14.2599992752075f0),
        ]
        for (dbh, ht, proot, rslop, gold) in cyc1
            @test FVSjl.rd_root(dbh, ht, proot, rslop, sc1...) === gold
        end
        sc2 = (-0.329999998211861f-02, 2.21770000457764f0, 489.120605468750f0,
               1.10000002384186f0, 6.32937335968018f0, 106.872207641602f0)
        cyc2 = [
            (12.1517381668091f0,   78.5633239746094f0, 14.2600002288818f0, 1.0f0, 21.6604709625244f0),
            (0.107538998126984f0,   4.42365455627441f0, 14.2600002288818f0, 1.0f0,  2.00382804870605f0),  # sub-3.5"
            (7.74255752563477f0,   36.2947998046875f0, 14.5f0,             1.0f0, 14.0333852767944f0),
            (8.80230426788330f0,   82.8553771972656f0, 14.2600002288818f0, 1.0f0, 15.6901063919067f0),
        ]
        for (dbh, ht, proot, rslop, gold) in cyc2
            @test FVSjl.rd_root(dbh, ht, proot, rslop, sc2...) === gold
        end
    end

end
