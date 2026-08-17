# Douglas-fir Beetle (DFB) impact model — beachhead (keyword reader + the four
# DETERMINISTIC routines DFBIND / DFBDBH / DFBER / DFBPRB), plus the inert-seam
# guarantee.
#
# DOCTRINE: the deterministic routines are validated BIT-EXACT (Float32) against
# an instrumented g16 build of the PRISTINE dfb/*.f Fortran. A standalone driver
# (scratchpad/dfb/driver_dfb.f) links the real dfber.o/dfbprb.o/dfbdbh.o/dfbind.o
# (+ dfblkdbc.o block data), populates the ARRAYS/CONTRL/PLOT/DFBCOM commons with
# the synthetic all-Douglas-fir stand below, and dumps Float32 bit patterns
# (scratchpad/dfb/dfb_golden.txt). The golden hex constants embedded here are that
# dump verbatim; the Julia port must reproduce every bit.
#
# Golden stand (12 DF records, IND1 order = tree order):
#   DBH  3.0  5.0  8.0  8.7  9.0 10.0 12.4 15.0 18.0 20.0 24.0 30.0
#   TPA 40.0 30.0 25.0 20.0 18.0 15.0 12.0 10.0  8.0  6.0  4.0  2.0
#   CFV = DBH²·0.30 ,  BA = Σ 0.005454154·DBH²·TPA  (Float32, in order)

using Test
using FVSjl
const _F = FVSjl

_hex(x::Float32) = uppercase(string(reinterpret(UInt32, x); base = 16, pad = 8))

const _DFB_DBH = Float32[3.0, 5.0, 8.0, 8.7, 9.0, 10.0, 12.4, 15.0, 18.0, 20.0, 24.0, 30.0]
const _DFB_TPA = Float32[40.0, 30.0, 25.0, 20.0, 18.0, 15.0, 12.0, 10.0, 8.0, 6.0, 4.0, 2.0]

# --- golden Float32 bit patterns from the instrumented FVSbc/dfb g16 oracle ---
const _G_BA     = "42DE3BEC"   # stand BA
const _G_BA9    = "42B028DA"   # BADF9==BA9 (all DF)
const _G_A45DBH = "4125D868"
const _G_PBADF4 = "3F7B79F3"
const _G_PROTBK = "3F666666"   # BADF9/BA9 = 0.982 > 0.9 ⇒ capped
const _G_START = ("00000000","42200000","41F00000","42340000","42040000",
                  "41400000","00000000","41200000","41000000","40C00000",
                  "00000000","40800000","00000000","00000000","40000000",
                  "00000000","00000000","00000000","00000000","00000000")
# DFBIND: (dbh-hex, expected size class) for a spread of diameters
const _G_DFBIND = (("3ECCCCCD",1),("3F800000",1),("3FF33333",1),("40000000",1),
                   ("410B3333",4),("41100000",5),("41200000",5),("41466666",6),
                   ("41980000",10),("41A00000",10),("41C00000",12),("41F00000",15),
                   ("42200000",20),("42340000",20))

# ---- keyword-reader fixtures (top level; const inside @testset is unsupported) ----
const _DFB_TRE = """
   1      248112       0101   011LP 11510   0734   00111     0  0
   2      248112       0101   031DF 001     0026   00222     0  0
   3      248112       0102   011WH 06523   0308   00111     0  0
   4      248112       0102   011WL 07906   0753   00111     0  0
   6      248112       0103   011WL 08007   0633   96222     0 56
   7      248112       0103   011GF 06220   0385   00111     0  0
   9      248112       0103   011LP 09511   0603   00111     0  0
  10      248112       0104   011DF 040     0203   00111    50  0
  17      248112       0106   011DF 10010   0654   00111     0  0
  19      248112       0106   011DF 12716   0674   00111     0  0
  22      248112       0108   011DF 10409   0555   97222     0 49
  25      248112       0109   011DF 09418   0604   00111     0  0
"""
_dfb_head(title) = """
SCREEN
NOAUTOES
NOTRIPLE
STATS
STDIDENT
S248112  $title
DESIGN                                        11.0       1.0
STDINFO     11406001     570.0      60.0     315.0      30.0      34.0
INVYEAR       1990.0
NUMCYCLE         5.0
TREEFMT
(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,
T52,I2,T66,5I1,T54,7I1,T75,F3.0)
TREEDATA
"""
const _DFB_BLOCK = """
DFBEETLE
MANSTART
MANSCHED           1
OLENGTH            6
EXYRMORT         8.0       3.0
MORTDIS
STOPROB          0.4
END
"""

@testset "Douglas-fir Beetle (DFB) — beachhead" begin

    @testset "DFBIND size-class mapping (bit-exact vs g16)" begin
        for (dh, cls) in _G_DFBIND
            d = reinterpret(Float32, parse(UInt32, dh; base = 16))
            @test _F.dfb_ind(d) == cls
        end
        # explicit boundary cases from dfbind.f
        @test _F.dfb_ind(0.4f0) == 1     # < 1"  ⇒ class 1
        @test _F.dfb_ind(9.0f0) == 5     # IFIX 9 → 9/2+1 = 5
        @test _F.dfb_ind(8.7f0) == 4     # IFIX 8 → 8/2+0 = 4
        @test _F.dfb_ind(45.0f0) == 20   # 23 → clamped to 20
    end

    # BA accumulated exactly as the Fortran driver (Float32, in tree order)
    ba = 0.0f0
    for i in eachindex(_DFB_DBH)
        ba += _F.DFB_BAF * _DFB_DBH[i] * _DFB_DBH[i] * _DFB_TPA[i]
    end

    @testset "DFBER stand statistics (bit-exact vs g16)" begin
        @test _hex(ba) == _G_BA
        r = _F.dfb_er(_DFB_DBH, _DFB_TPA, _DFB_DBH, _DFB_TPA, ba)
        @test _hex(r.ba9)    == _G_BA9
        @test _hex(r.badf9)  == _G_BA9      # all DF ⇒ BADF9 == BA9
        @test _hex(r.a45dbh) == _G_A45DBH
        @test _hex(r.pbadf4) == _G_PBADF4
        @test r.lmin == true                # ≥1 TPA of DF ≥4.5"
    end

    @testset "DFBPRB outbreak probability (bit-exact vs g16)" begin
        r = _F.dfb_er(_DFB_DBH, _DFB_TPA, _DFB_DBH, _DFB_TPA, ba)
        p = _F.dfb_prb(r.a45dbh, r.pbadf4, r.badf9, r.ba9)
        @test _hex(p) == _G_PROTBK
        # gate branches (dfbprb.f)
        @test _F.dfb_prb(8.9f0, 0.5f0, 10.0f0, 20.0f0) == 0.0f0   # A45DBH < 9
        @test _F.dfb_prb(9.0f0, 0.2f0, 10.0f0, 20.0f0) == 0.0f0   # PBADF4 < 0.25
        @test _F.dfb_prb(9.0f0, 0.5f0, 10.0f0, 20.0f0) == 0.5f0   # BADF9/BA9, uncapped
    end

    @testset "DFBDBH START array (bit-exact vs g16)" begin
        st = _F.dfb_dbh_start(_DFB_DBH, _DFB_TPA)
        @test length(st) == 20
        for i in 1:20
            @test _hex(st[i]) == _G_START[i]
        end
    end

    # ---- keyword reader (dfbin.f) via the real keyword dispatch ----
    dir = mktempdir()
    write(joinpath(dir, "shared.tre"), _DFB_TRE)
    ctrl_key = joinpath(dir, "ctrl.key")
    dfb_key  = joinpath(dir, "dfb.key")
    write(ctrl_key, _dfb_head("DFB CONTROL") * "ECHOSUM\nPROCESS\nSTOP\n")
    write(dfb_key,  _dfb_head("DFB PARSED ") * _DFB_BLOCK * "ECHOSUM\nPROCESS\nSTOP\n")
    cp(joinpath(dir, "shared.tre"), joinpath(dir, "ctrl.tre"))
    cp(joinpath(dir, "shared.tre"), joinpath(dir, "dfb.tre"))
    v = FVSjl.Kootenai()   # engine-of-convenience: the DFB seam is variant-independent & inert

    @testset "DFBEETLE keyword reader populates s.dfb (dfbin.f)" begin
        d = nothing
        for s in FVSjl.each_stand(dfb_key; variant = v)
            d = s.dfb
            break
        end
        @test d !== nothing
        @test d.active == true
        @test d.ismeth == 1               # MANSTART
        @test d.idbsch == 1               # MANSCHED
        @test d.mansched_years == Int32[1]
        @test d.ilenth == 6               # OLENGTH 6
        @test d.expctd == 8.0f0           # EXYRMORT
        @test d.exstdv == 3.0f0
        @test d.lbamod == true            # MORTDIS ⇒ BA method
        @test d.lepi == true              # STOPROB
        @test d.epiprb == 0.4f0
    end

    @testset "reader defaults (unset sub-keywords keep DFBINT defaults)" begin
        d = _F.dfb_defaults!()
        @test d.ismeth == 2 && d.idbsch == 2 && d.ilenth == 4
        @test d.expctd == 6.0f0 && d.exstdv == 2.0f0
        @test d.orseed == 55329.0f0 && d.dbevnt == 0.05f0
        @test d.prpwin == 0.8f0 && d.mwinht == 20.0f0
        @test d.lbamod == false && d.lepi == false && d.linprg == false
    end

    @testset "INERT seam — DFBEETLE block leaves the .sum byte-identical" begin
        # No per-cycle DFB engine seam is wired yet, so a stand carrying an active
        # DfbState must project byte-identically to the same stand without DFB.
        sumrows(key) = filter(l -> !startswith(l, "-999"),
                              split(strip(FVSjl.run_keyfile(key; variant = v, output = :sum)), '\n'))
        ctrl = sumrows(ctrl_key)
        dfb  = sumrows(dfb_key)
        @test !isempty(ctrl)
        @test ctrl == dfb
    end
end
