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
DFB
MANSTART
MANSCHED           1
OLENGTH            6
EXYRMORT         8.0       3.0
MORTDIS
STOPROB          0.4
END
"""
# DFB block that fires an outbreak in cycle 2 (MANSCHED 2 → FVS ICYC 2), MANSTART deterministic gate,
# defaults for EXYRMORT (6.0/2.0) — the exact configuration of the dump-replay golden below.
const _DFB_FIRE = """
DFB
MANSTART
MANSCHED           2
END
"""
# DFB block scheduled beyond the 5-cycle run (MANSCHED 99) — no outbreak fires, so the seam is inert.
const _DFB_NEVER = """
DFB
MANSTART
MANSCHED          99
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

    @testset "DFB keyword reader populates s.dfb (dfbin.f)" begin
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

    # -------------------------------------------------------------------------
    # DFBRAN + BACHLO + DFBMOD + DFBMRT — dump-replay BIT-EXACT vs a relinked
    # FVSie_dfb g16 oracle (dfb/*.o swapped for the exdfb.o no-op; iet01 stand,
    # DFB MANSTART + MANSCHED 2, cycle-2 outbreak).  The golden hex constants are
    # that oracle's stderr dump verbatim (scratchpad/dfb): the DFBRAN stream from
    # the DFBSCH-reset seed ORSEED+1128 = 56457, DFBMOD's DFKILL, and DFBMRT's
    # per-record TAMORT.  Every bit must reproduce.
    # -------------------------------------------------------------------------
    @testset "DFBRAN Lehmer stream (bit-exact vs g16 oracle, seed 56457)" begin
        d = _F.dfb_defaults!(); d.active = true    # orseed 55329 ⇒ lazy-seed 56457
        for g in ("3EE23A99", "3E6A6620", "3E56B516")
            @test _hex(_F.dfb_rand!(d)) == g
        end
    end

    @testset "DFBMOD DFKILL (bit-exact vs g16 oracle)" begin
        d = _F.dfb_defaults!(); d.active = true     # EXPCTD 6.0, EXSTDV 2.0 (defaults)
        badf9 = reinterpret(Float32, 0x41B78CB8)    # DF ≥9″ basal area, cycle 2
        ba9   = reinterpret(Float32, 0x42824F2E)    # stand ≥9″ basal area
        dfkill = _F.dfb_mod(d, badf9, ba9, 4, 2)    # NUMYRS = min(ILENTH 4, IFINT 10); ICYC 2 (new-outbreak)
        @test _hex(dfkill) == "412518D2"            # = 10.31856; BACHLO uniform branch (no log)
    end

    @testset "DFBMRT per-record TAMORT + WK2 (bit-exact vs g16 oracle)" begin
        # The 8 iet01 Douglas-fir records at cycle 2 (IND1 order): (DBHhex, PROBhex, TAMORThex|nothing).
        recs = (("3DDB85ED","424CABC9",nothing), ("40A12BB6","41D2B8B1",nothing),
                ("401B2F71","41CF0B9F",nothing), ("405E1D1E","41DB99B5",nothing),
                ("4147E132","40E3E549","4031CDC0"), ("416A2885","408CA304","40008AB8"),
                ("41414764","40CF9482","401C9AF1"), ("4143886C","410157C8","40456FE1"))
        df_dbh = Float32[reinterpret(Float32, parse(UInt32, r[1]; base=16)) for r in recs]
        df_tpa = Float32[reinterpret(Float32, parse(UInt32, r[2]; base=16)) for r in recs]
        dfkill = reinterpret(Float32, 0x412518D2)
        badf9  = reinterpret(Float32, 0x41B78CB8)
        n = length(recs)
        old_tpa = copy(df_tpa)
        t = (tpa = copy(df_tpa),)          # WK2_before = 0 (t.tpa == old_tpa) ⇒ WK2_after = TAMORT
        _F.dfb_mrt!(t, old_tpa, collect(1:n), df_dbh, df_tpa, dfkill, badf9, false, 0.0f0)
        # Compare the SURVIVING TPA the engine stores (old_tpa − TAMORT). Recovering WK2 as
        # old_tpa − t.tpa would lose a ULP (Float32 catastrophic cancellation), so assert the
        # stored t.tpa against old_tpa − TAMORT_gold — how DFBMRT applies WK2 to PROB.
        for k in 1:n
            if recs[k][3] === nothing
                @test t.tpa[k] == old_tpa[k]   # DBH < 9″ (DCLAS < 5): untouched
            else
                tamort = reinterpret(Float32, parse(UInt32, recs[k][3]; base=16))
                @test t.tpa[k] == old_tpa[k] - tamort   # WK2 == TAMORT, bit-exact
            end
        end
    end

    # -------------------------------------------------------------------------
    # End-to-end DFB engine seam (InlandEmpire, IDFSPC=3): a scheduled outbreak
    # kills large Douglas-fir (.sum-visible extra MORT at the outbreak cycle),
    # while an out-of-range schedule leaves the projection byte-identical.
    # (The absolute .sum-DELTA vs the oracle is CORNERED by the documented IE
    #  cycle-1 growth/mortality straddle — DFKILL scales with the stand's own
    #  cycle-start BADF9/BA9 — so this asserts the seam's behaviour, not a match.)
    # -------------------------------------------------------------------------
    vie = FVSjl.InlandEmpire()
    ie_ctrl = joinpath(dir, "ie_ctrl.key")
    ie_fire = joinpath(dir, "ie_fire.key")
    ie_never = joinpath(dir, "ie_never.key")
    write(ie_ctrl,  _dfb_head("IE CONTROL ") * "ECHOSUM\nPROCESS\nSTOP\n")
    write(ie_fire,  _dfb_head("IE DFB FIRE") * _DFB_FIRE  * "ECHOSUM\nPROCESS\nSTOP\n")
    write(ie_never, _dfb_head("IE DFB NEVR") * _DFB_NEVER * "ECHOSUM\nPROCESS\nSTOP\n")
    for k in ("ie_ctrl", "ie_fire", "ie_never")
        cp(joinpath(dir, "shared.tre"), joinpath(dir, "$k.tre"))
    end
    ierows(key) = filter(l -> !startswith(l, "-999"),
                         split(strip(FVSjl.run_keyfile(key; variant = vie, output = :sum)), '\n'))

    @testset "DFB seam INERT when no outbreak is scheduled in range" begin
        @test ierows(ie_ctrl) == ierows(ie_never)   # MANSCHED 99 never fires ⇒ byte-identical
    end

    @testset "DFB outbreak fires — extra Douglas-fir mortality is .sum-visible" begin
        ctrl = ierows(ie_ctrl); fire = ierows(ie_fire)
        @test ctrl != fire                          # the outbreak changes the projection
        # cycle-1 (1990 row) is unchanged (outbreak is scheduled for cycle 2 = the 2000 row);
        # the 2000-row mortality (field 25, the MOR column) must be strictly greater with DFB on.
        mort(row) = parse(Int, split(row)[25])
        r2000_ctrl = ctrl[findfirst(r -> startswith(r, "2000"), ctrl)]
        r2000_fire = fire[findfirst(r -> startswith(r, "2000"), fire)]
        @test ctrl[1] == fire[1]                    # 1990 inventory/first-cycle row identical
        @test mort(r2000_fire) > mort(r2000_ctrl)   # DFB adds Douglas-fir mortality at cycle 2
    end

    # =========================================================================
    # ACTIVATION MODES — dump-replay BIT-EXACT vs relinked FVSie_dfb g16 oracle
    # (scratchpad/dfb/mrun; instrumented DFBGO/DFBSCH/DFBMOD/DFBMRT/DFBWIN, each
    #  instrumented .sum verified byte-identical to the clean relink first).
    # =========================================================================
    _fh(h) = reinterpret(Float32, parse(UInt32, h; base = 16))
    _fd(i) = reinterpret(Float32, UInt32(i & 0xffffffff))

    @testset "RANSTART (ISMETH=2) inclusion draw precedes BACHLO (bit-exact)" begin
        # RANSTART + MANSCHED 2 + STOPROB 0.9: the DFBGO inclusion draw consumes the FIRST 56457
        # uniform (3EE23A99 < PROTBK 0.9 ⇒ fires), so DFBMOD's BACHLO starts at the SECOND uniform
        # ⇒ DFKILL 4116B301, distinct from the MANSTART value 412518D2.
        d = _F.dfb_defaults!(); d.active = true; d.ismeth = Int32(2); d.lepi = true; d.epiprb = 0.9f0
        rand1 = _F.dfb_rand!(d)
        @test _hex(rand1) == "3EE23A99"             # inclusion draw = first uniform
        @test rand1 < d.epiprb                       # < PROTBK ⇒ stand included
        dfkill = _F.dfb_mod(d, _fh("41B78CB8"), _fh("42824F2E"), 4, 2)
        @test _hex(dfkill) == "4116B301"            # BACHLO from the SECOND uniform (order proven)
    end

    @testset "RANSCHED (IDBSCH=2) DFBSCH auto-scheduler (bit-exact schedule + stream)" begin
        # DFBSCH walks the DFBRAN stream seeded at ORSEED (55329, NOT +1128), IWAIT 10, DBEVNT 0.2,
        # IPAST 1980, a 10-cycle 1990-start run ⇒ regional outbreaks in cycles 2,3,4,6,7,8,9.
        rs = replace(_dfb_head("IE RANSCHED"), "NUMCYCLE         5.0" => "NUMCYCLE        10.0") *
             "DFB\nMANSTART\nRANSCHED          10       0.2    1980\nEND\n" *
             "ECHOSUM\nPROCESS\nSTOP\n"
        rs_key = joinpath(dir, "ie_ransched.key"); write(rs_key, rs)
        cp(joinpath(dir, "shared.tre"), joinpath(dir, "ie_ransched.tre"))
        sched = Int[]
        for s in FVSjl.each_stand(rs_key; variant = vie)
            FVSjl.setup_growth!(s)
            sched = copy(s.dfb.scheduled_cycles); break
        end
        @test sched == [2, 3, 4, 6, 7, 8, 9]
        # cycle-loop DFKILL (MANSTART, no inclusion draw) — the continuous 56457 BACHLO stream across
        # all 7 outbreaks, each on the oracle's per-cycle BADF9/BA9.
        d = _F.dfb_defaults!(); d.active = true
        stream = (("41B78CB8","42824F2E","412518D2"), ("418F5603","42F894F8","405C41FF"),
                  ("418DA1B2","431B769B","3FE122C0"), ("419FEA8A","434C315E","405D8B99"),
                  ("41C64A72","435F4FF0","40124B29"), ("41C22339","4368AA76","40480C6D"),
                  ("41AF7072","43734D1E","4001A775"))
        for (b9, ba9, gk) in stream
            @test _hex(_F.dfb_mod(d, _fh(b9), _fh(ba9), 4, 2)) == gk
        end
    end

    @testset "CUROUTBK (LINPRG, ICYC=1) in-progress outbreak (bit-exact both branches)" begin
        badf9 = _fh("41800000"); ba9 = _fh("42000000")   # DFBER stats at ICYC 1 (BADF9=16, BA9=32)
        # User-entered kill (PREKLL 30, IYOUT 2): DFKILL = 30/PERDD[2] − 30 = 30/0.8 − 30 = 7.5 (no RNG)
        dp = _F.dfb_defaults!(); dp.active = true; dp.linprg = true; dp.iyout = Int32(2); dp.prekll = 30.0f0
        @test _hex(_F.dfb_mod(dp, badf9, ba9, 4, 1)) == "40F00000"    # 7.5
        # From-treelist (PREKLL 0 ⇒ BACHLO branch): DFKILL = BACHLO·(BADF9/BA9)·4·(1−PERDD[2])
        dl = _F.dfb_defaults!(); dl.active = true; dl.linprg = true; dl.iyout = Int32(2); dl.prekll = 0.0f0
        @test _hex(_F.dfb_mod(dl, badf9, ba9, 4, 1)) == "403B88CA"    # 2.930224
        # LINPRG only fires on ICYC==1 — a later cycle takes the new-outbreak branch.
        dn = _F.dfb_defaults!(); dn.active = true; dn.linprg = true; dn.iyout = Int32(2); dn.prekll = 30.0f0
        @test _hex(_F.dfb_mod(dn, badf9, ba9, 4, 2)) != "40F00000"
        # IYOUT > 4 ⇒ outbreak assumed over ⇒ DFKILL 0.
        do_ = _F.dfb_defaults!(); do_.active = true; do_.linprg = true; do_.iyout = Int32(5); do_.prekll = 30.0f0
        @test _F.dfb_mod(do_, badf9, ba9, 4, 1) == 0.0f0
    end

    @testset "DFBINV pre-killed-DF count (dfbinv.f)" begin
        d = _F.dfb_defaults!(); d.active = true; d.linv = true
        _F.dfb_inv!(d, Float32[])            # iet01: no DFB damage codes ⇒ PREKLL unchanged (0)
        @test d.prekll == 0.0f0
        _F.dfb_inv!(d, Float32[3.0f0, 4.5f0, 2.5f0])   # Σ PROB of killed DF records
        @test d.prekll == 10.0f0
        d2 = _F.dfb_defaults!(); d2.active = true      # LINV false (user gave PREKLL) ⇒ no-op
        _F.dfb_inv!(d2, Float32[3.0f0])
        @test d2.prekll == 0.0f0
    end

    @testset "DFBWIN windthrow (dfbwin.f) — kernel + OKILL feed (bit-exact)" begin
        # iet01 cyc2 full stand (27 live records) from the g16 DFBWIN dump: I,ISPI,PCT,HT,DBH,PROB,WK2b.
        WIN = [(1,7,1119210546,1117877413,1095866914,1084716091,1053093848),
               (2,3,992607519,1082526874,1037796845,1112320969,1102856148),
               (3,5,1108081999,1108731464,1090940957,1099257460,1063429521),
               (4,2,1111305449,1118653018,1091240733,1092373894,1073282560),
               (5,2,1110452554,1116795530,1091219158,1091890631,1074073580),
               (6,10,1102348489,1113970676,1089803754,1099802463,1069448681),
               (7,2,1112240853,1094780366,1092024129,1091656165,1068281880),
               (8,7,1116317357,1116082629,1093168176,1088932695,1064544642),
               (9,3,1090223644,1107611434,1084304310,1104328881,1080096000),
               (10,10,1115853198,1118608378,1092367331,1092909410,1060967812),
               (11,3,1065181732,1100647408,1075523441,1104087967,1080141181),
               (12,3,1075623039,1102701166,1079909662,1104910773,1075730756),
               (13,10,1050871516,1093665869,1065589502,1115918648,1099112584),
               (14,4,1106145794,1108698973,1090575191,1103788637,1066268547),
               (15,3,1118691207,1117363663,1095229746,1088677193,1050562142),
               (16,4,1109643741,1111079147,1090973995,1100471752,1064440242),
               (17,3,1120403456,1117349834,1097476229,1082958596,1046349065),
               (18,7,1116802856,1116125866,1093434164,1088892147,1062009165),
               (19,3,1117350273,1115834110,1094797156,1087345794,1053009369),
               (20,7,1114888355,1116024156,1092130717,1091515150,1067464313),
               (21,4,1119859102,1117466134,1096447846,1086357050,1047992817),
               (22,3,1118048466,1116931406,1094944876,1090607048,1051132307),
               (23,10,1084041142,1108467784,1083304300,1104593652,1077544873),
               (24,8,1010687048,1091342318,1050562592,1101698042,1087550763),
               (25,8,1099945859,1109481731,1088702288,1101384742,1067844392),
               (26,8,1096104728,1108316722,1087670175,1105248895,1069920153),
               (27,4,1113946836,1109600721,1092083438,1099102310,1059924626)]
        species = Int[r[2] for r in WIN]
        pct  = Float32[_fd(r[3]) for r in WIN]; ht  = Float32[_fd(r[4]) for r in WIN]
        dbh  = Float32[_fd(r[5]) for r in WIN]; prob = Float32[_fd(r[6]) for r in WIN]
        wk2  = Float32[_fd(r[7]) for r in WIN]
        df9kil, telig = _F.dfb_win_kernel!(species, dbh, ht, prob, pct, wk2,
                                           _F._DFB_IFVSSP_IE, 80.0f0, 20.0f0, 0.8f0, 0.0f0, 3, 23)
        @test _hex(telig)  == "41F6CA95"           # TELIG 30.848917
        @test _hex(df9kil) == "412E5589"           # DF9KIL = OKILL 10.895883
        # eligible (PCT≥80, HT>20) WK2 after windthrow: DF recs 15/17/22, GF 21, LP 1.
        gold = Dict(1 => 1071438053, 15 => 1081956374, 17 => 1075708295,
                    21 => 1086232353, 22 => 1083113028)
        for (I, gw) in gold
            k = findfirst(r -> r[1] == I, WIN)
            @test wk2[k] == _fd(gw)
        end
        # OKILL feeds DFBMOD: DFKILL = BACHLO·(BADF9/BA9)·4 + OKILL = 10.31856 + 10.895883 = 21.214.
        dm = _F.dfb_defaults!(); dm.active = true; dm.okill = df9kil
        @test _hex(_F.dfb_mod(dm, _fh("41B78CB8"), _fh("42824F2E"), 4, 2)) == "41A9B72E"
    end

    @testset "DFBMRT windthrow-add path (OKILL>0, bit-exact)" begin
        # The 8 iet01 cyc2 Douglas-fir records: (DBH, PROB, incoming-WK2, final-WK2). Eligible windthrow
        # WK2 for 15/17/22 (hex), background for 2/9/11/12/19 (decimal). DFKILL 21.214 (with OKILL).
        recs = (("3DDB85ED","424CABC9",1102856148,"41BC3FD4"), ("40A12BB6","41D2B8B1",1080096000,"4060F500"),
                ("401B2F71","41CF0B9F",1080141181,"4061A57D"), ("405E1D1E","41DB99B5",1075730756,"401E5944"),
                ("4147E132","40E3E549",0x407D5816,"40E3E549"), ("416A2885","408CA304",0x401E0187,"408CA304"),
                ("41414764","40CF9482",1053009369,"40AD36DB"), ("4143886C","410157C8",0x408EFE44,"410157C8"))
        n = length(recs)
        df_dbh = Float32[_fh(r[1]) for r in recs]; df_tpa = Float32[_fh(r[2]) for r in recs]
        wk2in  = Float32[_fd(r[3]) for r in recs]
        old_tpa = copy(df_tpa)
        t = (tpa = Float32[old_tpa[i] - wk2in[i] for i in 1:n],)
        _F.dfb_mrt!(t, old_tpa, collect(1:n), df_dbh, df_tpa, _fh("41A9B72E"), _fh("41B78CB8"), false, 10.895883f0)
        for i in 1:n
            @test t.tpa[i] == old_tpa[i] - _fh(recs[i][4])   # surviving TPA = PROB − WK2(windthrow+DFB)
        end
    end

    @testset "DFBWIN/RANSCHED seams INERT when nothing is scheduled in range" begin
        # WINDTHR 99 (out of the 5-cycle run) ⇒ no windthrow; MANSCHED 99 disables the default
        # RANSCHED auto-scheduler (IDBSCH defaults to 2) and its outbreak is out of range ⇒ inert.
        win_never = _dfb_head("IE WIN NEVER") *
                    "DFB\nMANSTART\nMANSCHED          99\nWINDTHR           99\nEND\n" *
                    "ECHOSUM\nPROCESS\nSTOP\n"
        wn_key = joinpath(dir, "ie_winnever.key"); write(wn_key, win_never)
        cp(joinpath(dir, "shared.tre"), joinpath(dir, "ie_winnever.tre"))
        @test ierows(ie_ctrl) == ierows(wn_key)
        # RANSCHED with DBEVNT 0 ⇒ no outbreak ever scheduled ⇒ byte-identical to control.
        rs_never = _dfb_head("IE RSCH NONE") *
                   "DFB\nMANSTART\nRANSCHED          10       0.0    1980\nEND\n" * "ECHOSUM\nPROCESS\nSTOP\n"
        rn_key = joinpath(dir, "ie_rschnone.key"); write(rn_key, rs_never)
        cp(joinpath(dir, "shared.tre"), joinpath(dir, "ie_rschnone.tre"))
        @test ierows(ie_ctrl) == ierows(rn_key)
    end
end
