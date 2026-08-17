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
        dfkill = _F.dfb_mod(d, badf9, ba9, 4)       # NUMYRS = min(ILENTH 4, IFINT 10)
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
end
