# Douglas-fir Tussock Moth (DFTM) defoliator model — chunk 0 (keyword reader +
# the TMRANN RNG + the deterministic TMOTPR outbreak probability), plus the
# inert-seam guarantee.
#
# DOCTRINE: DFTM ships nowhere (every FVS*_buildDir links the base/exdftm.f
# no-op stub), so a numeric oracle is a RELINK (dftm/*.o swapped for exdftm.o —
# recipe in scratchpad/dftm/build_ie_dftm.sh).  For this chunk the validated
# numeric surface is small and self-contained, so the goldens come from tiny
# standalone gfortran-16 drivers over the PRISTINE Fortran / its exact formulas
# (scratchpad/dftm/driver_tmrann.f, driver_tmotpr.f):
#   * TMRANN  — the full pristine tmrann.f, seed 55329 → 6-draw stream (bit-exact).
#   * TMOTPR  — the three IPRBMT logistic formulas verbatim, fixed scalar inputs.
#     These carry transcendentals (exp/cos/sin/sqrt/log); Julia's Float32 libm and
#     gfortran's may differ by ≤1 ULP, so TMOTPR is asserted to ≈ (a few ULP), the
#     documented straddle class — the ALGEBRA is a faithful transcription.  Full
#     routine-level dump-replay (RELDEN/RELDSP/TPROB plumbing) is a handoff item.
# The keyword reader is validated by parse-correctness + the INERT-seam A/B
# (a DFTM-present stand projects byte-identically to one without — no engine seam).

using Test
using FVSjl
const _F = FVSjl

_hex(x::Float32) = uppercase(string(reinterpret(UInt32, x); base = 16, pad = 8))
_fromhex(h) = reinterpret(Float32, parse(UInt32, h; base = 16))
# ULP distance between two Float32 (same sign, finite)
_ulps(a::Float32, b::Float32) = abs(Int64(reinterpret(Int32, a)) - Int64(reinterpret(Int32, b)))

# --- goldens from scratchpad/dftm drivers over pristine dftm/tmrann.f + tmotpr formulas ---
const _G_TMRANN = ("3EDDB57A", "3F5AB21B", "3F630A07", "3F2734C9", "3EF4FB2C", "3F4AF646")
const _G_OTPR_M1 = "3F2DE146"
const _G_OTPR_M2 = "3DDC07CC"
const _G_OTPR_M3 = "3DF4118C"

const _DFTM_TRE = """
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
_dftm_head(title) = """
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
# DFTM block exercising a spread of sub-keywords across the TMCOM1/UPPER/LOWER commons.
const _DFTM_BLOCK = """
DFTM
REPORT             1
NODFRUN
MANSTART
RANSCHED          25       0.2      1970
PROBMETH           2       0.9
TOPO               3
ASHDEPTH        12.5
NUMCLASS          15        10       0.4
WEIGHT           2.0       3.0
BIOMASS            2
NPV2
REDIST           0.7
SALVAGE         40.0
END
"""
# DFTM block that only names the extension (all defaults) — for the inert-seam A/B.
const _DFTM_MIN = """
DFTM
END
"""

@testset "Douglas-fir Tussock Moth (DFTM) — chunk 0" begin

    # -------------------------------------------------------------------------
    # TMRANN (tmrann.f) — the model's own Lehmer/MINSTD LCG, seed 55329.
    # -------------------------------------------------------------------------
    @testset "TMRANN Lehmer stream (bit-exact vs pristine tmrann.f, seed 55329)" begin
        d = _F.dftm_defaults!(_F.InlandEmpire())
        @test d.rng_s0 == 55329.0 && d.rng_ss == 55329.0f0
        for g in _G_TMRANN
            @test _hex(_F.dftm_rand!(d)) == g
        end
    end

    @testset "TMRNSD reseed (odd forcing + default reset)" begin
        d = _F.dftm_defaults!(_F.InlandEmpire())
        _F.dftm_seed!(d, 100.0f0, true)        # even → forced odd 101
        @test d.rng_ss == 101.0f0 && d.rng_s0 == 101.0
        _F.dftm_seed!(d, 55.0f0, true)         # already odd
        @test d.rng_ss == 55.0f0 && d.rng_s0 == 55.0
        _F.dftm_rand!(d)                       # advance S0
        @test d.rng_s0 != 55.0
        _F.dftm_seed!(d, 0.0f0, false)         # LSET=false: reset S0 to SS (55)
        @test d.rng_s0 == 55.0
    end

    # -------------------------------------------------------------------------
    # TMOTPR (tmotpr.f) — deterministic stand-outbreak probability, 3 methods.
    # Transcendental straddle: assert ≈ within a few ULP of the gfortran golden.
    # -------------------------------------------------------------------------
    @testset "TMOTPR outbreak probability (≈ gfortran formula golden, all 3 methods)" begin
        # method 1 (Heller): elev 35 (hundred-ft), slope 0.30, aspect 1.2 rad, topo 2,
        # relden 120, reldsp3 40, reldsp4 25, tprob 200.
        p1 = _F.dftm_otpr(1, 2.0f0, 15.93f0; elev = 35.0f0, slope = 0.30f0, aspect = 1.2f0,
                          relden = 120.0f0, reldsp3 = 40.0f0, reldsp4 = 25.0f0, tprob = 200.0f0,
                          has_df = true)
        @test _ulps(p1, _fromhex(_G_OTPR_M1)) <= 4
        # method 2 (Mika-Moore + ash): topo 3, ash 15.93, pgfba 0.4, ba 180.
        p2 = _F.dftm_otpr(2, 3.0f0, 15.93f0; pgfba = 0.4f0, ba = 180.0f0)
        @test _ulps(p2, _fromhex(_G_OTPR_M2)) <= 4
        # method 3 (Mika-Moore no ash): topo 3, pgfba 0.4, ba 180.
        p3 = _F.dftm_otpr(3, 3.0f0, 15.93f0; pgfba = 0.4f0, ba = 180.0f0)
        @test _ulps(p3, _fromhex(_G_OTPR_M3)) <= 4
        # no-DF guard (method 1 returns 0)
        @test _F.dftm_otpr(1, 2.0f0, 15.93f0; has_df = false) == 0.0f0
    end

    # -------------------------------------------------------------------------
    # DFTM chunk 1 — LIVE g16 dump-replay vs the relinked FVSie_dftm oracle on a
    # DENSE DF/GF host stand (10 DF + 8 GF, DFTMGO L=T, INSCYC forced 5-yr cycle).
    # Goldens = Float32 hex from stderr TRANSFER dumps added to the pristine
    # dftm/{tmotpr,dftmgo,tmbmas}.o for the relink (the instrumented .sum is
    # byte-identical to the clean relink; scratchpad/dftm/ovr/). These UPGRADE the
    # TMOTPR formula golden to a live oracle golden and add the DFTMGO gate + the
    # TMBMAS method-2 biomass goldens. All DETERMINISTIC (no TMRANN draw).
    # -------------------------------------------------------------------------
    @testset "TMOTPR method 1 — LIVE PROTBK dump-replay (FVSie_dftm, bit-exact)" begin
        # Two DFTMGO/TMCOUP calls on the dense stand; identical topo/tmashd/elev/
        # slope/aspect, differing BA/RELDEN/RELDSP/TPROB. IN: elev slope aspect
        # topo tmashd ba relden reldsp3 reldsp4 tprob ; OUT PROTBK.
        for (inhex, gold) in (
            (("42080000","3E99999A","40AFEDE4","3F800000","417EE148",
              "42900001","42A6272F","422F4FD0","421CFE8E","42C2FBF0"), "3F464F97"),
            (("42080000","3E99999A","40AFEDE4","3F800000","417EE148",
              "42AB60CA","42BA5334","42419DA3","423308C6","42B57FF4"), "3F4DB2B7"),
        )
            elev, slope, aspect, topo, tmashd, ba, relden, reldsp3, reldsp4, tprob =
                _fromhex.(inhex)
            p = _F.dftm_otpr(1, topo, tmashd; elev = elev, slope = slope, aspect = aspect,
                             relden = relden, reldsp3 = reldsp3, reldsp4 = reldsp4,
                             tprob = tprob, ba = ba, has_df = true)
            @test _hex(p) == gold          # bit-exact vs live oracle PROTBK
        end
    end

    @testset "DFTMGO host-threshold gate — CNTDF/CNTGF dump-replay (bit-exact)" begin
        # Per-record host PROB (trees/acre) in IND1 order, from DBGGODF/DBGGOGF.
        dfp = _fromhex.(("4093DA38","4056CD2A","402E9B5F","40DFFF5C","410A80DF",
                         "40B75025","408210B3","404868AC","411D11A4","401A4C03"))
        gfp = _fromhex.(("40BD0B95","4077AC52","4101CD8C","40416294","409C4546",
                         "40E17C01","4065EB69","40335832"))
        g = _F.dftm_go_gate(collect(dfp), collect(gfp); ldf = true, lgf = true, nclas = (20, 20))
        @test _hex(g.cntdf) == "424E0DC2"   # Σ DF PROB = 51.5134, serial add
        @test _hex(g.cntgf) == "421CF227"   # Σ GF PROB = 39.2365
        @test g.idf == 10 && g.igf == 8
        @test g.naclas == (10, 8)           # MIN(count, NCLAS)
        @test g.ldf && g.lgf && g.l         # both hosts go
        # NODFRUN drops DF; a no-GF stand with DF-only still goes on DF.
        gdf = _F.dftm_go_gate(Float32[], collect(dfp); ldf = false, lgf = true, nclas = (20, 20))
        @test !gdf.ldf && gdf.lgf && gdf.l && gdf.naclas == (0, 10)
        # below-threshold host (Σ<0.01) drops out; empty ⇒ no go.
        @test _F.dftm_go_gate(Float32[], Float32[]).l == false
    end

    @testset "TMBMAS method 2 — FBIOMS/PCNEWF dump-replay (bit-exact / ≤1 ULP)" begin
        # DBGBMAS2 rows: IS ISP | slope aspect ba tprob relden dbh ht dgi cr pct | FBIOMS PCNEWF.
        rows = split(strip("""
        1 3 3E99999A 40AFEDE4 42AB60CA 42B57FF4 42BA5334 414FFE05 4291D1B4 3F5DD840 3EA3D70A 421F878D 4364E18B 41F2DB66
        2 3 3E99999A 40AFEDE4 42AB60CA 42B57FF4 42BA5334 416CB1A5 429994A9 3F301710 3EA8F5C3 428689B2 438276F4 41EADF32
        3 3 3E99999A 40AFEDE4 42AB60CA 42B57FF4 42BA5334 418BB9F6 42A5B651 3FA2AB30 3EA8F5C3 42B188B2 439970E9 420145DD
        4 3 3E99999A 40AFEDE4 42AB60CA 42B57FF4 42BA5334 413CC8C5 428D1ECB 3FC7A5C8 3EAE147B 41B60A30 43510607 420A0550
        5 3 3E99999A 40AFEDE4 42AB60CA 42B57FF4 42BA5334 412E2166 42859F16 3FD0FC24 3EB33333 41846AE0 4345A521 420DBEEB
        6 3 3E99999A 40AFEDE4 42AB60CA 42B57FF4 42BA5334 41491916 429144A2 3FAE1468 3EA8F5C3 420BA686 435CC33F 42070412
        7 3 3E99999A 40AFEDE4 42AB60CA 42B57FF4 42BA5334 41654016 4298D001 3F936470 3EAE147B 427A600F 437C4945 4203414E
        8 3 3E99999A 40AFEDE4 42AB60CA 42B57FF4 42BA5334 41859D62 42A316EE 3FBCDD38 3EAE147B 42A6F071 4392D15A 4208E148
        9 3 3E99999A 40AFEDE4 42AB60CA 42B57FF4 42BA5334 410C7930 426FF218 3F2D07E8 3EA3D70A 409A1FC4 43237A47 41D14398
        10 3 3E99999A 40AFEDE4 42AB60CA 42B57FF4 42BA5334 4192FC29 42A8E5B8 3F986248 3EA8F5C3 42BBE4F3 43A1BDCB 41F9DE38
        11 4 3E99999A 40AFEDE4 42AB60CA 42B57FF4 42BA5334 415712B7 42903299 400F01F0 3EB851EC 423AB501 43C80000 41700000
        12 4 3E99999A 40AFEDE4 42AB60CA 42B57FF4 42BA5334 415CB289 4293F527 3F39E3D0 3EAE147B 424D740E 43C17A71 41700000
        13 4 3E99999A 40AFEDE4 42AB60CA 42B57FF4 42BA5334 4120B5E6 427960C2 3F74A470 3EA8F5C3 41206AF2 438372C3 41700000
        14 4 3E99999A 40AFEDE4 42AB60CA 42B57FF4 42BA5334 4180BA9C 429F7867 3F7F9560 3EAE147B 429BD1C4 43C80000 418352E6
        15 4 3E99999A 40AFEDE4 42AB60CA 42B57FF4 42BA5334 415CCB84 42937F0A 3FD2C780 3EB851EC 42652161 43C80000 41700000
        16 4 3E99999A 40AFEDE4 42AB60CA 42B57FF4 42BA5334 413F17FB 4288877B 3FE39B10 3EB33333 41E93899 43BCD74A 41700000
        17 4 3E99999A 40AFEDE4 42AB60CA 42B57FF4 42BA5334 417BBF1A 429CCB54 3FCB1AB8 3EB33333 4291DC13 43C80000 418411E4
        18 4 3E99999A 40AFEDE4 42AB60CA 42B57FF4 42BA5334 41935F28 42A6F58E 400DCD04 3EB33333 42C80000 43C80000 41B7B81C
        """), '\n')
        maxfb = 0; maxpn = 0
        for row in rows
            t = split(strip(row))
            isp = parse(Int, t[2])
            slope, aspect, ba, tprob, relden, dbh, ht, dgi, cr, pct = _fromhex.(Tuple(t[3:12]))
            gfb, gpn = _fromhex(t[13]), _fromhex(t[14])
            fb, pn = isp == 3 ?
                _F.dftm_bmas2_df(slope, aspect, ba, tprob, relden, dbh, ht, dgi, cr, pct) :
                _F.dftm_bmas2_gf(slope, aspect, ba, tprob, relden, dbh, ht, dgi, cr, pct)
            maxfb = max(maxfb, _ulps(fb, gfb)); maxpn = max(maxpn, _ulps(pn, gpn))
        end
        # DF FBIOMS + all PCNEWF bit-exact; GF FBIOMS `exp` carries ≤1-ULP straddle.
        @test maxfb <= 1
        @test maxpn == 0
    end

    # -------------------------------------------------------------------------
    # TMINIT (tminit{,ec,em,so,tt}.f) defaults + variant IGFCOD crosswalk.
    # -------------------------------------------------------------------------
    @testset "TMINIT defaults + IGFCOD crosswalk" begin
        d = _F.dftm_defaults!(_F.InlandEmpire())
        @test d.idfcod == 3 && d.igfcod == 4          # base/IE
        @test d.itmeth == 1 && d.itmsch == 1 && d.iprbmt == 1
        @test d.ibmtyp == 4 && d.nclas == (Int32(20), Int32(20))
        @test d.tmevnt == 0.1f0 && d.tmwait == 30 && d.tmpast == 1492
        @test d.topo == 1.0f0 && d.tmashd == 15.93f0 && d.tmdefl == 50.0f0
        @test d.prbscl == 1.0f0 && d.tmpn1 == 0.5f0
        @test d.dfegg == (11.0f0, 9.0f0, 7.0f0) && d.gfegg == (15.0f0, 10.0f0, 7.0f0)
        @test d.dffbio == (213.8f0, 64.2f0) && d.gffbio == (227.0f0, 63.7f0)
        @test length(d.b0) == 66 && length(d.r0) == 18 && length(d.b1) == 36
        @test d.b0[52] == 0.80f0 && d.b0[62] == 0.25f0   # B0 spot checks (DATA order)
        @test d.b1[26] == 52.0626838f0                    # B1(26) spot check
        # variant IGFCOD differences
        @test _F.dftm_igfcod(_F.EastCascades()) == 6
        @test _F.dftm_igfcod(_F.Teton()) == 9
        @test _F.dftm_igfcod(_F.Utah()) == 4              # base default
    end

    # -------------------------------------------------------------------------
    # Keyword reader (dftmin.f) via the real keyword dispatch.
    # -------------------------------------------------------------------------
    dir = mktempdir()
    write(joinpath(dir, "shared.tre"), _DFTM_TRE)
    dftm_key = joinpath(dir, "dftm.key")
    write(dftm_key, _dftm_head("DFTM PARSED") * _DFTM_BLOCK * "ECHOSUM\nPROCESS\nSTOP\n")
    cp(joinpath(dir, "shared.tre"), joinpath(dir, "dftm.tre"))
    v = FVSjl.InlandEmpire()

    @testset "DFTM keyword reader populates s.dftm (dftmin.f)" begin
        d = nothing
        for s in FVSjl.each_stand(dftm_key; variant = v)
            d = s.dftm
            break
        end
        @test d !== nothing
        @test d.active == true
        @test d.itmrep == 1                # REPORT 1
        @test d.ldf == false               # NODFRUN
        @test d.lgf == true
        @test d.itmeth == 1                # MANSTART
        @test d.itmsch == 2                # RANSCHED
        @test d.tmwait == 25 && d.tmevnt == 0.2f0 && d.tmpast == 1970
        @test d.iprbmt == 2 && d.prbscl == 0.9f0   # PROBMETH 2 0.9
        @test d.topo == 3.0f0              # TOPO 3
        @test d.tmashd == 12.5f0           # ASHDEPTH
        @test d.nclas == (Int32(15), Int32(10)) && d.tmpn1 == 0.4f0
        @test d.weight == (2.0f0, 3.0f0)
        @test d.ibmtyp == 2                # BIOMASS 2
        @test d.lnpv2 == true && d.b0[9] == 0.036f0 && d.b0[12] == 0.072f0
        @test d.b0[62] == 0.7f0            # REDIST 0.7
        @test d.itmslv == 1 && d.tmdefl == 40.0f0  # SALVAGE 40
    end

    # -------------------------------------------------------------------------
    # INERT seam: a DFTM-present stand must project byte-identically to a stand
    # with no DFTM keyword (no per-cycle engine seam is wired yet).
    # -------------------------------------------------------------------------
    ctrl_key = joinpath(dir, "ctrl.key")
    min_key  = joinpath(dir, "min.key")
    write(ctrl_key, _dftm_head("DFTM CTRL  ") * "ECHOSUM\nPROCESS\nSTOP\n")
    write(min_key,  _dftm_head("DFTM MIN   ") * _DFTM_MIN * "ECHOSUM\nPROCESS\nSTOP\n")
    write(dftm_key, _dftm_head("DFTM FULL  ") * _DFTM_BLOCK * "ECHOSUM\nPROCESS\nSTOP\n")
    for k in ("ctrl", "min", "dftm")
        cp(joinpath(dir, "shared.tre"), joinpath(dir, "$k.tre"); force = true)
    end
    rows(key) = filter(l -> !startswith(l, "-999"),
                       split(strip(FVSjl.run_keyfile(key; variant = v, output = :sum)), '\n'))

    @testset "DFTM seam is INERT (.sum byte-identical with/without a DFTM block)" begin
        base = rows(ctrl_key)
        @test rows(min_key)  == base       # DFTM/END only ⇒ byte-identical
        @test rows(dftm_key) == base       # full DFTM block ⇒ still byte-identical (inert seam)
    end
end
