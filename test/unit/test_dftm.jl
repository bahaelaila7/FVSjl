# Douglas-fir Tussock Moth (DFTM) defoliator model — keyword reader + the TMRANN
# RNG + the deterministic TMOTPR outbreak probability + DFTMGO gate + TMBMAS
# biomass (chunks 0-1), PLUS the DFTMOD population integrator (→ DPCENT) and the
# TMCOUP per-tree damage functions, both LIVE dump-replayed bit-exact against the
# instrumented FVSie_dftm oracle on the dense DF/GF host stand.  Inert-seam kept.
#
# VALIDATED HERE (dump-replay, 0 ULP): DFTMOD DPCENT (17 non-empty classes + the
# empty-GF-class NaN reset); TMCOUP T/ITAB, WK2, DG, HTG (no-top-kill), top-kill
# class-1 leader-loss + class-2 crown PCKILL/HTGLOS; GARBEL/GRCLAS classification
# (the RDPSRT-sorted pointer + ISC sector pointers incl. the empty class 11,9 and
# the cross-block underflow 10,11 + Z4/Z2/Z3, all 18 classes); and the full
# TMRANN-stream ordering (TMBCHL vs pristine driver; the RANLARVA egg X7 in JCLAS2
# order; the DO-380 per-tree PRTOPK + K≥2 RANDOM).  STILL engine-inert (no
# simulate.jl seam wires DFTM): the INSCYC hook + the gated coupling seam remain.
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

@testset "Douglas-fir Tussock Moth (DFTM) — integrator + coupling" begin

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
    # DFTMOD (dftmod.f) integrator → DPCENT — LIVE dump-replay, BIT-EXACT.
    # Goldens = instrumented FVSie_dftm (dftmod.f DBGIN/DBGDP hex, dense.key).
    # Per-class ICON inputs (IZ6,Z2,Z3,X5,X6,X7,Z4,Z5) → DPCENT1/DPCENT2.
    # Empty GF class 11 (PROB 0 ⇒ Z3=0) drives the X(1)=NaN "RESET TO 0.0"
    # guard; DPCENT there = NaN and must reproduce bit-exact.
    # -------------------------------------------------------------------------
    @testset "DFTMOD integrator DPCENT (dftmod.f) — dump-replay bit-exact" begin
        # I IZ6  Z2 Z3 X5 X6 X7 Z4 Z5  (hex)
        din = split(strip("""
        1 1 4208E148 4392D15A 42C8F6D4 4341274A 411037CC 404868AC 3F800000
        2 1 420145DD 439970E9 42C65B93 434FB408 40BE097C 402E9B5F 3F800000
        3 1 41F9DE39 43A1BDCB 42CA11ED 435E72A0 40E14C34 401A4C03 3F800000
        4 1 4203414E 437C4945 42A591BB 43298068 412B422C 408210B3 3F800000
        5 1 420DBEEB 4345A521 428C13A5 42FF369D 411C728E 410A80DF 3F800000
        6 1 420A0550 43510607 42903F78 4308E64B 41081FCC 40DFFF5C 3F800000
        7 1 42070412 435CC33F 42950849 43123F1A 4116C770 40B75025 3F800000
        8 1 41EADF32 438276F4 4299365F 433852B8 412579AE 4056CD2A 3F800000
        9 1 41F2DB67 4364E18B 428AF6A3 431F663A 4160FDDB 4093DA38 3F800000
        10 1 41D14398 43237A47 422B0CCA 42F16E29 40C7D5C6 411D11A4 3F800000
        11 2 00000000 00000000 00000000 00000000 410B5200 00000000 3F800000
        12 2 41CB9768 43580048 425BE13F 432107F8 415F715D 4149E7B0 3F800000
        13 2 418411E4 43C80000 428411E4 43A6FB87 417C2A91 4065EB69 3F800000
        14 2 418352E6 43C80000 428352E6 43A72B46 412A6EF0 40416294 3F800000
        15 2 41700000 43C7FFFF 426FFFFF 43A9FFFF 41116C52 412CA86E 3F800000
        16 2 41700000 43C17A71 42682C88 43A474E0 40C8B8BC 4077AC52 3F800000
        17 2 41700000 43BCD74A 42629BF2 43A083CC 41514E52 40E17C01 3F800000
        18 2 41700000 438372C3 421DBCEA 435F764C 4147426F 4101CD8C 3F800000
        """), '\n')
        # I DPCENT1 DPCENT2 (hex); class 11 = NaN (FFC00000)
        ddp = split(strip("""
        1 41E91A00 42ADC58D
        2 41A76DB5 4297BB7F
        3 41B18900 4295CAE6
        4 42251962 42C77519
        5 424D1F4E 42C80000
        6 42262DC9 42C80000
        7 422B9C93 42C80000
        8 421DFC30 42BD4838
        9 427ADA59 42C80000
        10 423CF3D2 42C80000
        11 FFC00000 FFC00000
        12 427F2043 429B5B65
        13 4218E48D 4281375C
        14 41DC0513 424F85A1
        15 41C1A1EC 423A67B8
        16 4190320A 422697A6
        17 420B2F51 426F7E89
        18 422C5937 424D16F3
        """), '\n')
        iz6 = Int[]; z2 = Float32[]; z3 = Float32[]; x5 = Float32[]
        x6 = Float32[]; x7 = Float32[]; z4 = Float32[]; z5 = Float32[]
        for l in din
            t = split(strip(l))
            push!(iz6, parse(Int, t[2])); push!(z2, _fromhex(t[3])); push!(z3, _fromhex(t[4]))
            push!(x5, _fromhex(t[5])); push!(x6, _fromhex(t[6])); push!(x7, _fromhex(t[7]))
            push!(z4, _fromhex(t[8])); push!(z5, _fromhex(t[9]))
        end
        gd1 = Float32[]; gd2 = Float32[]
        for l in ddp
            t = split(strip(l)); push!(gd1, _fromhex(t[2])); push!(gd2, _fromhex(t[3]))
        end
        g = _F.dftm_integ(iz6, z2, z3, x5, x6, x7, z4, z5,
                          _F.DFTM_B0_DEFAULT, _F.DFTM_R0_DEFAULT, _F.DFTM_B1_DEFAULT)
        _F.dftmod!(g)
        # NaN reset guard: empty GF class 11 → DPCENT = NaN, exactly as the oracle.
        @test isnan(g.dpcent[11, 1]) && isnan(g.dpcent[11, 2])
        maxu = 0
        for i in 1:18
            i == 11 && continue
            maxu = max(maxu, _ulps(g.dpcent[i, 1], gd1[i]), _ulps(g.dpcent[i, 2], gd2[i]))
        end
        @test maxu == 0                      # DPCENT bit-exact, all 17 non-empty classes
    end

    # -------------------------------------------------------------------------
    # TMCOUP (tmcoup.f) per-tree damage — LIVE dump-replay, BIT-EXACT.
    # Goldens = instrumented FVSie_dftm (tmcoup.f DBGMORT/DBGPRE/DBGPOST/DBGTK
    # hex, dense.key). Deterministic damage: T + ITAB (tree defoliation class),
    # WK2 (mortality), DG (diameter-growth loss), HTG (no-top-kill height loss),
    # and the top-kill class-1 (leader) + class-2 (crown-truncation) HTG paths
    # (fed the oracle's exact PRTOPK/RANDOM draws — the TMRANN stream itself is
    # never FFI'd; full-stream ordering is a later engine-seam item).
    # -------------------------------------------------------------------------
    @testset "TMCOUP damage functions (tmcoup.f) — dump-replay bit-exact" begin
        # IP ITAB DPMAX PROB PRNORM WK2out G19MAX T (hex)
        mort = split(strip("""
        8   5 42ADC58D 404868AC 3CA21380 3F957706 00000000 42B614BA
        3   3 4297BB7F 402E9B5F 3CC7A420 3E6D7724 00000000 4249E68F
        10  3 4295CAE6 401A4C03 3CD2E820 3E51D803 00000000 42346F02
        7   6 42C77519 408210B3 3D0490A0 4019D756 00000000 42C3883E
        5   7 42C80000 410A80DF 3CAFD3A0 41022A5E 41880000 42C80000
        4   7 42C80000 40DFFF5C 3CB7B2C0 40D28349 41000000 42C80000
        6   7 42C80000 40B75025 3CE0E260 40AC4712 41100000 42C80000
        2   6 42BD4838 4056CD2A 3D5FB941 3FFE1130 00000000 42C1BC93
        1   7 42C80000 4093DA38 3D45B1E0 408AF3A2 41B80000 42C80000
        9   7 42C80000 411D11A4 3D99E3A0 41139D05 41200000 42C80000
        9  10 429B5B65 411D11A4 3F7096BB 41139D05 00000000 42714636
        18 10 429B5B65 40335832 3C368A40 3E1B870E 00000000 42714636
        17  8 4281375C 4065EB69 3CA3B500 3E1E1E43 00000000 40F180A5
        14  8 424F85A1 40416294 3D1BB710 3E04384F 00000000 3EEB4F02
        15  8 423A67B8 409C4546 3CA71B20 3E3789D5 00000000 3E0F594C
        11  8 423A67B8 40BD0B95 3C4EB540 3E4F33EC 00000000 3E0F594C
        12  8 422697A6 4077AC52 3D6A8C80 3E62EB71 00000000 3D3B905D
        16  8 426F7E89 40E17C01 3C9DCEC1 3E94AD25 00000000 402EA5DB
        13  8 424D16F3 4101CD8C 3D5B1690 3EDE2C8F 00000000 3ECD3E22
        """), '\n')
        # DBGPRE: IP DG HTG DPMAX FINT  (row-aligned with DBGPOST below)
        pre = split(strip("""
        8 3F30F390 407809BE 42ADC58D 40A00000
        3 3F1A8830 40678DEB 4297BB7F 40A00000
        10 3EE98C00 404D241D 4295CAE6 40A00000
        7 3F074990 4074C17E 42C77519 40A00000
        5 3F5244F0 40A2AA2B 42C80000 40A00000
        4 3F1BCCD0 408D15DA 42C80000 40A00000
        6 3F097B10 40831D82 42C80000 40A00000
        2 3F16F7F0 407A7153 42BD4838 40A00000
        1 3F07AEA0 40809BFC 42C80000 40A00000
        9 3EDAFB30 4098F5D7 42C80000 40A00000
        9 3E29EDEE 4074BC8B 429B5B65 40A00000
        18 3F59F7C0 409F2177 429B5B65 40A00000
        17 3F477BB0 409EB718 4281375C 40A00000
        14 3F04CFA0 40893653 424F85A1 40A00000
        15 3F849640 40B2EA9A 423A67B8 40A00000
        11 3F17C850 4094809F 423A67B8 40A00000
        12 3F0A8110 408F0AA0 422697A6 40A00000
        16 3F26FCB0 409C0239 426F7E89 40A00000
        13 3EE18F20 408B46AC 424D16F3 40A00000
        """), '\n')
        # DBGPOST: IP DGout HTGout PRTOPK K
        post = split(strip("""
        8 3EC95EE7 404CEFA2 3F5296BB 0
        3 3EC84601 40446B89 3F25B9DE 0
        10 3E9756A0 4010587A 3E9B5800 2
        7 3E51F74A 4043EFFE 3F71ACDA 0
        5 3EA32B42 408221BC 3F0B199E 0
        4 3E71CD3C 4061BC90 3E7B7239 0
        6 3E555EBE 4051C8D0 3BFE8F63 1
        2 3E6A4D9A 404B09D2 3F10EF5A 0
        1 3E529424 404DC660 3EA41399 0
        9 3E29EDEE 4074BC8B 3BA367BE 1
        9 3DD56E6D 404EB708 3F4FDFD7 0
        18 3F08E23E 4086689D 3EE0FCD1 0
        17 3F1C9823 408A34E5 3EF61041 0
        14 3ED0836A 4075F32A 3EA5337D 0
        15 3F502956 40A23DAE 3F6BB11E 0
        11 3EEE4C6A 4086A960 3F393140 0
        12 3ED97397 4083206F 3EA0A271 0
        16 3F0315B8 40799D28 3C8AE2A7 1
        13 3EB1105B 4079FCF1 3F713212 0
        """), '\n')
        drows(s) = [split(strip(l)) for l in s]

        # --- dftm_tree_defol: percent-tree-defoliation T + AMORT row ITAB ---
        maxT = 0; badI = 0
        for r in drows(mort)
            itab = parse(Int, r[2]); dpmax = _fromhex(r[3]); g19 = _fromhex(r[7]); Tg = _fromhex(r[8])
            iz6 = itab > 7 ? 2 : 1
            (t, it) = _F.dftm_tree_defol(dpmax, g19, iz6)
            maxT = max(maxT, _ulps(t, Tg)); (it == itab) || (badI += 1)
        end
        @test maxT == 0 && badI == 0

        # --- dftm_mortality: DFTM WK2 = PROB·PRMORT ---
        maxW = 0
        for r in drows(mort)
            itab = parse(Int, r[2]); dpmax = _fromhex(r[3]); prob = _fromhex(r[4])
            prnorm = _fromhex(r[5]); wk2g = _fromhex(r[6])
            maxW = max(maxW, _ulps(_F.dftm_mortality(itab, dpmax, prob, prnorm * prob), wk2g))
        end
        @test maxW == 0

        # --- dftm_dgloss (DG) + dftm_htgloss_notopkill (HTG, K==0) ---
        pr = drows(pre); po = drows(post)
        itab_ip(ip, k) = begin
            ms = [m for m in drows(mort) if m[1] == ip]
            length(ms) == 1 ? parse(Int, ms[1][2]) : (k == 10 ? 7 : 10)   # IP=9 dup: row10 DF, row11 GF
        end
        maxDG = 0; maxH = 0; nH = 0
        for k in 1:length(pr)
            ip = pr[k][1]; itab = itab_ip(ip, k)
            (dgnew, _) = _F.dftm_dgloss(itab, _fromhex(pr[k][2]))
            maxDG = max(maxDG, _ulps(dgnew, _fromhex(po[k][2])))
            if parse(Int, po[k][5]) == 0     # no-top-kill height-loss path
                (htgnew, _) = _F.dftm_htgloss_notopkill(_fromhex(pr[k][4]), _fromhex(pr[k][3]), _fromhex(pr[k][5]))
                maxH = max(maxH, _ulps(htgnew, _fromhex(po[k][3]))); nH += 1
            end
        end
        @test maxDG == 0
        @test maxH == 0 && nH == 15

        # --- top-kill class-1 (leader): HTGLOS = HTG/FINT ---
        maxTK1 = 0; nTK1 = 0
        for k in 1:length(pr)
            parse(Int, po[k][5]) == 1 || continue
            htgpre = _fromhex(pr[k][3]); fint = _fromhex(pr[k][5])
            htgnew = htgpre - htgpre * (1.0f0 / fint)
            maxTK1 = max(maxTK1, _ulps(htgnew, _fromhex(po[k][3]))); nTK1 += 1
        end
        @test maxTK1 == 0 && nTK1 == 3

        # --- top-kill class-2 (crown truncation): PCKILL/HTGLOS from oracle RANDOM ---
        # DBGTK: IP=10 K=2 RANDOM=3EAE860F PCKILL=3D0B9E73 CROWN=41DEF1C8 HTGLOS=3F732E8A
        rnd = _fromhex("3EAE860F"); crown = _fromhex("41DEF1C8")
        pckill = _F.DFTM_TKBOT[2] + rnd * (_F.DFTM_TKTOP[2] - _F.DFTM_TKBOT[2])
        @test _ulps(pckill, _fromhex("3D0B9E73")) == 0
        @test _ulps(crown * pckill, _fromhex("3F732E8A")) == 0
    end

    # -------------------------------------------------------------------------
    # GARBEL/GRCLAS/GRPSUM/IQRSRT classification (garbel.f/grclas.f) — LIVE
    # dump-replay, BIT-EXACT.  Goldens = instrumented FVSie_dftm (tmcoup.f DBGA_*
    # hex, dense.key): the DF then GF GARBEL calls over the global IPT pointer,
    # each on its species sub-block, producing the sorted pointer, the class
    # sector pointers ISC, and the per-class Z4=ΣPROB / Z2 (weighted PCNEWF) /
    # Z3 (weighted FBIOMS).  The GF block reproduces the EMPTY class (ISC 11,9 —
    # start>end) and the cross-block sector UNDERFLOW (ISC 10,11, reaching into
    # the DF block's sorted tail) bit-exactly via the shared global gipt.
    # -------------------------------------------------------------------------
    @testset "GARBEL/GRCLAS classification (garbel.f) — dump-replay bit-exact" begin
        # per-record PROB PCNEWF FBIOMS (record index 1..18), from DBGA_REC.
        recs = split(strip("""
        1 4093DA38 41F2DB66 4364E18B
        2 4056CD2A 41EADF32 438276F4
        3 402E9B5F 420145DD 439970E9
        4 40DFFF5C 420A0550 43510607
        5 410A80DF 420DBEEB 4345A521
        6 40B75025 42070412 435CC33F
        7 408210B3 4203414E 437C4945
        8 404868AC 4208E148 4392D15A
        9 411D11A4 41D14398 43237A47
        10 401A4C03 41F9DE38 43A1BDCB
        11 40BD0B95 41700000 43C80000
        12 4077AC52 41700000 43C17A71
        13 4101CD8C 41700000 438372C3
        14 40416294 418352E6 43C80000
        15 409C4546 41700000 43C80000
        16 40E17C01 41700000 43BCD74A
        17 4065EB69 418411E4 43C80000
        18 40335832 41B7B81C 43C80000
        """), '\n')
        prob = zeros(Float32, 18); pcnewf = zeros(Float32, 18); fbioms = zeros(Float32, 18)
        for l in recs
            t = split(strip(l)); i = parse(Int, t[1])
            prob[i] = _fromhex(t[2]); pcnewf[i] = _fromhex(t[3]); fbioms[i] = _fromhex(t[4])
        end
        # golden sorted IPT (DBGA_IPT), ISC (DBGA_ISC, global), Z4/Z2/Z3 (DBGB_ICN).
        gipt_gold = Int[8,3,10,7,5,4,6,2,1,9, 18,17,14,15,11,12,16,13]
        isc_gold = [(1,1),(2,2),(3,3),(4,4),(5,5),(6,6),(7,7),(8,8),(9,9),(10,10),
                    (11,9),(10,11),(12,12),(13,13),(14,15),(16,16),(17,17),(18,18)]
        zg = split(strip("""
        1 404868AC 4208E148 4392D15A
        2 402E9B5F 420145DD 439970E9
        3 401A4C03 41F9DE39 43A1BDCB
        4 408210B3 4203414E 437C4945
        5 410A80DF 420DBEEB 4345A521
        6 40DFFF5C 420A0550 43510607
        7 40B75025 42070412 435CC33F
        8 4056CD2A 41EADF32 438276F4
        9 4093DA38 41F2DB67 4364E18B
        10 411D11A4 41D14398 43237A47
        11 00000000 00000000 00000000
        12 4149E7B0 41CB9768 43580048
        13 4065EB69 418411E4 43C80000
        14 40416294 418352E6 43C80000
        15 412CA86E 41700000 43C7FFFF
        16 4077AC52 41700000 43C17A71
        17 40E17C01 41700000 43BCD74A
        18 4101CD8C 41700000 438372C3
        """), '\n')
        z4g = Float32[]; z2g = Float32[]; z3g = Float32[]
        for l in zg
            t = split(strip(l)); push!(z4g, _fromhex(t[2])); push!(z2g, _fromhex(t[3])); push!(z3g, _fromhex(t[4]))
        end
        # NACLAS=(10,8) ⇒ nrecs==nclas per block ⇒ method-1 identity (NCL2=0), plus the gap-bubble
        # empty/underflow artifact.  WEIGHT default (1,1); TMPN1 0.5.  IND1 = identity 1..18.
        gipt = Int32.(collect(1:18))
        z4d, z2d, z3d, isc1d, isc2d, kd = _F.dftm_garbel(gipt, 1, 10, prob, pcnewf, fbioms;
            w1 = 1.0f0, w2 = 1.0f0, nclas = 10, pn1 = 0.5f0)
        z4f, z2f, z3f, isc1f, isc2f, kf = _F.dftm_garbel(gipt, 11, 8, prob, pcnewf, fbioms;
            w1 = 1.0f0, w2 = 1.0f0, nclas = 8, pn1 = 0.5f0)
        @test kd == 0 && kf == 0
        @test Int.(gipt) == gipt_gold                 # RDPSRT-sorted pointer, both blocks
        isc_out = Vector{Tuple{Int,Int}}(undef, 18)
        for i in 1:10; isc_out[i] = (isc1d[i], isc2d[i]); end               # DF offset 0
        for i in 1:8;  isc_out[10+i] = (isc1f[i] + 10, isc2f[i] + 10); end  # GF offset ISCT(IGFCOD,1)-1=10
        @test isc_out == isc_gold                      # sector pointers incl. empty (11,9) + underflow (10,11)
        z4o = vcat(z4d[1:10], z4f[1:8]); z2o = vcat(z2d[1:10], z2f[1:8]); z3o = vcat(z3d[1:10], z3f[1:8])
        maxz = 0
        for i in 1:18
            maxz = max(maxz, _ulps(z4o[i], z4g[i]), _ulps(z2o[i], z2g[i]), _ulps(z3o[i], z3g[i]))
        end
        @test maxz == 0                                # Z4/Z2/Z3 bit-exact, all 18 classes
    end

    # -------------------------------------------------------------------------
    # Full TMRANN-stream ordering (tmcoup.f) — LIVE dump-replay, BIT-EXACT.
    # Goldens = instrumented FVSie_dftm: a per-draw tmrann.f dump (DBGRNG, 98
    # draws) + a standalone gfortran driver over pristine dftm/tmbchl.f.  The DFTM
    # stream on the dense stand is: 18 RANLARVA egg TMBCHL draws (IEGTYP=1) in
    # descending-DBH JCLAS2 order (draws 1-78, w/ rejections) BEFORE DFTMOD, then
    # the DO-380 per-tree PRTOPK (+ the K≥2 crown-truncation RANDOM) in class×IPT
    # order (draws 79-98).  Never FFI'd — Julia's own dftm_rand!/dftm_tmbchl!.
    # -------------------------------------------------------------------------
    @testset "TMRANN stream — RANLARVA eggs + per-tree PRTOPK (dump-replay bit-exact)" begin
        # (1) TMBCHL standalone vs pristine dftm/tmbchl.f driver (seed 55329):
        #     TMBCHL(9,2)×3 then TMBCHL(11,3)×3.
        d = _F.dftm_defaults!(_F.InlandEmpire())
        bchl = vcat([_F.dftm_tmbchl!(d, 9.0f0, 2.0f0) for _ in 1:3],
                    [_F.dftm_tmbchl!(d, 11.0f0, 3.0f0) for _ in 1:3])
        bg = ("40E14C34","40BE097C","411037CC","41503684","4158E341","4194BE64")
        @test all(_hex(bchl[i]) == bg[i] for i in 1:6)

        # (2) Whole-stand RANLARVA egg allocation X7 (dftm_alloc_eggs! → dftm_tmbchl!),
        #     drawn in JCLAS2 (descending class DBH) order from a fresh seed 55329.
        jclas2 = Int[3,2,1,8,4,9,7,6,5,10, 14,13,16,15,11,17,12,18]
        iz6    = Int[1,1,1,1,1,1,1,1,1,1, 2,2,2,2,2,2,2,2]
        x7gold = ("411037CC","40BE097C","40E14C34","412B422C","411C728E","41081FCC",
                  "4116C770","412579AE","4160FDDB","40C7D5C6","410B5200","415F715D",
                  "417C2A91","412A6EF0","41116C52","40C8B8BC","41514E52","4147426F")
        d2 = _F.dftm_defaults!(_F.InlandEmpire())
        x7 = _F.dftm_alloc_eggs!(d2, iz6, jclas2, 9.0f0, 2.0f0, 11.0f0, 3.0f0)
        @test all(_hex(x7[i]) == x7gold[i] for i in 1:18)   # egg counts bit-exact, all classes

        # (3) The eggs consume EXACTLY draws 1-78; the next 20 raw draws must be the
        #     DO-380 topkill stream (draws 79-98, incl. the K=2 RANDOM at draw 82).
        stream = ("3F5296BB","3F25B9DE","3E9B5800","3EAE860F","3F71ACDA","3F0B199E",
                  "3E7B7239","3BFE8F63","3F10EF5A","3EA41399","3BA367BE","3F4FDFD7",
                  "3EE0FCD1","3EF61041","3EA5337D","3F6BB11E","3F393140","3EA0A271",
                  "3C8AE2A7","3F713212")
        after = ntuple(_ -> _hex(_F.dftm_rand!(d2)), 20)
        @test after == stream                                # stream ordering + count bit-exact
    end

    # -------------------------------------------------------------------------
    # TMCOUP DRIVER SEAM (dftm_couple!) — the WHOLE coupler composition, LIVE
    # cycle-2 dump-replay BIT-EXACT.  Goldens = instrumented FVSie_dftm (dense.key):
    #   * DBGIN_TRE — the pre-DFTM per-record entry state (WK2 background / HT / HTG /
    #     DBH / DG / PCT / PROB in hex, + IMC/ICR/NORMHT/ITRUNC), captured at TMCOUP
    #     entry (gradd, after MORTS/DGDRIV/HTGF, before UPDATE);
    #   * DBGA_REC — the TMBMAS per-record FBIOMS/PCNEWF (biomass predict);
    #   * DBGC_TRE — the FINAL per-record treelist after TMCOUP (IMC/ICR/NORMHT/
    #     ITRUNC + WK2/HT/HTG/DBH/DG in hex).
    # Driving dftm_couple! with the oracle's exact cycle-2 entry state must reproduce
    # DBGC_TRE for all 18 records — proving the driver ORDER/PLUMBING (GARBEL DF→GF,
    # JCLAS2 avg-DBH RDPSRT incl. the NaN empty class, ICOND fill + RANLARVA egg
    # allocation, DFTMOD, DO-320 mortality, DO-430 top-kill/DG-loss) AND the
    # class-sector double-processing of the empty/cross-block-underflow classes is
    # exact.  This is the engine-seam validation the individual kernels don't cover.
    # -------------------------------------------------------------------------
    @testset "TMCOUP driver seam (dftm_couple!) — cycle-2 dump-replay bit-exact" begin
        # DBGIN_TRE: I IMC ICR NORMHT ITRUNC WK2 HT HTG DBH DG PCT PROB (7 hex reals)
        din = split(strip("""
        1   1  32   0       0 3E645B51 4291D1B4 40809BFC 414FFE05 3F07AEA0 421F878D 4093DA38
        2   1  33   0       0 3E3BB828 429994A9 407A7153 416CB1A5 3F16F7F0 428689B2 4056CD2A
        3   1  33   0       0 3D882AB8 42A5B651 40678DEB 418BB9F6 3F1A8830 42B188B2 402E9B5F
        4   1  34   0       0 3E20BBF2 428D1ECB 408D15DA 413CC8C5 3F1BCCD0 41B60A30 40DFFF5C
        5   1  35   0       0 3E3E412E 42859F16 40A2AA2B 412E2166 3F5244F0 41846AE0 410A80DF
        6   1  33   0       0 3E21083A 429144A2 40831D82 41491916 3F097B10 420BA686 40B75025
        7   1  34   0       0 3E06B42E 4298D001 4074C17E 41654016 3F074990 427A600F 408210B3
        8   1  34   0       0 3D7DC302 42A316EE 407809BE 41859D62 3F30F390 42A6F071 404868AC
        9   1  32   0       0 3F3CD668 426FF218 4098F5D7 410C7930 3EDAFB30 409A1FC4 411D11A4
        10  1  33   0       0 3D7E3C85 42A8E5B8 404D241D 4192FC29 3EE98C00 42BBE4F3 401A4C03
        11  1  36   0       0 3D98A52A 42903299 4094809F 415712B7 3F17C850 423AB501 40BD0B95
        12  1  34   0       0 3E62EB71 4293F527 408F0AA0 415CB289 3F0A8110 424D740E 4077AC52
        13  1  33   0       0 3EDE2C8F 427960C2 408B46AC 4120B5E6 3EE18F20 41206AF2 4101CD8C
        14  1  34   0       0 3DEB41F2 429F7867 40893653 4180BA9C 3F04CFA0 429BD1C4 40416294
        15  1  36   0       0 3DCC037F 42937F0A 40B2EA9A 415CCB84 3F849640 42652161 409C4546
        16  1  35   0       0 3E0AFF28 4288877B 409C0239 413F17FB 3F26FCB0 41E93899 40E17C01
        17  1  35   0       0 3D930773 429CCB54 409EB718 417BBF1A 3F477BB0 4291DC13 4065EB69
        18  1  35   0       0 3CFFC31C 42A6F58E 409F2177 41935F28 3F59F7C0 42C80000 40335832
        """), '\n')
        # DBGA_REC: I PROB PCNEWF FBIOMS
        arec = split(strip("""
        1 4093DA38 41F2DB66 4364E18B
        2 4056CD2A 41EADF32 438276F4
        3 402E9B5F 420145DD 439970E9
        4 40DFFF5C 420A0550 43510607
        5 410A80DF 420DBEEB 4345A521
        6 40B75025 42070412 435CC33F
        7 408210B3 4203414E 437C4945
        8 404868AC 4208E148 4392D15A
        9 411D11A4 41D14398 43237A47
        10 401A4C03 41F9DE38 43A1BDCB
        11 40BD0B95 41700000 43C80000
        12 4077AC52 41700000 43C17A71
        13 4101CD8C 41700000 438372C3
        14 40416294 418352E6 43C80000
        15 409C4546 41700000 43C80000
        16 40E17C01 41700000 43BCD74A
        17 4065EB69 418411E4 43C80000
        18 40335832 41B7B81C 43C80000
        """), '\n')
        # DBGC_TRE expected: I IMC ICR NORMHT ITRUNC WK2 HT HTG DBH DG
        dexp = split(strip("""
        8   1  34   0       0 3F957706 42A316EE 404CEFA2 41859D62 3EC95EE7
        3   1  33   0       0 3E6D7724 42A5B651 40446B89 418BB9F6 3EC84601
        10  1  32   0       0 3E51D803 42A6FF5B 4010587A 4192FC29 3E9756A0
        7   1  34   0       0 4019D756 4298D001 4043EFFE 41654016 3E51F74A
        5   1  35   0       0 41022A5E 42859F16 408221BC 412E2166 3EA32B42
        4   1  34   0       0 40D28349 428D1ECB 4061BC90 413CC8C5 3E71CD3C
        6   1  33   0       0 40AC4712 429144A2 4051C8D0 41491916 3E555EBE
        2   1  33   0       0 3FFE1130 429994A9 404B09D2 416CB1A5 3E6A4D9A
        1   1  32   0       0 408AF3A2 4291D1B4 404DC660 414FFE05 3E529424
        9   1  32   0       0 41139D05 426FF218 404EB708 410C7930 3DD56E6D
        18  1  35   0       0 3E1B870E 42A6F58E 4086689D 41935F28 3F08E23E
        17  1  35   0       0 3E1E1E43 429CCB54 408A34E5 417BBF1A 3F1C9823
        14  1  34   0       0 3E04384F 429F7867 4075F32A 4180BA9C 3ED0836A
        15  1  36   0       0 3E3789D5 42937F0A 40A23DAE 415CCB84 3F502956
        11  1  36   0       0 3E4F33EC 42903299 4086A960 415712B7 3EEE4C6A
        12  1  34   0       0 3E62EB71 4293F527 4083206F 415CB289 3ED97397
        16  1  35   0       0 3E94AD25 4288877B 40799D28 413F17FB 3F0315B8
        13  1  33   0       0 3EDE2C8F 427960C2 4079FCF1 4120B5E6 3EB1105B
        """), '\n')
        n = 18
        prob=zeros(Float32,n); wk2=zeros(Float32,n); dbh=zeros(Float32,n); ht=zeros(Float32,n)
        dg=zeros(Float32,n); htg=zeros(Float32,n); pct=zeros(Float32,n)
        icr=zeros(Int32,n); imc=ones(Int32,n); normht=zeros(Int32,n); itrunc=zeros(Int32,n)
        kutkod=zeros(Int32,n); fbioms=zeros(Float32,n); pcnewf=zeros(Float32,n)
        for l in din
            t = split(strip(l)); i = parse(Int, t[1])
            imc[i]=parse(Int32,t[2]); icr[i]=parse(Int32,t[3]); normht[i]=parse(Int32,t[4]); itrunc[i]=parse(Int32,t[5])
            wk2[i]=_fromhex(t[6]); ht[i]=_fromhex(t[7]); htg[i]=_fromhex(t[8]); dbh[i]=_fromhex(t[9])
            dg[i]=_fromhex(t[10]); pct[i]=_fromhex(t[11]); prob[i]=_fromhex(t[12])
        end
        for l in arec
            t = split(strip(l)); i = parse(Int, t[1]); pcnewf[i]=_fromhex(t[3]); fbioms[i]=_fromhex(t[4])
        end
        d = _F.dftm_defaults!(_F.InlandEmpire())     # seed 55329, IBMTYP def; set method-2 + RANLARVA
        d.active = true; d.ibmtyp = Int32(2); d.iegtyp = Int32(1)
        gipt = Int32.(collect(1:18))                 # dense IND1 = identity
        isct = zeros(Int, _F.MAXSP, 2)
        isct[3,1]=1; isct[3,2]=10; isct[4,1]=11; isct[4,2]=18   # ISCT: DF (1,10), GF (11,18)
        _F.dftm_couple!(d, gipt, isct, (10,8), prob, wk2, dbh, ht, dg, htg, icr, pct, normht,
            itrunc, imc, kutkod, fbioms, pcnewf; fint=5.0f0, weight=(1.0f0,1.0f0), tmpn1=0.5f0,
            iegtyp=1, ldf=true, lgf=true, itmslv=0, tmdefl=50.0f0)
        maxu = 0; badint = 0
        for l in dexp
            t = split(strip(l)); i = parse(Int, t[1])
            (imc[i]==parse(Int,t[2]) && icr[i]==parse(Int,t[3]) &&
             normht[i]==parse(Int,t[4]) && itrunc[i]==parse(Int,t[5])) || (badint += 1)
            maxu = max(maxu, _ulps(wk2[i], _fromhex(t[6])), _ulps(ht[i], _fromhex(t[7])),
                       _ulps(htg[i], _fromhex(t[8])), _ulps(dbh[i], _fromhex(t[9])),
                       _ulps(dg[i], _fromhex(t[10])))
        end
        @test badint == 0     # IMC/ICR/NORMHT/ITRUNC exact, all 18 records
        @test maxu == 0       # WK2/HT/HTG/DBH/DG bit-exact, all 18 records — whole TMCOUP composition
    end

    # -------------------------------------------------------------------------
    # DFTM engine seam end-to-end (simulate.jl): a scheduled MANSTART/MANSCHED
    # tussock-moth outbreak fires through the LIVE cycle loop — DFTMGO+INSCYC forces
    # the TMBASE 5-yr outbreak cycle (a NEW boundary appears) and TMCOUP raises the
    # host mortality + growth loss.  The .sum-DELTA vs the relinked FVSie_dftm oracle
    # is CORNERED by the documented IE #206 OLDRN growth straddle (the pre-DFTM cyc-2
    # stand differs; the coupling math itself is the bit-exact dump-replay above), so
    # here we assert the qualitative seam behaviour: the outbreak cycle is inserted
    # and the host TPA collapses.  A DFTM block WITHOUT a MANSCHED outbreak stays
    # byte-identical to no DFTM (the inert-seam guarantee, also covered below).
    # -------------------------------------------------------------------------
    @testset "DFTM engine seam — MANSCHED outbreak fires + INSCYC (live cycle loop)" begin
        v = FVSjl.InlandEmpire()
        dir2 = mktempdir()
        tre = """
   1      248112       0101   011DF 12016   0654   00111     0  0
   2      248112       0101   011DF 14018   0704   00111     0  0
   3      248112       0102   011DF 16020   0754   00111     0  0
   4      248112       0102   011DF 10012   0604   00111     0  0
   5      248112       0103   011DF 09011   0554   00111     0  0
   6      248112       0103   011DF 11014   0634   00111     0  0
   7      248112       0104   011DF 13017   0684   00111     0  0
   8      248112       0104   011DF 15019   0734   00111     0  0
   9      248112       0105   011DF 08010   0504   00111     0  0
  10      248112       0105   011DF 17021   0774   00111     0  0
  11      248112       0106   011GF 11013   0604   00111     0  0
  12      248112       0106   011GF 13016   0664   00111     0  0
  13      248112       0107   011GF 09011   0524   00111     0  0
  14      248112       0107   011GF 15018   0714   00111     0  0
  15      248112       0108   011GF 12014   0634   00111     0  0
  16      248112       0108   011GF 10012   0574   00111     0  0
  17      248112       0109   011GF 14017   0684   00111     0  0
  18      248112       0109   011GF 16019   0724   00111     0  0
"""
        write(joinpath(dir2, "dn.tre"), tre)
        head = _dftm_head("DFTM SEAM  ")   # SCREEN/NOAUTOES/NOTRIPLE/… INVYEAR 1990 NUMCYCLE 5
        onblk  = "DFTM\nMANSTART\nMANSCHED           2\nBIOMASS            2\nEND\n"
        write(joinpath(dir2, "dn_on.key"),  head * onblk * "ECHOSUM\nPROCESS\nSTOP\n")
        write(joinpath(dir2, "dn_off.key"), head *          "ECHOSUM\nPROCESS\nSTOP\n")
        cp(joinpath(dir2, "dn.tre"), joinpath(dir2, "dn_on.tre"))
        cp(joinpath(dir2, "dn.tre"), joinpath(dir2, "dn_off.tre"))
        yr_tpa(key) = begin
            rows = filter(l -> !startswith(l, "-999") && !isempty(strip(l)),
                          split(strip(FVSjl.run_keyfile(key; variant = v, output = :sum)), '\n'))
            [(parse(Int, split(r)[1]), parse(Int, split(r)[3])) for r in rows]  # (year, TPA)
        end
        off = yr_tpa(joinpath(dir2, "dn_off.key"))
        on  = yr_tpa(joinpath(dir2, "dn_on.key"))
        offyears = first.(off); onyears = first.(on)
        # INSCYC forces the 5-yr outbreak cycle ⇒ a 2005 boundary the DFTM-off run lacks.
        @test !(2005 in offyears)
        @test 2005 in onyears
        # The outbreak collapses the host TPA in the inserted cycle (>30% drop 2000→2005).
        tpa2000 = only(t for (y, t) in on if y == 2000)
        tpa2005 = only(t for (y, t) in on if y == 2005)
        @test tpa2005 < 0.7 * tpa2000
        # Inert seam: a DFTM block that schedules NO MANSCHED outbreak is byte-identical to no DFTM.
        write(joinpath(dir2, "dn_ns.key"), head * "DFTM\nMANSTART\nBIOMASS 2\nEND\n" * "ECHOSUM\nPROCESS\nSTOP\n")
        cp(joinpath(dir2, "dn.tre"), joinpath(dir2, "dn_ns.tre"); force = true)
        strip999(key) = filter(l -> !startswith(l, "-999"),
                               split(strip(FVSjl.run_keyfile(key; variant = v, output = :sum)), '\n'))
        @test strip999(joinpath(dir2, "dn_ns.key")) == strip999(joinpath(dir2, "dn_off.key"))
    end

    # -------------------------------------------------------------------------
    # INSCYC (inscyc.f) — force the TMBASE (5-yr) outbreak cycle into the
    # schedule.  Integer cycle-year math ⇒ bit-exact.  Golden = the dense.key
    # DFTMGO run ("INSCYC:IY= 1990 2000 2005 2010 2020 2030 2040", ISPOT=3,
    # IFINT=5) with ICYC=2, IBOUND=TMBASE=5.
    # -------------------------------------------------------------------------
    @testset "INSCYC cycle-forcing (inscyc.f) — golden IY schedule" begin
        iy = [1990, 2000, 2010, 2020, 2030, 2040, 0]   # ncyc=5 (+1 growth slot)
        (nc, fi, sp) = _F.dftm_inscyc!(iy, 5, 10, 2, 5)
        @test iy[1:nc+1] == [1990, 2000, 2005, 2010, 2020, 2030, 2040]
        @test nc == 6 && fi == 5 && sp == 3            # NCYC+1, IFINT=5, inserted at subscript 3
        # a boundary already at the target year ⇒ no insertion (ISPOT 0, unchanged).
        iy2 = [1990, 2000, 2005, 2010, 2020, 0]
        (nc2, _, sp2) = _F.dftm_inscyc!(iy2, 5, 5, 2, 5)
        @test sp2 == 0 && nc2 == 5 && iy2[1:6] == [1990, 2000, 2005, 2010, 2020, 0]
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
