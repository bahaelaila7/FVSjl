# Douglas-fir Tussock Moth (DFTM) defoliator model — keyword reader + the TMRANN
# RNG + the deterministic TMOTPR outbreak probability + DFTMGO gate + TMBMAS
# biomass (chunks 0-1), PLUS the DFTMOD population integrator (→ DPCENT) and the
# TMCOUP per-tree damage functions, both LIVE dump-replayed bit-exact against the
# instrumented FVSie_dftm oracle on the dense DF/GF host stand.  Inert-seam kept.
#
# VALIDATED HERE (dump-replay, 0 ULP): DFTMOD DPCENT (17 non-empty classes + the
# empty-GF-class NaN reset); TMCOUP T/ITAB, WK2, DG, HTG (no-top-kill), top-kill
# class-1 leader-loss + class-2 crown PCKILL/HTGLOS.  NOT yet dump-replayed (draft
# code present in dftm.jl, engine-inert): GARBEL/GRCLAS classification (Z2/Z3/Z4 +
# ISC sector pointers) and the full TMRANN-stream ordering (RANLARVA egg draws +
# per-tree PRTOPK), plus the INSCYC engine hook + simulate.jl coupling seam.
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
