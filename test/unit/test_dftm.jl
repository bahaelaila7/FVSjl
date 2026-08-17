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
