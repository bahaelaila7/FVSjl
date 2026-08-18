# Western Spruce Budworm (WSBWE) beachhead — the `WSBW … END` stand-level keyword
# block reader (kw_wsbwe!) + the BWERAN MINSTD RNG + BWEINT/BLOCK-DATA defaults,
# with the INERT-seam A/B.  WSBWE is a STAND-LEVEL defoliation model (BWEGO from
# base GRINCR/GRADD, BWEINT from INITRE — no PPE dependency), so a stand oracle IS
# relinkable and runs (scratchpad/wsbwe/FVSem_wsbwe = FVSem .o − base/exbudl.o +
# wsbwe/*.o; WSBW-off ≡ stock is byte-identical, cksum-verified).  This beachhead
# ports the reachable, self-contained numeric surface — the RNG — plus the block
# reader; the defoliation → growth-loss/mortality effect path (BWEGO→BWEDR) is
# DEFERRED, so the reader is faithfully INERT (no engine seam wired yet).
#
# VALIDATED HERE:
#   * BWERAN  — pristine wsbwe/bweran.f, seed 55329 → 8-draw Float32 hex stream,
#     BIT-EXACT vs the standalone gfortran-16 driver over the PRISTINE Fortran
#     (scratchpad/wsbwe/driver_bweran.f).  IDENTICAL to LPMPB(MPRANN)/DFTM(TMRANN)/
#     WWPB(BMRANN) — same MINSTD LCG (a=16807, m=2^31-1, ÷2^31), same seed 55329.
#   * BWERSD  — reseed odd-forcing + LSET=false default-reset (bweran.f entries).
#   * kw_wsbwe! — DEFOL/GENDEFOL/DAMAGE/OBSCHED parse into s.wsbwe.
#   * INERT seam — a stand carrying a WsbweState projects .sum-byte-identically to
#     one with no WSBW block.

using Test
using FVSjl
const _FW = FVSjl

_hexw(x::Float32) = uppercase(string(reinterpret(UInt32, x); base = 16, pad = 8))

# --- golden from scratchpad/wsbwe/driver_bweran.f over pristine wsbwe/bweran.f ---
const _G_BWERAN = ("3EDDB57A", "3F5AB21B", "3F630A07", "3F2734C9",
                   "3EF4FB2C", "3F4AF646", "3F6E7B68", "3F67EA59")

_kwrec(kw, fields...) = rpad(kw, 10) * join(lpad(string(f), 10) for f in fields)

# EM host stand (S248112): DF (sp 3 → budworm class 2) is a host.
const _WSBWE_TRE = """
   1      248112       0101   011LP 11510   0734   00111     0  0
   2      248112       0101   031DF 001     0026   00222     0  0
   3      248112       0102   011PP 06523   0308   00111     0  0
  17      248112       0106   031DF 10010   0654   00111     0  0
"""

_wsbwe_head(title) = """
SCREEN
NOAUTOES
STATS
STDIDENT
S248112  $title
DESIGN                                        11.0       1.0
STDINFO        112.0     260.0      60.0     315.0      30.0      54.0
INVYEAR       1990.0
NUMCYCLE        10.0
TREEFMT
(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,
T52,I2,T66,5I1,T54,7I1,T75,F3.0)
TREEDATA
"""

# A WSBW block exercising manual DEFOL (deterministic activation), DAMAGE output,
# and OBSCHED (user outbreak window).
const _WSBWE_BLOCK = "WSBW\n" *
    _kwrec("DEFOL", 1990, 3, 0, "0.90", "0.90", "0.90", "0.90") * "\n" *
    "DAMAGE\n" *
    _kwrec("OBSCHED", 3, 1995, 2005) * "\n" *
    "END\n"

const _WSBWE_MIN = "WSBW\nDAMAGE\nEND\n"

@testset "Western Spruce Budworm (WSBWE) — WSBW block + BWERAN RNG (inert seam)" begin

    @testset "BWERAN Lehmer stream (bit-exact vs pristine bweran.f, seed 55329)" begin
        w = _FW.wsbwe_defaults!()
        @test w.dseed == 55329.0f0
        for g in _G_BWERAN
            @test _hexw(_FW.wsbwe_rand!(w)) == g
        end
    end

    @testset "BWERSD reseed (odd forcing + default reset)" begin
        w = _FW.wsbwe_defaults!()
        _FW.wsbwe_seed!(w, 100.0f0; lset=true)   # even → forced odd 101
        @test w.dseed == 101.0f0 && w.rng_s0 == 101.0
        _FW.wsbwe_seed!(w, 55.0f0; lset=true)    # already odd
        @test w.dseed == 55.0f0 && w.rng_s0 == 55.0
        _FW.wsbwe_rand!(w)                        # advance S0
        @test w.rng_s0 != 55.0
        _FW.wsbwe_seed!(w, 0.0f0; lset=false)     # LSET=false: reset S0 to stored seed (55)
        @test w.rng_s0 == 55.0
    end

    @testset "BWEINT / BLOCK-DATA defaults" begin
        w = _FW.wsbwe_defaults!()
        @test w.active == false
        @test w.ldefol == false && w.lbudl == false
        @test w.iobloc == 2 && w.iwopt == 1
        @test w.dseed == 55329.0f0
        @test isempty(w.defol_sched) && isempty(w.obsched)
    end

    dir = mktempdir()
    write(joinpath(dir, "shared.tre"), _WSBWE_TRE)
    v = FVSjl.EasternMontana()

    wsbwe_key = joinpath(dir, "wsbwe.key")
    write(wsbwe_key, _wsbwe_head("WSBWE PARSED") * _WSBWE_BLOCK * "ECHOSUM\nPROCESS\nSTOP\n")
    cp(joinpath(dir, "shared.tre"), joinpath(dir, "wsbwe.tre"))

    @testset "WSBW keyword reader populates s.wsbwe (bwein.f)" begin
        w = nothing
        for s in FVSjl.each_stand(wsbwe_key; variant = v)
            w = s.wsbwe
            break
        end
        @test w !== nothing
        @test w.active == true
        @test w.ldefol == true          # DEFOL activated the manual path
        @test w.lbwdam == true          # DAMAGE output flag
        @test length(w.defol_sched) == 1
        @test w.defol_sched[1][1] == 1990.0f0      # idt
        @test w.defol_sched[1][2] == 3.0f0         # species field (DF)
        @test w.iobopt == 3
        @test w.obsched == [(Int32(1995), Int32(2005))]
    end

    ctrl_key = joinpath(dir, "ctrl.key")
    min_key  = joinpath(dir, "min.key")
    write(ctrl_key, _wsbwe_head("WSBWE CTRL ") * "ECHOSUM\nPROCESS\nSTOP\n")
    write(min_key,  _wsbwe_head("WSBWE MIN  ") * _WSBWE_MIN * "ECHOSUM\nPROCESS\nSTOP\n")
    for k in ("ctrl", "min", "wsbwe")
        cp(joinpath(dir, "shared.tre"), joinpath(dir, "$k.tre"); force = true)
    end
    rows(key) = filter(l -> !startswith(l, "-999"),
                       split(strip(FVSjl.run_keyfile(key; variant = v, output = :sum)), '\n'))

    @testset "WSBWE seam gating (no-DEFOL block ⇒ byte-identical; wsbwe_go off)" begin
        base = rows(ctrl_key)
        # A WSBW block WITHOUT a DEFOL activity never enters wsbwe_apply! (wsbwe_go=false),
        # so it is byte-identical. (A firing DEFOL on a budworm-host stand now DOES change
        # the .sum — the effect path is LIVE; that path is validated end-to-end vs the live
        # relinked oracle FVSem_wsbwe in scratchpad/wsbwe/ — see GOLDENS.md.)
        @test rows(min_key) == base        # WSBW/DAMAGE/END only (no DEFOL) ⇒ byte-identical
    end

    @testset "Effect kernels + BWERNP/BWEBET damage RNG (glibc, bit-exact goldens)" begin
        # BWEDAM growth-loss (bwedam.f:184/188) — dump-replay goldens (GOLDENS.md)
        @test _hexw(_FW.wsbwe_rdds(reinterpret(Float32,0x3E199998), 1.0f0))       == "3EBA93DE"
        @test _hexw(_FW.wsbwe_rhtg(reinterpret(Float32,0x3E199998), 1.0f0))       == "3F030E28"
        # BWEPDM Marsden mortality logistic (bwepdm.f:640) — DF, ELEV=54, MFT=MFM=10, KTK=0
        @test _hexw(_FW.wsbwe_mort_pr(2, 54.0f0, reinterpret(Float32,0x400A809A),
                    reinterpret(Float32,0x400A809A), 10.0f0, 10.0f0, 0.0f0))      == "3E402809"
        # BWERNP(0.85,.06) then the topkill BWERAN — the exact per-tree draw order of the
        # first BWEPDM host tree (seed DSEEDD=55329). Validates BWEBET rejection-loop order.
        w = _FW.wsbwe_defaults!()
        w.rng_s0 = Float64(w.dseed)                     # BWERPT(DSEEDD=55329)
        @test _hexw(_FW.wsbwe_bernp(w, 0.85f0, 0.06f0)) == "3F609607"   # AVDEF (2 BWERAN draws)
        @test _hexw(_FW.wsbwe_rand!(w))                 == "3F630A07"   # topkill BWERAN (draw #3)
    end

    @testset "Per-variant host/biomass dispatch (EM + TT)" begin
        # TT (Teton, MAXSP=18): IBWSPM/IBIOMP from bwebktt.f/bwebmstt.f. Hosts (IBWSPM<7)
        # at sp3=DF(2)/sp8=ES(5)/sp9=AF(4) — same host classes as EM. Host-class tables +
        # ICVOPT=2 biomass coeffs are identical to EM (only these two maps differ).
        @test _FW.wsbwe_ibwspm_for(_FW.EasternMontana()) === _FW.WSBWE_IBWSPM_EM
        @test _FW.wsbwe_ibwspm_for(_FW.Teton())          === _FW.WSBWE_IBWSPM_TT
        @test _FW.wsbwe_ibiomp_for(_FW.Teton())          === _FW.WSBWE_IBIOMP_TT
        @test _FW.WSBWE_IBWSPM_TT == Int[7,7,2,7,7,7,7,5,4,7,7,7,7,7,7,7,7,7]
        @test _FW.WSBWE_IBIOMP_TT == Int[1,4,3,11,8,11,7,8,9,10,11,11,11,11,11,11,11,11]
        @test length(_FW.WSBWE_IBWSPM_TT) == 18                       # TT MAXSP
        @test (_FW.WSBWE_IBWSPM_TT[3], _FW.WSBWE_IBWSPM_TT[8], _FW.WSBWE_IBWSPM_TT[9]) == (2, 5, 4)  # DF/ES/AF
        @test _FW.wsbwe_ibwspm_for(_FW.InlandEmpire()) === nothing    # unsupported variant → inert
        # end-to-end .sum-DELTA validated in the main loop (2000 ΔTPA -478 = FVStt_wsbwe -478,
        # bit-exact; multi-cycle cornered by TT #206) — the oracle isn't available to the unit suite.
    end
end
