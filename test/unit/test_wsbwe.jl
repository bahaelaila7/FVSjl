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

    @testset "WSBWE is INERT (.sum byte-identical with/without a WSBW block)" begin
        base = rows(ctrl_key)
        @test rows(min_key)   == base       # WSBW/DAMAGE/END only ⇒ byte-identical
        @test rows(wsbwe_key)  == base       # full WSBW block ⇒ still byte-identical (inert)
    end

    # --- DEFOLIATION EFFECT KERNEL goldens (dump-replay validated bit-exact vs the
    # instrumented live oracle FVSem_wsbwe on a heavy DF/ES/AF host stand under
    # sustained 85/80/70/60% DEFOL; see scratchpad/wsbwe/GOLDENS.md). Each golden is
    # (inputs Float32-hex → output Float32-hex) captured from bwedam/bwepdm dumps. ---
    _f(h) = reinterpret(Float32, parse(UInt32, h; base = 16))

    @testset "BWEDAM Nichols growth-loss RDDS/RHTG (bit-exact, 90/90 full-run)" begin
        # (avprbo_top, avprbo_whole, rddsm1, rhtgm1) → (rdds, rhtg)
        for (a1, a2, r1, h1, rd, rh) in (
            ("3E199998","3E199998","3F800000","3F800000","3EBA93DE","3F030E28"), # yr1 carry=1
            ("3E199998","3E199997","3F800000","3F800000","3EBA93DC","3F030E28"),
            ("3E199998","3E199998","3EBA93DE","3F030E28","3E6E7E17","3EBA67EF"), # yr2 carry≠1
            ("3E199998","3E199998","3EBA93DC","3F030E28","3E6E7E15","3EBA67EF"))
            @test _hexw(_FW.wsbwe_rdds(_f(a2), _f(r1))) == rd
            @test _hexw(_FW.wsbwe_rhtg(_f(a1), _f(h1))) == rh
        end
    end

    @testset "BWEPDM Marsden mortality-prob logistic (bit-exact, 81/81 full-run)" begin
        # (ih, elev, pntba, pnthba, mft, mfm, ktk) → pr
        for (ih, el, pb, hb, mft, mfm, ktk, pr) in (
            (2,"42580000","400A809A","400A809A","41200000","41200000","00000000","3E402809"), # DF
            (2,"42580000","3FAE5663","3FAE5663","41200000","41200000","00000000","3E3F2A59"),
            (5,"42580000","400A809A","400A809A","41200000","41200000","00000000","3E402809"), # ES
            (4,"42580000","3FDB330C","3FDB330C","41200000","41200000","00000000","3E879871")) # AF
            @test _hexw(_FW.wsbwe_mort_pr(ih, _f(el), _f(pb), _f(hb),
                                          _f(mft), _f(mfm), _f(ktk))) == pr
        end
        # bwepdm.f:645 guard — zero missing-foliage ⇒ PR=0
        @test _FW.wsbwe_mort_pr(2, _f("42580000"), _f("400A809A"), _f("400A809A"),
                                0.0f0, 5.0f0, 0.0f0) == 0.0f0
    end

    @testset "wsbwe_go gate (BWEGO manual-DEFOL branch) + apply! inert" begin
        w = _FW.wsbwe_defaults!()
        @test _FW.wsbwe_go(w) == false                 # not active
        w.active = true; w.ldefol = true
        @test _FW.wsbwe_go(w) == false                 # active+ldefol but no schedule
        push!(w.defol_sched, (1990f0,3f0,0f0,85f0,80f0,70f0,60f0))
        @test _FW.wsbwe_go(w) == true                  # LDEFOL + scheduled ⇒ fires
        w.ldefol = false; w.lbudl = true
        @test _FW.wsbwe_go(w) == false                 # BUDLITE branch deferred
    end
end
