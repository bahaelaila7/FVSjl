# Westwide Pine Beetle (WWPB) model — the reachable BMIN stand-level output-block
# reader (kw_wwpbin!) + the BMRANN MINSTD RNG + BLOCK-DATA/BMINIT defaults, with
# the INERT-seam A/B.  WWPB is a landscape Parallel-Processing-Extension model
# whose outbreak driver (BMDRV from ALSTD2) + FVS mortality hand-back (BMKILL/WK2
# from PPMAIN) live in a PPE spatial harness that is ABSENT from this FVS tree —
# every FVS*_buildDir links the base/exbm.f no-op stub and no sourceList compiles
# wwpb/*.f.  So the outbreak is unreachable via any single-stand run; the only
# reachable entry is the stand-level BMIN output block (keywds.f option 126,
# initre.f `CALL BMIN`), which merely schedules .bm* report activities and applies
# no mortality — hence faithfully inert.
#
# VALIDATED HERE:
#   * BMRANN  — the full pristine wwpb/bmrann.f, seed 55329 → 8-draw Float32 hex
#     stream, BIT-EXACT vs the standalone gfortran-16 driver over the PRISTINE
#     Fortran (scratchpad/wwpb/driver_bmrann.f).  The first 6 draws coincide with
#     DFTM's TMRANN — same MINSTD LCG (a=16807, m=2^31-1, ÷2^31), same seed.
#   * BMRNSD  — reseed odd-forcing + LSET=false default-reset.
#   * kw_wwpbin! — MAINOUT/TREEOUT/BKPOUT/VOLOUT/END parse into s.wwpb (flags +
#     the OPNEW output requests with faithful IDT/nyears/incr field decode).
#   * INERT seam — a stand carrying a WwpbState projects .sum-byte-identically to
#     one with no BMIN block.
#
# DOCTRINE: like DFB/DFTM/WPBR, WWPB ships nowhere; a numeric oracle would be a
# RELINK, but here the reachable numeric surface (the RNG) is small + self-
# contained, so its golden comes from a tiny gfortran-16 driver over the pristine
# Fortran.  The full beetle dynamics await a port of the absent PPE harness.

using Test
using FVSjl
const _FW = FVSjl

_hexw(x::Float32) = uppercase(string(reinterpret(UInt32, x); base = 16, pad = 8))

# --- golden from scratchpad/wwpb/driver_bmrann.f over pristine wwpb/bmrann.f ---
const _G_BMRANN = ("3EDDB57A", "3F5AB21B", "3F630A07", "3F2734C9",
                   "3EF4FB2C", "3F4AF646", "3F6E7B68", "3F67EA59")

# Column-safe keyword-record builder: keyword left-justified in cols 1-10, then
# each numeric field right-justified in a 10-col block (cols 11-20, 21-30, …),
# exactly the FVS fixed-format the KeywordReader parses.
_kwrec(kw, fields...) = rpad(kw, 10) * join(lpad(string(f), 10) for f in fields)

const _WWPB_TRE = """
   1      248112       0101   011LP 11510   0734   00111     0  0
   2      248112       0101   031DF 001     0026   00222     0  0
   9      248112       0103   011LP 09511   0603   00111     0  0
  17      248112       0106   011DF 10010   0654   00111     0  0
"""

_wwpb_head(title) = """
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

# A BMIN block exercising all four output keywords: two with explicit fields, two
# defaulted.  MAINOUT idt=1 nyears=100 incr=5 ; TREEOUT defaults (1/100/5) ;
# BKPOUT idt=2 nyears=50 incr=10 ; VOLOUT defaults.
const _WWPB_BLOCK = "BMIN\n" *
    _kwrec("MAINOUT", 1, 100, 5) * "\n" *
    "TREEOUT\n" *
    _kwrec("BKPOUT", 2, 50, 10) * "\n" *
    "VOLOUT\n" *
    "END\n"

# A BMIN block that only names the extension (all defaults) — for the inert A/B.
const _WWPB_MIN = "BMIN\nEND\n"

@testset "Westwide Pine Beetle (WWPB) — BMIN block + BMRANN RNG (inert seam)" begin

    @testset "BMRANN Lehmer stream (bit-exact vs pristine bmrann.f, seed 55329)" begin
        w = _FW.wwpb_defaults!(_FW.InlandEmpire())
        @test w.rng_s0 == 55329.0 && w.rng_ss == 55329.0f0
        for g in _G_BMRANN
            @test _hexw(_FW.wwpb_rand!(w)) == g
        end
    end

    @testset "BMRNSD reseed (odd forcing + default reset)" begin
        w = _FW.wwpb_defaults!(_FW.InlandEmpire())
        _FW.wwpb_seed!(w, 100.0f0, true)        # even → forced odd 101
        @test w.rng_ss == 101.0f0 && w.rng_s0 == 101.0
        _FW.wwpb_seed!(w, 55.0f0, true)         # already odd
        @test w.rng_ss == 55.0f0 && w.rng_s0 == 55.0
        _FW.wwpb_rand!(w)                       # advance S0
        @test w.rng_s0 != 55.0
        _FW.wwpb_seed!(w, 0.0f0, false)         # LSET=false: reset S0 to SS (55)
        @test w.rng_s0 == 55.0
    end

    @testset "BMINIT / BLOCK-DATA defaults" begin
        w = _FW.wwpb_defaults!(_FW.InlandEmpire())
        @test w.active == false
        @test w.pbspec == 1                     # MPB
        @test w.upsiz == Float32[3, 6, 9, 12, 15, 18, 21, 25, 30, 50]
        @test w.iscmin == Int32[3, 3, 2]
        @test w.rng_s0 == 55329.0
        @test isempty(w.outreqs)
    end

    # -------------------------------------------------------------------------
    # Keyword reader (bmin.f) via the real keyword dispatch.
    # -------------------------------------------------------------------------
    dir = mktempdir()
    write(joinpath(dir, "shared.tre"), _WWPB_TRE)
    v = FVSjl.InlandEmpire()

    wwpb_key = joinpath(dir, "wwpb.key")
    write(wwpb_key, _wwpb_head("WWPB PARSED") * _WWPB_BLOCK * "ECHOSUM\nPROCESS\nSTOP\n")
    cp(joinpath(dir, "shared.tre"), joinpath(dir, "wwpb.tre"))

    @testset "BMIN keyword reader populates s.wwpb (bmin.f)" begin
        w = nothing
        for s in FVSjl.each_stand(wwpb_key; variant = v)
            w = s.wwpb
            break
        end
        @test w !== nothing
        @test w.active == true
        @test w.lbmain && w.lbmtre && w.lbmbkp && w.lbmvol
        @test length(w.outreqs) == 4
        # MAINOUT explicit (1/100/5)
        @test w.outreqs[1].myact == 2701
        @test (w.outreqs[1].idt, w.outreqs[1].nyears, w.outreqs[1].incr) == (Int32(1), Int32(100), Int32(5))
        # TREEOUT defaulted (1/100/5)
        @test w.outreqs[2].myact == 2702
        @test (w.outreqs[2].idt, w.outreqs[2].nyears, w.outreqs[2].incr) == (Int32(1), Int32(100), Int32(5))
        # BKPOUT explicit (2/50/10)
        @test w.outreqs[3].myact == 2703
        @test (w.outreqs[3].idt, w.outreqs[3].nyears, w.outreqs[3].incr) == (Int32(2), Int32(50), Int32(10))
        # VOLOUT defaulted
        @test w.outreqs[4].myact == 2704
        @test (w.outreqs[4].idt, w.outreqs[4].nyears, w.outreqs[4].incr) == (Int32(1), Int32(100), Int32(5))
    end

    # -------------------------------------------------------------------------
    # INERT seam: a WWPB-present stand must project .sum-byte-identically to a
    # stand with no BMIN keyword (WWPB applies no mortality — output-only, and
    # FVSjl runs no PPE loop).
    # -------------------------------------------------------------------------
    ctrl_key = joinpath(dir, "ctrl.key")
    min_key  = joinpath(dir, "min.key")
    write(ctrl_key, _wwpb_head("WWPB CTRL  ") * "ECHOSUM\nPROCESS\nSTOP\n")
    write(min_key,  _wwpb_head("WWPB MIN   ") * _WWPB_MIN * "ECHOSUM\nPROCESS\nSTOP\n")
    for k in ("ctrl", "min", "wwpb")
        cp(joinpath(dir, "shared.tre"), joinpath(dir, "$k.tre"); force = true)
    end
    rows(key) = filter(l -> !startswith(l, "-999"),
                       split(strip(FVSjl.run_keyfile(key; variant = v, output = :sum)), '\n'))

    @testset "WWPB is INERT (.sum byte-identical with/without a BMIN block)" begin
        base = rows(ctrl_key)
        @test rows(min_key)  == base       # BMIN/END only ⇒ byte-identical
        @test rows(wwpb_key) == base       # full BMIN block ⇒ still byte-identical (inert)
    end

    @testset "WWPB DISPERSE outbreak fires end-to-end (kills lodgepole)" begin
        # The DISPERSE keyword activates the reconstructed single-stand outbreak
        # (wwpb_apply! → wwpb_outbreak_cycle! → bmkill!). On this LP-host stand the
        # MPB kills lodgepole, so the .sum diverges from the inert BMIN-only run and
        # shows lower live TPA. End-to-end composition of the ten bit-exact kernels;
        # the outbreak magnitude is reconstruction, not bit-exact.
        disp_key = joinpath(dir, "disp.key")
        write(disp_key, _wwpb_head("WWPB DISP  ") *
              "BMIN\nDISPERSE         1.0       5.0       3.0      20.0\nEND\n" *
              "ECHOSUM\nPROCESS\nSTOP\n")
        cp(joinpath(dir, "shared.tre"), joinpath(dir, "disp.tre"); force = true)
        base = rows(ctrl_key); disp = rows(disp_key)
        @test disp != base                                    # the outbreak alters the projection
        @test parse(Float64, split(disp[end])[3]) < parse(Float64, split(base[end])[3])  # LP killed
    end

    # -------------------------------------------------------------------------
    # LANDSCAPE model (synthetic PPE harness, USER-approved) — chunk 1-2.
    # -------------------------------------------------------------------------
    @testset "BMDBHC — DBH → size class (bmdbhc.f vs UPSIZ breakpoints)" begin
        up = copy(FVSjl.WWPB_UPSIZ_DEFAULT)     # 3,6,9,12,15,18,21,25,30,50
        # class k covers [UPSIZ(k-1), UPSIZ(k)); class 1 = <3"; class 10 = ≥30"
        @test FVSjl.wwpb_dbh_class(up, 0.0)  == 1
        @test FVSjl.wwpb_dbh_class(up, 2.99) == 1
        @test FVSjl.wwpb_dbh_class(up, 3.0)  == 2     # DBH<UPSIZ(1)=3 false ⇒ advance
        @test FVSjl.wwpb_dbh_class(up, 5.99) == 2
        @test FVSjl.wwpb_dbh_class(up, 6.0)  == 3
        @test FVSjl.wwpb_dbh_class(up, 29.9) == 9
        @test FVSjl.wwpb_dbh_class(up, 30.0) == 10    # ≥UPSIZ(9)=30 ⇒ last class
        @test FVSjl.wwpb_dbh_class(up, 99.0) == 10
    end

    @testset "BMSDIT — FVS→BM tree bridge (binning + host split + averages)" begin
        w = FVSjl.wwpb_defaults!(v)             # PBSPEC=1 (MPB) ⇒ host = LP
        t = FVSjl.TreeList(100)
        # two LP (host) trees in size class 3 (6-9"): d=7, and one DF (nonhost) d=13 (class 5)
        function put!(i, sp, d, h, cr, htg, cfv, tpa)
            t.species[i]=Int32(sp); t.dbh[i]=Float32(d); t.height[i]=Float32(h)
            t.crown_pct[i]=Int32(cr); t.ht_growth[i]=Float32(htg); t.cuft_vol[i]=Float32(cfv)
            t.tpa[i]=Float32(tpa)
        end
        put!(1, 7,  7.0, 40.0, 50, 1.0,  8.0, 20.0)   # LP host, class 3
        put!(2, 7,  8.0, 44.0, 40, 1.2, 10.0, 10.0)   # LP host, class 3
        put!(3, 3, 13.0, 60.0, 55, 0.9, 25.0, 15.0)   # DF nonhost, class 5
        t.n = 3
        alpha(sp) = sp == 7 ? "LP" : (sp == 3 ? "DF" : "OT")
        st = FVSjl.WwpbStand()
        FVSjl.bmsdit!(st, t, w, alpha)
        # host TPA lands in class 3 (=30), nonhost in class 5 (=15)
        @test st.tree[3, 1] ≈ 30.0f0
        @test st.tree[5, 2] ≈ 15.0f0
        @test st.tree[FVSjl.WWPB_NSCL+1, 1] ≈ 30.0f0   # host summary slot
        @test st.tree[FVSjl.WWPB_NSCL+1, 2] ≈ 15.0f0   # nonhost summary slot
        # TPA-weighted mean height in class 3 = (40·20+44·10)/30
        @test st.hts[3, 1] ≈ Float32((40.0f0*20 + 44.0f0*10) / 30)
        # dominant host species = LP(7), nonhost = DF(3)
        @test st.isph[1] == 7
        @test st.isph[2] == 3
        # BA host summary = Σ d²·(π/576)·tpa for the two LP
        baexp = (7.0f0^2*20 + 8.0f0^2*10) * FVSjl.WWPB_PI24
        @test st.bah[FVSjl.WWPB_NSCL+1] ≈ baexp
        @test st.otpa[3, 1] ≈ 30.0f0                   # initial snapshot
    end

    @testset "BMMORT — fast/slow tree decrement + ledgers (bmmort.f)" begin
        st = FVSjl.WwpbStand()
        st.tree[3, 1] = 100.0f0                         # 100 host TPA in class 3
        st.bah[3]     = 40.0f0
        st.oakill[3, 1] = 0.10f0                        # 10% fast (windthrow/fire) kill, as proportion
        st.pbkill[3]  = 5.0f0                           # 5 TPA beetle kill (already TPA in BMISTD)
        FVSjl.bmmort!(st, false)                        # FAST pass: remove OAKILL
        @test st.tree[3, 1] ≈ 90.0f0                    # 100 − 0.10·100
        @test st.tpbk[3, 1, 1] ≈ 10.0f0                 # fast ledger
        @test st.fastk[1] ≈ 10.0f0
        @test st.oakill[3, 1] ≈ 0.0f0                   # zeroed after use
        FVSjl.bmmort!(st, true)                         # SLOW pass: remove PBKILL (host only)
        @test st.tree[3, 1] ≈ 85.0f0                    # 90·(1 − 5/90) = 85
        @test st.tpbk[3, 1, 3] ≈ 5.0f0                  # beetle ledger
        @test st.tree[FVSjl.WWPB_NSCL+1, 1] ≈ 85.0f0    # summary recomputed
    end

    @testset "BMCGRF — GRF/GRFSTD/RVDNST BIT-EXACT vs gfortran-16 driver (bmcgrf.f)" begin
        # golden: scratchpad/wwpb/driver_bmcgrf.f over the pristine wwpb/bmcgrf.f,
        # gfortran-16. Stand: BAH(icls)=2+1.5·icls, BANH=0.5, RVDSC=0.90, RVDFOL=0
        # (→1), all other stressors 0, PBSPEC=1, LCDENS=T. Float32 hex (Z8.8).
        hx(x) = reinterpret(UInt32, x)
        st = FVSjl.WwpbStand(); w = FVSjl.wwpb_defaults!(v)
        basum = 0.0f0
        for icls in 1:10
            st.bah[icls] = 2.0f0 + icls*1.5f0; st.banh[icls] = 0.5f0
            st.grf[icls] = 0.0f0; st.rvdsc[icls] = 0.90f0; st.rvdfol[icls] = 0.0f0
            basum += st.bah[icls] + st.banh[icls]
        end
        st.bah[11] = basum - 0.5f0*10; st.banh[11] = 0.5f0*10
        oldgrf = zeros(Float32, 10)
        FVSjl.bmcgrf!(st, w, oldgrf; lcdens = true)
        @test hx(st.grfstd) == 0x3F71870A     # GRFSTD (= 0.90 · RVDNST)
        @test hx(st.rvdnst) == 0x3F862E94     # stand-density RV (exp+powf via glibc)
        @test all(hx(st.grf[i]) == 0x3F666666 for i in 1:10)   # GRF = 0.90 (drought only)
        @test all(oldgrf[i] == 1.0f0 for i in 1:10)            # OLDGRF from GRF≤0 → 1
    end

    @testset "BMINIT coeffs — MSBA/UPBA/INC BIT-EXACT vs gfortran-16 (bminit.f)" begin
        # golden: scratchpad/wwpb/driver_bminit.f (verbatim bminit.f MSBA/UPBA/INC
        # lines), gfortran-16, UPSIZ=3,6,9,12,15,18,21,25,30,50. Float32 hex.
        hx(x) = reinterpret(UInt32, x)
        c = FVSjl.wwpb_init_coeffs(copy(FVSjl.WWPB_UPSIZ_DEFAULT))
        msba = UInt32[0x3C490FD0,0x3DE231CA,0x3E9D145B,0x3F19F01B,0x3F7E7803,
                      0x3FBE10F3,0x4004BB70,0x4038A7DB,0x4083FD8C,0x410BA051]
        upba = UInt32[0x3D490FD0,0x3E490FD0,0x3EE231CA,0x3F490FD0,0x3F9D145B,
                      0x3FE231CA,0x4019F01B,0x405A2A7E,0x409D145B,0x415A2A7E]
        inc1 = UInt32[0x3F0CCCCC,0x3F59999A,0x3F933333,0x3FB9999A,0x3FE00000,
                      0x40033333,0x40166666,0x402CCCCC,0x4049999A,0x40800000]
        @test all(hx(c.msba[i]) == msba[i] for i in 1:10)
        @test all(hx(c.upba[i]) == upba[i] for i in 1:10)
        @test all(hx(c.inc[1,i]) == inc1[i] for i in 1:10)
        @test all(hx(c.inc[2,i]) == inc1[i] for i in 1:10)   # INC(2,·) = INC(1,·)
        @test hx(c.inc[3,1]) == 0x3D6147AD                    # INC(3,·) = INC(1,1)·0.1
        @test hx(c.inc[3,10]) == 0x3D6147AD                   # Ips constant across classes
    end

    @testset "BMCBKP — BKP brood BIT-EXACT vs gfortran-16 driver (bmcbkp.f)" begin
        # golden: scratchpad/wwpb/driver_bmcbkp.f over the pristine wwpb/bmcbkp.f
        # (with GPGET2/GPADD stubbed to LOK=false ⇒ the normal, non-bad-year path),
        # gfortran-16. PBSPEC=1, NBGEN=1, NIBGEN=2, PBKILL(3)=10, PBKILL(5)=5, all
        # else 0 ⇒ BKP = Σ MSBA·PBKILL·INC(1,·). Float32 hex (Z8.8).
        hx(x) = reinterpret(UInt32, x)
        st = FVSjl.WwpbStand(); w = FVSjl.wwpb_defaults!(v)
        c = FVSjl.wwpb_init_coeffs(copy(FVSjl.WWPB_UPSIZ_DEFAULT))
        st.pbkill[3] = 10.0f0; st.pbkill[5] = 5.0f0
        FVSjl.bmcbkp!(st, w, c; nbgen = 1, nibgen = 2, ipson = false, ipsmin = 2)
        @test hx(st.bkp)    == 0x41439CF3     # BKP brood from PBKILL
        @test hx(st.bkpips) == 0x00000000     # no Ips input
        @test hx(st.oldbkp) == 0x41439CF3     # OLDBKP = BKP (PBSPEC≠3)
        @test all(st.pbkill[i] == 0.0f0 for i in 1:10)    # PBKILL consumed
    end

    @testset "BMCSPT — special-tree proportion BIT-EXACT vs gfortran-16 (bmcspt.f)" begin
        # golden: scratchpad/wwpb/driver_bmcspt.f over the pristine wwpb/bmcspt.f,
        # gfortran-16, PBSPEC=1. MEASURED: the nested overlap correction is dead
        # code (uninitialized JK loop never fires) ⇒ SPCLT = clamp(Σ SP, 0, 1) with
        # ==0→0 / ==1→1 short-circuits. Float32 hex.
        hx(x) = reinterpret(UInt32, x)
        st = FVSjl.WwpbStand(); w = FVSjl.wwpb_defaults!(v)
        for i in 1:10; st.tree[i,1] = 100.0f0; end
        @test hx(FVSjl.bmcspt!(st, w, 1, 1)) == 0x00000000            # all SP=0 ⇒ 0
        st.pitch[2] = 0.30f0
        @test hx(FVSjl.bmcspt!(st, w, 2, 1)) == 0x3E99999A            # one SP=0.3 ⇒ 0.3
        st.pitch[3] = 0.30f0; st.strike[3] = 0.40f0
        @test hx(FVSjl.bmcspt!(st, w, 3, 1)) == 0x3F333334            # 0.3+0.4=0.7 (NO overlap)
        st.pitch[4] = 0.20f0; st.strike[4] = 0.30f0; st.scorch[4] = 0.40f0
        @test hx(FVSjl.bmcspt!(st, w, 4, 1)) == 0x3F666666            # 0.2+0.3+0.4=0.9
        st.scorch[5] = 1.0f0
        @test hx(FVSjl.bmcspt!(st, w, 5, 1)) == 0x3F800000            # any SP==1 ⇒ 1.0
    end

    @testset "BMCNUM — attractiveness numerator BIT-EXACT vs gfortran-16 (bmcnum.f)" begin
        # golden: scratchpad/wwpb/driver_bmcnum.f over pristine bmcnum.f+bmcspt.f,
        # gfortran-16. PBSPEC=1, no special trees, USERA=1, REPPHE→1, SSBATK=0,
        # GRFSTD=0.8, BAH(i)=2+1.5i, BANH=0.5, TREE=50, ISCMIN(1)=3. Float32 hex.
        hx(x) = reinterpret(UInt32, x)
        st = FVSjl.WwpbStand(); w = FVSjl.wwpb_defaults!(v)
        c = FVSjl.wwpb_init_coeffs(copy(FVSjl.WWPB_UPSIZ_DEFAULT))
        basum = 0.0f0
        for i in 1:10
            st.bah[i] = 2.0f0 + i*1.5f0; st.banh[i] = 0.5f0; st.tree[i,1] = 50.0f0
            basum += st.bah[i] + st.banh[i]
        end
        st.bah[11] = basum - 0.5f0*10; st.banh[11] = 0.5f0*10
        st.grfstd = 0.80f0; st.repphe = 0.0f0; st.ssbatk = 0.0f0
        FVSjl.bmcnum!(st, w, c; ipson = false)
        @test hx(st.numer[1]) == 0x46455D00   # (1)·BAIS·(BAHG)/GRFSTD, BAHG = Σ BAH[3:10]
        @test hx(st.tfood[1]) == 0x42A46405   # BAHG·(BAHG/BAIS)
        @test st.repphe == 1.0f0              # 0 → 1
    end

    @testset "BMATCT — single-stand BKP saturation BIT-EXACT vs gfortran-16 (bmatct.f)" begin
        # golden: scratchpad/wwpb/driver_bmatct.f over pristine bmatct.f (SPLAAR/
        # SPLALO/SPLADS/GPGET2 stubbed), gfortran-16, BMSTND=1, OUTOFF=T. The
        # redistribution collapses to BKP = TFOOD·(1−exp(−ALPHA·BKP/(TFOOD+1e-6))),
        # ALPHA=0.365 (SDD=0). BKP=100, TFOOD=50 ⇒ 50·(1−exp(−0.73)) ≈ 25.9.
        hx(x) = reinterpret(UInt32, x)
        st = FVSjl.WwpbStand(); w = FVSjl.wwpb_defaults!(v)
        st.bkp = 100.0f0; st.bkpips = 0.0f0; st.tfood[1] = 50.0f0; st.tfood[2] = 0.0f0
        FVSjl.bmatct_single!(st, w; sdd = 0.0f0, ipson = false)
        @test hx(st.bkp)    == 0x41CF3C85     # saturated BKP into the stand
        @test hx(st.bkpips) == 0x00000000     # NEWBKP(2)=0 ⇒ BKPIPS 0 (no Ips DV)
    end

    @testset "BMCBET — beta-dist weights BIT-EXACT vs gfortran-16 (bmcbet.f AS245+AS63)" begin
        # golden: scratchpad/wwpb/driver_bmcbet.f over the pristine wwpb/bmcbet.f
        # (self-contained: ALNGAM AS245 + BETAIN AS63 are defined IN the file),
        # gfortran-16. Two ABETA/range cases. Float32 hex (Z8.8).
        hx(x) = reinterpret(UInt32, x)
        b1 = FVSjl.bmcbet!(15.0f0, 1, 5)      # BKP>6 branch
        g1 = UInt32[0x2FEA2FEC,0x37342309,0x3B56FEA4,0x3E0CBEB3,0x3F5BF0D6,
                    0x3727C5AC,0x3727C5AC,0x3727C5AC,0x3727C5AC,0x3727C5AC]
        @test all(hx(b1[i]) == g1[i] for i in 1:10)
        b2 = FVSjl.bmcbet!(2.5f0, 1, 8)       # BKP≈3.6 branch, wider range
        g2 = UInt32[0x3C903FF1,0x3D93F002,0x3E05F726,0x3E3553F9,0x3E4D65DE,
                    0x3E45E894,0x3E17F108,0x3D75D110,0x3727C5AC,0x3727C5AC]
        @test all(hx(b2[i]) == g2[i] for i in 1:10)
        # AS245 ALNGAM sanity: log Γ(5) = log 24 (exercises the 4≤x<12 branch)
        @test isapprox(FVSjl._wwpb_alngam(5.0f0), log(24.0f0); atol = 1.0f-4)
    end

    @testset "BMISTD — stochastic kill kernel BIT-EXACT vs gfortran-16 (bmistd.f)" begin
        # golden: scratchpad/wwpb/driver_bmistd.f over pristine bmistd.f+bmcbet.f+
        # bmrann.f (SPLAAR stubbed), gfortran-16, seed 55329, SAREA=1, GRF=1,
        # TVOL=10, no special trees. The exact BMRANN stream + bmcbet! + kill
        # phases must all match. Float32 hex (Z8.8).
        hx(x) = reinterpret(UInt32, x)
        mk() = (st = FVSjl.WwpbStand(); w = FVSjl.wwpb_defaults!(v);
                c = FVSjl.wwpb_init_coeffs(copy(FVSjl.WWPB_UPSIZ_DEFAULT));
                for i in 1:10; st.grf[i] = 1.0f0; st.tvol[i,1] = 10.0f0; end;
                st.iqptyp[1] = 1; (st, w, c))
        # case A — individual-kill path (BKP=3, trees in 3 & 5) → 3 TPA in class 5
        stA, wA, cA = mk(); stA.tree[3,1] = 50.0f0; stA.tree[5,1] = 30.0f0; stA.bkp = 3.0f0
        FVSjl.bmistd!(stA, wA, cA; sarea = 1.0f0)
        @test hx(stA.pbkill[5]) == 0x40400001    # ≈3.0 TPA killed
        @test hx(stA.pitch[5])  == 0x3D088889
        @test hx(stA.bkp)       == 0x00000000
        @test hx(stA.final[2])  == 0x40A00000    # last kill in size class 5
        # case B — group-kill path (BKP=40, trees in 3/5/7) → ~19 TPA in class 7
        stB, wB, cB = mk(); stB.tree[3,1] = 50.0f0; stB.tree[5,1] = 30.0f0; stB.tree[7,1] = 40.0f0
        stB.bkp = 40.0f0
        FVSjl.bmistd!(stB, wB, cB; sarea = 1.0f0)
        @test hx(stB.pbkill[7]) == 0x41980000    # ≈19.0 TPA
        @test hx(stB.pitch[7])  == 0x3CCCCCCD
        @test hx(stB.bkp)       == 0x00000000
        @test hx(stB.final[2])  == 0x40E00000    # last kill in size class 7
    end

    @testset "BMKILL — BM→FVS mortality handback BIT-EXACT vs gfortran-16 (bmkill.f)" begin
        # golden: scratchpad/wwpb/driver_bmkill.f over pristine bmkill.f+bmdbhc.f
        # (OPBISR/SVMORT stubbed), gfortran-16. 2 LP trees, class 3, PROB=25, WK2=2,
        # OTPA(3,host)=50. Reconciled WK2 = MFAST+MBTL. Float32 hex.
        hx(x) = reinterpret(UInt32, x)
        alpha(sp) = sp == 7 ? "LP" : "OT"
        mktree() = (t = FVSjl.TreeList(10); t.species[1] = 7; t.species[2] = 7;
                    t.dbh[1] = 7.0f0; t.dbh[2] = 7.0f0; t.tpa[1] = 25.0f0; t.tpa[2] = 25.0f0;
                    t.n = 2; t)
        # case A — beetle kill 30/50 of class 3 ⇒ WK2 = min(30/50,1)·25 ≈ 14.99999
        stA = FVSjl.WwpbStand(); wA = FVSjl.wwpb_defaults!(v)
        stA.otpa[3,1] = 50.0f0; stA.tpbk[3,1,3] = 30.0f0
        wk2A = Float32[2.0, 2.0]; FVSjl.bmkill!(stA, wA, mktree(), wk2A, alpha)
        @test hx(wk2A[1]) == 0x416FFFFF && hx(wk2A[2]) == 0x416FFFFF
        # case B — TPBK 60 > OTPA 50 ⇒ MBTL clamps to PROB, then PROB−WK2<1e-6 bound
        stB = FVSjl.WwpbStand(); wB = FVSjl.wwpb_defaults!(v)
        stB.otpa[3,1] = 50.0f0; stB.tpbk[3,1,3] = 60.0f0
        wk2B = Float32[2.0, 2.0]; FVSjl.bmkill!(stB, wB, mktree(), wk2B, alpha)
        @test hx(wk2B[1]) == 0x41C7FFFF     # 25 − 1e-6 (PROB bound)
    end

    @testset "WWPB outbreak orchestration — end-to-end composition (behavioral)" begin
        # The reconstructed single-stand PPMAIN/BMDRV loop composing the ten
        # bit-exact kernels. NOT bit-exact (the outer harness source is absent) —
        # asserts the composition is COHERENT: seeded inventory damage drives a
        # multi-year outbreak that produces beetle mortality scaling with the seed.
        c = FVSjl.wwpb_init_coeffs(copy(FVSjl.WWPB_UPSIZ_DEFAULT))
        alpha(sp) = sp == 7 ? "LP" : "OT"   # LP = lodgepole, the MPB host
        mktree() = begin
            t = FVSjl.TreeList(100)
            p(i, sp, d, h, cr, htg, cfv, tp) = (t.species[i] = Int32(sp); t.dbh[i] = Float32(d);
                t.height[i] = Float32(h); t.crown_pct[i] = Int32(cr); t.ht_growth[i] = Float32(htg);
                t.cuft_vol[i] = Float32(cfv); t.tpa[i] = Float32(tp))
            p(1, 7, 7.0, 45, 50, 1.0, 8.0, 60); p(2, 7, 13.0, 65, 45, 0.9, 25.0, 40)
            p(3, 7, 19.0, 80, 40, 0.7, 60.0, 25); t.n = 3; t
        end
        runkill(seedtpa) = begin
            st = FVSjl.WwpbStand(); w = FVSjl.wwpb_defaults!(v); t = mktree()
            seed = zeros(Float32, 10); seed[3] = seedtpa
            FVSjl.wwpb_outbreak_cycle!(st, w, c, t, alpha; sarea = 5.0f0, iyr1 = 1, iyr2 = 5,
                                       seed_pbkill = seed)
            wk2 = zeros(Float32, 3); FVSjl.bmkill!(st, w, t, wk2, alpha); sum(wk2)
        end
        m0, m8, m20 = runkill(0.0f0), runkill(8.0f0), runkill(20.0f0)
        @test m0 == 0.0f0            # no seed ⇒ no outbreak (BKP starts at 0)
        @test m8 > 0.0f0             # seeded damage drives real beetle mortality
        @test m20 > m8               # bigger seed ⇒ bigger outbreak (monotone)
        @test m20 < 125.0f0          # bounded by the stand (never over-kills)
    end
end
