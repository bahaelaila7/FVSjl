# White Pine Blister Rust (WPBR) canker model — keyword reader (brin.f, "BRUST"
# keyword, keywds.f option 75) + the BRANN Lehmer/MINSTD RNG + the BRINIT/BRBLKD
# defaults (chunk 0) PLUS the canker generation/growth/status/mortality MATH
# kernels (dynamics chunk). Engine-inert: no simulate.jl seam wires WPBR into a
# projection, so a BRUST-present stand projects BYTE-IDENTICALLY to one without.
#
# VALIDATED HERE (bit-exact, 0 ULP unless noted):
#   * BRANN  — the pristine brann.f MINSTD LCG (seed 55329) → 6-draw stream, and
#     the BRNSED reseed semantics (even→forced-odd, LSET=false reset-to-SS).
#     Goldens from scratchpad/wpbr/driver_brann.f (gfortran-16 over the exact
#     brann.f algorithm). The stream is IDENTICAL to DFB's DFBRAN and DFTM's
#     TMRANN (same 16807/2^31-1 generator, same 55329 seed) — measured, not
#     assumed. Never FFI'd — Julia's own wpbr_rand!/wpbr_seed!.
#   * BRINIT/BRBLKD defaults + the variant BRSPM host-species map (IE/SO/CR/base).
#   * The brin.f keyword reader populating s.wpbr (parse-correctness).
#   * The INERT-seam A/B: a full BRUST block ⇒ .sum byte-identical to no BRUST.
#   * DYNAMICS kernels — dump-replay vs the relinked+g16-instrumented FVSie_wpbr
#     oracle (HEX goldens, so 0 ULP): BRSETP BRGD/BRHTBC init, BRGI/BRSTAR GI+
#     TSTARG (GI ≤1 ULP on one libm-powf edge, report-only), BRTARG RI + BRIBA,
#     BRECAN RITEM/TNEWC/PLI/NUMTIM + canker placement TOUT/PLETH, BRCGRO bole
#     girdle growth + status→kill. Fed the oracle's own hex inputs (equal-inputs
#     dump-replay); the end-to-end .sum-DELTA is cornered by the IE #206 straddle.
#
# DOCTRINE: WPBR ships nowhere (every FVS*_buildDir links the base/exbrus.f no-op
# stub), so a numeric oracle is a RELINK (wpbr/*.o swapped for exbrus.o — recipe
# in scratchpad/wpbr/build_ie_wpbr.sh, FVSie_wpbr). For chunk 0 the validated
# numeric surface (the RNG) is self-contained, so the goldens come from a tiny
# standalone gfortran-16 driver over the PRISTINE brann.f.

using Test
using FVSjl
const _FW = FVSjl

_hexw(x::Float32) = uppercase(string(reinterpret(UInt32, x); base = 16, pad = 8))

# --- goldens from scratchpad/wpbr/driver_brann.f over pristine brann.f ---
const _G_BRANN   = ("3EDDB57A", "3F5AB21B", "3F630A07", "3F2734C9", "3EF4FB2C", "3F4AF646")
const _G_SEED101 = "3A4F3718"      # even 100 → forced odd 101, first draw
const _G_SEED55A = "39E1AE10"      # odd 55, first draw
const _G_SEED55B = "3E70354E"      # odd 55, second draw

# Right-justify each numeric field in a 10-col slot after the (≤10-col) keyword —
# the FVS fixed-field keyword layout the reader expects.
_kw(name, vals...) = rpad(name, 10) * join(lpad(string(v), 10) for v in vals)

# IE host stand: a few western white pine (WP) records + DF, for the inert-seam A/B.
const _WPBR_TRE = """
   1      248112       0101   011WP 12014   0654   00111     0  0
   2      248112       0101   011WP 08010   0504   00111     0  0
   3      248112       0102   011DF 11014   0634   00111     0  0
   4      248112       0102   011WP 15019   0734   00111     0  0
   5      248112       0103   011DF 09011   0524   00111     0  0
"""

_wpbr_head(title) = """
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

# A BRUST block exercising a spread of sub-keywords across the BRCOM commons.
const _WPBR_BLOCK = join([
    "BRUST",
    _kw("BROUT", 1, 1),                        # lbrsum on, lbrdbh on
    _kw("INACT", 0, 0.08, 0.03),               # ratinv (branch,bole)
    _kw("BRSEED", 101),                         # reseed 101 (odd)
    _kw("RUSTINDX", 0.025, "", 2, 0.2, 0.3, 0.5),   # ridef, rimeth=2, ribprp
    _kw("GROWRATE", 1, 1, 6.0, 5.5),           # WP stock1: bogrth,brgrth
    _kw("DEVFACT", 0, 1, 0.5, 0.5, 0.5, 0.5),  # WP dev factors
    _kw("STOCK", 0, 1, 2, 0.3, 0.4),           # WP stock-type 2 mix
    _kw("PRNSPECS", 0, 0.60, 7.0, 5.0, 20.0),  # prune thresholds
    _kw("EXSPECS", 0, 4.0, 5.0, 40.0, 95.0, 2.0),   # excise thresholds
    "END",
], "\n") * "\n"

# A BRUST block that only names the extension (all defaults) — inert-seam A/B.
const _WPBR_MIN = "BRUST\nEND\n"

@testset "White Pine Blister Rust (WPBR) — keyword reader + BRANN RNG + defaults" begin

    @testset "BRANN Lehmer stream (bit-exact vs pristine brann.f, seed 55329)" begin
        w = _FW.wpbr_defaults!(_FW.InlandEmpire())
        for g in _G_BRANN
            @test _hexw(_FW.wpbr_rand!(w)) == g
        end
        # Identical generator to DFB/DFTM — the shared 16807/2^31-1 MINSTD, seed 55329.
        w2 = _FW.wpbr_defaults!(_FW.InlandEmpire())
        d2 = _FW.dftm_defaults!(_FW.InlandEmpire())
        @test all(_hexw(_FW.wpbr_rand!(w2)) == _hexw(_FW.dftm_rand!(d2)) for _ in 1:6)
    end

    @testset "BRNSED reseed (odd forcing + default reset)" begin
        w = _FW.wpbr_defaults!(_FW.InlandEmpire())
        _FW.wpbr_seed!(w, 100.0f0, true)          # even → forced odd 101
        @test _hexw(_FW.wpbr_rand!(w)) == _G_SEED101
        _FW.wpbr_seed!(w, 55.0f0, true)           # already odd
        @test _hexw(_FW.wpbr_rand!(w)) == _G_SEED55A
        @test _hexw(_FW.wpbr_rand!(w)) == _G_SEED55B
        _FW.wpbr_seed!(w, 0.0f0, false)           # LSET=false: reset BRS0 to SS (55)
        @test _hexw(_FW.wpbr_rand!(w)) == _G_SEED55A
    end

    @testset "BRINIT/BRBLKD defaults + variant BRSPM crosswalk" begin
        w = _FW.wpbr_defaults!(_FW.InlandEmpire())
        @test w.ibrdam == 36 && w.brpi == 3.14159f0
        @test w.brspc == ("WP", "LM") && w.brspm[1] == 1 && w.brspm[13] == 2
        @test w.srate == (0.9f0, 0.5f0) && w.htprpr == 0.50f0
        @test w.exdmin == 3.0f0 * 2.54f0 && w.htmin == 3.0f0 * 2.54f0
        @test w.htmax == (8.0f0 * 30.48f0, 6.0f0 * 30.48f0)
        @test w.girmax == 50.0f0 && w.girmrt == 100.0f0
        @test w.outdst == 6.0f0 * 2.54f0 && w.outnld == 24.0f0 * 2.54f0
        @test w.ratinv == (0.05f0, 0.01f0)
        @test w.gidef == 15.0f0 && w.ridef == 0.015f0 && w.rimeth == 0
        @test w.ribprp == (0.0f0, 0.5f0, 0.5f0)
        @test w.ribus[1, 2] == 200.0f0 && w.ribus[2, 2] == 75.0f0
        @test w.rsf == (2.3f0, 1.0f0, 0.64f0)
        @test all(w.dfact .== 0.33f0) && all(w.riaf .== 1.0f0)
        @test all(w.brgrth .== 5.0f0) && all(w.bogrth .== 4.5f0)
        @test w.resist[1, 1] == 1.0f0 && w.resist[1, 4] == 0.11f0
        @test w.prpstk[1, 1] == 1.0f0 && w.prpstk[1, 2] == 0.0f0
        @test w.icin == 55 && w.idtout == 56 && w.idcout == 57
        @test w.lbrsum == true && w.lbrdbh == false && !w.active
        # variant host maps
        so = _FW.wpbr_defaults!(_FW.SouthCentralOregon())
        @test so.brspc == ("WP", "SP") && so.brspm[1] == 1 && so.brspm[2] == 2
        cr = _FW.wpbr_defaults!(_FW.CentralRockies())
        @test cr.brspc == ("LP", "PP") && cr.brspm[11] == 1 && cr.brspm[13] == 2
    end

    dir = mktempdir()
    write(joinpath(dir, "shared.tre"), _WPBR_TRE)
    v = FVSjl.InlandEmpire()
    wpbr_key = joinpath(dir, "wpbr.key")
    write(wpbr_key, _wpbr_head("WPBR PARSE ") * _WPBR_BLOCK * "ECHOSUM\nPROCESS\nSTOP\n")
    cp(joinpath(dir, "shared.tre"), joinpath(dir, "wpbr.tre"))

    @testset "BRUST keyword reader populates s.wpbr (brin.f)" begin
        w = nothing
        for s in FVSjl.each_stand(wpbr_key; variant = v)
            w = s.wpbr; break
        end
        @test w !== nothing
        @test w.active == true                    # BRYES
        @test w.lbrdbh == true                    # BROUT field 2 = 1
        @test w.ratinv == (0.08f0, 0.03f0)        # INACT
        @test w.ridef == 0.025f0 && w.rimeth == 2 # RUSTINDX
        @test w.ribprp == (0.2f0, 0.3f0, 0.5f0)   # ribes proportions (sum 1.0)
        @test w.bogrth[1, 1] == 6.0f0 && w.brgrth[1, 1] == 5.5f0   # GROWRATE (WP,stock1)
        @test w.ldfact == true && all(w.dfact[1, :] .== 0.5f0)     # DEVFACT (BR host 1)
        @test w.prpstk[1, 2] == 0.3f0 && w.resist[1, 2] == 0.4f0   # STOCK
        @test w.htprpr == 0.60f0 && w.outdst == 5.0f0 * 2.54f0     # PRNSPECS
        @test w.exdmin == 4.0f0 * 2.54f0 && w.girmax == 40.0f0     # EXSPECS
        @test w.girmrt == 95.0f0 && w.htmin == 2.0f0 * 2.54f0
        # reseed via BRSEED 101 leaves BRS0=BRSS=101 (odd) — first draw matches golden
        @test _hexw(_FW.wpbr_rand!(w)) == _G_SEED101
    end

    # -------------------------------------------------------------------------
    # INERT seam: a BRUST-present stand must project byte-identically to a stand
    # with no BRUST keyword (no per-cycle engine seam is wired yet).
    # -------------------------------------------------------------------------
    ctrl_key = joinpath(dir, "ctrl.key")
    min_key  = joinpath(dir, "min.key")
    write(ctrl_key, _wpbr_head("WPBR CTRL  ") * "ECHOSUM\nPROCESS\nSTOP\n")
    write(min_key,  _wpbr_head("WPBR MIN   ") * _WPBR_MIN * "ECHOSUM\nPROCESS\nSTOP\n")
    write(wpbr_key, _wpbr_head("WPBR FULL  ") * _WPBR_BLOCK * "ECHOSUM\nPROCESS\nSTOP\n")
    for k in ("ctrl", "min", "wpbr")
        cp(joinpath(dir, "shared.tre"), joinpath(dir, "$k.tre"); force = true)
    end
    rows(key) = filter(l -> !startswith(l, "-999"),
                       split(strip(FVSjl.run_keyfile(key; variant = v, output = :sum)), '\n'))

    @testset "WPBR seam is INERT (.sum byte-identical with/without a BRUST block)" begin
        base = rows(ctrl_key)
        @test rows(min_key)  == base       # BRUST/END only ⇒ byte-identical
        @test rows(wpbr_key) == base       # full BRUST block ⇒ still byte-identical (inert seam)
    end
end

# =============================================================================
# WPBR dynamics kernels — dump-replay BIT-EXACT vs the relinked (+g16-instrumented)
# FVSie_wpbr oracle on the IE goldens stand (S248112, RUSTINDX 0.05, 3 WP hosts +
# 2 DF). Goldens are the oracle's REAL*4 values printed in HEX (Z8 edit descriptor,
# scratchpad/wpbr/instr2.py → unit 66; the instrumented .sum stays byte-identical
# to the clean relink), so every golden round-trips to the exact IEEE Float32 bit
# pattern (an F-format decimal is under-precise for large magnitudes → 1-ULP
# artifacts). Inputs are fed to the Julia kernels from the SAME oracle hex, so this
# validates the deterministic init + index + canker generation/growth/status math
# on EQUAL INPUTS, 0 ULP. The end-to-end .sum-DELTA is CORNERED by the IE #206
# growth straddle (MEASURED — FVSjl-off already ≠ FVSie_wpbr-off on this stand), so
# the dump-replay here is the doctrine-valid check; no engine seam is wired (inert).
# =============================================================================
_fhw(h) = reinterpret(Float32, parse(UInt32, h; base = 16))   # hex string → Float32
@testset "WPBR dynamics kernels — dump-replay bit-exact HEX (FVSie_wpbr, 0 ULP)" begin

    # ULP distance between two Float32 (for the one report-only libm-powf edge).
    _ulps(a, b) = abs(Int(reinterpret(Int32, Float32(a))) - Int(reinterpret(Int32, Float32(b))))

    @testset "BRSETP per-tree init (BRGD ground-diam + BRHTBC crown-base, cm)" begin
        # HSETP goldens: (HT, DBH, ICR) → hex(BRGD), hex(BRHTBC)
        @test _hexw(_FW.wpbr_brgd(65, 12)) == "42015D1B"
        @test _hexw(_FW.wpbr_brgd(50,  8)) == "41AFB401"
        @test _hexw(_FW.wpbr_brgd(73, 15)) == "4220A13F"
        @test _hexw(_FW.wpbr_brhtbc(65, 35)) == "44A0F8F6"
        @test _hexw(_FW.wpbr_brhtbc(50, 35)) == "4477A667"
        @test _hexw(_FW.wpbr_brhtbc(73, 35)) == "44B4C8D6"
        # BRGD floored at BRDBH only when BRHT<1.14 m (HT<3.74 ft) makes the ratio blow up
        @test _FW.wpbr_brgd(3, 20) == Float32(20 * 2.54f0)
    end

    @testset "BRGI growth index + BRSTAR sum-target (GI, TSTARG)" begin
        # HTARG cyc0 goldens: iiag=60 (stand age), hht = HT·0.3048.
        # TSTARG is bit-exact; GI matches ≤1 ULP — the McDonald `base**(-2.071822)`
        # is an openlibm-vs-glibc powf 1-ULP edge (r4). GI is REPORT-ONLY (AVGGI),
        # NOT in the BRECAN/BRCGRO mortality path, so it never reaches the .sum.
        g1, t1 = _FW.wpbr_brgi(60, Float32(65 * 0.3048f0))
        @test _hexw(g1) == "41820533" && _hexw(t1) == "463995D7"
        g2, t2 = _FW.wpbr_brgi(60, Float32(50 * 0.3048f0))
        @test _hexw(g2) == "4173D70A" && _hexw(t2) == "462A0637"
        g4, t4 = _FW.wpbr_brgi(60, Float32(73 * 0.3048f0))
        @test _ulps(g4, _fhw("41921034")) <= 1 && _hexw(t4) == "4659855B"
        @test _FW.wpbr_brstar(19.812f0) > 0.0f0
        # GI clamps: old + short → 15.24 floor; young + tall → 38.10 cap
        gs, _ = _FW.wpbr_brgi(400, 0.1f0); @test gs == 15.24f0
        gt, _ = _FW.wpbr_brgi(2, 60f0);    @test gt == 38.10f0
    end

    @testset "BRTARG per-tree RI + BRIBA basal-area rust index" begin
        # RI = RIDEF·RESIST(sp,stock)·RIAF; golden path RI = 0.05·1·1
        @test _FW.wpbr_ri(0.05f0, 1.0f0, 1.0f0) == 0.05f0
        @test _FW.wpbr_ri(0.05f0, 0.33f0, 1.0f0) == Float32(0.05f0 * 0.33f0)
        # HIBA cyc0 golden: BA=0x41A00001, RIBPRP=(0.2,0.3,0.5), RSF=(2.3,1,0.64) → RIDEF
        @test _hexw(_FW.wpbr_briba(_fhw("41A00001"), (0.2f0, 0.3f0, 0.5f0),
                                   (2.3f0, 1.0f0, 0.64f0))) == "3BEF9CBE"
    end

    @testset "BRECAN expected-canker probabilities (RITEM, TNEWC, PLI, NUMTIM)" begin
        # HECAN cyc1 trees 1,2,4 (inputs HITE/SSTAR from oracle hex), DFACT=0.33, RI=0.05
        r1 = _FW.wpbr_brecan_probs(_fhw("41BF2B6E"), 0.05f0, _fhw("44042CF7"), 0.33f0)
        @test _hexw(r1[1]) == "3C234CED" && _hexw(r1[2]) == "40A8A0B2"
        @test _hexw(r1[3]) == "3F5A9DBB" && r1[4] == 6
        r2 = _FW.wpbr_brecan_probs(_fhw("41A2A13D"), 0.05f0, _fhw("43BBFBC6"), 0.33f0)
        @test _hexw(r2[1]) == "3CD52967" && _hexw(r2[2]) == "411C86E3"
        @test _hexw(r2[3]) == "3F66AE78" && r2[4] == 10
        r4 = _FW.wpbr_brecan_probs(_fhw("41D2B193"), 0.05f0, _fhw("441A61AF"), 0.33f0)  # >25 m ⇒ RI·0.1
        @test _hexw(r4[1]) == "3BA3D70B" && _hexw(r4[2]) == "40459BC2"
        @test _hexw(r4[3]) == "3F4887E2" && r4[4] == 4
    end

    @testset "BRECAN canker placement (TOUT distance-out + PLETH lethality)" begin
        # HPLAC cyc1 tree1 canker1: SSTHT/CRLEN/TUP from oracle hex → TOUT, PLETH
        tout, pleth = _FW._wpbr_tout_pleth(_fhw("41A1C36C"), _fhw("4430EDC0"), _fhw("44EFDD40"))
        @test _hexw(tout) == "41B77909" && _hexw(pleth) == "3F1B8E61"
        # far-out canker (2nd HPLAC line, power-law PLETH branch) stays ≥0
        _, pl2 = _FW._wpbr_tout_pleth(_fhw("41A1C36C"), _fhw("4430EDC0"), _fhw("44CE9E86"))
        @test _hexw(pl2) == "3D0C601D" && pl2 >= 0.0f0
    end

    @testset "BRCGRO bole-canker girdle growth + status→kill" begin
        # HCGBO record2 cyc1 first grown year: GIRAMT/GROBOL from oracle hex, GIRD 0→
        @test _hexw(_FW._wpbr_girdle_update(0.0f0, _fhw("40CC3238"), _fhw("400B9771"))) == "4208B917"
        # cap at 100
        @test _FW._wpbr_girdle_update(90.0f0, 6.38113f0, 6.38113f0) == 100.0f0
        # status classification (shared BRCGRO/BRCSTA): a fully-girdled bole canker
        # below the base of crown kills the tree (code 7 — the CGKIL golden event).
        th = (exht = 6.0f0 * 30.48f0, htmin = 3.0f0 * 2.54f0, exdmin = 3.0f0 * 2.54f0,
              girmax = 50.0f0, girmrt = 100.0f0, outnld = 24.0f0 * 2.54f0, outdst = 6.0f0 * 2.54f0)
        # gird=100 (≥GIRMRT), up below crown (≤htbcr) ⇒ 7 kill
        @test _FW.wpbr_canker_status(0.0f0, 200.0f0, 100.0f0, 30.0f0;
              exht = th.exht, htmin = th.htmin, exdmin = th.exdmin, girmax = th.girmax,
              girmrt = th.girmrt, htbcr = 400.0f0, phtst = 300.0f0, outnld = th.outnld, outdst = th.outdst) == 7
        # gird=100 but up above crown ⇒ 5 top-kill
        @test _FW.wpbr_canker_status(0.0f0, 500.0f0, 100.0f0, 30.0f0;
              exht = th.exht, htmin = th.htmin, exdmin = th.exdmin, girmax = th.girmax,
              girmrt = th.girmrt, htbcr = 400.0f0, phtst = 300.0f0, outnld = th.outnld, outdst = th.outdst) == 5
        # small bole canker (gird≤GIRMAX, in excise window) ⇒ 3 excisable
        @test _FW.wpbr_canker_status(0.0f0, 100.0f0, 10.0f0, 30.0f0;
              exht = th.exht, htmin = th.htmin, exdmin = th.exdmin, girmax = th.girmax,
              girmrt = th.girmrt, htbcr = 400.0f0, phtst = 300.0f0, outnld = th.outnld, outdst = th.outdst) == 3
        # far-out branch canker ⇒ 1 non-lethal
        @test _FW.wpbr_canker_status(200.0f0, 100.0f0, 0.0f0, 30.0f0;
              exht = th.exht, htmin = th.htmin, exdmin = th.exdmin, girmax = th.girmax,
              girmrt = th.girmrt, htbcr = 400.0f0, phtst = 300.0f0, outnld = th.outnld, outdst = th.outdst) == 1
    end

    @testset "BRANN stream starts clean at the first BRECAN draw (seed 55329)" begin
        # BRINIT resets BRS0→BRSS=55329; BRSTYP (default stock) + BRCREM (no prune/
        # excise) draw NOTHING, so the run's first BRANN is BRECAN tree1 year1.
        # Oracle ECAND cyc1 tree1 J=1 = 0.433025181 = first draw from 55329 = 3EDDB57A.
        w = _FW.wpbr_defaults!(_FW.InlandEmpire())
        @test _hexw(_FW.wpbr_rand!(w)) == _G_BRANN[1]
        @test _FW.wpbr_rand!(w) == _fhw("3F5AB21B")   # 2nd draw = ECAND canker up-position
    end
end
