# test_tt_regent_dg_blend.jl — TT (Teton) regent small/large-tree coupling regression guard.
#
# Guards the coupled M331D woodland over-growth fix (2026-09-03):
#  (1) TT_RG_XMIN/XMAX must be the VERBATIM tt/regent.f:161-171 DATA (1.5/3.0 for the TTVAR conifers),
#      NOT the old mis-read 2.0/4.0.
#  (2) The DIAMETER increment is NOT XWT-blended: FVS regent.f blends only HEIGHT (regent.f:799); for
#      D<BKPT the regent small-tree DG fully replaces the large-tree dgf DG. The invariant that makes the
#      default-conifer DG pure-regent is BKPT == XMAX == 3.0 (so within the blend loop, d<XMAX ⇒ d<BKPT).
#      A regression that re-introduced the DG blend, or moved XMAX off BKPT, over-grew the 1.5–3" woodland
#      trees (measured 40-stand M331D total-rel: density-fix-only 9.84 → XMIN-alone 17.44 → combined 5.79).
using Test
using FVSjl

@testset "TT regent XMIN/XMAX DATA + DG-not-blended coupling" begin
    XMIN = FVSjl.TT_RG_XMIN
    XMAX = FVSjl.TT_RG_XMAX
    BREAK = FVSjl.TT_RG_BREAK

    # (1) Verbatim tt/regent.f DATA XMIN / XMAX.
    @test XMIN == Float32[1.5, 1.5, 1.5, 90.0, 1.5, 1.5, 1.5, 1.5, 1.5, 2.0,
                          90.0, 90.0, 90.0, 2.0, 0.5, 90.0, 1.5, 0.5]
    @test XMAX == Float32[3.0, 3.0, 3.0, 99.0, 3.0, 3.0, 3.0, 3.0, 3.0, 5.0,
                          99.0, 99.0, 99.0, 4.0, 2.0, 99.0, 3.0, 2.0]

    # The TTVAR default conifers must NOT carry the old 2.0/4.0 window.
    for sp in (1, 2, 3, 5, 7, 8, 9, 17)
        @test XMIN[sp] == 1.5f0
        @test XMAX[sp] == 3.0f0
    end

    # (2) DG-purity invariant. The blend loop `continue`s at d≥XMAX and the DG is set from the regent value
    # only when d<BKPT. For the default conifers BKPT == XMAX == 3.0, so every tree the loop sees (d<XMAX)
    # also has d<BKPT ⇒ pure regent DG, never a blend of the large-tree dgf DG.
    for sp in (1, 2, 3, 5, 6, 7, 8, 9, 17)
        @test BREAK[sp] == 3.0f0
        @test XMAX[sp] == 3.0f0
    end
    # MM (14) is the one default species whose window is wider than its break (XMAX=4 > BKPT=3): a
    # 3–4" MM tree enters the loop but keeps the large-tree dgf DG (d≥BKPT), so the `d < BKPT` gate is
    # load-bearing — it must not be simplified back to an unconditional overwrite.
    @test BREAK[14] == 3.0f0
    @test XMAX[14] == 4.0f0


    # ── Firing test: BI/MC (13,16) inventory Curtis-Arney HT-DBH (regent.f:874-901) ──
    # jl formerly used the 0.1·HTG rule-of-thumb; but LHTDRG(13)=LHTDRG(16)=.FALSE. means FVS ALWAYS
    # takes the inventory equation `DG=(DK−DKK)·bark`. _tt_bimc_dk maps a grown height to the predicted
    # DBH. Pinned to the FVStt_dbg-measured branch (dense Gambel-oak stand 51031230020004 cyc1): a BI
    # tree grown to HK≈4.84 gets DK≈0.335 (oracle 0.3324 at its ZZRAN'd HK=4.816), NOT 0.1·HTG≈0.383.
    # Deleting the inventory branch (reverting to 0.1·HTG) changes this value ⇒ the test is not vacuous.
    let
        # BI(13): P2/P3/P4 = 76.5170 / 2.2107 / −0.6365
        p2, p3, p4 = 76.5170f0, 2.2107f0, -0.6365f0
        hat3 = 4.5f0 + p2 * FVSjl.fexp(-p3 * FVSjl.fpow(3.0f0, p4))
        @test hat3 == 30.005816f0
        @test FVSjl._tt_bimc_dk(4.8407f0, hat3, p2, p3, p4) == 0.33502105f0   # HK<HAT3 rational form
        # MC(16): P2/P3/P4 = 1709.7229 / 5.8887 / −0.2286
        q2, q3, q4 = 1709.7229f0, 5.8887f0, -0.2286f0
        hat3m = 4.5f0 + q2 * FVSjl.fexp(-q3 * FVSjl.fpow(3.0f0, q4))
        @test hat3m == 22.017498f0
        @test FVSjl._tt_bimc_dk(10.0f0, hat3m, q2, q3, q4) == 1.1466658f0
        # DG must be (DK−DKK)·bark (DKK=D since original H≤4.5), NOT 0.1·HTG. For DK=0.335, DKK=0.1,
        # bark≈0.99: DG ≈ 0.233 — distinctly below the old rule-of-thumb 0.1·3.83 ≈ 0.383.
        dg_inventory = (0.33502105f0 - 0.1f0) * 0.99f0
        @test 0.22f0 < dg_inventory < 0.24f0
        @test dg_inventory < 0.1f0 * 3.83f0   # strictly less than the removed 0.1·HTG over-grow
    end

    # ── Firing test: TT aspen(6)/MM(14) REGENT small-tree HEIGHT self-calibration CON = exp(HCOR) ──
    # tt_regent_hcor_aspen_init! (tt/regent.f:1097-1362): CORNEW = Σ(HTG·SCALE3·P)/Σ(EDH·P), EDH = NPER·SMHTGF,
    # HCOR = ln(CORNEW). small_tree_growth! then multiplies the aspen/MM height increment by CON = exp(HCOR).
    # Pinned to FVStt_g16 stand 325585226489998 cyc1 (M331D pure-aspen woodland): the 5 sub-5" aspen records
    # carrying measured HTG are HT=[13,12,17,14,16], HTG=[4,2,1,2,2], SITEAR(6)=34 (⇒ RSIMOD=0.52857),
    # FINTH=10 ⇒ SCALE3=10/FINTH=1, NPER=2. Oracle: SNX=6.5179 (mean EDH), CORNEW=0.337531, HCOR=−1.08610,
    # and after the dgdriv WCI attenuation (WCI=0.0797, CORMLT=exp(−0.2773)=0.75797) the growth-cycle
    # CON = 0.4476. Without this the aspen HEIGHT increment ran CON=1 ⇒ ~3.6× over-growth (i19 htg 1.69 vs
    # 0.756=1.69·0.4476). jl takes ZRAND=0 in the calibration (the SMHTGF ZRAND draw is a separate RNG-stream
    # item), landing CORNEW/CON within ~0.7% of the oracle.
    let
        HT  = Float32[13, 12, 17, 14, 16]     # SMHTGF reads HT(I) (current height)
        HTG = Float32[4, 2, 1, 2, 2]          # measured height increment
        si6 = 34.0f0
        nper = 2f0; scale3 = 1.0f0            # 10/FINTH with FINTH=10
        snx = 0f0; sny = 0f0
        for k in 1:5
            # aspen SMHTGF (sp6) folds in the ·0.75 and RSIMOD; TPCCF/CR are unused for the aspen closed form.
            edh = nper * FVSjl._tt_smhtgf(6, HT[k], 25f0, 100f0, 0f0, si6)
            snx += edh; sny += HTG[k] * scale3
        end
        # per-subcycle EDH for the first tree (HTGR·0.75·RSIMOD), the load-bearing SMHTGF kernel.
        @test FVSjl._tt_smhtgf(6, 13.0f0, 25f0, 100f0, 0f0, si6) ≈ 3.1993902f0 rtol=1f-5
        cornew = sny / snx
        hcor_init = log(cornew)
        @test cornew ≈ 0.3398311f0 rtol = 1f-5           # jl CORNEW (regression pin)
        @test abs(cornew - 0.337531f0) < 0.004f0          # within ~0.7% of the FVStt_g16 oracle CORNEW
        # dgdriv WCI attenuation → the applied growth-cycle CON (dgdriv.f:213). With the oracle WCI the cyc1
        # CON lands at the oracle 0.4476 to within the ZRAND=0 residual.
        cormlt_h = exp(-0.02773f0 * 10f0)                 # exp(−0.2773) = 0.75797
        con_cyc1 = exp(0.0797f0 + cormlt_h * (hcor_init - 0.0797f0))
        @test abs(con_cyc1 - 0.4476f0) < 0.004f0          # CON matches the oracle 0.4476 (±0.7%)
        # CON must be the dominant reduction, not ≈1 (the pre-fix over-growth): a >2× height brake.
        @test con_cyc1 < 0.5f0
    end

    # ── Firing test: MM(14) UTVAR FINDAG height + Wykoff DBH + XWT blend (M331D sparse-woodland fix) ──
    # MM (Rocky Mtn maple, FIA 321) is UTVAR (regent.f:396 CASE 4,11:16,18), NOT the aspen/TTVAR small-tree
    # path. jl formerly mis-routed it through the default+aspen path (SMHTGF CASE(6) +5-yr Sheppard SUBCYCLED
    # twice, no ·0.75, no RSIMOD) which over-grew MM height ~1.5× → compounded 2.2–2.5× into stand BA on the
    # sparse-woodland stands (2766510 O:87 vs jl:153; 2764247 O:159 vs 246). FVS MM = FINDAG aspen-height
    # (AG2=SITAGE+10, ·RSIMOD·0.75, applied ONCE) + HTDBH Wykoff H-D DBH (acd/htdbh.f:458 D=HT2/(ln(H-4.5)-HT1)-1).
    # After the fix + the XWT height blend (regent.f:793-799, XMIN(14)=2/XMAX(14)=4): 2766510→90, 2764247→160.
    @testset "TT MM(14) UTVAR FINDAG + XWT-blend routing" begin
        # (1) MM(14) routes through the UTVAR pass, NOT the default/aspen path (reverting step-1 flips these).
        @test FVSjl._tt_rg_utvar(14) == true
        @test FVSjl._tt_rg_default(14) == false

        # (2) FINDAG height applied ONCE + XWT blend. Args:
        #   (sp, h, d, cr, sitear, pctred, con, bark, dgmax, diam, scale2, htg_large)
        # con=1, sitear=30 ⇒ RSIMOD=0.5·(1+(30−5)/25)=1.0, h=10 ⇒ HTGRL=(H(sitage+10)−H(sitage))/(30.48)·0.75.
        bark = FVSjl.tt_bratio(14, 1.0f0)
        @test bark == 0.95f0
        call(d, htgL) = FVSjl._tt_utvar_regent(14, 10.0f0, d, 25.0f0, 30.0f0, 1.0f0, 1.0f0,
                                               bark, 2.5f0, 0.1f0, 1.0f0, htgL)

        # FINDAG height at d=1 (XWT=0 ⇒ pure small-tree, blend inert). Regression pin: 12.1616 ft — NOT the
        # subcycled default+aspen value the pre-fix produced (~1.5× larger). Applied ONCE (·0.75, RSIMOD=1).
        h_d1, dg_d1 = call(1.0f0, 0.0f0)
        @test h_d1 ≈ 12.161556f0 rtol = 1f-5
        # Wykoff DBH increment is non-zero (guards the HTDBH branch vs a DG=0 seedling stub).
        @test dg_d1 ≈ 1.6859541f0 rtol = 1f-5

        # d≤2 ⇒ XWT=0 ⇒ height increment is INDEPENDENT of the large-tree HTG (seedling cycles stay bit-exact).
        @test call(1.0f0, 0.0f0)[1] == call(1.0f0, 99.0f0)[1]

        # d=3 ⇒ XWT=(3−2)/(4−2)=0.5 ⇒ height increment blends HALF-way toward the large-tree HTG.
        # (a) With htg_large=0 the blend halves the pure small-tree value (load-bearing: the pre-fix omitted
        #     the blend, so MM height over-grew once D crossed XMIN=2).
        @test call(3.0f0, 0.0f0)[1] ≈ 0.5f0 * h_d1 rtol = 1f-5
        # (b) A Δhtg_large of 4 shifts the blended increment by exactly XWT·4 = 2.0.
        @test (call(3.0f0, 6.0f0)[1] - call(3.0f0, 2.0f0)[1]) ≈ 2.0f0 rtol = 1f-4
    end
end