# =============================================================================
# simulate.jl — the projection cycle loop (FVS / TREGRO orchestration)
#
# Ported from: base/fvs.f (cycle loop) + base/tregro.f (per-cycle driver).
#
# After initialization, FVS projects the stand cycle by cycle: recompute stand
# density, grow diameters and heights, apply mortality and regeneration, then
# recompute statistics. This is the skeleton; mortality (C4), regeneration (C4)
# and volume/output (C5) are wired in as those chunks land.
# =============================================================================

"""
    setup_growth!(state)

One-time growth setup after initialization (the FVS LSTART pass): compute the
site-dependent diameter-growth constants (DGCONS) and run the diameter-growth
calibration against the input measured growth (COR). Needs density set first.
"""
function setup_growth!(s::StandState)
    build_cycle_schedule!(s)             # CYCLEAT/TIMEINT → cycle-boundary year array (IY)
    # OC ORGANON: the CRATET ORGANON section (oc/cratet.f:155-401) runs BEFORE the FVS-native
    # missing-value dubbing — it dubs valid-ORGANON trees' missing HT/CR via ORGANON PREPARE
    # (PRDHT/PRDCR) and computes ACALIB. A blank-height ORGANON tree gets its ORGANON dub here;
    # dub_missing_heights! then only dubs the NON-ORGANON records (their HT is still 0). No-op
    # for other variants and for OC stands without a big-6 tree.
    s.variant isa OregonCoast && oc_organon_prepare!(s)
    s.variant isa Olympic && op_organon_prepare!(s)   # OP ORGANON NWO PREPARE (op/cratet.f) — dub valid-ORGANON HT/CR + ACALIB
    s.variant isa Utah && ut_cratet_site_adjust!(s)   # ut/cratet.f:99-150 50-yr-base SITEAR (TEMCCF on NOTRE-expanded PROB)
    dub_missing_heights!(s)              # CRATET — dub HT=0 / resolve broken-top NORMHT
    apply_growth_input_types!(s)         # GROWTH IDG/IHTG=1/3 — past DBH/HT field ⇒ increment
    setup_volume_equations!(s)           # VOLEQDEF — per-species NVEL equation ids
    isempty(s.control.sp_bf_vol_eq) && (s.control.sp_bf_vol_eq = copy(s.species.vol_eq)) # VEQNNB = default
    apply_voleqnum_overrides!(s)         # VOLEQNUM — user overrides of those equation ids (cubic only)
    compute_forest_type!(s)              # FORTYP — needed by dgf!'s forest-type term
    compute_density!(s)
    snapshot_esb_inputs!(s)              # ESFLTR (fvs.f:201): freeze the AUTOES ESB inventory-calibration inputs
    ak_esfltr!(s)                        # AK ESFLTR (fvs.f:201): inventory per-point overstory BAAINV/TPAAINV
                                         # (small-tree TPACRE + per-point overstory BAAINV) BEFORE any growth
    root_disease_setup!(s)               # WRD fvs.f RDMN1 init seam — inert unless an RDIN block is active
    dfb_setup!(s)                        # DFB fvs.f DFBSCH init seam — RANSCHED auto-schedule; inert unless a DFB block is active
    dftm_schedule!(s)                    # DFTM DFTMGO→INSCYC seam — force the outbreak cycle to TMBASE=5yr; inert unless a DFTM MANSCHED outbreak is due
    wpbr_setup!(s)                       # WPBR fvs.f BRSETP init seam — per-tree canker init; inert unless a BRUST block is active with a host pine
    s.variant isa SoutheastAlaska || sdi_max_check!(s)   # SDICHK — reset species SDImax if over-dense (AK: at the
                                         # END of CRATET, after the crown dub + calibration — see the AK branch)
    # The DG-constant + calibration pass is variant-specific. NE's DGCONS is trivial
    # (ne/dgf.f:188 zeros DGCON/ATTEN/SMCON; the DG model reads B1/B2/B3 + SITEAR directly),
    # and an uncalibrated NE stand (no measured-DG input) has COR=0 — so the SN LSTART
    # calibration (which needs SN-only coefficient columns) is skipped for NE.
    # NOTE: the LSTART calibration is the SHARED dgdriv.f framework — BOTH variants calibrate DG/HT to
    # the stand's measured past growth (ne/dgdriv.f:8-11), only the `dgf!(s, s.variant)` DDS prediction
    # differs — so it runs for NE too (it was wrongly skipped: net01 HAS measured DG ⇒ COR≠0). NE's DGCONS
    # is just the bark copy (ne/dgf.f:188 zeros DGCON/ATTEN); init_crown_ratios! is SN CRATET (NE DG uses BAL).
    dfint = s.control.growth_fint
    # DG-calibration TERM scale = YR/FINT (dgdriv.f:325). YR = the variant's NATIVE measurement period
    # (htg_period: 5 SN, 10 NE/CS), NOT a hardcoded 5 — the old `5f0/dfint` under/over-scaled NE/CS explicit
    # GROWTH FINT (growth_fint10: NE 14300→17746 bit-exact). The GROWTH FINT struct-default is 5 ("native"), so
    # treat dfint==5 (or unset) as YR ⇒ scale 1 for ALL variants (keeps NE/CS DEFAULT stands bit-exact —
    # net01/cst01 have no GROWTH keyword); only an EXPLICIT non-native FINT scales. Mirrors the `meas_fint`
    # convention in diameter_growth!. `growth_dg_set` (set in kw_growth! when the GROWTH keyword gives FINT)
    # distinguishes an EXPLICIT NE/CS `GROWTH FINT=5` (non-native ⇒ scale 10/5=2, growth_idg1) from the struct
    # default 5 — a plain default of 0 was tried and regressed ~1637 tests (the 5 default is relied on elsewhere).
    yr = htg_period(s.variant)
    dgscale = (s.control.growth_dg_set && dfint > 0f0) ? yr / dfint : 1f0
    if s.variant isa Southern
        dgcons!(s)                        # sets bark_a/bark_b + the SN DGCON
        init_crown_ratios!(s)             # CRATET — dub inventory crown (DENSE backdated-dbh CCF) before calibrate
        calibrate_diameter_growth!(s; scale = dgscale)
    elseif s.variant isa Northeast
        ne_dgcons!(s)                     # bark copy (BKRAT); DGCON/ATTEN = 0
        calibrate_diameter_growth!(s; scale = dgscale)
        # D39: NE DG uses BAL not crown, so the inventory crown was never dubbed ⇒ the cycle-0 FVS_TreeList
        # reported PctCr=0 (live's NE CRATET dubs it). Dub it here (backdated-CCF CRATET, shared init_crown_ratios!)
        # for the report; NE .sum is bit-exact with crown=0, so this must NOT regress it (verified).
        init_crown_ratios!(s)
    elseif s.variant isa CentralStates
        cs_dgcons!(s)                     # DGCON=0, ATTEN=OBSERV, bark copy (BKRAT)
        _cs_init_crowns!(s)               # CRATET: dub missing crowns (backdated-dbh BA) before calibrate — cs/dgf.f reads CR
        calibrate_diameter_growth!(s; scale = dgscale)
    elseif s.variant isa LakeStates
        ls_dgcons!(s)                     # DGCON=0, ATTEN=OBSERV, bark copy (BKRAT) — same as CS
        _ls_init_crowns!(s)               # CRATET: dub missing crowns (ls/dgf.f reads CR); no-op when inventory crowns present
        calibrate_diameter_growth!(s; scale = dgscale)
    elseif s.variant isa CentralRockies
        cr_dgcons!(s)                     # DGCON=0, ATTEN, bark inert; enables c.sigma=SIGMAR for DG serial-corr
        compute_density!(s)               # current-stand density before the CRATET DENSE/crown dub
        cr_misscr = cr_any_missing_crown(s)   # cratet.f:503-522 MISSCR, before the dub fills the crowns
        crown_init_lstart_dead_inclusive!(s)  # cratet.f (== bm core) backdated dead-inclusive DENSE → CROWN. CRATET dub of MISSING (ICR=0) inventory crowns (cr/crown.f);
                                          # eastern variants call init_crown_ratios! here. Without it, 0.1" seedlings keep
                                          # crown_pct=0 ⇒ VARMRT CRI=0 ⇒ EFFTR (100−CRI)/100 = 20× too high ⇒ seedling
                                          # over-kill cascades to the whole stand's mortality distribution.
        _cr_dub_ages!(s; misscr = cr_misscr)  # CRATET age dub (cratet.f:535-552 FINDAG, AFTER :522 CROWN): ABIRTH from
                                          # height for un-aged trees, reading the :175 DENSE's BA/RELDEN snapshot and
                                          # CROWN's BADIST. Before calibration (FVS CRATET→DGDRIV order). Without it
                                          # htgf's AP floors to 1 ⇒ tall trees over-grow height 2-3× (the TopHt drift).
        calibrate_diameter_growth!(s; scale = dgscale)
    elseif s.variant isa Kootenai
        kt_dgcons!(s)                     # KT DGCON (DGHAB+DGFOR+elev/slope-aspect), ATTEN=OBSERV, bark=BKRAT
        compute_density!(s)               # current-stand density (BA/AVH/PCCF/PCT) + RELDEN before the dub
        crown_init_lstart_dead_inclusive!(s)  # kt/cratet.f:598 `IF(MISSCR)CALL CROWN` — KT scans the live AND the
                                          # cycle-0 dead records for a missing crown and dubs them, against the
                                          # SAME backdated dead-inclusive DENSE as every other variant
                                          # (kt/cratet.f:182-184 `LBKDEN = IDG.LT.2; CALL DENSE`). jl had NO LSTART
                                          # crown call for KT at all — the same defect OC had — so every
                                          # missing-crown inventory record kept crown_pct=0 and KT's PCR/DUBSCR
                                          # crown model ran off it. MEASURED on the KT WRD control fixture with the
                                          # live crowns blanked: worst |jl-live| was 31 TPA / 51 BA, the largest
                                          # LSTART-dub residual of any variant.
        calibrate_diameter_growth!(s; scale = dgscale)
    elseif s.variant isa InlandEmpire
        ie_cratet_site_adjust!(s)         # ie/cratet.f:80-118 LM/PY 50-yr-base SITEAR (TEMCCF from the tree list), once
        ie_dgcons!(s)                     # IE DGCON (DGHAB+DGFOR+MAPDSQ/MAPCCF+elev/slope-aspect+site adj), ATTEN=OBSERV
        compute_density!(s)               # current-stand density for the crown dub
        crown_init_lstart_dead_inclusive!(s)  # cratet.f (== bm core) backdated dead-inclusive DENSE → CROWN. CRATET dub of MISSING (ICR=0) inventory crowns (ie/crown.f).
                                          # Was MISSING (like EM; #137 sibling) ⇒ missing-CR seedlings kept crown_pct=0
                                          # ⇒ ie regent HTG1=beta1+beta2·cr loses the crown term ⇒ never cross 4.5'
                                          # ⇒ DBH skipped ⇒ QMD frozen. IE's crown model already dubs d<3 at lstart.
        calibrate_diameter_growth!(s; scale = dgscale)
    elseif s.variant isa EasternMontana
        em_dgcons!(s)                     # EM DGCON/DGDSQ/DGCCF (DGHAB+DGFOR+MAPDSQ+elev/slope-aspect+site adj), ATTEN=OBSERV
        crown_init_lstart_dead_inclusive!(s)  # em/cratet.f (== bm core) backdated dead-inclusive DENSE → CROWN:
                                          # dub of MISSING (ICR=0) inventory crowns + the cycle-0 DEAD records
                                          # (em/crown.f DO 79). The shared helper also supplies the fvs.f:196
                                          # SDICLS SDIAC, so the EM-specific crown_sdi line master carried while
                                          # EM was still calling crown_ratio_update! directly is no longer needed.
                                          # Was live-only current density ⇒ 0.1" seedlings kept
                                          # crown_pct=0 ⇒ _em_smhtgf beta2·cr term = 0 ⇒ HTGR under-predicts ⇒
                                          # never cross 4.5' ⇒ DBH growth skipped ⇒ QMD frozen ⇒ dense self-thin
                                          # holds at the tiny-QMD target (#137). Live dubs these to CR 51-79%.
        calibrate_diameter_growth!(s; scale = dgscale)
        em_regent_aspen_calib!(s)         # aspen/PB (12,17): dub ABIRTH from height (cratet→pothtg) + seed the
                                          # REGENT small-tree HEIGHT self-calib HCOR (regent.f:1298-1365). Without
                                          # both, tiny aspen ran ABIRTH=1 + CON=1 ⇒ HEIGHT over-grew ~2× ⇒ the
                                          # inflated hk fed the aspen inverse-Wykoff DK ⇒ BA cyc1 ~2× (stand
                                          # 373796912489998: oracle 68 → jl 123). Runs after calibrate (dbh restored,
                                          # ht_growth = measured increment). The shared per-cycle attenuation decays CON.
    elseif s.variant isa Teton
        tt_dgcons!(s)                     # TT DGCON (DGSIC·XSITE + DGFOR + aspect/slope/elev), DGDSQ, DGCCF, ATTEN, bark
        _tt_dub_ages!(s)                  # NC/OH (sp15,18) GENGYM height needs ABIRTH dubbed from height (cratet FINDAG,
                                          # IMODTY=4); no-op unless the stand has NC/OH. Other TT species use SBB (no age).
        compute_density!(s)               # density for the crown dub (fresh scalars for the ndead=0 path + dub_ages)
        crown_init_lstart_dead_inclusive!(s)  # tt/cratet.f (== bm core: :242 `LBKDEN=IDG.LT.2; CALL DENSE` over live+dead
                                          # → :639 CROWN) — the shared backdated dead-inclusive DENSE, so DUBSCR reads the
                                          # BACKDATED point CCF / BA and the pre-dub AVHT40. TT's own dead-inclusive-only
                                          # init skipped the LBKDEN backdating: FIA 2783239010690 seedling DUBSCR TPCCF
                                          # 353.42 vs live 355.27 ⇒ CR .52695 vs .525 ⇒ ICR 53 vs 52 ⇒ SMHTGF drift.
        calibrate_diameter_growth!(s; scale = dgscale)
    elseif s.variant isa Utah
        # NB: the CRATET age-50 site-curve conversion (ut/cratet.f) is ALREADY applied once in site_setup!
        # (ut_cratet_site_adjust!, site_index.jl). Calling ut_cratet_siteconv! here too DOUBLE-CONVERTED the
        # site index for the RM-29/RM-32/Meyer site species (1,2,4,5,7,8,9,10,23) — e.g. LP-site FIA stand
        # 3626079010690: raw 46 → 30 (correct=live) → 20.8 (wrong), dropping DGCON ~0.178 ⇒ conifer DG/BA
        # under-grew ~8%. utt01 masked it (DF site species, not in the conversion groups). Removed the dup call.
        ut_dgcons!(s)                     # UT DGCON (DGSIC·XSITE + DGFOR + aspect/slope/elev), DGDSQ, DGCCF, ATTEN, bark
        compute_density!(s)               # density for the crown dub
        ut_misscr = ut_any_missing_crown(s)   # ut/cratet.f:622-644 MISSCR, before the dub fills the crowns
        crown_init_lstart_dead_inclusive!(s)  # cratet.f (== bm core) backdated dead-inclusive DENSE → CROWN. CRATET dub of MISSING crowns (ut/crown.f) — was MISSING
                                          # (EM #137 sibling): missing-CR seedlings kept crown_pct=0 ⇒ ut regent VIGOR(CR)
                                          # lost the crown term ⇒ QMD freeze. UT crown model already dubs missing at lstart.
        _ut_dub_ages!(s; misscr = ut_misscr)  # ut/cratet.f:677 FINDAG after CROWN, on the CRATET-DENSE BA + backdated BAU
                                          # (CR-surrogate 17:19,22 htgf + aspen/oak/MC ABIRTH); no-op without aged species.
        calibrate_diameter_growth!(s; scale = dgscale)
    elseif s.variant isa BlueMountains
        bm_dgcons!(s)                     # BM DGCON + SMCON (habitat-group SMHAB) + DGDSQ/DGCCF/ATTEN, POWER bark
        crown_init_lstart_dead_inclusive!(s)          # CRATET (before DGDRIV): dead-inclusive DENSE → DUBSCR dub of missing-CR
                                          # inventory crowns. Without it dense read 0.1" seedlings keep crown_pct=0 ⇒
                                          # VIGOR floors at 0.30 ⇒ HTGR under-predicts ⇒ never cross 4.5' ⇒ DBH growth
                                          # skipped ⇒ small-tree DG/BA ~2× low (#149). bm/crown.f reads BA/AVH/TPCCF/RMAI.
        _bm_dub_ages!(s)                  # CRATET age dub (bm/cratet.f:656 FINDAG): ABIRTH from height for un-aged
                                          # trees. BM height growth recomputes SITAGE fresh each cycle (doesn't read
                                          # ABIRTH), so this feeds ONLY the Climate-FVS DMORT BIRTHYR (clmorts.f:170);
                                          # inert for non-climate runs. Aged per-cycle at line ~861 (gradd.f:205).
        calibrate_diameter_growth!(s; scale = dgscale)
    elseif s.variant isa CentralIdaho
        ci_dgcons!(s)                     # CI DGCON (DGHAB via ICHBCL + DGFOR + elev/slope-aspect + site adj), ATTEN
        compute_density!(s)               # current-stand density for the crown dub
        crown_init_lstart_dead_inclusive!(s)  # cratet.f (== bm core) backdated dead-inclusive DENSE → CROWN. CRATET dub of MISSING (ICR=0) inventory crowns (ci/crown.f).
                                          # Without it, 0.1" seedlings keep crown_pct=0 ⇒ rcr=1 starves the REGENT
                                          # small-tree height-growth crown term (CI_RG_CRSQ·rcr²) ⇒ seedlings never
                                          # reach breast height (4.5') ⇒ DBH growth skipped ⇒ small-tree DG ~35× low
                                          # (0.024 vs live 0.833) ⇒ mortality g-term too low ⇒ rip too HIGH (logistic
                                          # decreasing in g) ⇒ multi-cycle small-tree OVER-KILL. Mirrors the CR branch.
        calibrate_diameter_growth!(s; scale = dgscale)
    elseif s.variant isa BritishColumbia
        bc_dgcons!(s)                     # BC V3 DGCON (ZNKONST/SSKONST via BEC PrettyName match) — chunk 3, V3 zones only
        compute_density!(s)               # current-stand density before the dub
        crown_init_lstart_dead_inclusive!(s)  # canada/bc/cratet.f:529+ `IF(MISSCR)CALL CROWN` over the live AND the
                                          # cycle-0 dead records, against the same backdated dead-inclusive DENSE
                                          # (bc/cratet.f:206 `LBKDEN = IDG.LT.2`). jl had NO LSTART crown call for BC
                                          # at all — the same defect KT and OC had. The sub-2cm DUBSCR route inside
                                          # bc/crown.f stays deferred (documented), so this dubs the >=2cm
                                          # missing-crown records via the V3 PCR model.
        calibrate_diameter_growth!(s; scale = dgscale)
        dm_init!(s)                       # #196 C1: seed NEWSPRED per-tree initial DMR from damage codes (inert until C6)
    elseif s.variant isa Klamath
        nc_dgcons!(s)                     # NC DGCON (DGFOR/MAPLOC default + DGLAT2 site sp2/6/9 + redwood ln SITEAR) — chunk 3
        # NC calibration SCALE = YR/FINT_meas: NC's DG MODEL basis is 5-yr (blkdat YR=5.0) but the measured past-DG
        # period defaults to IFINT=10 (nct01 "TALLY 2 AT 10 YEARS") ⇒ SCALE=5/10=0.5 (dgdriv.f:419,328 TERM·SCALE).
        # When a GROWTH keyword sets the remeasurement FINT, dgscale=yr/dfint IS that YR/FINT_meas already (so use
        # it directly); only the NO-GROWTH default (dgscale=1) needs the 10-yr-measurement 0.5. Other western
        # variants have YR=IFINT ⇒ scale 1; NC is the unique YR=5-with-10yr-default-measurement case.
        crown_init_lstart_dead_inclusive!(s)  # CRATET DENSE (DEAD-INCLUSIVE) → DUBSCR/Weibull dub of MISSING (ICR=0)
                                          # inventory crowns (nc/crown.f). Was MISSING (like IE #137/EM) ⇒ 0.1"
                                          # seedlings kept crown_pct=0 ⇒ htgr5 CR² term=0 ⇒ QMD frozen. Now dubs
                                          # them; the dead-inclusive AVH (standing-dead heights enter AVHT40)
                                          # feeds DUBSCR's BCR8·AVH term ⇒ correct seedling crown (else AVH≈
                                          # seedling height ⇒ crown ~4pts high ⇒ htgr5 over-growth ⇒ over-mort).
        calibrate_diameter_growth!(s; scale = s.control.growth_dg_set ? dgscale : 0.5f0)
    elseif s.variant isa SoutheastAlaska
        ak_dgcons!(s)                     # AK DGCON (0 + ln COR2 if READCORD), ATTEN=OBSERV; AK bark via ak_bratio in the driver
        compute_density!(s)               # current-stand density for the crown dub (point BA/BAL/TPA + point-Zeide inputs)
        crown_init_lstart_dead_inclusive!(s)  # cratet.f (== bm core) backdated dead-inclusive DENSE → CROWN. CRATET/DUBSCR dub of MISSING (ICR=0) inventory crowns (ak/crown.f);
                                          # D<1 seedlings draw a bounded-normal crown (ak/dubscr.f, RNG-aligned via bachlo).
        calibrate_diameter_growth!(s; scale = dgscale)
        sdi_max_check!(s)                 # ak/cratet.f:664 CALL SDICHK is the LAST step of CRATET — after the :522 CROWN dub
                                          # (whose PRD = ZRD/XMAXPT reads the UNRESET SDIDEF) and the :600 DGDRIV + REGENT
                                          # calibration. Resetting first fed an over-dense stand's dub the reset SDImax:
                                          # FIA 10708179010497 XMAXPT 711.63 vs live 614.10 ⇒ PRD 0.8497 vs 0.9846 ⇒ crowns
                                          # 1-3 pts high ⇒ diverged from 2007.
    elseif s.variant isa WestCascades
        wc_dgcons!(s)                     # WC DGCON (DGFOR/MAPLOC + elev/aspect + King's-SI WO transform) — chunk 3
        compute_density!(s)               # current-stand density (RELDEN) for the crown dub SCALE
        crown_init_lstart_dead_inclusive!(s)  # cratet.f (== bm core) backdated dead-inclusive DENSE → CROWN. CRATET dub of MISSING (ICR=0) inventory crowns
                                          # (wc/crown.f). Missing-CR seedlings would keep crown_pct=0 ⇒ the regent/
                                          # mortality crown term starves (EM #137 class). wct01's inventory crowns are
                                          # all present ⇒ this lstart pass bypasses every tree (verified vs the live
                                          # LSTART CROWN dump: ICR array unchanged), but it is on the setup path and
                                          # gates real-FIA seedling stands. d<1" missing crowns → wc/dubscr.f (bachlo).
        calibrate_diameter_growth!(s; scale = dgscale)
    elseif s.variant isa PacificNorthwest
        pn_dgcons!(s)                     # PN DGCON (20-group; SS g18, WO King's-SI g19, no JFOR remap) — chunk 3
        crown_init_lstart_dead_inclusive!(s)  # pn/cratet.f:162-164 DENSE over live+dead (LBKDEN) → CROWN: Weibull
                                          # dub of missing live crowns + DUBSCR of D<1 AND of the cycle-0 DEAD records
                                          # (pn/crown.f DO 79). Was live-only density + no dead dub ⇒ the dead-record
                                          # BACHLO draws were never consumed ⇒ the whole DGSCOR stream ran 3 draws
                                          # behind live FVS (WRD fixture S248112: ctrl BA +62 by 2090).
        calibrate_diameter_growth!(s; scale = dgscale)
    elseif s.variant isa EastCascades
        ec_dgcons!(s)                     # EC DGCON (32-species uncompressed; WO King's-SI sp28, MH/OS ×3.281) — chunk 3
        compute_density!(s)               # current-stand density (RELDEN) for the crown dub SCALE
        crown_init_lstart_dead_inclusive!(s)  # cratet.f (== bm core) backdated dead-inclusive DENSE → CROWN. CRATET/DUBSCR dub of MISSING (ICR=0) inventory crowns (ec/crown.f)
        calibrate_diameter_growth!(s; scale = dgscale)
    elseif s.variant isa CentralCalifornia
        ca_dgcons!(s)                     # CA DGCON (50-species / 13-group-compressed; ln(SITEAR) form for GS/RW) — chunk 3
        compute_density!(s)               # current-stand density (RELDEN) for the crown dub SCALE
        crown_init_lstart_dead_inclusive!(s)  # cratet.f (== bm core) backdated dead-inclusive DENSE → CROWN. CRATET/DUBSCR dub of MISSING (ICR=0) inventory crowns (ca/crown.f)
        calibrate_diameter_growth!(s; scale = dgscale)
    elseif s.variant isa SouthCentralOregon
        so_dgcons!(s)                     # SO DGCON (33-species uncompressed; RMAI/maical + DGSIC/DGMAI Prognosis terms) — chunk 3
        crown_init_lstart_dead_inclusive!(s)  # so/cratet.f:164-172 (identical to bm) backdated dead-inclusive DENSE →
                                          # CROWN: DUBSCR d<1 seedlings + the cycle-0 DEAD records (so/crown.f DO 79)
        calibrate_diameter_growth!(s; scale = dgscale)
    elseif s.variant isa WestSierra
        ws_dgcons!(s)                     # WS DGCON (43-species uncompressed; ws/dgf.f ENTRY DGCONS) — chunk 3
        ws_htcons!(s)                     # WS HTCON site intercept (ws/htgf.f ENTRY HTCONS) — chunk 4a
        crown_init_lstart_dead_inclusive!(s)  # ws/cratet.f (== bm core) backdated dead-inclusive DENSE → CROWN: Weibull
                                          # dub of missing crowns, ws/dubscr.f for d<1 seedlings AND the cycle-0 DEAD
                                          # records (DO 79). RELDEN is set by compute_density! for every variant now.
        calibrate_diameter_growth!(s; scale = dgscale)
    elseif s.variant isa OregonCoast
        compute_density!(s)               # current-stand density (BA/AVH/PCCF/PCT) + RELDEN before the dub
        crown_init_lstart_dead_inclusive!(s)  # oc/cratet.f:463-499 (== bm core) backdated dead-inclusive DENSE →
                                          # :856 `IF(MISSCR) CALL CROWN`. oc/crown.f LSTART dubs the records whose
                                          # crown is STILL missing after oc_organon_prepare! (which only reloads
                                          # HT/CR for IORG=1 trees): ca/dubscr.f for d<1, the rank-Weibull /
                                          # RW-GS logistic for d>=1, plus the cycle-0 DEAD records (DO 79 :448).
                                          # Was entirely MISSING ⇒ every IORG=0 seedling and every dead record kept
                                          # crown_pct=0, and the DO 79 DUBSCR draws were never consumed.
    elseif s.variant isa Olympic
        op_dgcons!(s)                     # op/dgf.f ENTRY DGCONS — per-species site DGCON (FVS-native IORG=0 trees)
        compute_density!(s)               # current-stand density (BA/AVH/PCCF/PCT) + RELDEN (op/ccfcal.f CCF)
        crown_init_lstart_dead_inclusive!(s)  # cratet.f (== bm core) backdated dead-inclusive DENSE → CROWN. op/crown.f LSTART dub of MISSING inventory crowns
                                          # (Weibull rank d≥1; op/dubscr.f d<1); ORGANON-dubbed HT/CR already set by op_organon_prepare!.
        calibrate_diameter_growth!(s; scale = dgscale)     # op/dgdriv.f LSTART large-tree DG COR (SIGMAR/OBSERV/
                                          # PSIGSQ=0.0898, op_bratio in BOTH the backdating AND the TERM bark).
                                          # WF COR = +0.03379 bit-exact vs live FVSop_clean (op2c_dbg.out:344).
    elseif s.variant isa Ontario
        on_dgcons!(s)                     # canada/on/dgf.f ENTRY DGCONS: DGCON=0, SMCON=0, ATTEN=OBSERV; Penner
                                          # large-tree DG reads its coeffs directly, bark via on_bratio in the driver.
        compute_density!(s)               # current-stand density (BA) for the crown dub
        on_cratet_dead_snapshot!(s)       # cycle-0 dead records' BAPctile/PtBAL (cratet.f:128-160 DENSE)
        crown_ratio_update!(s, s.variant; lstart = true)  # CRATET dub of MISSING (ICR=0) inventory crowns
                                          # (canada/on/crown.f, shared TWIGS NC-125 kernel with NE/CS/LS): current
                                          # inventory BA + DBH (no backdating). ON_BCR1..4 from data/ontario CSV.
        calibrate_diameter_growth!(s; scale = dgscale)   # canada/on/cratet.f → DGDRIV (LSTART, after CROWN/AVHT40):
                                          # dgdriv.f DO 155-202 — SIGMA=SIGMAR (on/blkdat.f), COR fit for measured-DG
                                          # species (FN≥FNMIN), else DO 191 OLDRN=BACHLO(0,SIGMA) bounded by DGSD·SIGMA
                                          # (RNG draws, species order), DO 195 VARDG=(e^{σ²}−1)e^{σ²}/VMLT. Was MISSING ⇒
                                          # VARDG=0 ⇒ SSIGMA=0 ⇒ the tripled records all grew at the expected DG (no
                                          # FU/FM/FL spread; ont01 cyc1 PW DG .0658 cm ×3 vs live .0615/.0795/.0509) and
                                          # the whole RNG stream sat ahead of live by the skipped OLDRN draws.
    end
    cratet_findag_dub!(s)                 # cratet.f "ESTIMATE MISSING TOTAL TREE AGES" (FINDAG → ABIRTH) for the
                                          # variants whose own growth never reads ABIRTH (Climate-FVS BIRTHYR only)
    return s
end

"""
    build_cycle_schedule!(s)

Build the cycle-boundary year array (`control.cycle_year`, the FVS IY) from NUMCYCLE +
TIMEINT + CYCLEAT, mirroring base/fvs.f:106-135. Each cycle's length is its TIMEINT
per-cycle override (`cycle_lengths[k]`) or the uniform `control.year` (default 5); the
lengths are cumulated onto the inventory year (`cycle_year[1]`, INVYEAR) to give the
calendar year at each boundary. CYCLEAT then inserts each requested year as a NEW boundary
strictly inside the run (never extending the end or moving the start), bumping the effective
cycle count `ncycle_eff` (capped at MAXCYC). Idempotent: recomputes purely from the
immutable inputs (cycle_year[1], ncycle, cycle_lengths, cycleat_years), so the per-stand
double-call from setup is safe. For uniform cycles `cycle_year[k+1] == cycle_year[1] + k·per`
exactly, so routing the year derivations through this array is bit-exact for snt01.
"""
function build_cycle_schedule!(s::StandState)
    c = s.control
    ncyc = Int(c.ncycle); ncyc < 1 && (ncyc = 1); ncyc > MAXCYC && (ncyc = MAXCYC)
    per = round(Int, c.year); per < 1 && (per = 5)
    iy = c.cycle_year
    @inbounds for k in 2:MAXCY1                         # cumulate per-cycle lengths → boundary years
        len = c.cycle_lengths[k] > 0 ? Int(c.cycle_lengths[k]) : per
        iy[k] = iy[k-1] + Int32(len)
    end
    for yr in c.cycleat_years                           # CYCLEAT insert (fvs.f:116-135)
        (yr <= iy[1] || yr >= iy[ncyc+1]) && continue   # don't extend end / move start
        @inbounds for i in 1:ncyc
            if iy[i] < yr < iy[i+1]
                ncyc += 1; ncyc > MAXCYC && (ncyc = MAXCYC)
                for k in ncyc+1:-1:i+2; iy[k] = iy[k-1]; end
                iy[i+1] = Int32(yr)
                break
            end
        end
    end
    c.ncycle_eff = Int32(ncyc)
    return s
end

"""
Calendar year at the start of cycle `cyc` (0-based) from the IY schedule (build_cycle_schedule!).
Falls back to the uniform derivation `cycle_year[1] + cyc·per` when the schedule has not been
built yet (a boundary of 0), so direct callers that skip `setup_growth!` still get the right year.
"""
function cycle_year_at(c::Control, cyc::Integer)
    y = Int(c.cycle_year[cyc + 1])
    y > 0 && return y
    per = round(Int, c.year); per < 1 && (per = 5)
    return Int(c.cycle_year[1]) + Int(cyc) * per
end
"Calendar year at the start of the current cycle (`control.cycle`)."
current_cycle_year(s::StandState) = cycle_year_at(s.control, Int(s.control.cycle))
"Length in years of cycle `cyc` (0-based) = next boundary − this boundary."
cycle_period_at(c::Control, cyc::Integer) = Int(c.cycle_year[cyc + 2] - c.cycle_year[cyc + 1])

"""
    compute_density!(state)

Recompute the per-cycle stand density quantities the growth models read:
basal area, average dominant height (AVH), and per-point basal area (PTBAA).
"""
function compute_density!(s::StandState; cratet_ind::Bool = false)
    s.plot.basal_area = stand_ba(s)
    s.plot.avg_height = stand_top_height(s; cratet_ind = cratet_ind)
    # RMSQD (stand quadratic mean diameter, inches) — DENSE computes it into COMMON; ON's Penner
    # dgf! (ontario/diameter_growth.jl) reads it as `p.qmd*ON_INtoCM`. No other variant reads
    # p.qmd (summary QMD comes from stand_qmd() directly), so this is inert elsewhere; gate to
    # Ontario to keep the shared density path byte-identical for every other variant.
    s.variant isa Ontario && (s.plot.qmd = stand_qmd(s))
    point_basal_area!(s; cratet_ind = cratet_ind)
    point_density!(s)                  # PCCF/PTPA per point (regen crown ratio + TCONDMLT weights)
    stand_pct!(s; cratet_ind = cratet_ind)  # PCT = stand BA percentile (for DGF competition)
    # RELDEN = stand CCF, set by DENSE for EVERY variant (dense.f → CCFCAL sum). This was a per-variant whitelist
    # (KT IE EM TT UT BC BM CI EC SO WS OP); PN/WC/NC/CA read p.relative_density in crown.f SCALE (and NC in dgf!)
    # but were never set ⇒ RELDEN=0 ⇒ crown SCALE capped at 1.0 ⇒ every crown dubbed/updated high ⇒ one-directional
    # growth over-prediction (PN WRD fixture S248112: BA +59 by 2090). Engine consumers (LPMPB, COVER, DFTM,
    # establishment) likewise read 0 for any non-whitelisted variant. Set at whatever t.n is current: the backdated
    # calibration pass runs dead-inclusive (RELDM1), the growth-cycle pass live-only — FVS's DENSE→DGF/CROWN flow.
    s.plot.relative_density = s.variant isa BritishColumbia ? bc_stand_ccf(s) : stand_ccf(s)
    return s
end

# DENSE's per-point CCF (point_density!) over the CURRENT records, leaving density.point_ccf/point_tpa untouched.
function _fresh_point_ccf(s::StandState)
    sc = copy(s.density.point_ccf); st = copy(s.density.point_tpa)
    point_density!(s)
    out = copy(s.density.point_ccf)
    copyto!(s.density.point_ccf, sc); copyto!(s.density.point_tpa, st)
    return out
end

"Number of early cycles that use deterministic record tripling (FVS ICL4)."
const TRIPLE_CYCLE_LIMIT = 2

# Western variants whose establishment runs the SHARED `estb/` model AND run the recurring AUTOES natural
# tally: estab.f:1493 calls ESGENT (esgent.f:49 → SPESRT) whenever new regen records were added this cycle,
# rebuilding IND1 in ascending physical order and DISCARDING the post-TRIPLE REASS lineage interleave. Without
# the reset, a stand that keeps adding AUTOES natural regen across tripled cycles mis-orders the per-tree
# DGSCOR/REGENT RNG draws (the IE `none` +BA over-production, fixed for IE at 80d96174). Only IE and EM run the
# recurring AUTOES tally in jl (`ie_autoes_establish!`, simulate.jl above), so only they exercise this gap —
# a stratified 120-stand `none` A/B vs FVSem_g16 measured EM at 42/47 changed stands improved (many
# −1 BA rel−25% → bit-exact), meanFinalBA bias −0.72→−0.41. UT/BM/TT/CI do NOT run AUTOES: on `none` they add
# no records (reorder never fires) and on `plant` the one-shot bare-stand planting fires the gate only where
# sort_key is already ascending (measured `changed=false`), so the reorder is a measured no-op for them — they
# are deliberately EXCLUDED until a regime that exercises their tripled-lineage establishment is measured.
uses_estab_spesrt(v) = v isa InlandEmpire || v isa EasternMontana

"""
    fertilizer_growth!(s; fint)

FERTILIZE / FFERT (ffert.f): the 200-lb-N fertilizer response — a multiplicative boost to each
tree's squared-diameter change (DDS) and height growth for up to 10 years after application,
applied for the `iflen` of the cycle's years that fall in that window and scaled by the application
efficacy. Carries over across cycles via `ifert_date`/`ifert_eff`. No-op until a FERTILIZE keyword
fires (so default stands are unchanged). SN is outside the model's calibrated DF/GF range — Fortran
warns but still applies the (species-agnostic) factor, so we match.
"""
function fertilizer_growth!(s::StandState; fint::Float32 = 5f0)
    c = s.control
    (isempty(c.fertilize_events) && c.ifert_date < 0) && return s
    yr = cycle_year_at(c, Int(c.cycle))   # IY schedule (TIMEINT/CYCLEAT-aware) = cycle START
    # OPCYCL containing-cycle bucketing: a FERTILIZE at a mid-cycle date (1995 in a 10-yr NE cycle 1990→2000)
    # activates THIS cycle, and FVS sets IFFDAT=IY(ICYC) — the cycle-START year, NOT the keyword date (ffert.f:75)
    # — so `ifert_date = yr` here is correct (full-cycle effect from the cycle start). Boundary dates unchanged.
    ce = cycle_year_at(c, Int(c.cycle) + 1); ce <= yr && (ce = yr + 1)
    @inbounds for ev in c.fertilize_events           # a fert scheduled this cycle becomes active (OPDONE)
        (yr <= Int(ev.year) < ce) && (c.ifert_date = Int32(yr); c.ifert_eff = ev.params[1])
    end
    c.ifert_date < 0 && return s
    ifint  = round(Int, fint)
    ifstrt = yr - Int(c.ifert_date)
    ifstrt > 10 && return s                          # > 10 yr since application ⇒ effect gone
    iflen  = min(ifstrt + ifint, 10) - ifstrt        # years of fertilizer effect within this cycle
    iflen <= 0 && return s
    feff = c.ifert_eff
    t = s.trees; ba = s.plot.basal_area
    ba_a = s.calib.bark_a; ba_b = s.calib.bark_b
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        dib = d * bark_ratio(ba_a, ba_b, t.species[i], d)
        bal = (1f0 - t.crown_ratio[i] / 100f0) * ba   # basal area in larger trees (PCT in crown_ratio)
        rdds = exp(0.1108f0 * log(d) + 0.003004f0 * bal / log(d + 1f0))
        rdds > 2.6f0 && (rdds = 2.6f0)
        dg  = t.diam_growth[i]
        dds = 2f0 * dib * dg + dg * dg                # squared-diameter change this cycle
        ddsit = (dds / 5f0) * (rdds * iflen * feff + ifint - iflen)
        t.diam_growth[i] = sqrt(dib * dib + ddsit * 5f0 / fint) - dib
        htgit = (t.ht_growth[i] / 5f0) * (1.1626f0 * iflen * feff + ifint - iflen)
        t.ht_growth[i]   = htgit * 5f0 / fint
    end
    return s
end

"""
    _fire_due(s) -> Bool

A scheduled SIMFIRE is due THIS cycle iff the cycle's year range [cycle_start, cycle_end)
CONTAINS `fire.fire_year` — FVS OPCYCL (opcycl.f:58-64): an activity at date D is assigned to
the cycle with IY(i) ≤ D < IY(i+1), NOT only when D is a cycle boundary. So a SIMFIRE scheduled
at a non-boundary year fires in its containing cycle (on the cycle-start stand, the FMBURN basis).
For a boundary fire year this reduces EXACTLY to current_cycle_year == fire_year (boundary fires —
net01/snt01 — are unchanged), so it only ADDS mid-cycle support.
"""
function _fire_due(s::StandState)::Bool
    return _due_fire_index(s) != 0
end

"""
    _due_fire_index(s) -> Int

Index into `s.fire.fire_schedule` of the earliest scheduled SIMFIRE whose year falls in the current
cycle's range [cycle_start, cycle_end) (FVS OPCYCL bucketing), or 0 if none. Multiple SIMFIRE keywords
(fire_repeat) each schedule their own event; this picks the one due THIS cycle. Falls back to the legacy
scalar `fire_year` when no schedule list is present (defensive — every SIMFIRE now populates the list).
"""
function _due_fire_index(s::StandState)::Int
    (s.fire === nothing || !s.fire.active) && return 0
    cyc = Int(s.control.cycle)
    cs = current_cycle_year(s)
    ce = cycle_year_at(s.control, cyc + 1)
    ce <= cs && (ce = cs + 1)              # last/degenerate cycle ⇒ exact-match only
    sched = s.fire.fire_schedule
    if isempty(sched)
        fy = Int(s.fire.fire_year)
        return (fy != 0 && cs <= fy < ce) ? -1 : 0   # -1 = legacy scalar path
    end
    @inbounds for i in 1:length(sched)
        fy = Int(sched[i][1])
        cs <= fy < ce && return i
    end
    return 0
end

# Duff fuel moisture MOIS(1,5) at a fire in calendar `year`, mirroring fmburn!'s own moisture selection
# (MOISTURE-keyword override else the FMMOIS dryness-model table). Used to derive the FMCONS EXPOSR
# mineral-soil exposure for the AUTOES BurnPrep seam. Row 1 col 5 = duff (fmburn.f:373-380).
function _fire_duff_moisture(s::StandState, year::Int)::Float32
    (s.fire === nothing) && return 0f0
    movr = _active_moisture_override(s, year)
    mois = movr === nothing ? fuel_moisture(Int(s.fire.fmois), s.variant) : _moisture_matrix(movr)
    return Float32(mois[1, 5])
end

"""
    _maybe_burn!(s, fint) -> fire_mort

Run a scheduled SIMFIRE if it is due this cycle (`_fire_due`, FMMAIN). Operates on the current
(post-MORTS, post-TRIPLE) records at cycle-start dimensions, and returns the periodic mortality
VOLUME of the fire-killed TPA (each record's lost TPA × its cycle-start cubic volume), for the
caller to add to OMORT. No-op (returns 0) otherwise.
"""
function _maybe_burn!(s::StandState, fint::Float32)::Float32
    di = _due_fire_index(s)
    di == 0 && return 0f0
    if di > 0                                # load this scheduled event's conditions into the scalars
        ev = s.fire.fire_schedule[di]
        s.fire.fire_year = Int32(ev[1]); s.fire.swind = ev[2]; s.fire.fmois = Int32(ev[3])
        s.fire.atemp = ev[4]; s.fire.mortcode = Int32(ev[5]); s.fire.psburn = ev[6]
        s.fire.burnseas = Int32(ev[7])
    end
    yr = current_cycle_year(s)   # IY schedule (TIMEINT/CYCLEAT-aware); fire fires on this cycle's stand
    t = s.trees
    pre_tpa = Float32[t.tpa[i] for i in 1:t.n]
    pre_cfv = Float32[t.cuft_vol[i] for i in 1:t.n]
    fmburn!(s; atemp = s.fire.atemp, wind = s.fire.swind, fmois = Int(s.fire.fmois),
            psburn = s.fire.psburn, mortcode = Int(s.fire.mortcode),
            burnseas = Int(s.fire.burnseas), flmult = s.fire.flmult, crburn = s.fire.crburn,
            year = yr, cyclen = fint)
    fm = 0f0
    @inbounds for i in 1:length(pre_tpa)
        fm += (pre_tpa[i] - t.tpa[i]) * pre_cfv[i]
    end
    # AUTOES disturbance seam (fmcons.f:247-258): "a fire is a disturbance so the year of the fire is the year
    # of the disturbance." When automatic tallies are on (LAUTAL), FVS's FMCONS OPADDs a TALLY (activity 427,
    # PRMS(1)=FLOAT(IYR)) at the fire year — this forces ESNUTR to fire the disturbance regen tally (NTALLY=1)
    # INDEPENDENTLY of any harvest removal (ONTREM/OCVREM are cut-only, so a fire alone never trips the LAUTAL
    # removal path). Without it jl fell through to the small LINGRW ingrowth and heavily UNDER-regenerated after a
    # fire (measured 1143093321290487 simfire: post-fire TPA 285 vs oracle 1244). FMCONS also OPADDs a BurnPrep
    # (491, PRMS(1)=EXPOSR mineral-soil exposure) which shifts species COMPOSITION only (esetpr/burnprep); the
    # dominant TPA effect is the 427. Gated to the AUTOES variants (IE/EM, which read 427 via ie_autoes_establish!)
    # and to LAUTAL, mirroring FVS's IF(LAUTAL). Injected post-FMBURN, pre-ESTAB — the FMBURN→ESNUTR order FVS
    # uses; the existing scheduled-427 path in ie_autoes_establish! then consumes it (idt=fire year sits in
    # [year,next_year) at the fire cycle ⇒ fires exactly once, inert every other cycle). ESB1 uses est.inv_baaold
    # (the ESFLTR-frozen inventory BA) so the post-fire tally calibrates against inventory, not the depleted BA.
    if s.estab.lautal && (s.variant isa InlandEmpire || s.variant isa EasternMontana)
        push!(s.control.schedule, ScheduledActivity(Int32(yr), Int32(427),
              (Float32(yr), 0f0, 0f0, 0f0, 0f0, 0f0)))
        # EXPOSR (fmcons.f:186-208): PRDUF(%) = 83.7 − 0.426·m_duff%, floored 0; EXPOSR = (−8.98 + 0.899·PRDUF)·
        # PSBURN/100, zeroed when PRDUF<10. Added only when >0 so the (usual) EXPOSR=0 fire stays byte-identical.
        mduff = _fire_duff_moisture(s, Int(yr))
        prduf = 83.7f0 - 0.426f0 * mduff * 100f0
        prduf < 0f0 && (prduf = 0f0)
        exposr = prduf < 10f0 ? 0f0 : (-8.98f0 + 0.899f0 * prduf) * (Float32(s.fire.psburn) / 100f0)
        exposr > 0f0 && push!(s.control.schedule, ScheduledActivity(Int32(yr), Int32(491),
              (Float32(yr), exposr, 0f0, 0f0, 0f0, 0f0)))
    end
    # One-shot: drop the just-fired event (index `di`, captured before firing — no deletions in between, so
    # it is still valid) from the schedule; resync the scalars to the next pending fire so a later cycle's
    # SIMFIRE (fire_repeat 2020 after 2000) still fires. Empty schedule ⇒ clear fire_year.
    di > 0 && deleteat!(s.fire.fire_schedule, di)
    if isempty(s.fire.fire_schedule)
        s.fire.fire_year = Int32(0)
    else
        ev = s.fire.fire_schedule[1]
        s.fire.fire_year = Int32(ev[1]); s.fire.swind = ev[2]; s.fire.fmois = Int32(ev[3])
        s.fire.atemp = ev[4]; s.fire.mortcode = Int32(ev[5]); s.fire.psburn = ev[6]
        s.fire.burnseas = Int32(ev[7])
    end
    compute_density!(s)
    return fm
end

# FMKILL's ICR=-FMICR travels with the record in FVS (the sign is part of ICR). jl keeps the kept crown %
# in s.fire.crown_bypass, so a COMCUP swap-from-end must carry it too (composed with the RDTDEL hook).
function _record_move_hook(s::StandState)
    rd = rd_tdel_hook(s)
    # BRTDEL (tredel.f, alongside RDTDEL): WPBR per-slot state follows the moved record.
    br = (s.wpbr !== nothing && (s.wpbr isa WpbrState) && (s.wpbr::WpbrState).active) ? true : nothing
    byp = s.fire === nothing ? nothing : s.fire.crown_bypass
    if byp === nothing || isempty(byp)
        br === nothing && return rd
        return (iv, ir) -> (rd === nothing || rd(iv, ir); wpbr_tdel!(s, iv, ir))
    end
    return (iv, ir) -> begin
        rd === nothing || rd(iv, ir)
        br === nothing || wpbr_tdel!(s, iv, ir)
        if ir <= length(byp)
            iv <= length(byp) && (byp[iv] = byp[ir])
            byp[ir] = Int32(0)
        end
    end
end

# Drop bypass entries past the compacted record count, so regen appended later never inherits one.
function _trim_crown_bypass!(s::StandState)
    s.fire === nothing && return
    byp = s.fire.crown_bypass
    length(byp) > s.trees.n && resize!(byp, s.trees.n)
    return
end

"""
    crown_ratio_update_fvs!(s; kwargs...)

CROWN with FVS's negative-ICR bypass (every variant's crown.f: "IF ICR(I) IS NEGATIVE, CROWN RATIO CHANGE WAS
COMPUTED IN A PEST DYNAMICS EXTENSION. SWITCH THE SIGN ON ICR(I) AND BYPASS CHANGE CALCULATIONS"). The bypassed
records keep the value FMKILL set (s.fire.crown_bypass). They reach no random draw in any variant's update
(draws are only on the ICR=0 / DBH≤0 / LSTART dub paths), so running the update and restoring them is exact.
"""
function crown_ratio_update_fvs!(s::StandState; kwargs...)
    byp = s.fire === nothing ? nothing : s.fire.crown_bypass
    crown_ratio_update!(s, s.variant; kwargs...)
    (byp === nothing || isempty(byp)) && return
    t = s.trees
    @inbounds for j in 1:min(length(byp), t.n)
        byp[j] != 0 && (t.crown_pct[j] = byp[j])
    end
    empty!(byp)
    return
end

"""
    mortality_and_fire!(s; fint) -> StandState

Apply this cycle's MORTS background/density mortality and, if a SIMFIRE is scheduled for the
current year, the fire kill — on the cycle-start (pre-growth) stand, mirroring FVS GRINCR(MORTS→WK2)
→ GRADD/FMKILL(WK2=MAX(WK2,FIRKIL)). A tree dies from whichever is LARGER (density MORTS or fire),
not both summed (fmkill.f:86); the EXCESS density kill over the fire kill is booked to snags
(fmkill.f:135). Non-fire cycles are byte-identical to a plain `mortality!` (fk≡0 ⇒ kill = mk).

Factored out of `grow_cycle!` so the fire half can also be driven from the report path at the FVS
sample phase (FMBURN before FMCRBOUT/annual loop) — see docs/audit/BACKLOG.md item 8 (#28).
"""
function mortality_and_fire!(s::StandState; fint::Float32 = 5f0,
                             stash = nothing,
                             post_fire::Union{Nothing,Function} = nothing,
                             book_snags::Bool = true)
    t = s.trees
    fire_now = _fire_due(s)   # OPCYCL: fires in the cycle whose range contains fire_year (incl. mid-cycle)
    if !fire_now
        # Non-fire cycle with a deferred FFE annual loop (R6 variants, summary.jl): FMMAIN runs after GRINCR's
        # increment draws and before FMKILL/UPDATE apply the mortality (FMPROB = PROB), so run it here, ahead of the
        # mortality booking. jl's MORTS draws no main-stream RANN, so the state is FVS's FMMAIN state.
        post_fire === nothing || post_fire(s)
        mortality!(s, s.variant; fint = fint, book_snags = book_snags)   # MORTS (FVS GRINCR order)
        return (0f0, false)        # non-fire OMORT is computed by the caller (pre-TRIPLE originals)
    end
    # FIRE CYCLE — FVS order: MORTS (GRINCR, on the ORIGINAL ITRN records — VARMRT distributes a stand
    # total over records, so it MUST see the UN-tripled set) → TRIPLE (triple.f, splits TPA+WK2 .60/.25/
    # .15) → FMBURN (GRADD fire draws an independent XRAN per TRIPLED record) → FMKILL WK2=MAX(MORTS,FIRKIL).
    nrec = t.n
    pre  = Float32[t.tpa[i] for i in 1:nrec]
    mortality!(s, s.variant; fint = fint, book_snags = false)  # MORTS on the un-tripled stand
    mk   = _morts_wk2(s, pre, nrec)                            # per-ORIGINAL density+bkgd+cap kill (MORTS WK2 itself)
    @inbounds for i in 1:nrec; t.tpa[i] = pre[i]; end          # restore PROB for the fire pass
    tripled = stash !== nothing
    if tripled                                                 # split TPA + the MORTS kill onto the 3 records
        triple_records!(s, stash)                              # central=i, upper=nrec+2i-1, lower=nrec+2i
        mks = zeros(Float32, t.n)
        @inbounds for i in 1:nrec
            mks[i]         = mk[i] * 0.60f0
            mks[nrec+2i-1] = mk[i] * 0.25f0
            mks[nrec+2i]   = mk[i] * 0.15f0
        end
        mk = mks
    end
    n   = t.n
    pre = Float32[t.tpa[i] for i in 1:n]                       # cycle-start TPA on the (now tripled) set
    empty!(s.fire.fmicr)                                       # FMICR of THIS burn only (set by fmburn!)
    _maybe_burn!(s, fint)                                      # FMBURN/FIRKIL — independent XRAN per record
    # FVS FMMAIN order: FMBURN (just done) → FMCRBOUT carbon report → annual fuel loop (FMSNAG/FMCWD/
    # FMCADD) — all BEFORE FMKILL's WK2 combine below. `post_fire` runs the carbon sample + the FFE annual
    # fuel update here so they see the post-fire, start-of-cycle fuel pools + fresh fire snags (#28).
    post_fire === nothing || post_fire(s)
    extra = Vector{Float32}(undef, n)
    mort  = 0f0
    akfk = length(s.fire.firkil) >= min(n, t.n)
    @inbounds for j in 1:min(n, t.n)
        fk = pre[j] - t.tpa[j]                                 # fire kill (FIRKIL) on this record
        # FMKILL(1) reads FIRKIL itself (clamped ≤ PROB, fmkill.f:75) and WK2 = MAX(WK2,FIRKIL) is MortPA — not the
        # PROB−survivor difference, which rounds (akffe 2013 FVS_TreeList MortPA 1.6238010 live vs 1.6238011).
        akfk && (fk = min(s.fire.firkil[j], pre[j]))
        t.tpa[j] = pre[j] - max(mk[j], fk)                     # WK2 = MAX(MORTS, fire), per fmkill.f:86
        extra[j] = max(0f0, mk[j] - fk)                        # regular snags = WK2 − FIRKIL (fmkill.f:135)
        m = akfk ? max(mk[j], fk) : pre[j] - t.tpa[j]
        mort += m * t.cuft_vol[j]                              # OMORT on the cycle-start per-record CFV
        t.mort_pa[j] = m                                       # FVS_TreeList MortPA (post-TRIPLE)
    end
    # WPBR BRTREG (gradd.f:126) follows FMKILL(1) (:122): on a FIRE tripling cycle it runs here, on the tripled
    # records at full PROB (`pre`), reading WK2 = MAX(MORTS, fire) (a canker kill overrides WK2=PROB·0.99999).
    if tripled && s.wpbr !== nothing && (s.wpbr::WpbrState).active
        surv = Float32[t.tpa[j] for j in 1:n]
        wpbr_brtreg!(s, fint, pre; wk2_hint = Float32[pre[j] - surv[j] for j in 1:n])
        @inbounds for j in 1:n
            d = surv[j] - t.tpa[j]
            d != 0f0 && (mort += d * t.cuft_vol[j]; t.mort_pa[j] += d)
        end
    end
    book_mortality_snags!(s, extra, n, fint)                   # FMSDIT snags for the EXCESS MORTS only (FMKILL)
    # FMKILL crown hand-back (fmkill.f:92-94): IF(FMICR<1) FMICR=1; IF(FMICR<|ICR|) ICR=-FMICR. The negative
    # ICR makes the next CROWN keep FMICR instead of recomputing (crown.f "ICR(I) WAS CALCULATED ELSEWHERE");
    # jl stores the kept value in crown_bypass, applied by crown_ratio_update_fvs!.
    fm = s.fire.fmicr
    if length(fm) == n
        byp = s.fire.crown_bypass
        resize!(byp, n); fill!(byp, Int32(0))
        @inbounds for j in 1:n
            f = max(fm[j], Int32(1))
            if f < abs(t.crown_pct[j])
                t.crown_pct[j] = f
                byp[j] = f
            end
        end
    end
    empty!(fm)
    compute_density!(s)
    return (mort, tripled)
end

# grincr.f:74 LTRIP's ITRN (see Control.itrn_grincr): the record count at GRINCR entry of the current cycle, latched once.
function latch_itrn_grincr!(s::StandState)
    c = s.control
    if c.itrn_grincr_cycle != c.cycle
        c.itrn_grincr = Int32(s.trees.n); c.itrn_grincr_cycle = c.cycle
    end
    return c.itrn_grincr
end

"""
    grow_cycle!(state; fint=5f0) -> (; accretion, mortality)

Advance the stand by one growth cycle: recompute density, grow diameters/heights
(with record tripling), apply mortality, update dimensions, and recompute volumes.
Returns the period's per-acre cubic accretion and mortality (OACC/OMORT,
vols.f:190 / update.f:60): accretion = Σ(newCFV−oldCFV)·survivingTPA / fint,
mortality = Σ killedTPA·oldCFV / fint, both ÷ gross area. Requires the cycle-start
volumes to be present in `trees.cuft_vol` (run `compute_volumes!` once at setup).
"""
# The MORTS kill WK2 per original record, EXACT: the applied kill buffer (mortality! leaves `killed` in
# s.scratch.mort_killed and sets tpa = max(0, P − killed)) when it reproduces the observed survivor bit-for-bit, else
# the old−new difference. Recovering WK2 by subtraction alone is inexact (P − fl(P−W) ≠ W), and the tripling seams need
# the exact WK2 because FVS forms survivors as PROB·w − WK2·w (triple.f:30-32/75-76, then UPDATE).
function _morts_wk2(s::StandState, old_tpa::Vector{Float32}, nlive::Int)::Vector{Float32}
    t = s.trees; kb = s.scratch.mort_killed
    w = Vector{Float32}(undef, nlive)
    on = s.variant isa Ontario        # ON VARMRT can leave WK2 > PROB (varmrt.f exhaustion); WK2 is reported unclamped
    @inbounds for i in 1:nlive
        k = i <= length(kb) ? (on ? kb[i] : min(kb[i], old_tpa[i])) : -1f0
        w[i] = (k >= 0f0 && t.tpa[i] == max(0f0, old_tpa[i] - k)) ? k : old_tpa[i] - t.tpa[i]
    end
    return w
end

function grow_cycle!(s::StandState; fint::Float32 = 5f0,
                     carbon_hook::Union{Nothing,Function} = nothing,
                     fuel_period::Union{Nothing,Real} = nothing,
                     ffe_init_period::Union{Nothing,Real} = nothing,
                     wwpb_barrier::Union{Nothing,Function} = nothing,
                     fmmain_hook::Union{Nothing,Function} = nothing)
    # BM: the first grow cycle's DGDRIV reads the PCT that CRATET's DENSE (cratet.f:692) built over CRATET's IND
    # (IND1-seeded RDPSRT, see bm_cratet_ind!), not a fresh gradd.f:186-style sort; a thin re-sorts (cuts.f:302).
    compute_density!(s; cratet_ind = (_fvs_ind_lifecycle(s.variant) &&
                                      s.control.cycle == Int32(0)))   # CI: ci/cratet.f:230-233/:337 → :732 DENSE, same as BM
    # ECON: ECSETP (fvs.f:148, once before cycling — default STRTECON at IY(1), revenue-class sort) then
    # ECSTATUS(…,0) (grincr.f:273, cycle start before CUTS). Inert unless an ECON block is active.
    econ_cycle_start!(s)
    root_disease_mn2!(s, fint)           # WRD grincr.f RDMN2 seam (cycle start) — inert unless an RDIN block is active
    # Climate-FVS: realize the cycle-scheduled GrowMult/MortMult weights for this cycle (FVS ICYC = jl cycle+1)
    # BEFORE growth/mortality read growmult/mortmult. Inert unless a CLIMATE block parsed GrowMult/MortMult events.
    (s.climate !== nothing && s.climate.active) && apply_climate_schedule!(s, Int(s.control.cycle) + 1)
    # clgmult.f runs every cycle inside DGDRIV even with ITRN=0 (SPWTS=0 ⇒ SPGMULT=1); jl skips growth on a bare
    # stand, so start each cycle at 1 (climate_growth_wk4! overwrites it when it runs) and clear last cycle's report.
    (s.climate !== nothing && s.climate.active) &&
        (fill!(s.climate.spgmult, 1f0); s.climate.pending_report = nothing)
    # Climate SPCALIB (clmorts.f:57-75 ICYC==1): set at cycle 0 from INVENTORY presence, BEFORE establishment
    # adds regen — so an empty-at-cycle-1 establishment stand correctly gets SPCALIB=−1 (matches the oracle),
    # not a mis-calibration from a later cycle's established cohort. Inert unless CLIMATE is active.
    (s.climate !== nothing && s.climate.active && s.control.cycle == Int32(0)) && init_climate_spcalib!(s)
    # IE crown OLDPCT init (cratet.f:513): at the first grow cycle, seed OLDPCT from the BACKDATED percentile that
    # the initial-CRATET backdating DENSE (cratet.f:217-219) computes — NOT the plain inventory PCT. The backdating
    # runs over the un-deleted inventory (dead-inclusive) with each live diameter backdated to start-of-growth; see
    # ie_seed_backdated_oldpct!. Later cycles get OLDPCT from the post-crown snapshot (crown_ratio) below.
    # EM runs the same cratet.f:481 OLDPCT=PCT after the same backdating DENSE (em/dense.f ≡ ie/dense.f), with em/bratio.f.
    if s.variant isa EasternMontana && s.control.cycle == Int32(0)
        ie_seed_backdated_oldpct!(s; bratio = em_bratio)
        s.plot.old_ba = s.plot.basal_area                        # cycle-1 OBA/RDM1 = inventory density (as IE)
        s.plot.relative_density_prev = s.plot.relative_density
    end
    if s.variant isa InlandEmpire && s.control.cycle == Int32(0)
        # cratet.f:513 saves the PCT of the :219 backdating DENSE, which crown_init_lstart_dead_inclusive! already
        # built (c.cratet_pct) over CRATET's IND1-seeded RDPSRT(.FALSE.) order; ie_seed_backdated_oldpct!'s own
        # identity-seeded sort swapped every equal-DBH pair's OLDPCT (FIA 3027007010690: three tied pairs ⇒ crown
        # ICR ±1 ⇒ cycle-2 DDS ±0.0154 on 11 records).
        if length(s.calib.cratet_pct) == s.trees.n
            copyto!(s.trees.old_crown_pct, 1, s.calib.cratet_pct, 1, s.trees.n)
        else
            ie_seed_backdated_oldpct!(s)
        end
        ie_dub_aspen_birthage!(s)   # cratet.f:544-563 CALL FINDAG: dub ABIRTH=SITAGE for sp18/20/21 (AS/MM/PB)
        # IE crown DCR backdates against OLDBA/RELDM1 = the PREVIOUS cycle's stand BA/RELDEN (dense.f:239-240,
        # threaded start-of-cycle; crown.f:277-281 reads them for DCRCON). Seed cycle-1's pair from the
        # inventory (pre-growth) density — analog of the OLDPCT seed above (matches oracle OBA[1]=inventory BA).
        s.plot.old_ba = s.plot.basal_area
        s.plot.relative_density_prev = s.plot.relative_density
    end
    # KT: kt/cratet.f:578 OLDPCT=PCT right before the LSTART CROWN — the PCT of the backdating DENSE, which the shared
    # crown_init_lstart_dead_inclusive! snapshots as cratet_pct — and cycle-1 OBA/RDM1 = inventory density (kt/crown.f
    # is IE's crown model; same seeds as IE above).
    if s.variant isa Kootenai && s.control.cycle == Int32(0)
        length(s.calib.cratet_pct) == s.trees.n && copyto!(s.trees.old_crown_pct, 1, s.calib.cratet_pct, 1, s.trees.n)
        s.plot.old_ba = s.plot.basal_area
        s.plot.relative_density_prev = s.plot.relative_density
    end
    apply_setsite!(s)                                      # SETSITE (act 120): mid-run site change (RCON), before growth
    # FVS latches LTRIP (grincr.f:74) at cycle start from the CURRENT NOTRIP, BEFORE COMCUP (:391) may set
    # NOTRIP=.TRUE. So capture NOTRIP here: a COMPRESS this cycle suppresses tripling only from NEXT cycle.
    notrip_start = s.control.no_tripling
    # grincr.f:74 evaluates LTRIP's ITRN.LE.(MAXTRE/3) at GRINCR entry — BEFORE CUTS (grincr.f:292) TREDELs the zero-PROB
    # records the previous cycle's mortality left (cuts.f:255-275) and before COMCUP. Latched by the .sum driver ahead
    # of its own CUTS call (latch_itrn_grincr!), else here.
    itrn_grincr = Int(latch_itrn_grincr!(s))
    compressed = apply_compress!(s)                        # COMPRESS (act 250): cluster records → NCLAS (sets NOTRIP for later cycles)
    # ECON: zero the cycle's harvest accumulators; cuts!/_log_cut! values each removed tree.
    econ_on = s.econ !== nothing && s.econ.active
    econ_on && (s.econ.cycle_cost = 0f0; s.econ.cycle_rev = 0f0)
    # Zero-PROB record deletion happens at the START of the next cycle in FVS, not at the end of the cycle that killed
    # them: cuts.f:255-275 ("THEY GET HERE WITH ZERO PROB FROM PREVIOUS CYCLE MORTALITY") TREDELs PROB≤1E-10 on CUTS
    # entry, and COMCUP (grincr.f:391, after CUTS, before growth) TREDELs PROB≤1E-5. A record killed outright this
    # cycle (MORTS/FMKILL via UPDATE in GRADD) therefore survives THROUGH this cycle's ESNUTR, so new regen/sprouts are
    # appended AFTER it and the next cycle's swap-from-end TREDEL moves them into its slot. jl formerly ran COMCUP at
    # the end of growth (before ESNUTR), appending regen to an already-compacted list ⇒ a different physical record
    # order ⇒ different per-record REGENT/DGSCOR ZZRAN assignment. MEASURED FVSie_g16 24829032010900 post-SIMFIRE:
    # fire-killed ingrowth stubs (TPA 0) persist in the 2024 treelist; 193/193 regen height increments desynced.
    cuts_entry_tredel!(s)                                   # cuts.f:255-275 CUTS-entry zero-PROB TREDEL (+RDTDEL, +FMKILL crown carry)
    rem = cuts!(s; fint = fint)                             # CUTS — thin (accrues econ per cut tree; stashes AUTOES XTES)
    # comcup.f:103-140: when COMCUP deletes records (NDEL>0) it re-runs SPESRT … DENSE, so the growth DGF reads a fresh
    # PTBALT for the moved records (TREMOV does not carry PTBALT). AK only here (base code; other variants untested):
    # FIA 644916319126144 cycle 4 — 25 PROB≤1E-5 records deleted, ES020605 moved into slot 511, live PBAL 255.53
    # (its own) vs jl 351.94 (the deleted slot-511 record's) ⇒ WK2 −1.930 vs −2.002 ⇒ DG/LTHG/HTG ⇒ mortality.
    ncomcup = s.variant isa SoutheastAlaska ? count(i -> s.trees.tpa[i] <= 1f-5, 1:s.trees.n) : 0
    comcup!(s.trees; onmove = _record_move_hook(s))         # COMCUP (grincr.f:391): PROB≤1E-5, after CUTS, before growth
    ncomcup > 0 && compute_density!(s)
    _trim_crown_bypass!(s)
    rem.tpa > 0f0 && compute_density!(s)                    # recompute post-thin density
    if s.fire !== nothing && s.fire.active
        apply_salvage!(s)                                  # SALVAGE (act 2520) — remove snags (FMSALV from CUTS)
        apply_fuelmove!(s)                                 # FUELMOVE (act 2530) — transfer fuel between pools (FMTRET)
        apply_pileburn!(s)                                 # PILEBURN (act 2523) — pile/jackpot burn (FMTRET)
        # FVS runs the cut-phase FFE activities (FMSALV/FMTRET, in cuts.f) BEFORE FMMAIN's FMBURN, so the
        # salvaged snags' released CWD2B crown debris (and any fuel moved) is part of the down-wood the fire
        # samples. The summary driver stashes `fire_smlg` at the cycle START (pre-salvage); nothing between
        # then and here touches `cwd`, so on a fire cycle re-stash it to reflect the post-salvage down wood.
        fuel_period !== nothing && (s.fire.fire_smlg = _small_large_fuel(s.fire))
        # FFE-init year (non-fire), DEFERRED from the pre-grow driver: FVS FMMAIN loads the initial dead-fuel
        # pools (FMCBA) AFTER the cut phase, so the one-time load reads the POST-THIN stand (matches live PERCOV).
        # No-op for any stand without an init-year thin (pre==post state) ⇒ eastern FFE unaffected.
        ffe_init_period !== nothing && ffe_fuel_update!(s, Int(ffe_init_period))
    end
    econ_on && econ_status!(s, Int(s.control.cycle) + 1, 1)   # ECSTATUS(…,1) after CUTS (grincr.f:370)
    if econ_on
        yr = current_cycle_year(s)   # IY schedule (TIMEINT/CYCLEAT-aware)
        s.econ.base_year < 0 && (s.econ.base_year = Int32(yr))
        (s.econ.cycle_cost > 0f0 || s.econ.cycle_rev > 0f0) &&
            push!(s.econ.harvests, (Float32(yr), s.econ.cycle_cost, s.econ.cycle_rev))
    end
    # FFE: the standing-snag falldown (update_snags!) now runs in ffe_fuel_update!'s annual loop
    # INTERLEAVED with the FMCWD decay (FVS FMMAIN order FMSNAG→FMCWD→FMCADD), so the freshly-fallen
    # bole decays in the year it falls. ffe_fuel_update! is called every cycle before grow_cycle! with
    # the same current_cycle_year, so the snag year-accounting is unchanged. (Was here as a single
    # nyrs-at-once call, which let the fresh bole skip its cycle's decay → the DDW size-4/5 overshoot.)
    apply_volume_overrides!(s; fint = fint)  # VOLUME/BFVOLUME merch-standard overrides (volkey.f)
    t = s.trees
    nlive = t.n                              # ORIGINAL live records (pre-tripling)
    # Cycle-start volume + TPA of the originals, for the period accounting.
    old_cfv = Float32[t.cuft_vol[i] for i in 1:nlive]
    old_tpa = Float32[t.tpa[i]      for i in 1:nlive]
    # Tripling is active only for the first ICL4 cycles (s.control.icl4; default 2, set to 0
    # by NOTRIPLE / to n by NUMTRIP); afterwards growth is the stochastic serial-correlation path.
    # grincr.f:31 LTRIP = (ICYC.LE.ICL4 .AND. ITRN.LE.(MAXTRE/3) .AND. .NOT.NOTRIP): the ITRN≤MAXTRE/3 guard is
    # LOAD-BEARING — without it a stand whose (tripled) record count would exceed MAXTRE overflows the MAXTRE-sized
    # tree arrays (intermittent SIGSEGV in the volume loop `1:(t.n+t.ndead)`, e.g. dense ttt01 realizations at ~1000+
    # live records × 3). jl stores the dead block UPWARD (t.n+1 … t.n+ndead), unlike FVS's downward IREC2…MAXTRE, so
    # the tripled live block (3·nlive) plus the dead block must fit MAXTRE: nlive ≤ (MAXTRE−ndead)/3. Reduces to
    # FVS's MAXTRE/3 when ndead=0 (the common case); tighter only when inventory dead records are present.
    trip = !notrip_start && Int(s.control.cycle) < Int(s.control.icl4) && max(nlive, itrn_grincr) <= (variant_maxtre(s.variant) - Int(t.ndead)) ÷ 3   # (ON MAXTRE=6000) NOTRIP (prior-cycle COMPRESS) suppresses tripling
    crown_sdi = stand_sdi_reineke(s)   # pre-growth Reineke SDI for CROWN's RELSDI (SDIBC, grincr.f:241)
    # grincr.f:240/322 SDICAL(0,…) sets the common BAMAX = XMAX·0.5454154·PMSDIU every cycle (sdical.f:203-204, unless the
    # user BAMAX); MORTS's SDICAL overwrites it later, but a stand with no records at MORTS (bare-ground PLANT, cycle 1)
    # keeps this value for the cycle-end CROWN's PP RELSDI=BA/BAMAX (tt/crown.f:179-181). There CRATET skipped SDICHK
    # (cratet.f:731 ITRN≤0 ⇒ GO TO 500), so PMSDIU is still the grinit PERCENT (85) until cycle-1 MORTS divides it
    # (morts.f:194-195): XMAX=1 (no BA) ⇒ BAMAX=46.36. MEASURED FVStt_g16 bare PLANT cycle 1: CROWN sp10 RELSDI 0.2016 =
    # BA 9.3441/46.36; jl left BAMAX 0 ⇒ the SDI form ⇒ RELSDI 0 ⇒ top-rank PP crowns 90 vs live 84.
    if s.variant isa Teton && s.control.ba_max <= 0f0
        pm = s.plot.pct_sdimax_mort_hi > 0f0 ? s.plot.pct_sdimax_mort_hi : 0.85f0
        (Int(s.control.cycle) == 0 && s.trees.n == 0) && (pm = pm * 100f0)
        s.control.sdical_bamax = stand_sdimax(s) * 0.5454154f0 * pm
    end
    # BC/EM V2 mortality WK1 (morts.f) = DG(I) at the START of dgdriv (dgdriv.f:141/144 WK1=DG), i.e. the PRE-prediction
    # DG: the measured input increment at cycle 1, or the previous cycle's DG later. Snapshot it before
    # diameter_growth! overwrites diam_growth. Without it WK1=0 ⇒ the Hamilton G collapses to the DGT floor ⇒
    # RIP over-predicts ⇒ over-kill. (KT/IE/TT/CI use the post-update snapshot at :791, which misses cycle-1's measured DG.)
    # EM: the .tre carries a measured past DG (intree.f:151 reads DG(I) into diam_growth); FVS's cycle-1 WK1 is that
    # value verbatim (oracle FVSem_g16 WK1={1.0,2.3,0.6,0.7} == em_LM.tre DG field). jl's post-cycle snapshot left
    # dg_prev=0 at cycle 1 ⇒ the LM/added-species Hamilton G collapsed ⇒ ~2.5× first-cycle mortality over-kill.
    # ON: the same dgdriv.f DO 5 `WK1(I)=DG(I)` — ON's MORTS does not read it, but a record MORTS empties keeps it: VOLS
    # skips P≤0 (vols.f:125), so its FVS_TreeList MCuM (WK1·FT3toM3) is this DG (compute_volumes_on!).
    (s.variant isa BritishColumbia || s.variant isa EasternMontana || s.variant isa Ontario) &&
        (@inbounds for i in 1:t.n; t.dg_prev[i] = t.diam_growth[i]; end)
    (s.variant isa Ontario && Int(s.control.cycle) == 0 && length(s.calib.dub_wk2) == t.n) &&
        (@inbounds for i in 1:t.n; t.dg_prev[i] = on_do220_dg(s, i); end)   # cycle 1's WK1 = the LSTART DO-220 DG
    (s.variant isa EasternMontana && Int(s.control.cycle) == 0) && em_cycle0_wk1!(s)   # dgdriv.f DO 220 precedence
    # IE: same ie/morts.f Hamilton path (WK1=DG at dgdriv.f:142). ie/morts.f:273 override
    # (ICYC.EQ.1 .OR. WK1==0) .AND. DG>0.5 ⇒ G=DG/(BARK·10) MASKS WK1=0 for every measured-DG tree whose
    # PREDICTED cycle-1 DG>0.5 (verified on iet01 STDINFO: all 27 measured-DG trees fire the override). The
    # residual gap is the tree carrying a measured past DG>0 whose PREDICTED cycle-1 DG≤0.5 (override does NOT
    # fire): FVS's G reads WK1=measured-DG, jl read dg_prev=0 ⇒ Hamilton G collapse ⇒ cycle-1 over-kill. This
    # bites FIA stands too — FVS_TREEINIT carries a measured DG (oracle FVSie_g16 stand 373781950489998: 38/44
    # trees WK1>0, ~10 with predDG≤0.5). Fix: snapshot the measured DG (intree.f DG field, held in diam_growth
    # pre-growth). CYCLE-0 ONLY: at cycle≥1 the post-update snapshot at :793 already gives WK1=DG (this cycle's
    # applied DG → next cycle's vigor) exactly; an every-cycle snapshot here instead OVER-KILLS the shelterwood
    # auto-regen path (iet01 THN3) by racing establishment/tripling churn — measured as ~2.8× cyc-2040 mort.
    # A/B (373781950489998, no regen): moves TPA toward oracle every cycle (2025 2328→2333 vs 2337, …).
    # The faithful cycle-1 WK1 is ie/dgdriv.f's DO-220 result (ie_cycle0_wk1!, from the calibration's DGF(WK3) re-call
    # and the calibration OLDRN — so BEFORE diameter_growth! advances OLDRN). The snapshot + post-DGDRIV stand-in below
    # is kept only as the fallback when that calibration stash is absent.
    _ie_wk1_do220 = s.variant isa InlandEmpire && Int(s.control.cycle) == 0 && length(s.calib.dub_wk2) == t.n
    _ie_wk1_do220 && ie_cycle0_wk1!(s)
    (s.variant isa InlandEmpire && Int(s.control.cycle) == 0 && !_ie_wk1_do220) &&
        (@inbounds for i in 1:t.n; t.dg_prev[i] = t.diam_growth[i]; end)
    # CI: ci/dgdriv.f:169-172 WK1(I)=DG(I) at the top of DGDRIV — at cycle 1 that is the DO-220 calibration DG (the
    # measured increment, capped at the inside-bark DBH when IDG<2; 0 at HT≤4.5; else the DGF dub), which ci/morts.f
    # reads as the vigor term G=WK1/(BARK·OLDFNT). jl fed its own cycle-1 prediction: FIA 753188889290487 LP (past DBH
    # 6.8→8.0, DG_MEASURE 10) WK1 0.47 vs live ~1.13 ⇒ G halved ⇒ LP cycle-1 kill 1.152 vs live 0.568 of 6.
    (s.variant isa CentralIdaho && Int(s.control.cycle) == 0 && length(s.calib.dub_wk2) == t.n) &&
        (@inbounds for i in 1:t.n; t.dg_prev[i] = ci_do220_dg(s, i, t.dbh[i]); end)
    # DFTM DFTMGO+TMBMAS predict seam (grincr.f:402/424, BEFORE DGDRIV): on a scheduled tussock-moth
    # outbreak this cycle, gate on host presence and compute the IBMTYP=2 foliage biomass/percent-new
    # from the PRIOR-cycle DG (t.diam_growth still holds it here) for the gradd TMCOUP coupler. Inert
    # (no-op, byte-identical) unless a DFTM block is active and an outbreak is due.
    dftm_predict!(s)
    # MPSVDG (mpgr.f): save the pre-growth DG for the LPOPDY MPGR resistance, BEFORE diameter_growth!
    # overwrites diam_growth. Inert unless an LPOPDY MPB block is active.
    s.mpb !== nothing && mpb_svdg!(s)
    rd_cycle_start!(s)                     # WRD: size driver (RDESTB for last cycle's regen) + WK1=DG snapshot (dgdriv.f)
    stash = diameter_growth!(s, s.variant; tripling = trip, sfint = fint)  # DGs only; no records yet
    # IE cycle-1 WK1 dub (dgdriv.f:755-795 LSTART "DUB IN DBH INCREMENT FOR TREES ON WHICH IT WAS NOT
    # MEASURED"): the calibration pass sets DG(I) per dgdriv.f:774-795, and that value becomes cycle-1 WK1
    # (dgdriv.f:142 WK1(I)=DG(I) before DGF recomputes DG). The FVS precedence, EXACTLY:
    #   • HT≤4.5 (seedling)          ⇒ DG(I)=0  (dgdriv.f:784-785) ⇒ WK1=0   [unconditional — the keep-measured
    #                                   branch at :774 requires HT>4.5, so a sub-breast-height tree is ALWAYS zeroed]
    #   • HT>4.5, measured DG>0       ⇒ WK1 = measured DG           (dgdriv.f:774-782)
    #   • HT>4.5, no measured DG      ⇒ WK1 = DGF dub                (dgdriv.f:786-792)
    # The snapshot above (665-666) copied `diam_growth` = the measured input DG, which for a FIA/no-remeasurement
    # stand is the −1 MISSING sentinel (apply_growth_input_types!, IDG=1). jl's own DGF prediction (diam_growth
    # just filled by diameter_growth!) is the faithful analog of the FVS dub for the HT>4.5 unmeasured case.
    # ★ The HT≤4.5 ⇒ 0 branch is LOAD-BEARING for dense seedling stands: a sub-breast-height FIA seedling carries a
    # nominal DBH (~0.1") but HT<4.5, so FVS gives it WK1=0 ⇒ Hamilton G collapses to the DGT floor (0.05/BARK) ⇒
    # HIGH cycle-1 mortality. Without this branch jl fed WK1=predicted-DG (often >0.5 for a vigorous DF seedling),
    # inflating 11.2007·G + 6.07129·G/D ⇒ RIP collapses ⇒ massive under-kill (oracle 39607788010690 cyc-1 drop
    # 1614 vs jl 344). Verified vs FVSie_g16: all 3 cyc-1 seedling records HT=1.01 ⇒ WK1=0 ⇒ mortG≈0.079. CYCLE-0 only.
    (s.variant isa InlandEmpire && Int(s.control.cycle) == 0 && !_ie_wk1_do220) &&
        (@inbounds for i in 1:t.n
            if t.height[i] <= 4.5f0
                t.dg_prev[i] = 0f0                                   # dgdriv.f:784-785
            elseif t.dg_prev[i] <= 0f0
                t.dg_prev[i] = t.diam_growth[i]                      # dgdriv.f:786-792 (DGF dub analog)
            end
        end)
    # CR dwarf mistletoe diameter growth-loss (misdgf.f, dgdriv.f:230): DG·=DGPDMR(sp,DMR); applied to the
    # central + tripled DGs right after the DG driver, using START-of-cycle DMR (before cr_mistoe! spread).
    s.variant isa CentralRockies && cr_dm_growth_loss!(s, stash)
    _dm_effects_variant(s.variant) && ie_dm_growth_loss!(s, stash)   # western MISTOE DG-loss (misdgf.f) — IE+KT/EM/BM/UT/TT/CI + BC(NEWSPRED); START-of-cycle DMR
    # IE REGENT(LESTB) (ie_esgent!) reads the START-of-cycle (post-thin, pre-growth) stand density as TEMAHT/
    # TEMBA/TEMCCF (grincr.f:316-320 ATAVH/ATBA/ATCCF → regent.f:218-222; regent.f:326-328 AH=TEMAHT,R=TEMCCF;
    # :1082-1083 BANEXT=TEMBA,RDNEXT=TEMCCF) — NOT the post-growth density. Capture it here (the same values
    # small_tree_growth! reads) before growth+establishment recompute avg_height upward; ie_esgent! uses it for
    # the DADJ AH-vs-4.5 relh clamp (else on a stand whose trees grew past 4.5 ft, relh<0⇒0 ⇒ DADJ=0.65 instead
    # of the pre-growth all-seedling relh=1 ⇒ DADJ=0.65−DELMAX≈0.77: every crossing seedling's dubbed DBH ~0.12"
    # too small, halving the birth-cohort BA — the dominant plant-regime residual). Same #194 ATAVH class CI
    # already fixed. IE-only (only ie_esgent! consumes this trio; every other variant path is unchanged).
    es_at_avh = s.plot.avg_height
    es_at_ba = s.plot.basal_area
    es_at_relden = s.plot.relative_density
    height_growth!(s, s.variant; scale = fint / htg_period(s.variant))   # HTG scaled to cycle (YR: SN=5, NE=10)
    # IE htgf.f (317-347) recomputes each TRIPLED large-tree copy's HTG from the copy's spread DG (the
    # NI-section species); height_growth! only computed the central record's HTG, leaving the copies flat.
    # Deterministic (no RNG) ⇒ stream untouched. Restores the copy height spread the oracle produces.
    s.variant isa InlandEmpire && ie_triple_htg!(s, stash; scale = fint / htg_period(s.variant))
    s.variant isa EasternMontana && em_triple_htg!(s, stash; scale = fint / htg_period(s.variant))
    s.variant isa CentralIdaho && ci_triple_htg!(s, stash; scale = fint / htg_period(s.variant))
    s.variant isa Kootenai && kt_triple_htg!(s, stash; scale = fint / htg_period(s.variant))
    s.variant isa WestSierra && ws_triple_htg!(s, stash)   # ws/htgf.f:921-956 copy HTG = TEMHTG·DG(copy)/DG
    small_tree_growth!(s, stash, s.variant; fint = fint)  # REGENT overrides DG/HTG for small trees (SN <3", NE <5")
    apply_fix_scalers!(s, stash, :fixdg, fint)   # FIXDG/FIXHTG: one-shot DG/HTG scalers,
    apply_fix_scalers!(s, stash, :fixhtg, fint)  # after all growth, before MORTS (grincr.f:451)
    # The report driver's hook onto the FMMAIN point (gradd.f:118 — after REGENT's direct small-tree DBH, before UPDATE)
    fmmain_hook === nothing || fmmain_hook(s, stash)
    # CR dwarf mistletoe spread/intensification (mistoe.f MISTOE, gradd.f:96 — after growth+FIXHTG, before
    # UPDATE; uses HTG). Updates per-tree DMR, drawing rann! in ISCT order (RNG-aligned to FVS). No-op for
    # non-CR and for mistletoe-free stands (SMR=0 ⇒ zero draws). The DM mortality it enables is max-combined
    # in mortality! (below); the DM diameter growth-loss is applied in diameter_growth!.
    cr_mistoe!(s; fint = fint)
    # MISTOE TRIPLING SEAM (mirrors the RD rd_post_triple pattern below). FVS runs MISTOE at gradd.f:96 —
    # in GRADD, AFTER GRINCR's MORTS+TRIPLE (grincr.f:535/543) — so on a TRIPLING cycle the spread draws its
    # per-infected-tree rann! on the ALREADY-TRIPLED record list (ITRN×3), not the un-tripled ITRN. jl ran
    # ie_mistoe! here (pre-mortality/pre-triple), so on a tripling cycle it drew 1/3 of FVS's mistletoe draws
    # (measured stand 1143092700290487: jl 2 vs oracle 6 at cyc0, jl 6 vs oracle 18 at cyc1) — desyncing the
    # SHARED main rann! stream for every downstream DGSCOR/REGENT draw (the #206 straddle). DEFER the spread to
    # the post-triple block (right after mortality_and_fire!) when this cycle triples; the pre-mortality dmr is
    # unchanged for the DM growth-loss (start-of-cycle, applied at diameter_growth!). DM mortality is NOT part of
    # MORTS: MISMRT is called inside MISTOE after the spread and MISINF (mistoe.f:517-522), so on a tripling cycle
    # it too moves post-triple (mis_post below). Non-tripling cycles keep the existing seam (MORTS/TRIPLE draw no
    # rann!, and spread → MISINF → DM max-combine in mortality! is the FVS order) — byte-identical there.
    mis_defer = _ie_mis_variant(s.variant) && (stash !== nothing)
    # mis_post: on a NON-fire tripling cycle the WHOLE GRADD MISTOE call (mistoe.f: spread → MISINF :517 →
    # MISMRT :522) runs post-TRIPLE on the tripled records at full pre-UPDATE PROB (the post-triple block
    # below) — not just the spread. MISINF must follow the spread (else the spread intensifies the
    # freshly-forced DMR in the same cycle) and MISMRT must read the post-spread DMR on the tripled PROB
    # (else the DM kill uses the pre-spread DMR). dm_mrt_defer keeps MORTS from applying the DM kill.
    # Fire cycles triple inside mortality_and_fire! and keep the previous order.
    mis_post = mis_defer && !_fire_due(s)
    if !mis_defer
        _ie_mis_variant(s.variant) && ie_mistoe!(s; fint = fint)   # western MISTOE spread (mistoe.f) — shared across N-Rockies Wykoff
    end
    mis_post || dm_misinf!(s)   # MISTPINF forced initial DM infection (misinf.f MISINF, mistoe.f:517 — after spread, before DM mortality); inert w/o a card
    # BC NEWSPRED spatial dwarf-mistletoe spread (canada/newmist DMTREG) — updates per-tree DMR via the
    # spatial model, then publishes ms.dmr→t.dmr for the base misdgf/mismrt effects. Self-guards on the
    # NEWSPRED/MISTOE keyword (ms.active||newmod); inert on non-DM BC stands. lastyr = cycle length (yr).
    s.variant isa BritishColumbia && dm_tregro!(s, round(Int, fint))
    # WRD RDTREG seam ORDERING. FVS runs the ENTIRE root-disease chain (RDCNTL: RDINSD/RDSPRD/RDINF/RDMORT/
    # RDSTP, then RDEND, then RDGROW) in GRADD (gradd.f:131) AFTER GRINCR's TRIPLE (grincr.f:543). MORTS only
    # sets the WK2 kill array (PROB is untouched); TRIPLE (triple.f:67/129) splits PROB *and* WK2 proportionally;
    # RDTREG's RDINSD/RDMORT/RDSTP therefore see the ALREADY-tripled records at FULL pre-mortality PROB, and
    # UPDATE (gradd.f:180) subtracts WK2 only afterwards. Running the chain on the un-tripled list (ITRN, not
    # ITRN×3) coarsens RDINSD's NUMTRE=INT(IRINIT/ITRN) Monte-Carlo (rdinsd.f:290) → over-infection → an
    # undersized cyc-1 stump the faithful rd_inoc_decay! destroys → disease fizzles. So when this cycle BOTH
    # triples (stash) AND has no SIMFIRE, DEFER the whole chain to the post-triple block below (rd_post_triple);
    # otherwise (no tripling, or a fire cycle) FVS's RDTREG also runs un-tripled, so keep the pre-triple call.
    rd_stand = (s.root_disease !== nothing) && rd_active(s.root_disease) &&
               s.root_disease.iroot != 0 && s.root_disease.driver !== nothing
    rd_post_triple = rd_stand && (stash !== nothing) && !_fire_due(s)
    rd_post_triple || root_disease_treg!(s, fint)  # WRD gradd.f RDTREG seam — inert unless an RDIN block is active
    # FFE SIMFIRE this cycle? FVS computes MORTS (GRINCR) on the FULL pre-fire stand into WK2,
    # then GRADD's FMKILL sets WK2(I)=MAX(WK2(I),FIRKIL(I)) (fmkill.f:86) — a tree dies from
    # whichever is LARGER, density/background MORTS or fire, NOT both summed. The old code ran
    # fire FIRST (density mortality then saw the thinned post-fire stand → under-kill on a dense
    # burn, fire_fuel9 2010 TPA 155 vs FVSsn 143). Replicate the MAX combine, each kill measured
    # on the same pre-fire stand. Non-fire path is byte-identical (fk≡0 ⇒ kill = mk).
    # FFE fire cycle: the carbon report (carbon_hook, FVS FMCRBOUT) and the deferred FFE annual fuel update
    # (fuel_period, the FMSNAG/FMCWD/FMCADD loop) run inside mortality_and_fire! right after FMBURN and
    # before the WK2 combine — the FVS fmmain.f order (FMBURN:170 → FMCRBOUT:206 → annual loop:228 → FMKILL).
    # So the fire consumes/snags the START-of-cycle fuels (the report driver withholds the pre-grow
    # ffe_fuel_update! for the fire cycle and hands its period here). Non-fire cycles pass neither ⇒ no-op.
    pf = (carbon_hook !== nothing || fuel_period !== nothing) ?
         (st -> (carbon_hook === nothing || carbon_hook(st);
                 fuel_period === nothing || ffe_fuel_update!(st, fuel_period))) : nothing
    # FIRE cycle (FVS): MORTS on the originals → TRIPLE → fire on the tripled set → FMKILL MAX-combine.
    # mortality_and_fire! does that internally and returns its OMORT + `tripled` so we don't TRIPLE twice;
    # the NON-fire path keeps MORTS-then-TRIPLE here (VARMRT must see the un-tripled ITRN records).
    s.control.dm_mrt_defer = mis_post
    # MISTOE post-TRIPLE seam (mis_post, below): FVS books the cycle's snags at FMKILL(2) (gradd.f → FMSADD(YEAR,4))
    # from the FINAL post-TRIPLE WK2 — AFTER MISMRT has MAX-combined the dwarf-mistletoe kill (mistoe.f:522). Booking
    # inside mortality! (pre-TRIPLE, pre-MISMRT) dropped every DM-killed tree from the snag pools (EM 196378260020004
    # cycle 1: LP snags 11.77 vs live 15.38 TPA — the 3.6 TPA MISMRT added). Defer the booking to that seam.
    # Every NON-fire tripling cycle books its snags the same way (#277): FMKILL(2) (fmkill.f:129-143) hands FMSADD each
    # TRIPLED record's WK2 (the per-original MORTS kill ·.60/.25/.15, after MISMRT/BRTREG/RDEND), in record order — not
    # the un-tripled kill inside MORTS. FMSADD's zero-initialized class means (fmsadd.f:334-340) round differently on
    # the three parts (MEASURED SN carbon_jenkins: a 12" record's DBHS 11.999999 from the whole kill, 12.0 from the parts).
    post_book = stash !== nothing && !_fire_due(s)
    (mortf, tripled) = mortality_and_fire!(s; fint = fint, stash = stash, post_fire = pf, book_snags = !post_book)
    s.control.dm_mrt_defer = false
    # WRD rd/rdend.f: reconcile the RD infected-tree kill (RRKILL) with FVS's just-applied
    # MORTS WK2 (= old_tpa − t.tpa) and re-apply the RD-adjusted WK2 — FVS runs RDEND at
    # MORTS time (GRINCR MORTS → GRADD RDTREG/RDEND). Non-fire, non-tripled RD path only;
    # gated so a no-RD stand is byte-identical (root_disease === nothing ⇒ no-op).
    if !tripled && !rd_post_triple && (s.root_disease !== nothing) && rd_active(s.root_disease) &&
       s.root_disease.iroot != 0 && s.root_disease.driver !== nothing
        rd_end_apply!(s.root_disease, s, old_tpa)   # RDEND (un-tripled path only; tripled RD runs post-triple below)
    end
    # DFB dfb/dfbdrv.f (DFBDBH→DFBMOD→DFBMRT) gated by DFBGO: on a cycle with a scheduled Douglas-fir
    # Beetle outbreak, raise the large-DF WK2 mortality to MAX(background, DFKILL). FVS calls DFBDRV in
    # GRADD after MORTS has set WK2 (gradd.f:74), so this sits with the RD seam, on the non-tripled
    # cycle stand, reading cycle-start old_tpa/DBH. Inert (no-op, byte-identical) unless a DFB block is
    # active and an outbreak is due this cycle.
    if !tripled && s.dfb !== nothing && s.dfb.active
        dfb_win!(s, old_tpa)               # DFBWIN windthrow (gradd.f:72) — WK2 windthrow kill + OKILL feed, BEFORE DFBDRV
        dfb_apply!(s, old_tpa, fint)       # DFBDRV (gradd.f:74)
        s.dfb.okill = 0.0f0                # DFBMRT clears OKILL each cycle; guard leaks if DFBDRV early-returns
    end
    # DFTM tmcoup.f (GARBEL→DFTMOD→mortality+growth-loss+top-kill) gated by DFTMGO: on a scheduled
    # tussock-moth outbreak, raise host WK2 mortality and reduce DG/HTG + apply top-kill. FVS calls
    # TMCOUP in GRADD after MORTS/MISTOE (gradd.f:103), so this sits after the DFB block on the
    # non-tripled cycle stand, reading cycle-start old_tpa. Inert (byte-identical) unless a DFTM block
    # is active and dftm_predict! armed the coupling this cycle.
    if !tripled && s.dftm !== nothing && s.dftm.active
        dftm_apply!(s, old_tpa, fint)      # TMCOUP (gradd.f:103)
    end
    # WPBR wpbr/brtreg.f (BRTREG, gradd.f:126): after MORTS/MISTOE (and the DFB/DFTM
    # couplers), grow this cycle's blister-rust cankers and impose canker mortality —
    # a fully girdled bole canker sets WK2=PROB·0.99999 (t.tpa override), a top-kill
    # canker lowers ITRUNC/NORMHT and reduces the crown. Per-record canker state
    # persists across cycles in w.recs. Non-fire, non-tripled path only; inert
    # (byte-identical) unless a BRUST block is active with a host pine present.
    # On a TRIPLING cycle FVS's BRTREG (gradd.f:126) runs AFTER GRINCR's TRIPLE, on the tripled records at full
    # pre-UPDATE PROB (BRANN draws per tripled record) — see the post-triple blocks below (br_post).
    br_on = s.wpbr !== nothing && (s.wpbr::WpbrState).active
    br_post = br_on && !tripled && stash !== nothing
    if !tripled && br_on && !br_post
        wpbr_brtreg!(s, fint, old_tpa)
    end
    # WWPB (Westwide Pine Beetle, wwpb/*.f): the landscape Parallel-Processing-
    # Extension outbreak. The PPE spatial multi-stand harness (PPMAIN/ALSTD2/SPLAEX)
    # is ABSENT from this FVS tree, so the outer orchestration is a USER-approved
    # (2026-08-21) reconstructed single-stand harness (wwpb_outbreak_cycle!); every
    # beetle kernel it composes is bit-exact vs the pristine wwpb/*.f. INERT unless
    # a DISPERSE keyword activated the outbreak (w.outbreak) — a BMIN-only stand
    # still projects byte-identically. wwpb_apply! runs one cycle's outbreak and
    # reconciles the beetle kills into the FVS mortality (BMKILL/WK2 handback).
    # PPE mode-2 LIVE landscape coupling (ppe_landscape.jl ppe_run_landscape_live!): when a
    # `wwpb_barrier` closure is threaded in, it REPLACES the single-stand wwpb_apply! at this exact
    # seam (post-MORTS/pre-GRADD, records un-tripled) — the barrier rendezvouses every landscape
    # stand here, runs the bmdrv_multi! dispersal ONCE across the whole landscape, and hands each
    # stand's beetle kill back. `wwpb_barrier === nothing` (every ordinary run) is byte-identical to
    # the prior single-stand path (the guard test_multicycle stays 339/11).
    if !tripled
        if wwpb_barrier !== nothing
            wwpb_barrier(s, old_tpa, fint)
        elseif s.wwpb !== nothing && (s.wwpb::WwpbState).outbreak
            wwpb_apply!(s, old_tpa, fint)
        end
    end
    # LPMPB (Mountain Pine Beetle, lpmpb/*.f): stand-level Cole rate-of-loss
    # mortality (MPBGO→MPBCUP→COLDRV, gradd.f:63). Inert (byte-identical) unless an
    # MPB block is active, a scheduled outbreak (OPFIND 555) is due this cycle, the
    # variant has a lodgepole host, and the stand meets the MPBER minimum condition.
    if !tripled && s.mpb !== nothing && (s.mpb::MpbState).active
        mpb_apply!(s, old_tpa, fint)
    end
    # WSBWE (Western Spruce Budworm, wsbwe/*.f): stand-level DEFOL defoliation
    # (BWEGO→BWECUP→BWEDR→BWEDAM/BWEPDM, grincr.f:414 / gradd.f:108). The reader is
    # ported and the deterministic payload kernels (wsbwe_rdds/wsbwe_rhtg/wsbwe_mort_pr)
    # are dump-replay validated bit-exact vs FVSem_wsbwe. INERT: wsbwe_apply! early-
    # returns until the BWESIT→BWEAGE→BWEDEF→BWEDAM foliage/PRBIO feeder + the BWEPDM
    # per-tree apply (draws the damage RNG) are ported (scratchpad/wsbwe/HANDOFF.md
    # "NEXT CHUNKS"), so this projects byte-identically even when the gate fires.
    if !tripled && s.wsbwe !== nothing && (s.wsbwe::WsbweState).active &&
       wsbwe_go(s.wsbwe::WsbweState)
        wsbwe_apply!(s, old_tpa, fint)
    end
    # tripled MORTS WK2 per record (triple.f WEIGHT split of the per-original kill) — BRTREG's WK2 input
    _wk2_trip(wk2_u, nl) = (v = zeros(Float32, 3nl);
                            @inbounds(for i in 1:nl; v[i] = wk2_u[i] * 0.60f0; v[nl+2i-1] = wk2_u[i] * 0.25f0;
                                                     v[nl+2i] = wk2_u[i] * 0.15f0; end); v)
    g = s.plot.gross_space
    # Mortality volume (OMORT): MORTS deaths AND the fire kill (the MAX per record), reduced t.tpa from
    # the cycle-start old_tpa at the same cycle-start CFV. Fire cycle: computed inside (on the tripled set).
    mort = mortf
    if !tripled
        mort = 0f0
        # FVS_TreeList MortPA (dbstrls.f DP=WK2/GROSPC) and OMORT (Σ WK2·CFV) read MORTS's WK2 itself, not the
        # PROB−(PROB−WK2) difference, which rounds to the survivor's ULP (±10-20 ULP of WK2 on small kills).
        wk2_0 = _morts_wk2(s, old_tpa, nlive)
        _on_wk = s.variant isa Ontario
        @inbounds for i in 1:nlive
            m = wk2_0[i]
            mort += (_on_wk ? min(m, old_tpa[i]) : m) * old_cfv[i]   # update.f:72-74 WK6 = MIN(WK2,PROB)·CFV
            t.mort_pa[i] = m                   # per-record period mortality (FVS_TreeList MortPA), pre-TRIPLE
        end
        if mis_post && !rd_post_triple
            # ==== FVS-faithful MISTOE seam on a tripling cycle (gradd.f:96, after MORTS+TRIPLE, before UPDATE):
            # recover the MORTS kill (WK2), restore full PROB, TRIPLE (splitting PROB), run the spread and MISINF on
            # the tripled full-PROB records, then survivors = PROB − WK2·weight and MISMRT MAX-combines the DM kill.
            wk2_u = _morts_wk2(s, old_tpa, nlive)
            @inbounds for i in 1:nlive; t.tpa[i] = old_tpa[i]; end
            triple_records!(s, stash)
            n2 = t.n
            full_prob = Float32[t.tpa[i] for i in 1:n2]
            ie_mistoe!(s; fint = fint)         # mistoe.f spread (rann! over ITRN×3)
            dm_misinf!(s)                      # mistoe.f:517 MISINF
            @inbounds for i in 1:nlive
                t.tpa[i]          = full_prob[i]          - wk2_u[i] * 0.60f0
                t.tpa[nlive+2i-1] = full_prob[nlive+2i-1] - wk2_u[i] * 0.25f0
                t.tpa[nlive+2i]   = full_prob[nlive+2i]   - wk2_u[i] * 0.15f0
            end
            wk2t = _wk2_trip(wk2_u, nlive)            # triple.f WK2·WEIGHT per record
            ie_dm_mismrt_post!(s, full_prob, fint; wk2 = wk2t)   # mistoe.f:522 MISMRT → WK2=MAX(WK2,PROB·rate)
            br_post && wpbr_brtreg!(s, fint, full_prob; wk2_hint = _wk2_trip(wk2_u, nlive))   # gradd.f:126 BRTREG
            mort = 0f0
            @inbounds for c in 1:n2
                # MortPA/OMORT = WK2; a record BRTREG re-killed falls back to the survivor difference
                m = t.tpa[c] == full_prob[c] - wk2t[c] ? wk2t[c] : full_prob[c] - t.tpa[c]
                mort += m * t.cuft_vol[c]
                t.mort_pa[c] = m
            end
        elseif rd_post_triple
            # ==== FVS-faithful WRD seam: the whole RD chain on the TRIPLED, FULL pre-mortality PROB list ====
            # Mirror gradd.f: MORTS set WK2 (jl applied it eagerly → t.tpa are survivors); TRIPLE splits FULL
            # PROB and WK2 proportionally; RDTREG (RDCNTL RDINSD/RDMORT/RDSTP) runs on the tripled full-PROB
            # records; RDEND combines RRKILL into WK2; UPDATE subtracts WK2. Reconstruct here: recover the
            # MORTS kill, restore t.tpa→full PROB, triple, RDTRIP the driver, run RDCNTL on full PROB, then
            # re-apply the (proportionally-tripled) WK2 and let RDEND fold in the RD kill.
            wk2_u = _morts_wk2(s, old_tpa, nlive)                          # per-original MORTS kill (WK2)
            @inbounds for i in 1:nlive; t.tpa[i] = old_tpa[i]; end       # restore full pre-mort PROB
            triple_records!(s, stash)                                    # splits FULL PROB .60/.25/.15
            rd_triple_driver!(s.root_disease, stash.nlive)               # RDTRIP: split RD per-record arrays
            n2 = t.n
            full_prob = Float32[t.tpa[i] for i in 1:n2]                  # tripled pre-mort PROB (= RDTREG input)
            if mis_post                                                  # gradd.f:96 MISTOE precedes :131 RDTREG
                ie_mistoe!(s; fint = fint)
                dm_misinf!(s)
            end
            root_disease_treg!(s, fint)                                  # RDCNTL RDINSD/RDMORT/RDSTP on full PROB
            @inbounds for i in 1:nlive                                   # survivors = PROB − WK2 (triple.f WEIGHT split)
                t.tpa[i]          = full_prob[i]          - wk2_u[i] * 0.60f0
                t.tpa[nlive+2i-1] = full_prob[nlive+2i-1] - wk2_u[i] * 0.25f0
                t.tpa[nlive+2i]   = full_prob[nlive+2i]   - wk2_u[i] * 0.15f0
            end
            wk2t = _wk2_trip(wk2_u, nlive)
            mis_post && ie_dm_mismrt_post!(s, full_prob, fint; wk2 = wk2t)   # MISMRT into WK2 before RDEND
            br_post && wpbr_brtreg!(s, fint, full_prob; wk2_hint = _wk2_trip(wk2_u, nlive))   # gradd.f:126 BRTREG (before RDTREG's RDEND)
            rd_end_apply!(s.root_disease, s, full_prob)                  # RDEND: fold RRKILL into WK2, re-apply
            mort = 0f0                                                   # OMORT + MortPA from the final tripled kill
            @inbounds for c in 1:n2
                # WK2 where BRTREG/RDEND left the record's kill untouched; else the survivor difference
                m = t.tpa[c] == full_prob[c] - wk2t[c] ? wk2t[c] : full_prob[c] - t.tpa[c]
                mort += m * t.cuft_vol[c]
                t.mort_pa[c] = m
            end
        elseif stash !== nothing
            # gradd.f order: MORTS sets WK2 (jl applied it eagerly), TRIPLE splits PROB and WK2 SEPARATELY
            # (triple.f:30-32,75-76 PROB·w, WK2·w), UPDATE then subtracts ⇒ survivor = PROB·w − WK2·w. Splitting the
            # already-reduced TPA, (PROB−WK2)·w, rounds differently (1 ULP on ~40% of records — MEASURED vs live DENSE
            # PROB bits, bare-PLANT IE stand cycle 2). WK2 = the applied MORTS kill (scratch buffer) when it matches the
            # observed reduction, else recovered from old−new TPA.
            wk2_u = _morts_wk2(s, old_tpa, nlive)
            @inbounds for i in 1:nlive; t.tpa[i] = old_tpa[i]; end
            triple_records!(s, stash)          # TRIPLE splits the FULL pre-mortality PROB
            full_prob_p = br_post ? Float32[t.tpa[i] for i in 1:t.n] : Float32[]
            @inbounds for i in 1:nlive
                t.tpa[i]          = max(0f0, t.tpa[i]          - wk2_u[i] * 0.60f0)
                t.tpa[nlive+2i-1] = max(0f0, t.tpa[nlive+2i-1] - wk2_u[i] * 0.25f0)
                t.tpa[nlive+2i]   = max(0f0, t.tpa[nlive+2i]   - wk2_u[i] * 0.15f0)
                t.mort_pa[i] = wk2_u[i] * 0.60f0; t.mort_pa[nlive+2i-1] = wk2_u[i] * 0.25f0; t.mort_pa[nlive+2i] = wk2_u[i] * 0.15f0
            end
            if (s.root_disease !== nothing) && rd_active(s.root_disease) &&
               s.root_disease.iroot != 0 && s.root_disease.driver !== nothing
                rd_triple_driver!(s.root_disease, stash.nlive)
            end
            if br_post                         # gradd.f:126 BRTREG on the tripled full-PROB records
                surv = Float32[t.tpa[i] for i in 1:t.n]
                wpbr_brtreg!(s, fint, full_prob_p; wk2_hint = _wk2_trip(wk2_u, nlive))
                @inbounds for c in 1:t.n       # a BR canker kill (WK2=PROB·0.99999) adds to OMORT/MortPA
                    d = surv[c] - t.tpa[c]
                    d != 0f0 && (mort += d * t.cuft_vol[c]; t.mort_pa[c] += d)
                end
            end
        else
            triple_records!(s, stash)          # TRIPLE after mortality (splits surviving TPA)
            # RD driver must be tripled in lockstep (rd/triple.f RDTRIP): FVS RDGROW runs on the
            # already-tripled list, so the RD per-record arrays are split .60/.25/.15 to match.
            if stash !== nothing && (s.root_disease !== nothing) && rd_active(s.root_disease) &&
               s.root_disease.iroot != 0 && s.root_disease.driver !== nothing
                rd_triple_driver!(s.root_disease, stash.nlive)
            end
        end
    end
    # FMKILL(2) (fmkill.f:135-143): the non-fire tripling cycle's snags from the final post-TRIPLE WK2 per record — the
    # WK2 array itself (t.mort_pa), not PROB−survivor, which rounds to the survivor's ULP (FMSADD's density-weighted
    # HTDEAD then drifts: akffe YC 8.4"×5' ⇒ 4.9999995 live).
    post_book && book_mortality_snags!(s, Float32[t.mort_pa[c] for c in 1:t.n], t.n, fint)
    fertilizer_growth!(s; fint = fint)     # FFERT fertilizer DG/HTG boost (grincr.f:564, after TRIPLE)
    s.variant isa Ontario && on_gradd_dg_scale!(s, fint)   # gradd.f:79-90 DG → FINT basis (after GRINCR, before MISTOE)
    # MISTOE post-triple seam. FVS runs MISTOE at gradd.f:96 — in GRADD, AFTER GRINCR's MORTS+TRIPLE — so on a
    # TRIPLING cycle the spread draws its per-host-tree rann! on the ALREADY-TRIPLED record list (ITRN×3). jl
    # ran ie_mistoe! before TRIPLE (the pre-mortality seam above), so on a tripling cycle it drew only 1/3 of
    # FVS's mistletoe draws — measured stand 1143092700290487 (2 DMR-6 western-larch trees): jl drew 2 vs the
    # oracle's 6 at cyc0, 6 vs 18 at cyc1 — desyncing the SHARED main rann! stream for every downstream
    # DGSCOR/REGENT draw (the #206 straddle). `mis_defer` (set at the pre-mortality seam) deferred the spread on
    # a tripling cycle. On a NON-fire tripling cycle (mis_post) the whole MISTOE call — spread, MISINF and MISMRT —
    # already ran in the post-triple mortality block above (MISMRT is called INSIDE MISTOE at mistoe.f:522, so it
    # reads the POST-spread DMR). Only a FIRE tripling cycle still runs the spread here. The DM growth-loss
    # (start-of-cycle DMR, applied at diameter_growth!) is unchanged. Non-tripling cycles keep the pre-mortality
    # seam (MORTS/TRIPLE draw no rann! ⇒ identical RNG position) — byte-identical there.
    (mis_defer && !mis_post) && ie_mistoe!(s; fint = fint)   # fire tripling cycle only (non-fire: the seam above)
    htgstp!(s; fint = fint)                # HTGSTOP/TOPKILL top damage (gradd.f:158, before UPDATE)
    # WRD rd/rdgrow.f (+ tail rd/rdinoc.f decay): reduce the per-record DG/HTG by the infected-
    # root proportion, on the PRE-DBH-update increments (FVS RDGROW runs in RDTREG before UPDATE,
    # after RDEND). Non-tripled RD path only; no-RD stand is byte-identical (no-op).
    if !tripled && (s.root_disease !== nothing) && rd_active(s.root_disease) &&
       s.root_disease.iroot != 0 && s.root_disease.driver !== nothing
        rd_grow_apply!(s.root_disease, s, fint)
    end
    # Per-record cycle-start CFV (tripled records inherit the originals' cycle-0 vol).
    n = t.n
    old_cfv2 = Float32[t.cuft_vol[i] for i in 1:n]
    sd = s.coef.species
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    _cr_up = s.variant isa CentralRockies; _cr_up_imod = _cr_up ? Int(s.plot.model_type) : 0
    _tt_up = s.variant isa Teton   # TT bark = tt_bratio (PP sp10 IMAP=4 power model)
    _bm_up = s.variant isa BlueMountains   # BM bark = bm_bratio (POWER model, per-species groups)
    # gradd.f:205 ABIRTH(I)=ABIRTH(I)+FINT is shared by every variant; jl ages it where something reads ABIRTH:
    # the CR/TT/UT/IE/EM/BM growth models, Climate-FVS BIRTHYR (clgmult/clmorts) in every climate-wired variant, and
    # ON's Mowraski NMV (volont.f:405). It runs AFTER UPDATE/VOLS (see the single aging loop below compute_volumes!).
    _age_up = _cr_up || _tt_up || _bm_up || s.variant isa Utah || s.variant isa InlandEmpire || s.variant isa Ontario ||
              s.variant isa EasternMontana || s.variant isa CentralIdaho || s.variant isa Kootenai ||
              s.variant isa Klamath || s.variant isa PacificNorthwest || s.variant isa WestCascades ||
              s.variant isa EastCascades || s.variant isa SouthCentralOregon || s.variant isa CentralCalifornia ||
              s.variant isa WestSierra || s.variant isa OregonCoast || s.variant isa Olympic ||
              s.variant isa BritishColumbia
    _wc_up = s.variant isa WestCascades   # WC bark = wc_bratio (POWER a·Dᵇ for bark_imap=1; linear cannot express it)
    _pn_up = s.variant isa PacificNorthwest   # PN bark = wc_bratio (POWER, all imap=1) — same as WC
    _ec_up = s.variant isa EastCascades   # EC bark = wc_bratio (per-species bark_imap POWER/linear)
    _op_up = s.variant isa Olympic        # OP bark = op_bratio (op/bratio.f 3-path); the DDS→DG conversion used it,
                                          # so the DBH += DG/bark apply MUST use the same bark (else the round-trip breaks)
    _ak_up = s.variant isa SoutheastAlaska # AK bark = ak_bratio (3-type: power/linear/power)
    _ut_up = s.variant isa Utah    # UT ages ABIRTH (gradd.f:205); CR-surrogate (17:19,22) htgf reads it
    _ie_up = s.variant isa InlandEmpire   # IE ages ABIRTH (gradd.f:205) — needed by Climate-FVS BIRTHYR; IE reads
                                          # birth_age nowhere else ⇒ inert for climate-off IE runs (bit-exact).
    _on_up = s.variant isa Ontario        # ON update.f: HT updated before the DBH bark (see below)
    _em_up = s.variant isa EasternMontana # EM ages ABIRTH (gradd.f:205, ALL trees): aspen/PB (12,17) regent reads
                                          # it (HITE1/HITE2); em/dgf.f age-range term is NOT ported ⇒ conifer-inert.
    @inbounds for i in 1:n
        # DG is the INSIDE-bark increment; outside-bark DBH grows by DG/bark, with
        # bark evaluated at the pre-growth DBH (update.f:115 / update.jl:75). CR uses the GENGYM
        # BRATIO (cr/bratio.f = cr_bratio): its bark_a/bark_b are 0, so bark_ratio would floor to 0.80
        # and over-apply DG/bark (~0.89→0.80 ⇒ ~11% too much outside-bark DBH per cycle).
        # ON: canada/on/update.f DO 90 adds HTG to HT for every record BEFORE DO 110 `DBH=DBH+DG/BRATIO(IS,DBH,HT)`,
        # so ON's metric H/D bark (bratio.f: maple group/black spruce/cedar) reads the END-of-cycle height (vols.f:131
        # likewise). HT is ignored by every other variant's BRATIO, so the pre-growth height stays for them.
        bark = variant_bratio(s, t.species[i], t.dbh[i],
                              _on_up ? t.height[i] + t.ht_growth[i] : t.height[i])   # update.f:115 — the SAME BRATIO as dgdriv's DDS→DG
        # OC stashes its own oc_bratio(D_start) in the ORGANON hook (this generic bark_ratio floors to
        # 0.80 for OC's unset bark_a/bark_b → wrong CFTOPK/BFTOPK truncation on broken-top trees).
        s.variant isa OregonCoast || (t.vol_bark[i] = bark)   # stash BRATIO(D_start) for CFTOPK/BFTOPK (vols.f:150)
        (s.variant isa Kootenai || s.variant isa InlandEmpire || s.variant isa Teton ||
         s.variant isa CentralIdaho || s.variant isa EasternMontana) &&
            (t.dg_prev[i] = t.diam_growth[i])   # KT/IE/TT/CI/EM mortality WK1 (this cycle's applied DG → next cycle's
                                                # vigor). dgdriv.f:161 WK1(I)=DG(I). CI/EM were OMITTED ⇒ their
                                                # mortality read WK1=0 forever ⇒ wrong vigor G ⇒ mis-calibrated
                                                # mortality on mature real-FIA stands (the multi-cycle blow-up).
        # OC applied DBH/HT/NORMHT inline in organon_apply_growth! (ORGANON does growth+mort+crown in one
        # EXECUTE) and KEEPS diam_growth/ht_growth for the FVS_TreeList DG/HtG report — so this shared
        # apply-loop must SKIP OC (else it double-applies). Provably .sum-inert: OC previously zeroed
        # diam_growth/ht_growth so these lines already added 0 (the norm_ht trunc(N+0.5)=N was a no-op).
        s.variant isa OregonCoast && continue   # (OC's ABIRTH is aged with the others after VOLS, gradd.f:205)
        t.dbh[i]    += t.diam_growth[i] / bark
        t.height[i] += t.ht_growth[i]
        # Broken-top trees: the full (NORMHT) height grows by the same increment as the standing
        # height. MATCH FVS update.f:67 op order EXACTLY — `INT(REAL(NORMHT)+(HTG*100.+.5))`: the
        # (HTG*100+0.5) is grouped and evaluated in Float32 FIRST, then added to NORMHT. The old
        # `(norm_ht + htg*100) + 0.5` grouping is non-associatively different in Float32 and accretes
        # a ~0.05-ft NORMHT error over cycles for topkill trees (the broken-top cuft ±1 residual).
        t.norm_ht[i] > 0 &&
            (t.norm_ht[i] = trunc(Int32, Float32(t.norm_ht[i]) + (t.ht_growth[i] * 100f0 + 0.5f0)))
    end
    compute_volumes!(s)                     # end-of-period volumes
    # gradd.f:205 `ABIRTH(I)=ABIRTH(I)+FINT` — AFTER UPDATE and its VOLS (above), so an age-dependent volume (ON's
    # Mowraski NMV) uses the pre-increment age this cycle. The one ABIRTH aging for every variant that reads it.
    _age_up && @inbounds(for i in 1:n; t.birth_age[i] += fint; end)
    accr = 0f0
    @inbounds for i in 1:n
        d = t.cuft_vol[i] - old_cfv2[i]     # OACC over the tripled set; FVS clamps
        d > 0f0 && (accr += d * t.tpa[i])   # negative growth to 0 (vols.f: CFV>tcf ⇒ WK5=0)
    end
    # GRADD order (gradd.f): UPDATE → DENSE → ESNUTR → DENSE → CROWN → VOLS. Establish
    # scheduled regen AFTER growth+mortality (fresh, full TPA this period) but BEFORE
    # CROWN, so the new trees' crown ratio (ICR) is computed this cycle (not carried
    # bogus into next cycle's DGF/mortality).
    # esnutr.f:59 CALL ESADDT(1) is the FIRST thing ESNUTR does — before the sprout logic (:121) and before every
    # OPFIND 430/431 (:63/:163/:218/:320/:348 NPNATS) and the ESTAB tally (:401). The ADDTREES bridge must therefore
    # schedule its .es2 PLANT/NATURAL here, so the AUTOES tally (ie_autoes_establish!, NPNATS rule 6 + the NBEST
    # plant pool) sees them exactly like a keyword PLANT, and the .es1 summary excludes this cycle's sprouts.
    # (establish!'s own call is then a no-op: each 432 carries a `fired` guard.)
    isempty(s.estab.addtrees) || addtrees_bridge!(s, Int32(current_cycle_year(s)), round(Int, fint))
    # IE REGENT(LESTB) height growth reads RDNEXT/BANEXT = RELDEN/BA of the gradd.f:192 DENSE (regent.f:288-295 with
    # NTYR=5 ⇒ RDNEXT(1)=RELDEN): post-growth, BEFORE ESNUTR adds sprouts/AUTOES/PLANT (post-regen DENSE is :244).
    # establish! recomputes density WITH the new cohort, so snapshot it here (same pattern as BM below).
    es_ie_relden_pre, es_ie_ba_pre = s.variant isa InlandEmpire ? (stand_ccf(s), stand_ba(s)) : (-1f0, -1f0)
    # WC: gradd.f:192 DENSE (post-UPDATE, pre-ESNUTR). ESTAB's ESSUBH reads AVH from it and ESGENT (called inside
    # ESTAB, before gradd.f:244's DENSE) reads its PCCF/PTBAA/AVH — the post-growth PRE-regen values.
    local es_wc_ptba::Vector{Float32}, es_wc_pccf::Vector{Float32}, es_wc_avh::Float32
    if s.variant isa WestCascades || s.variant isa PacificNorthwest   # PN compiles the same gradd/estab/esgent
        compute_density!(s)
        es_wc_ptba = copy(s.density.point_ba); es_wc_pccf = copy(s.density.point_ccf)
        es_wc_avh = s.plot.avg_height
    end
    # TT/UT REGENT(LESTB) likewise reads GRADD's post-growth, PRE-regen DENSE: RELDEN/BA/AVH and the per-point PCCF
    # (tt/regent.f:160-180,274-330; ut/regent.f:162-171,233). establish! re-DENSEs with the seedlings and jl's
    # density.point_ccf is still the start-of-cycle one here, so recompute the point CCF over the current records.
    es_tu_relden_pre, es_tu_ba_pre, es_tu_avh_pre, es_tu_pccf_pre = (s.variant isa Teton || s.variant isa Utah) ?
        (stand_ccf(s), stand_ba(s), stand_top_height(s), _fresh_point_ccf(s)) : (-1f0, -1f0, -1f0, Float32[])
    # EM likewise: em/esgent.f → REGENT(LESTB) runs inside ESTAB, before gradd.f:244's post-regen DENSE, so its RELDEN/
    # BA (RDNEXT/BANEXT, PPCCF) and the per-point PCCF (TPCCF for SMHTGF/SMDGF, the seedling crown dub) are the
    # gradd.f:192 DENSE's — post-growth, PRE-ESNUTR. jl's AUTOES booking re-DENSEs with the new cohort first
    # (MEASURED FVSem_g16 196378260020004 @2031: live RELDEN 109.1242 / point-1 PCCF 88.5849 vs jl post-regen
    # 109.1673 / 88.6849 ⇒ every birth-cycle HTGRR ~2e-4 low).
    es_em_relden_pre, es_em_ba_pre, es_em_pccf_pre = s.variant isa EasternMontana ?
        (stand_ccf(s), stand_ba(s), _fresh_point_ccf(s)) : (-1f0, -1f0, Float32[])
    esuckr!(s; fint = fint)                 # ESNUTR — stump/root sprouts (LSPRUT; before ESTAB)
    es_nstart = s.trees.n                    # records before ESTAB (CR grows the new regen in its birth cycle)
    es_avh_pre = s.plot.avg_height           # #194: ci/regent.f ATAVH = PRE-regen avg height (0 on bare) for the
                                             # CIVAR esgent RELHT; establish! recomputes avg_height WITH the new regen.
    # Climate-FVS AutoEstb (clauestb.f) runs BEFORE ESNUTR (gradd.f:223 < :229): schedule NATURAL regen off the
    # current (post-growth) density + species viability, for a year WITHIN this cycle (IY(ICYC+1)-1) so the
    # establish! immediately below picks it up SAME cycle. Inert unless a CLIMATE block parsed AutoEstb.
    # FVS_Climate rows are written INSIDE CLAUESTB (clauestb.f:196-216) — post-growth/mortality/aging, pre-ESTAB — so
    # snapshot them here (the summary writer emits this snapshot; building the row after grow_cycle! counted this
    # cycle's regen and the post-aging ABIRTH).
    (s.climate !== nothing && s.climate.active) &&
        (s.climate.pending_report = climate_report(s; report_year = Int(current_cycle_year(s)), fint = fint))
    (s.climate !== nothing && s.climate.active) && clim_autoestb!(s, Int(s.control.cycle) + 1, fint)
    # AUTOES (IE/EM): the AUTOMATIC natural tally (esnutr.f scheduler → estab.f DO-99, indices 1..ITPP) must run
    # BEFORE establish! appends any scheduled PLANT/NATURAL trees (estab.f appends those at ITPP+1..ITPP+ITODO, AFTER
    # the natural tally). estab.f is ONE call — the natural tally's ESB small-tree stocking (estab.f:301-322 TPACRE)
    # is therefore computed from the PRE-plant treelist; running establish! first would let the just-planted
    # seedlings inflate ESB and suppress the natural cohort (the plant-regime ~½-ingrowth bug). esnutr.f:345-359
    # ALSO fires this tally in a PLANT cycle via the rule-6 PLANT/NATURAL-in-case catch-all (ie_autoes_schedule!).
    # For a NON-plant cycle establish! finds no due activity and returns early (engine/establishment.jl:340) with no
    # side effects or RNG draws, so this reordering is INERT outside the plant/natural regime — the certified `none`
    # floor and the ie_autoes RNG stream are unchanged.
    (s.variant isa InlandEmpire || s.variant isa EasternMontana) && ie_autoes_establish!(s; fint = fint)
    # ESTAB site-prep status for the non-AUTOES variants' PLANT/NATURAL catch-all ESTAB call (esnutr.f) — bookkeeping
    # only (ECON MECHCST/BURNCST); IE/EM record it inside ie_autoes_establish!.
    (s.variant isa InlandEmpire || s.variant isa EasternMontana || s.variant isa SoutheastAlaska) || estab_prep_esnutr!(s)
    # BM REGENT(LESTB) reads RELDEN/AVH from the GRADD DENSE that precedes ESNUTR (gradd.f UPDATE→DENSE→ESNUTR):
    # post-growth, PRE-regen. establish! recomputes density WITH the new seedlings, so snapshot it here.
    es_bm_relden_pre, es_bm_avh_pre = (s.variant isa BlueMountains || s.variant isa EastCascades) ?
        (stand_ccf(s), stand_top_height(s)) : (0f0, 0f0)
    establish!(s; fint = fint)              # ESNUTR — adds scheduled PLANT/NATURAL regen (ICR=0), recomputes density
    # AK: estb/esnutr.f + ak/estab.f (natural tally + PLANT/NATURAL) + ak/esgent.f, one ESNUTR call per cycle.
    s.variant isa SoutheastAlaska && ak_esnutr!(s; fint = fint)
    # WPBR BRESTB (estab.f, IE/EM): seed this cycle's new host records at their birth state before ESGENT grows them.
    (s.variant isa InlandEmpire || s.variant isa EasternMontana) && s.wpbr !== nothing && wpbr_brestb_new!(s, es_nstart)
    # CR-only: esgent.f grows the just-established regen IN their creation cycle via REGENT (eastern leaves them
    # ungrown per GRADD order — bit-exact). Fixes the ESTAB 1-cycle-offset (TopHt lag) on cr_estab.
    s.variant isa CentralRockies && cr_esgent!(s, es_nstart; fint = fint)
    s.variant isa Teton && tt_esgent!(s, es_nstart; fint = fint,
        atavh = es_at_avh, atba = es_at_ba, atrelden = es_at_relden,
        relden_pre = es_tu_relden_pre, ba_pre = es_tu_ba_pre, pccf_pre = es_tu_pccf_pre)   # TT: tt/esgent.f → REGENT(LESTB)
    (s.variant isa WestCascades || s.variant isa PacificNorthwest) && wc_esgent!(s, es_nstart; fint = fint, atavh = es_at_avh,
        avh_pre = es_wc_avh, ptba_pre = es_wc_ptba, pccf_pre = es_wc_pccf)   # WC: wc/esgent.f → REGENT(LESTB)
    s.variant isa EastCascades && ec_esgent!(s, es_nstart; fint = fint,
        atavh = es_at_avh, atrelden = es_at_relden,
        relden_pre = es_bm_relden_pre, avh_pre = es_bm_avh_pre)   # EC western: grow birth-cycle regen (ec/esgent.f)
    s.variant isa EasternMontana && em_esgent!(s, es_nstart; fint = fint,
        atba = es_at_ba, atccf = es_at_relden, atavh = es_at_avh,
        relden_pre = es_em_relden_pre, ba_pre = es_em_ba_pre, pccf_pre = es_em_pccf_pre)   # EM: em/esgent.f -> REGENT(LESTB) (#137)
    s.variant isa Utah && ut_esgent!(s, es_nstart; fint = fint,
        atavh = es_at_avh, atrelden = es_at_relden,
        relden_pre = es_tu_relden_pre, avh_pre = es_tu_avh_pre, pccf_pre = es_tu_pccf_pre)   # UT western: grow birth-cycle regen (ut/esgent.f, #184); #194-class start-of-cycle ATAVH/ATCCF blend for PCTRED
    s.variant isa CentralIdaho && ci_esgent!(s, es_nstart; fint = fint, avh_pre = es_avh_pre)   # CI western: grow birth-cycle regen (ci/esgent.f, #185); #194 pass pre-regen ATAVH
    s.variant isa BlueMountains && bm_esgent!(s, es_nstart; fint = fint,
        atavh = es_at_avh, atrelden = es_at_relden,
        relden_pre = es_bm_relden_pre, avh_pre = es_bm_avh_pre)   # BM western: grow birth-cycle regen (bm/esgent.f, #185); #194-class start-of-cycle ATAVH/ATCCF blend for PCTRED
    s.variant isa InlandEmpire && ie_esgent!(s, es_nstart; fint = fint,
        atavh = es_at_avh, atba = es_at_ba, atrelden = es_at_relden,
        relden_pre = es_ie_relden_pre, ba_pre = es_ie_ba_pre)   # IE western: grow birth-cycle regen (ie/esgent.f, #186; NIVAR). #194-class: start-of-cycle TEMAHT/TEMBA/TEMCCF for DADJ
    # dgdriv.f:142 WK1(I)=DG(I) runs at the START of the next cycle's DGDRIV over ALL records, so a tree born this
    # cycle enters next cycle's MORTS with WK1 = its birth-cycle DG (regent.f:941 LESTB NIVAR: DG(K)=DK). jl copies
    # dg_prev in the growth-apply loop above, which runs BEFORE establishment — so the new records kept WK1=0 and
    # took the morts.f `WK1.EQ.0 ⇒ G=DG/(BARK·10)` vigor branch (e.g. bare PLANT stand: a seedling past 4.5 ft
    # got RIPP 0.00024 vs live 0.0104 ⇒ ~40× under-kill). Copy for the new IE records now.
    # TT morts reads WK1 too (teton/mortality.jl) ⇒ the same copy for its birth-cycle DG (tt_esgent!).
    if s.variant isa InlandEmpire || s.variant isa Teton
        @inbounds for i in (es_nstart + 1):s.trees.n
            s.trees.dg_prev[i] = s.trees.diam_growth[i]
        end
    end
    # estab.f:1490-1493 — "IF NEW TREES HAVE BEEN ADDED TO THE TREELIST" the establishment
    # model calls ESGENT, which calls SPESRT (esgent.f:49), rebuilding IND1 in ASCENDING
    # physical order and DISCARDING the post-TRIPLE REASS (U,C,L) lineage interleave
    # (grincr.f:553). So a cycle that ADDS regen enters the next cycle's DGSCOR in ascending
    # order, not lineage order; a cycle with NO new trees leaves the REASS lineage intact.
    # jl only reset `sort_key` inside the TREDEL removal path (tredel_compact!); it never
    # modeled this establishment SPESRT, so the stale REASS lineage key mis-ordered the
    # per-tree DGSCOR/REGENT RNG draws from cyc3 on (IE none over-grows BA +3/+4 vs oracle).
    # Gate on "records were added this cycle" (t.n > es_nstart) — mirrors estab.f:1493 exactly;
    # NOT the coarse AUTOES-active flag (which over-fires on no-establishment cycles and moves
    # otherwise-bit-exact stands off the oracle). Applied to the AUTOES western variants IE + EM
    # (both measured bit-exact-or-improving vs FVS{ie,em}_g16); SN/eastern lineage sort unchanged.
    # See uses_estab_spesrt above for the per-variant measurement (UT/BM/TT/CI measured-inert, excluded).
    uses_estab_spesrt(s.variant) && s.trees.n > es_nstart && spesrt_reorder!(s.trees)
    # ECCALC (gradd.f:239): after ESNUTR, before DENSE — this cycle's ECON cash flows + FVS_EconSummary row.
    (s.econ !== nothing && s.econ.active) && econ_calc!(s, Int(s.control.cycle) + 1)
    compute_density!(s)                     # gradd.f DENSE-before-CROWN: refresh the POST-growth stand BA the
                                            # NE/CS crown model reads (was stale pre-growth ⇒ CS crown/DG drift).
                                            # SN's crown uses the pre-growth crown_sdi captured above, so unaffected.
    crown_ratio_update_fvs!(s; fint = fint, crown_sdi = crown_sdi)  # CROWN — pre-growth Reineke RELSDI
    # gradd.f:267 — snapshot PCT into OLDPCT AFTER crown, so next cycle's crown DCR reads this cycle's PCT.
    # (IE crown uses OLDPCT in the backdated DCR term; other variants approximate OLDPCT≈PCT so this is inert.)
    if s.variant isa InlandEmpire || s.variant isa BritishColumbia || s.variant isa EasternMontana ||
       s.variant isa Kootenai
        @inbounds for i in 1:s.trees.n; s.trees.old_crown_pct[i] = s.trees.crown_ratio[i]; end
    end
    # IE: snapshot the crown-time (current) stand BA/RELDEN into OLDBA/RELDM1 so NEXT cycle's CROWN backdates
    # DCRCON against them (dense.f threads RELDM1/OLDBA as the prior-cycle density; oracle OBA[N]=BA[N-1],
    # RDM1[N]=RELDEN[N-1]). Previously never assigned ⇒ OBA==BA, RDM1==RELDEN ⇒ DCRCON==XCRCON ⇒ EDCR too low
    # ⇒ CHG (=EXPPCR−EXPDCR) too high ⇒ ICR +1..3 too high every cycle. Verified vs FVSie_g16 on 3307603010690:
    # per-tree ICR at CROWN goes from 33/39 one-directional +diffs to ~5 mixed ±1 (residual = a small stand-BA gap).
    if s.variant isa InlandEmpire || s.variant isa EasternMontana || s.variant isa Kootenai   # EM/KT: same crown.f OBA/RDM1
        s.plot.old_ba = s.plot.basal_area
        s.plot.relative_density_prev = s.plot.relative_density
    end
    # NOTE: newly-established trees get NO volume in their birth cycle. The oracle's
    # VOLS in the establishment cycle runs before the records are inserted, so a planted
    # stand reports CFV=0 at cyc1 (verified: bare_plant 1997 cuft=0) and the regen first
    # gets volume from the next cycle's VOLS (next grow_cycle!'s compute_volumes!). The
    # crown pass above DOES set the new trees' ICR this cycle (DGF/mortality read it next).
    # WPBR BRPR (fvs.f:408, after TREGRO/DISPLY): BRTSTA tree statuses + BRSTAT stand statistics that the
    # next cycle's BRCREM/BRECAN read. Inert (no-op) unless a BRUST block is active with host pines.
    s.wpbr !== nothing && wpbr_brpr!(s)
    s.control.total_removal = 0f0            # fvs.f:432 ONTREM(7)=0 for the next cycle
    # RDSUM: FVS_RD_Sum row (rdpr.f at fvs.f:404, after TREGRO — so after GRADD's ESNUTR). estab.f:1247/1336/1426 call
    # RDESTB for every record ESTAB books, entering it into the disease area (PROBIU=PROB·PAREA, FPROB=PROB): size the
    # driver over this cycle's regen HERE, then report — post DBH-UPDATE (grown DBH for Live_BA), post the end-of-period
    # VOLS, post rd_grow_apply!→rdinoc decay. jl reported before establishment and entered the regen only at the next
    # cycle start (MEASURED FVSem_g16 196378260020004 rootdis 2032: UnInf_TPA live 341.25, jl 284.12 — the 2031 cohort).
    if !tripled && s.root_disease !== nothing && rd_active(s.root_disease) && s.root_disease.iroot != 0 &&
       s.root_disease.driver !== nothing
        let rd = s.root_disease, d = rd.driver::RDDriver
            if d.n != s.trees.n
                rd.wk1_nold = d.n                   # rd_cycle_start!: WK1 of these records is still 0 (estab.f WK1=0)
                rd.driver = _rd_resize_driver!(rd, d, s.trees.n, s)
            end
        end
        if s.control.dbs_rd_sum || s.control.dbs_rd_detail
            yr = cycle_year_at(s.control, Int(s.control.cycle) + 1)
            iage = Int(s.plot.stand_age) + (yr - Int(s.control.cycle_year[1]))
            s.control.dbs_rd_sum    && push!(s.root_disease.sum_rows, (yr, rd_sum_report(s.root_disease, s, yr, iage)))
            s.control.dbs_rd_detail && push!(s.root_disease.det_rows, (yr, rd_det_report(s.root_disease, s, yr)))
        end
    end
    s.control.cycle += Int32(1)
    return (; accretion = accr / fint / g, mortality = mort / fint / g)
end

"""
    run_keyfile(keypath; variant=Southern(), faithful=true, period=5) -> String

Full multi-stand run: project EVERY stand in `keypath` (the FVS `main.f` stand loop)
and return the concatenated `.sum` text — one `-999` header + per-cycle rows per
stand. Each stand is set up (`notre!`/`setup_growth!`/`compute_volumes!`) and projected
by `write_sum_file`, which runs the per-cycle loop including scheduled management
(CUTS/ESTAB/fire). Stands are independent (`each_stand` gives each a fresh state with
the tree format carried across), so this is also the unit of thread-parallelism.
"""
function run_keyfile(keypath::AbstractString;
                     variant::Union{AbstractVariant,Nothing} = nothing,
                     output::Union{Symbol,AbstractString,Nothing} = nothing,
                     faithful::Bool = true, period::Integer = 5,
                     date::AbstractString = "", time::AbstractString = "")
    # Output format: :sum (legacy fixed-column, the default) or :csv (named columns). Resolution
    # mirrors `variant` — an explicit `output=` wins; else the YAML's `output_format:`; else :sum.
    outfmt = _resolve_output(keypath, output)
    # FVS stamps the run date/time into the .sum -999 header (MM-DD-YYYY / HH:MM:SS); mirror
    # that when the caller didn't supply them, so the header carries a timestamp like FVSsn's
    # (its value is wall-clock and so won't byte-match a separate run — the data rows do).
    # `Base.Libc.strftime`/`Base.time` avoid a Dates stdlib dependency. (`time` is the kwarg.)
    isempty(date) && (date = Base.Libc.strftime("%m-%d-%Y", Base.time()))
    isempty(time) && (time = Base.Libc.strftime("%H:%M:%S", Base.time()))
    out = IOBuffer()
    csv_stands = outfmt === :csv ? Tuple[] : nothing   # (stand_id, mgmt_id, SummaryRows) per stand
    case = 0
    kt_ierrck = Int32(0)                          # kt/cratet.f IERRCK: a -fno-automatic static carried stand to stand
    ak_r10ra = false                              # r10tap.f saved ISP='RA' (AK DEM taper branch), carried the same way
    for s in each_stand(keypath; variant = variant, faithful = faithful)
        s.control.kt_cratet_ierrck = kt_ierrck
        s.control.ak_r10tap_ra = ak_r10ra
        notre!(s)
        setup_growth!(s)
        kt_ierrck = s.control.kt_cratet_ierrck
        compute_volumes!(s)
        # SVSTART seam (fvs.f:333, gated JSVOUT≠0): emit the cycle-0 inventory SVS picture at the
        # inventory state (post-setup, pre-growth). Only stands with an SVS keyword (svs_on) write files.
        if s.control.svs_on
            stem = isempty(s.control.svs_keystem) ?
                   (isempty(keypath) ? "svs" :
                    joinpath(dirname(keypath), first(splitext(basename(keypath))))) :
                   s.control.svs_keystem
            svs_write_cycle0_files(stem, s)
        end
        sid = strip(s.plot.stand_id)
        mid = strip(s.plot.mgmt_id); mid = isempty(mid) ? "NONE" : String(mid)
        # DBS output (DATABASE block): collect this stand's summary rows and/or per-cycle tree
        # snapshots and append them to the DSNOUT database, in addition to the text `.sum`.
        has_db = !isempty(s.control.dbs_out_file)
        sum_on = s.control.dbs_summary && has_db
        tl_on = s.control.dbs_treelist && has_db
        cp_on = s.control.dbs_compute && has_db && !isempty(s.control.compute_defs)
        cl_on = s.control.dbs_cutlist && has_db
        al_on = s.control.dbs_atrtlist && has_db
        rows = (sum_on || outfmt === :csv) ? SummaryRow[] : nothing   # also collected for CSV output
        tl_cycles = tl_on ? Tuple[] : nothing
        cp_rows = cp_on ? Tuple[] : nothing
        cl_cycles = cl_on ? Tuple[] : nothing
        al_cycles = al_on ? Tuple[] : nothing
        # FFE Stand Carbon Report (CARBREPT) / Potential Fire (POTFIRE): collect per cycle, same simulation.
        _fuels_db = s.control.dbs_fuels && s.control.ffe_fuelout      # FVS_Fuels gate (FUELSOUT + FUELOUT window)
        # fmmain.f runs FMSOUT/FMSSUM/FMPOFL/FMDOUT/FMCRBOUT/FMCHRVOUT every FFE year; each DBS table is written only
        # when its DATABASE toggle is set (and, for the snag/down-wood tables, its FMIN report keyword).
        ctl0 = s.control
        _ffe_tbls = has_db && (ctl0.dbs_carbrept || (ctl0.dbs_snagsum && ctl0.ffe_snagsum) ||
                               (ctl0.dbs_snagdet && ctl0.ffe_snagout) || (ctl0.dbs_dwdvol && ctl0.ffe_dwdvlout) ||
                               (ctl0.dbs_dwdcov && ctl0.ffe_dwdcvout))
        carb_rows = ((ctl0.carbon_report_on || _fuels_db || _ffe_tbls) && s.fire !== nothing && s.fire.active) ? Tuple[] : nothing
        pf_rows = (s.control.potfire_report_on && s.fire !== nothing && s.fire.active) ? Tuple[] : nothing
        hc_rows = (has_db && ctl0.dbs_carbrept && s.fire !== nothing && s.fire.active) ? Tuple[] : nothing
        clim_rows = (s.control.dbs_climate && s.climate !== nothing && s.climate.active) ? Tuple[] : nothing
        cprof_rows = (s.control.dbs_canprofile && s.fire !== nothing && s.fire.active) ? Tuple[] : nothing
        strcl_rows = s.control.dbs_strclass ? Tuple[] : nothing
        dm_rows = (s.control.dbs_mistoe && _dm_report_variant(s.variant)) ? Tuple[] : nothing
        dm_top4 = Int[]
        # PRTRLS(1) (fvs.f:328 pre-projection with the cycle-1 options, fvs.f:412 at each cycle end): one
        # FVS_TreeList block per TREELIST request accomplished this cycle (none without a TREELIST activity).
        # PrdLen = IFINT (dbstrls.f): per cycle grincr.f:65 sets it to the cycle just grown; the inventory list sees
        # the DB DG_MEASURE IFIX(FINT) (dbsstandin.f:702) or else grinit's IFINT (10; 5 in SN/OC/OP) — the GROWTH
        # keyword changes FINT only (initre.f:829 has the IFINT line commented out).
        hook = tl_on ? (st, yr, pl, cy) -> begin
            pl_eff = cy == 0 ? (st.control.dbs_ifint >= 0 ? Int(st.control.dbs_ifint) :
                                (st.variant isa Southern || st.variant isa OregonCoast || st.variant isa Olympic) ? 5 : 10) : pl
            for _ in prtrls_requests!(st, 1, cy == 0 ? 1 : cy; lstart = cy == 0)
                push!(tl_cycles, treelist_snapshot(st, yr, pl_eff; cycle = cy))
            end
        end : nothing
        write_sum_file(out, s; period = Int(period), stand_id = String(sid),
                       mgmt_id = mid, variant = variant_code(s.variant), date = date, time = time,
                       collect_rows = rows, cycle_hook = hook, compute_collect = cp_rows,
                       cutlist_collect = cl_cycles, atrtlist_collect = al_cycles, carbon_collect = carb_rows, potfire_collect = pf_rows,
                       hrvcarbon_collect = hc_rows, climate_collect = clim_rows,
                       dm_collect = dm_rows, dm_top4 = dm_top4,
                       canprof_collect = cprof_rows, strclass_collect = strcl_rows)
        carb_rows === nothing ||
            write_carbon_report_block(out, carb_rows; stand_id = String(sid), mgmt_id = mid)
        # COVER report (CVOUT): "CANOPY COVER STATISTICS" table, appended after the .sum
        # rows (report-only; never touches .sum/tree). Gated on the COVER activity 900.
        (s.cover !== nothing && s.cover.active) &&
            cover_report(s.cover, out, String(sid), mid, strip(s.control.title))
        csv_stands === nothing || push!(csv_stands, (String(sid), mid, strip(s.control.title), rows))
        if has_db
            case += 1
            caseid = string(sid, "-", case)
            # FVS_Cases registry + FVS_InvReference reference dump accompany any DBS output.
            kwfile = isempty(keypath) ? "" : first(splitext(basename(keypath)))
            write_dbs_cases!(s.control.dbs_out_file, caseid, String(sid);
                             mgmt_id = mid, variant = variant_code(s.variant),
                             keyword_file = kwfile, sampling_wt = s.plot.sample_weight,
                             run_datetime = strip(string(date, " ", time)))
            write_dbs_invref!(s.control.dbs_out_file, caseid, String(sid), s)
            write_dbs_error!(s.control.dbs_out_file, caseid, String(sid), s.control.error_msgs)   # DBSERROR rows
            # BC/ON link metric/dbsqlite: DBSSUMRY/DBSTRLS write the *_Metric tables (East naming for ON) instead.
            met = _metric_variant(s.variant); east = s.variant isa Ontario
            if sum_on
                met ? write_dbs_summary_metric!(s.control.dbs_out_file, caseid, String(sid), rows; east = east) :
                      write_dbs_summary!(s.control.dbs_out_file, caseid, String(sid), rows;
                                         mgmt_id = mid, variant = variant_code(s.variant))
            end
            # DBSTRLS/DBSCUTS create their table only when actually called for an accomplished list request
            if tl_on && !isempty(tl_cycles)
                met ? write_dbs_treelist_metric!(s.control.dbs_out_file, caseid, String(sid), tl_cycles; east = east) :
                      write_dbs_treelist!(s.control.dbs_out_file, caseid, String(sid), tl_cycles)
            end
            (cl_on && any(c -> !isempty(c[3]), cl_cycles)) &&
                write_dbs_cutlist!(s.control.dbs_out_file, caseid, String(sid), cl_cycles; metric = met, east = east)
            (al_on && any(c -> !isempty(c[3]), al_cycles)) &&
                write_dbs_atrtlist!(s.control.dbs_out_file, caseid, String(sid), al_cycles; metric = met, east = east)
            clim_rows === nothing ||
                write_dbs_climate!(s.control.dbs_out_file, caseid, String(sid), clim_rows, s.coef)
            cprof_rows === nothing ||
                write_dbs_canprofile!(s.control.dbs_out_file, caseid, String(sid), cprof_rows)
            strcl_rows === nothing ||
                write_dbs_strclass!(s.control.dbs_out_file, caseid, String(sid), strcl_rows, s.coef)
            if dm_rows !== nothing && !isempty(dm_rows)
                # FVS_DM_Stnd_Sum + FVS_DM_Spp_Sum (DBSMIS2/DBSMIS1); the by-DBH-class FVS_DM_Sz_Sum
                # (DBSMIS3) additionally needs the MISTPRT report keyword (PRTMIS).
                write_dbs_dm_stndsum!(s.control.dbs_out_file, caseid, String(sid), dm_rows)
                write_dbs_dm_sppsum!(s.control.dbs_out_file, caseid, String(sid), dm_rows, s.coef)
                s.control.mistprt_on &&
                    write_dbs_dm_szsum!(s.control.dbs_out_file, caseid, String(sid), dm_rows)
            end
            # FVS_BM_Main/Tree/Vol (WWPB MAINOUT/TREEOUT/VOLOUT): rows accumulated per outbreak cycle.
            if s.wwpb !== nothing && !isempty((s.wwpb::WwpbState).main_rows)
                bm_rows = (s.wwpb::WwpbState).main_rows
                s.control.dbs_bm_main &&
                    write_dbs_bm_main!(s.control.dbs_out_file, caseid, String(sid), bm_rows)
                s.control.dbs_bm_tree &&
                    write_dbs_bm_tree!(s.control.dbs_out_file, caseid, String(sid), bm_rows)
                s.control.dbs_bm_vol &&
                    write_dbs_bm_vol!(s.control.dbs_out_file, caseid, String(sid), bm_rows)
                s.control.dbs_bm_bkp &&
                    write_dbs_bm_bkp!(s.control.dbs_out_file, caseid, String(sid), bm_rows)
            end
            # FVS_RD_Sum (WRD root-disease summary): rows accumulated per cycle after rd_end_apply!.
            if s.root_disease !== nothing && s.control.dbs_rd_sum && !isempty(s.root_disease.sum_rows)
                write_dbs_rd_sum!(s.control.dbs_out_file, caseid, String(sid), s.root_disease.sum_rows)
            end
            # FVS_RD_Det (WRD per-species patch detail): rows accumulated per cycle (dbs/dbsrd.f DBSRD2).
            if s.root_disease !== nothing && s.control.dbs_rd_detail && !isempty(s.root_disease.det_rows)
                write_dbs_rd_det!(s.control.dbs_out_file, caseid, String(sid), s.root_disease.det_rows, s.coef)
            end
            s.control.dbs_calibstats &&
                write_dbs_calibstats!(s.control.dbs_out_file, caseid, String(sid), s.calib, s.coef)
            _fuels_db && carb_rows !== nothing &&
                write_dbs_fuels!(s.control.dbs_out_file, caseid, String(sid), carb_rows)
            if carb_rows !== nothing
                ctl1 = s.control
                ctl1.dbs_carbrept &&                                     # dbsfmcrpt.f ICMRPT (CARBREDB)
                    write_dbs_carbon!(ctl1.dbs_out_file, caseid, String(sid), carb_rows)
                (ctl1.dbs_snagsum && ctl1.ffe_snagsum) &&                # fmssum.f ISNGSM≠−1 + dbsfmssnag.f ISSUM
                    write_dbs_snagsum!(ctl1.dbs_out_file, caseid, String(sid), carb_rows)
                (ctl1.dbs_snagdet && ctl1.ffe_snagout) &&                # fmsout.f window + dbsfmdsnag.f ISDET
                    write_dbs_snagdet!(ctl1.dbs_out_file, caseid, String(sid), [(r[1], r[7]) for r in carb_rows], s.coef)
                (ctl1.dbs_dwdvol && ctl1.ffe_dwdvlout) &&                # fmdout.f LPRINT2 + dbsfmdwvol.f IDWDVOL
                    write_dbs_dwd_vol!(ctl1.dbs_out_file, caseid, String(sid), carb_rows)
                (ctl1.dbs_dwdcov && ctl1.ffe_dwdcvout) &&                # fmdout.f LPRINT3 + dbsfmdwcov.f IDWDCOV
                    write_dbs_dwd_cov!(ctl1.dbs_out_file, caseid, String(sid), carb_rows)
            end
            # Fire-EVENT DBS tables: one row per SIMFIRE event (captured by fmburn!), independent of CARBREPT
            if s.fire !== nothing && s.fire.active && !isempty(s.fire.burn_reports)
                br = s.fire.burn_reports
                # fmfout.f: each table needs its FMIN report window AND its DATABASE toggle (dbsfmburn/-mort/-fuel)
                ctl = s.control
                (ctl.ffe_burnrept && ctl.dbs_burnrept) &&
                    write_dbs_burnreport!(ctl.dbs_out_file, caseid, String(sid), br)
                (ctl.ffe_mortrept && ctl.dbs_mortrept) &&
                    write_dbs_mortality!(ctl.dbs_out_file, caseid, String(sid), br)
                (ctl.ffe_fuelrept && ctl.dbs_fuelcons) &&
                    write_dbs_consumption!(ctl.dbs_out_file, caseid, String(sid), br)
            end
            pf_rows === nothing ||
                write_dbs_potfire!(s.control.dbs_out_file, caseid, String(sid), pf_rows)
            # fmchrvout.f: ICHRVB defaults to 9999 (fminit.f:899), so the `ICHRVB .EQ. 0` exit never fires and DBSFMHRPT
            # writes a row every FFE year once CARBREDB set ICHRPT — zero rows included (jl required a removal).
            hc_rows === nothing || isempty(hc_rows) ||
                write_dbs_hrvcarbon!(s.control.dbs_out_file, caseid, String(sid), hc_rows)
            if cp_on
                var_names = String[nm for (_, nm, _) in s.control.compute_defs]
                write_dbs_compute!(s.control.dbs_out_file, caseid, String(sid), var_names, cp_rows)
            end
            # FVS_EconHarvestValue: log-graded HRVRVN harvest-value detail (echarv.f). Emit when an
            # active ECON run has accumulated log-graded (unit-4) revenue over the projection's cuts.
            # FVS_EconSummary (DBSECSUM, called from ECCALC each period): only when ECONRPTS set IDBSECON>0.
            if s.econ !== nothing && s.econ.calc !== nothing && s.econ.calc.dbs_econ > 0 && !isempty(s.econ.calc.rows)
                write_dbs_econsummary!(s.control.dbs_out_file, caseid, String(sid), s.econ.calc.rows)
            end
            # FVS_EconHarvestValue (DBSECHARV, eccalc.f:745-855): faithful per-cycle rows from revVolume (all
            # revenue units); written only when ECONRPTS set IDBSECON=2 (dbsecharv.f:16,86).
            if s.econ !== nothing && s.econ.calc !== nothing && s.econ.calc.dbs_econ == 2 &&
               !isempty(s.econ.calc.hv_rows)
                write_dbs_econharvest_rows!(s.control.dbs_out_file, caseid, s.econ.calc.hv_rows, s.coef)
            end
        end
        ak_r10ra = s.control.ak_r10tap_ra
    end
    if outfmt === :csv
        cio = IOBuffer(); write_sum_csv(cio, csv_stands); return String(take!(cio))
    end
    return String(take!(out))
end


