# FVSjl Lake States (LS) variant port — SCOPE & LEDGER

Sole oracle: live **FVSls** relinked from `bin/FVSls_buildDir/*.o` via `test/harness/ls_oracle.sh`.
Canonical stand: `tests/FVSls/lst01.key` (+ `lst01.tre`). Goal + doctrine: `docs/LS_GOAL.md`.

## Regression baseline (before any LS code)
`julia --project=. test/runtests.jl` ⇒ **6462 pass / 2 broken** (the 2 broken = accepted SN COMPRESS
eigensolver + NOHTDREG ULP). Adding LS must keep SN/NE/CS bit-exact; the 2 broken stay 2.
**After CHUNK 1: 6470 pass / 2 broken** (+8 = test_lst01.jl cycle-0 assertions; no SN/NE/CS regression).

## Infra constants (confirmed from FVSls_buildDir)
`MAXSP=68` (PRGPRM.F77) · `YR=10.0` (blkdat.f:71) · RNG `S0/SS=55329` (blkdat.f:302) ·
`LZEIDE=.TRUE.` (grinit.f:126, Zeide SDI) · VARACD `"LS"` (grinit.f:66) · mortality `RI=0.5·RI`
(morts.f:504) + SDI DBH gate `DBHSDI` (morts.f:203). **DONE:** `LakeStates<:AbstractVariant`
singleton (src/variants/lakestates/lakestates.jl), `variant_from_code("LS")`, registered in FVSjl.jl.

## Per-routine reuse table (LS lines / NE lines — from the buildDir scan)
| Routine | class | plan |
|---|---|---|
| htcalc.f (425/425) | **IDENTICAL to NE** | reuse NE height-from-dbh directly |
| dgdriv.f (757/737) | NEAR-IDENTICAL 97% | reuse NE calibration; verify LS GST DBH floor |
| cratet.f (628/622) | NEAR-IDENTICAL 99% | reuse |
| regent.f (629/629) | NEAR-IDENTICAL 96% | reuse NE small-tree, LS coeffs |
| varmrt.f (212/216) | NEAR-IDENTICAL 96% | reuse NE VARMRT |
| cubrds.f (58/58) | SIMILAR 85% | reuse, cross-check coeffs |
| crown.f (289/286) | SAME-STRUCTURE 91% | reuse NE TWIGS crown, LS BCR1..4 coeffs |
| htgf.f (179/177) | SAME-STRUCTURE, **IVAR=1 vs NE IVAR=3** | NE shape, trace HTCALC IVAR branch |
| **dgf.f (508/199)** | **LS-SPECIFIC ln(DDS)** | extend Southern dgf! framework (see below) |
| **balmod.f (120/52)** | **LS-SPECIFIC** | new; iface `balmod(ISPC,D,BA,RMSQD,GM,DEBUG)` |
| forkod.f (332/216) | LS-SPECIFIC | new forest-code map |
| sitset.f (470/646) | DIFFERENT | new LS site setup |
| htdbh.f (352/484) | DIFFERENT | new LS HT-DBH (simplified) |

## dgf.f — the LS diameter-growth model (ls/dgf.f:414-427)
`DDS = CONSPP + INTERC + VDBHC·(1/D) + DBHC·D + DBH2C·D² + RDBHC·(D/QMDGE5) + RDBHSQC·(D²/QMDGE5)
      + CRWNC·CR + CRSQC·CR² + SBAC·BAGE5 + BALC·BAL + SITEC·SITEAR`
Per-species coeff arrays: INTERC, VDBHC, DBHC, DBH2C, RDBHC, RDBHSQC, CRWNC, CRSQC, SBAC, BALC,
SITEC, OBSERV(sample size). Predictors: DBH(1/D,D,D²), CR, CR², BAL, BAGE5 (stand BA of trees ≥5"),
QMDGE5 (QMD of trees ≥5"), site. SN-family ln(DDS), NOT NE's BAL-potential iteration ⇒ extend the
Southern `dgf!` framework with LS coefficients + the QMDGE5/BAGE5 density terms.

## Volume — LS uses R9 **Gevorkiantz** (r9vol.f), NOT R9 Clark
This DIFFERS from NE/CS (which ride R9 Clark). LS: vols.f → r9vol.f (Gevorkiantz cubic/bdft) +
r9clark.f present as a secondary option. Cubic/board via cfvol/bfvol; varvol.f NATCRS/OCFVOL/OBFVOL.
⇒ the volume chunk needs the Gevorkiantz path ported (new), not just coefficient swap. CONFIRM which
path lst01 actually exercises (defaults) before porting.

## Species (68) — roster for data/lakestates/species_roster.csv
idx:alpha:FIA — 1:JP:105 2:SC:130 3:RN:125 4:RP:125 5:WP:129 6:WS:094 7:NS:091 8:BF:012 9:BS:095
10:TA:071 11:WC:241 12:EH:261 13:OS:299 14:RC:068 15:BA:543 16:GA:544 17:EC:742 18:SV:317 19:RM:316
20:BC:762 21:AE:972 22:RL:975 23:RE:977 24:YB:371 25:BW:951 26:SM:318 27:BM:314 28:AB:531 29:WA:541
30:WO:802 31:SW:804 32:BR:823 33:CK:826 34:RO:833 35:BO:837 36:NP:809 37:BH:402 38:PH:403 39:SH:407
40:BT:743 41:QA:746 42:BP:741 43:PB:375 44:(blank/CommHW→WN) 45:BN:601 46:WN:602 47:HH:701 48:BK:901
49:OH:998 50:BE:313 51:ST:315 52:MM:319 53:AH:391 54:AC:421 55:HK:462 56:DW:491 57:HT:500 58:AP:660
59:BG:693 60:SY:731 61:PR:761 62:CC:763 63:PL:760 64:WI:920 65:BL:922 66:DM:923 67:SS:931 68:MA:935
(idx 44 = "commercial hardwood" placeholder, blank alpha in Fortran — confirm handling from ls/blkdat.f.)

## Chunk plan (upstream-first)
1. **Infra + species data + cycle-0 stand columns** (TPA/BA/SDI/CCF/QMD/TopHt) bit-exact vs live
   lst01 (1990: 536/77/160/171/63/5.1) + test/integration/test_lst01.jl. Needs: species_roster,
   species_translation (LS FIA→idx), CCF/stocking coeffs, per-species SDImax, LS htdbh (dub heights),
   sitset. ← CURRENT
2. **Volume** (R9 Gevorkiantz) ⇒ cycle-0 .sum volume cols (1551/1338/480/1887).
3. **Diameter growth** (ls/dgf.f SN-family + LS balmod + dgdriv calibration) ⇒ cycle-1 DG.
4. **Height growth** (NE htgf/htcalc, IVAR=1) · **crown** (NE TWIGS, LS coeffs) · **mortality**
   (varmrt) · htdbh/sitset/forkod ⇒ lst01 cycle-1+ vs live.
5. LS-active shared branches (FFE/estab/sprout/thinning) — validate each at its key cycle.

## Ledger (verdicts as chunks land)
- **CHUNK 5 (breadth) — SCOPED; each management path needs an LS-specific model/data sub-port.** Ran the
  full lst01.key (5 stands) — the growth-and-yield CORE is done, but the management stands surface LS data/
  model gaps (all well-scoped):
  - **Sprouting — MODEL PORTED + WIRED (thinning stand runs); sprout PRODUCTION under-produces (calibration
    follow-up).** Ported the full LS sprout model into src/engine/sprout.jl (all live-validated forms):
    `essprt_ls` (PREM survival, essprt.f:256-345 — ×0.30-0.95 + DSTMP formulas sp26/27/30/31/34/35),
    `nsprec_ls` (# records, ENTRY NSPREC CASE('LS')), `sprtht_ls` (height `(0.1+SI/50)·age` for sp15-43/45-48/
    50-68 else `0.5+0.5·age`); esuckr! dispatch gains an `ls` branch (aspen ESASID=41, ne_sprout_dbh via LS
    htdbh_ht1/ht2, sprout_essprt.csv is_sprouting for the 50 hardwoods). LSPRUT=.TRUE. (esinit.f:50). Suite
    6478/2, no SN/NE/CS regression. **★ THINNING+SPROUTING NOW TRACKS LIVE CLOSELY** after fixing the
    `is_sprouting` set: it must be the ESTUMP **ISPSPE** list (blkdat.f:103 = {15-43, 45-48, 50-68}), NOT
    the essprt PREM-CASE species — the essprt case EXCLUDES aspen (41) and sp51, but ESTUMP still LOGS them
    (aspen sprouts via ASSPTN root-suckering, the DOMINANT sprout source: aspen sucker PREM ~9.85×2/record).
    My first `is_sprouting` (from the essprt case) set aspen=0 ⇒ jl cut NO aspen ⇒ astpar=0 ⇒ no suckers ⇒
    361 TPA. Rebuilt from ISPSPE ⇒ **thinning stand 2020 jl 512/107/197/199/67/6.2 vs live 517/110/200/201/
    67/6.2** (QMD/TopHt bit-exact); 2030 TPA 499/500; 2040 452/455; 2050 511/532 — within a few TPA all cycles
    (residual = the balmod-density class + sprout-RNG). Debug-verified: LS ASSPTN (essprt.f:1228) IDENTICAL to
    NE/CS (ISHAG=10, same 7th-order Crouch poly, ASBAR/ASTPAR via estump.f:110-111 = jl's accumulation exactly).
  - **FFE** (stand 4) — the one remaining LS subsystem (a distinct multi-turn sub-port, as for CS/NE). LS
    has NONE of the 5 fire data CSVs (fire_species_props / fire_biomass / fire_fuel_dead / fire_fuel_live /
    fire_fuel_models). SCOPE: the per-species fire props (V2T volume-to-weight, snag decay/fall/bark classes,
    biogrp) are set in **ls/fmvinit.f** (a 68-species SELECT CASE, e.g. V2T 24.9/25.6/21.2/… by species) →
    build fire_species_props.csv; ISPMAP (ls/fmcrow.f:50, LS sp→shared fire-index SPI, mostly identity) →
    crown-biomass FMCROWE; fuel-loading tables (fmcba/fmcfmd, Anderson-13 weighted models) → the 3 fuel CSVs;
    BIOGRP (fmcblk.f:27). Then wire the LS FFE like CS (fmeff/fmbrkt/fmcrowe/fmcblk). `:v2t` is ALSO read by
    cuts.jl (cut biomass) — so even non-FFE cut-biomass reporting needs fire_species_props. **Extraction
    feasibility CONFIRMED:** a Python parse of fmvinit.f's `DO I=1,MAXSP; SELECT CASE(I)` cleanly yields all 68
    species' V2T/TFALLCLS/LEAFLF/DKRCLS/SNAGCLS (sp1 V2T24.9/tfall4/leaf2/dkr4/snag2; sp41 V2T21.8). Remaining
    for the fire_species_props.csv: the derived snag_decayx/fallx/alldwn (fmvinit.f:796 TFALLCLS→ + :823
    SNAGCLS→ class maps), bark_eqnum, biogrp (fmcblk.f:27), ls_spi (fmcrow.f:50 ISPMAP). Then the 3 fuel CSVs +
    the FFE wiring (mirror CS: fmeff/fmbrkt/fmcrowe/fmcblk/fmcba/fmcfmd) + fire-stand validation. A scoped
    multi-turn sub-campaign.
    **PROGRESS (data side started; FFE stand now advances gap-by-gap):**
    - **fire_species_props.csv BUILT** — parsed ls/fmvinit.f SELECT CASE (V2T/TFALLCLS/LEAFLF/DKRCLS/SNAGCLS) +
      SNAGCLS→(decayx/fallx/alldwn) {1:(.4,1.66,10),2:(.8,1.33,30),3:(1,1.16,50),4:(1.2,1,50),5:(1.5,.83,50),6:(2.3,.53,50)}
      + BIOGRP (fmcblk.f, inline-comment-stripped) + ls_spi=ISPMAP (fmcrow.f, IDENTITY 1:68). bark_eqnum=0 placeholder
      (NOT read by jl FFE). Unblocked `:v2t` (also gates cuts.jl cut-biomass).
    - **fire_biomass.csv BUILT** — bio_group=biogrp (verified fmcbio.f:72 IGRP=BIOGRP; CS bio_group==biogrp ∀sp). Unblocked `:bio_group`.
    - **fire_fuel_models.csv** — copied from CS (Anderson-13; CS==NE identical).
    - **fire_fuel_dead.csv + fire_fuel_live.csv BUILT** — FUINI (26 dead types × 11 size classes) + FULIV (9 live
      types × herb/shrub) extracted from ls/fmcba.f DATA (column-major, inline-comments stripped). ⇒ **ALL 5 FFE
      DATA CSVs now present.**
    - **REMAINING = FFE CODE (a sub-campaign, mirrors the CS FFE port):** (1) `ffe_fuel_loading` LS branch — the
      forest-type→fuel-type mapping: FMLSFT (ls/fmlsft.f IFORTP→IFFEFT: {102-105,381}→1 wp/rp, {101}→2 jp, {121-127}
      →3 sf, {181}→4 rc, …) + ISZCL size class + FTLIVEFU/FTDEADFU (fmcba.f: IFFEFT{1,3,4}→(ISZCL)1-3, {2}→4-6,
      default→7-9; FTDEADFU similar at :253). (2) fire behavior (FMCFMD weighting), (3) effects/mortality (fmeff/
      fmbrkt bark groups), (4) crown biomass (fmcrowe via ls_spi), (5) snag dynamics, (6) carbon + fire-stand
      validation. Suite 6478/2 (all 5 CSVs are LS-only data).
    - **★★★ FFE FIRE STAND VALIDATED — fire behavior BIT-EXACT, fire mortality within ~3% of the kill.** Ported the
      LS fuel-loading mapping into fuel_loading.jl (`_ls_iffeft`=fmlsft IFORTP→IFFEFT; `_ls_dead/live_fuel_type` =fmcba
      FTDEADFU/FTLIVEFU with the redcedar IFFEFT=4→single-row-10 quirk; `ls_dead/live_fuel_loading`+fmcba.jl LS
      branches) + extracted the real `bark_eqnum`=EQNUM (ls/fmbrkt.f). THEN found & fixed the three real FFE-fire bugs
      (clean fire-ONLY key `ffe_fireonly.key`, no THINDBH confound; live-validated each step):
        1. **Fuel-model selection** — jl fell through to the SN forest-type path and picked model **6**; LS has its own
           cover-type × PERCOV × season selection (ls/fmcfmd.f). PORTED as `ls_select_fuel_models` (fuel_model.jl):
           cover-type metagroups (JPCT/NHCT/RPCT/MWCT/OACT/ABCT by BA) → per-type branch → post-selection natural-fuel
           candidates 10/12/13, resolved by FMDYN over the 22-class LS XPTS/IPTR (`_FMD_XPTS_LS`/`_FMD_IPTR_LS`,
           `_fmdyn` generalized to `length(eqwt)` + optional `iptr` class→model map). Our stand → ICT=MWCT →
           **model 10 @ 100%** == live. (Activity-fuel model-11 path deferred = LATFUEL false, same as SN/NE.)
        2. **Moisture table** — LS used the SN `_FM_MOIS`; the LS `fmmois.f` preset table is IDENTICAL to NE's
           `_FM_MOIS_NE` (verified all 4 conditions). Added `fm_mois_table(::LakeStates)=_FM_MOIS_NE`. Fixed jl 10hr
           0.07→0.08 and 3+″ 0.17→0.15.  ⇒ fire behavior now **flame 3.46/live 3.4, scorch 13.3/live 13.0**.
        3. **FMEFF mortality adjustments** — jl (a) leaked the SN Regelbrugge MORTGP groups into LS (fire_tree_mortality
           now forces group 6 for LS/NE, gated to SN/CS per fmeff.f:196) and (b) skipped the LS dormant-season
           reductions. PORTED the LS branch in `fire_mortality_adjust` (ls/fmeff.f:278-300): before-greenup conifers
           (sp≤14) ×½, balsam-fir (sp8) floor 0.7, maples {18,19,26,27,51,52}<4″ die, hardwoods (sp>14) ×0.8 (oaks
           30:36 ≥2.5″ ×½), hardwoods ≤1″ die.
      **RESULT vs live fire-only (SIMFIRE 2003, season 1):** 2003 **524/104 == live** (pre-fire, bit-exact); the fire
      kill lands 2003→2013 (the fire is at the start of that cycle): **jl 2013 188 vs live 177, 2023 184/173,
      2033 181/171** — was a 37+ over-kill (130) before these fixes; now Δ11 (~3% of the 347-TPA kill), jl slightly
      UNDER-kills. **No-fire baseline is BIT-EXACT** (jl 2013 502/134 vs live 504/136, ULP-scale) ⇒ density mortality
      is fine; the entire fire gap is FMEFF. Suite 6478/2, no regression (all changes LS/NE-gated).
      **★★★ 4th FFE BUG FIXED — the fire mortality is now BIT-EXACT.** The Δ11 was the WHITE-PINE fire bark: a
      per-tree live FMEFF stamp showed jack pine (sp1, bark-driven CSV=0) PMORT 0.265 == jl BIT-EXACT, but working
      back from live white-pine (sp5) PMORT 0.303 gave FMBRKT bark ≈ 0.27 while jl computed 0.5584 (2× thick → lower
      mortality → under-kill). ROOT: `fire_bark_thickness` applied the SN SHORTLEAF-PINE Harmon quadratic to sp5 for
      ALL non-CS variants (`slpine = CS ? 3 : 5`), but shortleaf pine is sp5 ONLY in SN (sn/fmbrkt.f:126) / sp3 in CS
      (cs/fmbrkt.f:133); NE and LS fmbrkt.f are plain `DBH·B1[EQNUM]` with NO special case, and **LS sp5 = eastern
      WHITE PINE** (EQNUM 24 → B1 .045). FIX: gate the Harmon to `Southern ? 5 : CS ? 3 : 0`. ⇒ white-pine bark
      0.5584→0.27==live ⇒ **fire-only 2013 177/89/159 + 2023 173/113/193 BIT-EXACT** (was 188/95). Also fixes a
      MASKED NE bug (NE sp5 was mis-barked; no NE-test regression). Suite 6493/2. The full lst01 FFE-fire stand
      (stand 3) is now maxΔTPA=0/ΔBA=1 vs live. Only-remaining fire residual = the 2003 pre-fire BA Δ1 (105/104,
      the SIGMAR tripling) + PERCOV ~3.4-low crown area (midflame 1.2/1.1 → flame 3.46/3.4 — cosmetic, doesn't move
      the now-bit-exact mortality).
      **★ PERCOV TRACED (cosmetic, fire mortality is bit-exact regardless):** live fmcba.f:193 uses the STORED
      `CRWDTH(I)` = the FOREST-grown crown width (cwidth.f CWCALC IWHO=0), NOT the open-grown. jl's fmcba! recomputes
      `crown_width(...,iwho=0)`. The mapping is correct (crown_width_species.csv JP 10501-forest/10503-open, verified
      vs ls/cwcalc.f), lat/long/elev match (47.38/94.6/10 → Hopkins 47.59) and the crown ratio is right (mortality
      bit-exact). KEY: the CCF that validated at cyc0 uses the OPEN-grown **ek** equation (10503, no climate term),
      but the fire percov uses the FOREST-grown **bechtold** equation (10501, a DIFFERENT never-validated path). So
      jl's percov ~4.5%-low totcra is the untested bechtold forest-grown crown width at the grown 2003 stand.
      ★ 10501 FORMULA VERIFIED IDENTICAL: live cwcalc.f CASE('10501') = `CW=0.7478+0.8712·D+0.0913·CR` (cap 25) ==
      jl's bechtold 10501 (a=0.7478,b=0.8712,cr_coef=0.0913,hi_coef=0). So the forest-grown crown-width FORMULA is
      correct — the ~4.5% percov residual is a subtle CROWN-RATIO-TIMING difference: live uses the STORED CRWDTH
      computed by the crown model at START-of-cycle CR, jl recomputes crown_width(iwho=0) at FIRE-time CR (crown
      recession shifts CR mid-cycle). Cosmetic (moves only the burn-report flame/scorch; fire mortality bit-exact).
  - **★★★ FFE CARBON — snag Stand-Dead pool VALIDATED (LS snag-dynamics port). Suite 6505/2.** The FFE Stand
    Carbon Report (CARBREPT) fire-year StandDead was jl 14.5 vs live 12.0 tons C. Diagnosis (live fmdout BIOSNAG
    stamp: `TOTSNG(1)+TOTSNG(2)` decomposed to snag-bole vs CWD2B-crown; + a single-tree SNAGINIT differential):
    the over-book is the snag **BOLE** (jl 9.99 vs live 7.27 C), not the crown. Two LS-specific snag mechanisms
    jl ran with SN defaults: **(1) snag HEIGHT LOSS** — LS has non-zero default HTX per snag class (fmvinit.f:
    823-875: class1=3.0/2=1.0/3-4=0/5=0.65/6=0.45, hemlock sp12=0), FMSNGHT shrinks HTIH yearly (SN/NE HTX=0).
    jl's `ffe_snag_height_loss!` was a no-op (empty `snag_htx`) + had a latent HTR bug (0.01 both regimes; FVS =
    HTR1=0.1 first-50%/HTR2=0.01 after, fmsnght.f:154-159). **(2) snag FALL** — LS FMSFALL BASE=−0.006·d+0.18
    (fmsfall.f:128, ~2.5× SN's −0.001679·d+0.064311), MODRATE-clamped [0.01,1], small-snag linear breakpoint 18″
    (12″ cedar/tamarack ksp 10,11,14). FIXES: `snag_htx` col in fire_species_props.csv (from snag_cls); seed
    `fs.params.snag_htx` for LS in `kw_fmin!` (SNAGBRK overrides); HTR1/HTR2 in `_snagbrk!`+`ffe_snag_height_loss!`
    (cancels for keyword path); LS branch in `snag_fall_density(...; variant)`; **wired the current-height bole
    truncation for LS/NE** in `snag_bole_carbon` — the CFTOPK block used `_fm_cuft` (=0 for LS ⇒ silently skipped),
    now r9clark at HTDEAD → Behre-truncate to htcur. RESULT bit-exact vs live: single-tree SNAGINIT 2003 den
    20.95/live 20.95, bole 3.702/live 3.70; ffe_carb 2003 den 415.27/live 415.23, StandDead 11.9/live 12.0 (Δ~1%
    = CFTOPK-form residual, r9-recompute vs exact NATCRS+CFTOPK); 2013 den 0.849/live 0.849. .sum cyc0/1 bit-exact
    (snags don't feed growth). Tests: +8 unit (test_snag.jl LS) +4 integration (test_lst01_ffe.jl, ffe_carb.key/
    .tre fixture). NE-note: NE FMSFALL uses an ALGSLP table {(1,5,12)→(0.20,0.0667,0.04)} — NOT ported (NE keeps
    SN form); latent NE snag-fall discrepancy IF NE runs SNAGINIT/snags (out of LS scope).
  - **★★★ Establishment** (stand 5 BARE-PLANT) — now **BIT-EXACT** vs live FVSls (was "tracks closely").
    The residual was a REAL BUG in the planted-seedling default-height RAN acceptance window (estab.f:483/490):
    SN=[0,1.5] but **NE, CS, AND LS all=[-2.5,2.5]**; jl's `Northeast ? (-2.5,2.5) : (0,1.5)` wrongly gave
    CS+LS the SN window, which rejects the low RAN tail ⇒ biased seedling heights HIGH ("BARE-PLANT
    over-sizing"), shifting BOTH the merch report AND the dense-stand self-thinning. FIX: `Southern ?
    (0,1.5) : (-2.5,2.5)` (establishment.jl:111). ⇒ stand5 @2002 seedling DBH mism 46/50→5/50; whole
    trajectory bit-exact (2022 mcuft 1332→1325, 2042 TPA 515→503 == live; 1-unit bdft ULP @2072 only). Also
    faithfully fixes CS (same latent bug, unexercised by CS tests). test_lst01_estab.jl (5 bit-exact asserts)
    added. ⇒ ALL 5 lst01 stands now bit-exact/cornered-ULP. (Original establishment port notes below.)
  - **★ Establishment** (stand 5 BARE-PLANT) — PORTED + VALIDATED (tracks live closely). Ported the LS ESSUBH: `ls_htcalc_height` (NC-128 forward via a MAPLS module const),
    `_LS_ES_HHTMAX` (ls/blkdat.f DATA HHTMAX) + `_LS_ESSUBH_REFAGE` (ls/essubh.f DATA MAPLS = CARAGE, distinct
    from the htcalc map), + LS branches in establish! (es_hhtmax / bc=nothing / the `(H/CARAGE)·min(5,TIME−DELAY)`
    base height). Result vs live BARE-PLANT: **2002 TPA 800 EXACT** (planting works), but seedling sizes LOW
    (2002 jl 2/9/7/8/0.6 vs live 10/38/37/14/1.5 — TopHt 8 vs 14, QMD 0.6 vs 1.5); converges by 2032 (610/142/
    306/286/54/6.5 vs 597/157/330/291/59/6.9). Suite 6478/2, no regression. NEXT: the planted-seedling height is
    too small. TRACED: the ESSUBH **base height is BIT-EXACT** (ZZEH stamp of live essubh.f: JP SI 61.6/CARAGE 20/
    H 28.26/HHT **7.066** == jl `ls_htcalc_height(1,61.6,20)`=28.265 → base HHT 7.066; RN 6.0 == live), AND jl's
    regent NC-128 increment for a 7.066-ft seedling is a healthy 14.94 ft/10yr. So both the base height and the
    growth-model coefficients are right — yet jl's 2002 seedling only reaches 8 ft (live 14). ⇒ it's a
    creation-cycle TIMING/wiring issue: jl isn't applying the full regent creation-cycle growth to the planted
    seedlings (they appear near their bare ESSUBH height instead of grown). NEXT: trace the establish!↔grow_cycle!
    ordering (establish! runs at GRADD/end-of-cycle; the planted seedlings must still take that cycle's REGENT
    growth — check whether jl creates them post-growth so they miss it, and how FVS ESTAB/REGENT sequence handles
    creation-cycle seedling growth). NOT a coefficient bug.
    **★ FIXED:** esgent.f:48 calls `REGENT(.TRUE.)` to GROW the newly-created seedlings from the 5-yr ESSUBH base
    to the cycle-end height. jl's establish! PHASE 2 had `ne_estab`/`cs_estab` growth branches but NO `ls_estab` —
    so LS seedlings entered at the 5-yr ESSUBH (7.066) ungrown. Added `ls_estab` (LS NC-128 increment via
    `ls_htcalc_incr/age/htmax` MAPLS + `ls_balmod` + `_LS_ES_HHTMAX` cap + htdbh⁻¹ DBH; scale=(per−5)/10=0.5; BARE
    stand ⇒ AVH=0). Result: **2002 TopHt 14 EXACT** (was 8), TPA 800. **★★★ THEN FIXED the seedling OVER-sizing
    (BA 12 vs 10) — now BIT-EXACT 2002-2032.** A live regent.f stamp showed the establishment ls_balmod GMOD=0.745
    for jack pine even at BA=0/RMSQD=0/AVH=0 (the rmsqd≤0→omega=b4 branch), but jl passed `rmsqd_e=stand_qmd(s)=0.626`
    — the POST-establishment QMD of the just-added seedlings — flipping ls_balmod to the else branch (GM 1.0) so the
    cohort skipped the 0.745 competition reduction and OVER-grew. FIX: snapshot `rmsqd_pre = stand_qmd(s)` BEFORE the
    seedlings are added (establishment.jl:128) and use it as rmsqd_e. ⇒ **BARE 2002 800/10/38, 2012 777/50/140,
    2022 708/114/266 BIT-EXACT** (was 12/43, 54/148, 699/115); 2032 Δ1. Suite 6493/2 (LS-only). Residual now: the
    terminal 2042 TPA Δ12 (jl 515/live 503) = the NOTRIPLE dense-mortality self-thinning phase (SDI~378), a separate
    small tail (same class as the Control terminal residual) the establishment fix UNMASKED.
  - Thinning SELECTIONS (the cut logic) are shared/work; only the post-cut sprout regen is gated on the above.
- **★ LS PORT STATUS: growth-and-yield CORE = validated faithful drop-in.** Cycle-0 all 10 columns BIT-EXACT;
  multi-cycle tracks within a few % for 40 yr; all CORE models ported + validated (dgf bit-exact, htgf/regent/
  crown/mortality, volume incl. the Scribner-board fix, sitset/htdbh/merch/forkod). Suite 6478/2, no SN/NE/CS
  regression. Remaining = CHUNK 5 breadth (sprout/FFE/estab data+models) + the balmod tripling-order density
  (compounds in dense cycles) + broken-top aspen cubic. NOT complete (don't touch docs/LS_COMPLETE) — the
  breadth sub-ports remain.

- **CHUNK 3 (diameter growth) — MODEL PORTED + WIRED (cycle-1 validation pending CHUNK 4 height growth).**
  LS large-tree DG is the SAME SN-family ln(DDS) regression as CS — `ls/dgf.f` DDS assembly is byte-for-byte
  the CS `dgf!` (DGCON+COR+INTERC+VDBHC/D+DBHC·D+DBH2C·D²+RDBHC·D/QMD5+RDBHSQC·D²/QMD5+CRWNC·CR+CRSQC·CR²+
  SBAC·BAGE5+BALC·BAL+SITEC·SI → DIAGRO=√(D²+eᴰᴰˢ)−D → bark → ln(IB DDS)). Ported `dgf!(::LakeStates)`
  (data/lakestates/dg_coeffs.csv already extracted; LS-specific only in the QMD/CR caps: QMD ISPC{1,2,10,37-39,59}→13
  {11}→15 {17,30-33}→25, CR {17}→60 {60,65}→85, ls/dgf.f:236-251), `ls_dgcons!` (DGCON=0/ATTEN=OBSERV/bark,
  == cs_dgcons!), `_ls_init_crowns!` (no-op — lst01 has inventory crowns). Wired into setup_growth!. Shared-engine
  changes (variant-gated): calibration GST DBH floor `gst_min` now 5 for LS too (ls/dgdriv.f:396 WK3<5.0, == CS);
  crown_ratio_update! Union += LakeStates (LS crown.f == NE TWIGS `10·(BCR1/(1+BCR2·BA)+BCR3·(1−exp(BCR4·D)))`,
  BCR4 used DIRECT like NE — no CS sign-flip; LS BCR1-4 added to species_coefficients.csv). Suite 6472/2, NO
  SN/NE/CS regression. **Blocker to cycle-1 .sum: height_growth!(::LakeStates) missing (CHUNK 4).**
- **★★ CHUNK 4 (full growth spine) — DONE; cycle-1 (2000) 3/6 BIT-EXACT + 3/6 within Δ1 vs live FVSls.**
  TPA 524 ✓ · TopHt 64 ✓ · QMD 6.0 ✓ (BIT-EXACT); BA 103/104 · SDI 202/203 · CCF 209/210 (Δ1 each). The
  ENTIRE LS growth spine now runs end-to-end. Ported this session:
  - **Height growth** (lakestates/height_growth.jl) — NE's NC-128 curve via `_ls_htcoef` (MAPLS into the
    SAME shared 127-row LTBHEC) + the LS-specific `ls_balmod` (ls/balmod.f: OMEGA/BETA/BAMAX1, GM floor 0.2,
    BATEMP=BA>1?BA:1). RMSQD=`p.qmd` (dense.f:250), BA=`p.basal_area`, AVH=`p.avg_height`. Structurally == NE
    htgf; only the map + balmod differ. Data: htgf_coeffs.csv (MAPLS + 8 balmod arrays).
  - **Small-tree growth** (lakestates/small_tree_growth.jl) — ls/regent.f == cs/regent.f (XMIN=3/XMAX=5/
    DGMAX=5/REGYR=10) with the two swaps: htcalc MAPLS + `ls_balmod(sp,D,BA,RMSQD)`.
  - **Mortality** — generic `mortality!` + LS background (Hamilton logistic `ri=1/(1+exp(PMSC+PMD·D))`,
    PMSC(4)/PMD(4) grouped by IMAPLS(68) → per-species mort_bkgd_intercept/mort_bkgd_dbh in species_coefficients.csv);
    VARMRT distribution — ls/varmrt.f efficiency formula is IDENTICAL to NE/CS ⇒ LS joined the `_varmrt_efftr!`
    Union. mort_ri_scale=0.5, SDI gate DBHSDI (set in the singleton).
  - Suite **6472/2** — NO SN/NE/CS regression (touched: mortality/crown/varmrt Unions, gst_min, coefficients
    loader +htgf_coeffs.csv, setup_growth!). test_lst01.jl now has cyc0 + cyc1 testsets.
  - **✓✓ cycle-1 Δ1 FULLY CLOSED (2026-07) — base stand now BIT-EXACT 1990-2040 (50 yrs).** The "balmod
    density basis" trace BELOW is a **MISREAD** (superseded): jl's balmod inputs (85.13/5.145) already match
    live's HTGF DENSE call exactly; the 67.49/114.53 the trace chased are the DGF-calibration and post-triple
    densities, NOT the HTGF one. The real cycle-1 Δ1 was TWO small-tree bugs: (1) the **HTDBH mode-1 DB floor**
    (htdbh.f:343, missing for NE/CS/LS in jl — a sugar-maple seedling got a negative Wykoff-inverse DBH →
    grew 0.68″ vs floored 0.20″; fixed via `db_floor` kwarg on `_htdbh_dbh`), and (2) the **crown-ratio BA
    basis** (crown.f uses RAW per-acre BA, jl used /gross_space; 2 un-crown-capped trees' CR came out 48 vs
    live 46, shifting their dgf DDS; fixed by gating crown BA to RAW for LS). Both live-stamp-verified
    BIT-EXACT (SM seedling 0.2064, the 2 trees' CR 46/dds 2.82518). NOTRIPLE bit-exact ALL 100 yrs. Remaining
    tail 2050-2090 (TPA Δ1-5, BA/SDI exact) = tripling-spread DGSCOR. Suite 6505/2. [Original stale trace kept below for history.]
  - **cycle-1 Δ1 — RMSQD BUG FOUND + FIXED (SDI now BIT-EXACT); remaining BA/CCF Δ1 traced to the balmod
    density BASIS.** Full trace (per doctrine):
    (1) ZZDG stamp of ls/dgf.f PROVED the **LS large-tree DG model is BIT-EXACT** (max|ΔDDS|=0 over all trees;
        BAGE5 80.00, QMDGE5 7.516, CONSPP=0 [no calibration, all FN<FNMIN], BAL, SI all match live).
    (2) Localized to the small-tree/regent path (d<5 = 7 trees/300 TPA = 56% of stand). ZZRG stamp: trees d<3
        (pure regent) BIT-EXACT; the d3-5 BLEND trees diverged (SM d4 jl HTG 3.88 vs live 4.95).
    (3) ROOT CAUSE: **`p.qmd` is NEVER populated (always 0)** — my ls_balmod/height_growth read `rmsqd=p.qmd=0`,
        forcing the `omega=b4` fallback. FIXED: `rmsqd = stand_qmd(s)` in height_growth.jl + small_tree_growth.jl.
        ⇒ SDI 202→**203 BIT-EXACT**; SM d4 HTG 3.88→**4.94** (live 4.95); suite 6478/2 (no regression).
    (4) REMAINING BA/CCF Δ1 (jl BA 103.28 vs live ~103.5): ZZBM stamp of ls/balmod.f shows live passes
        **BA=67.49, RMSQD=5.472** to balmod, but jl passes the CURRENT stand values **85.13, 5.145**. This is
        the precise remaining cause. **Basis UNRESOLVED** (candidates REFUTED by measurement): NOT current
        (85.13/5.145), NOT fully-backdated (jl `_backdate_dbh!` gives 59.49/4.301), NOT a clean DBH threshold
        (live RMSQD 5.472 is ABOVE current-all 5.145 but BELOW d≥0.5's 6.411 — no threshold yields both 67.49
        AND 5.472). Note dgf! uses the CURRENT BA and matched live BIT-EXACT, so balmod uses a DIFFERENT density
        than dgf in the same cycle. **ROOT NAILED (ls/dense.f ZZDN stamp):** dense's two-pass LREDO/LBKDEN logic
        (dense.f:242-267) sets `BA = OLDBA` = the **backdated BAT from pass 1 = 67.49**, and `RMSQD = SQRT(TSUMD2/
        TPROB)` from the **current pass = 5.472**; RAT=FINTH/FINT=0.5 only interpolates the NEXT-cycle OLDBA, not
        BA. BUT the current-pass **BAT=101.13** while jl's current stand_ba=**85.13** — a **1.1879× normalization
        gap** (jl stand_ba is per-STOCKABLE = 85.13/gross_space1.1=77.4 displayed; FVS dense BAT is higher-normalized),
        and dense **RMSQD=5.472 ≠ jl stand_qmd 5.145** (different tree set/weighting) — even though the DISPLAYED
        cycle-0 density matched BIT-EXACT (BA77/QMD5.1). So the fix is NOT a balmod-input swap; it's reconciling
        jl's internal density (stand_ba/stand_qmd) with FVS **dense's BAT/RMSQD basis** (the ×1.1879 + the RMSQD
        tree-set), then feeding ls_balmod the backdated-BAT + that RMSQD. Deep density-normalization trace = the
        ONLY open cycle-1 item (small: BA 103/104, CCF 209/210). Then multi-cycle + thinning/estab/sprout/FFE.
    (5) **DEEPEST ROOT NAILED (ls/dense.f ZZTR per-tree stamp) — the Δ1 is a TRIPLING-ORDER artifact.**
        Per-tree D and P (PROB) MATCH jl's dbh/tpa BIT-EXACT (tree1 P=5.545=jl tpa; tree13 D0.1 P=90.0=jl;
        tree23 D3.2 P=30.0=jl) — the tree DATA is right. But the dense list the balmod call iterates includes
        EXTRA high-index records (I=2999 D7.2 P28.29, I=3000 D34.6 P1.23, …) absent from jl's 27-record list
        ⇒ BAT 101.13 vs jl 85.13. These are FVS **tripling-expanded records** (GRINCR triples BEFORE the
        density/HTGF pass). jl DEFERS tripling (diameter_growth! returns a stash, "no records yet"), so at
        height_growth! time jl has the UNTRIPLED list ⇒ ls_balmod reads untripled BA/RMSQD. **This is the
        complete root of the cycle-1 Δ1.** NE/CS are UNAFFECTED (their htgf uses ne_badist! BAL-distribution,
        NOT RMSQD/BA) — why deferred-tripling was fine until LS. FIX (deferred, non-trivial): give ls_balmod
        the tripling-expanded density (materialize the triple before height_growth! for LS, or compute a
        tripled-list BA+RMSQD for balmod). Small residual (BA 103/104, CCF 209/210); root fully traced (NOT
        ULP) — documented + deferred like the accepted SN-COMPRESS class; cycle-1 test keeps a `≤1` tol on BA/CCF.
- **MULTI-CYCLE projection (lst01 stand-1, 10 cycles vs live) — ★★★ TERMINAL RESIDUAL ROOT-CAUSED & FIXED.**
  The 2040-2050 under-kill (was TPA 292/271 +21, 211/188 +23) traced NOT to balmod but to the **dg_resid_sd
  (SIGMAR) tripling-spread placeholder**. Decisive diagnostic: with **NOTRIPLE** jl==live BIT-EXACT at 2000
  (524/101/198) ⇒ the whole cycle-1 BA Δ1 + terminal under-kill is the TRIPLING DBH-spread, NOT balmod
  (balmod affects height not the DG-spread; everything matches with NOTRIPLE DESPITE the "wrong" balmod density
  — it is a REAL but INERT difference, do NOT chase). Chain: terminal Δ21 → low DIA0 (jl 8.483 vs live 8.54,
  SDIMAX bit-exact 455.9 ⇒ mortality FAITHFUL) → BA Δ1 → NOTRIPLE-exact → tripling spread. The .trl 2000 dump
  showed jl over-represented 5-9″ (BA 49 vs live 37) and under-represented 9-15″/15+ (44/0 vs 58/1.8) = too-NARROW
  spread. Spread magnitude = DG_FU·ssigma·rhocp, ssigma←vardg←**dg_resid_sd**, which jl had as a **uniform 0.6
  placeholder for all 68 species** (flagged CHUNK-3 TODO). FIX: extracted the real per-species SIGMAR (blkdat.f:253
  DATA SIGMAR, 0.552-0.913) into species_coefficients.csv. **RESULT: 1990 bit-exact; 2000 BA 103→104==live;
  2010-2030 Δ1-3; 2040 292→265 (live 271, Δ21→Δ6, BA 192/SDI 308 BIT-EXACT); 2050 211→183 (live 188, Δ23→Δ5,
  BA 193 exact).** Whole 60-yr trajectory now Δ1-6 (was Δ1-23). Suite 6493/2 (FFE test bounds updated — SIGMAR
  shifted the fire stand pre-fire BA 104→105 Δ1; fire behavior model-10/flame/scorch UNCHANGED). Remaining tail =
  2040-2050 TPA Δ5-6 (jl now slightly OVER-kills, BA/SDI exact) — a small residual, no longer the big lever.

- **CHUNK 2 (volume) — 3 of 4 columns BIT-EXACT; 1 small residual.** LS grinit defaults METHB=METHC=6
  ⇒ **900CLKE (R9 Clark)**, the SAME cubic model NE/CS use (DVEE/METHC=5 is a per-species opt-in, not
  hit by lst01). Routed LS through the shared `compute_volumes_ne!` (FIA-keyed/variant-agnostic; national
  r9clark_coef.csv). vs live lst01 1990 (1551/1338/480/1887): **SCuFt 480 ✓ BIT-EXACT · BdFt 1887 ✓
  BIT-EXACT**; TCuFt 1546 / MCuFt 1333 (Δ5, ~0.3% low = open).
  - **★ BOARD-FOOT FIXED (was 48% high 2786→1887).** Debug-stamp of live ls/vols.f (ZZBF per-tree)
    proved LS `.sum` BdFt is **Scribner** (vol2), NOT the International ¼" (vol10) NE/CS report — JP
    D11.5/H73 live=62 == jl now 62. r9clark.f R9BDFT computes BOTH: vol(2)=`nint(len·scrbnr(int(dib)))`
    (Scribner Decimal-C, 120-factor table) + vol(10)=International polynomial; r9cor's cf3(vol2)==cf4(vol10)
    and the DIB round INT(dib+0.499) are identical ⇒ Scribner differs ONLY in the per-log kernel. Ported:
    scribner_factor.csv (r9clark scrbnr(120)) + `_r9_scrib_log` + refactored `_r9_intlqtr_bf` into a shared
    `_r9_bucked_bf(...,logfn)` core (International/Scribner wrappers) + `board_scribner` kwarg on r9clark_cubic,
    set for LS in compute_volumes_ne!. NE/CS stay bit-exact (International path untouched). test_lst01.jl now
    asserts scuft+bdft; suite 6472/2.
  - **✓ CLOSED (was OPEN Δ5) — broken-top aspen cubic now BIT-EXACT.** The D8.0 top-killed aspen (`trunc=5600`,
    norm_ht=6300, .tre tree 6, jl tree idx 5): jl CFV now **8.962141 == live 8.96214** (was 8.889936). Re-verified
    2026-07 via a per-tree TreeList dump — the broken-top CFTOPK path reconciled somewhere in the later volume
    fixes; ALL cycle-0 `.sum` volume columns (TCuFt 1551 / MCuFt 1338 / SCuFt 480 / BdFt 1887) are now bit-exact,
    AND the per-tree top-killed aspen cubic is bit-exact. CHUNK 2 volume = COMPLETE at both stand and per-tree level.

- **★ CHUNK 1 COMPLETE — all 6 cycle-0 stand columns BIT-EXACT vs a FRESH live FVSls lst01 (1990):**
  TPA 536 · BA 77 · SDI 160 · CCF 171 · TopHt 63 · QMD 5.1. test/integration/test_lst01.jl written +
  green; suite 6470/2 (no SN/NE/CS regression). CCF closed by reusing the shared eastern crown-width
  equation library: LS crown_width_species.csv (alpha→CWEQ forest/open, parsed from ls/cwcalc.f SELECT
  CASE — all 68 species covered, 71 unique CWEQ codes all present in NE's crown_width_equations.csv,
  which was copied verbatim). `LakeStates` now exported. Full live lst01.sum all-cycle ground truth
  captured for later chunks: 2000 524/104/203/210/64/6.0 (vol 2243/2059/802/3313); 2010 505/135/251/252/
  63/7.0 (3122/2926/1600/6965); 2020 450/162/286/281/67/8.1; … forest-type col=401 (=FORTYP, not the NF).
  Landed earlier this session:
  - Infra: `LakeStates` singleton (MAXSP 68, YR 10, Zeide, RNG 55329, morts RI=0.5/DBHSDI gate);
    `variant_from_code("LS")`; FVSjl.jl includes.
  - Data (`data/lakestates/`): species_roster (68, verified vs blkdat JSP/FIAJSP/PLNJSP),
    species_translation (copied, col 5 target_ls, other=49 OH), species_coefficients (codes +
    bark_slope=BKRAT/bark_intercept=0 + varmrt_varadj=VARADJ + sdi_max_default=SDICON + dbh_min=5;
    **TODO** dg_resid_sd=0.6 placeholder [chunk3 SIGMA] + estab 0.5/20 [chunk5]), dg_coeffs (full LS
    DDS model: INTERC/VDBHC/DBHC/DBH2C/RDBHC/RDBHSQC/SBAC/BALC/CRWNC/CRSQC/SITEC/OBSERV),
    htdbh_coeffs (Wykoff HT1/HT2 + Curtis-Arney SNALL P2/P3/P4 + SNDBAL DB + IWYKCA), site_sicoef
    (two 68×68 SICOEF matrices, Fortran I-fastest flat order).
  - Code: species.jl (init_blockdata! + spctrn col 5 / other 49), site_index.jl (ls_site_index_setup!
    full SICOEF fan-out + aspen-41 2nd pass + SDIDEF=SDICON; ls_forkod_defaults! JFOR→IFOR + lat/long/elev),
    `_ls_merch` in r9clark_vol.jl + LS branch in init_merch_standards! (softwood≤14, aspen 40-42, TOPD=4).
  - **NEXT (finish CHUNK 1): CCF** — port LS crown-width (ls/cwcalc.f / ls/ccfcal.f) + stocking so
    stand_ccf ≠ 0 (target 171). Then CHUNK 2 volume (LS = R9 **DVEE Gevorkiantz** + Clark mix by
    per-species METHC, NOT pure Clark — new path).
  - **Watch (variant-hardening):** the engine `_htdbh_wykoff` ifor==3 Allegheny override is NE-specific;
    an LS stand mapping to IFOR=3 would wrongly trigger it (lst01 is IFOR=2, so inert now). Gate it NE-only
    before any LS IFOR-3 stand is validated.

- **★★★★★ LS PORT COMPLETE (LS_COMPLETE flipped).** Meets the SN/NE/CS validated-drop-in bar: bit-exact
  end-to-end barring only the accepted tripling-spread ULP class. Decisive proofs (vs live FVSls):
  (1) **NOTRIPLE ⇒ jl == live BIT-EXACT across the ENTIRE projection** (1990-2060+, all cols
  TPA/BA/SDI/TopHt/QMD) — the deterministic LS model is fully faithful. (2) **Default drop-in (tripling ON,
  lst01-first): BIT-EXACT 6 cycles / 50 yr (1990-2040 ΔTPA=ΔBA=ΔSDI=0)**, then Δ1-5 TPA / Δ0-1 BA late tail
  = accepted tripling-spread (cornered: per-species SIGMAR ~2-ULP spread × discrete self-thinning deletion).
  (3) cycle-0 all 10 cols bit-exact incl volumes (1551/1338/480/1887). (4) FFE fire, establishment, and
  post-fire + cut STUMP SPROUTING (fmkill.f:80 fire-kill→ESTUMP→ESUCKR) all live-validated; default-mode
  (no-NOAUTOES) surface probed across fire/cut/no-disturbance with NO new gaps. (5) No SN/NE/CS regression
  (suite 6513/2; the 2 broken = pre-existing accepted SN COMPRESS + NOHTDREG).
  - **Two prior "open" flags were STALE** (re-verified against the current live differential): the
    "balmod cycle-1 Δ1 basis UNRESOLVED" is a deferred-tripling artifact that is REAL but INERT (NOTRIPLE
    is bit-exact despite it) and its actual driver (SIGMAR) was fixed; the "TCuFt/MCuFt Δ5 open" volume
    residual is now bit-exact. Neither is a real gap.
  - **Residual class (accepted, = SN/NE/CS):** the tripling DBH-spread realization (Δ1-5 TPA in late
    self-thinning cycles, BA/SDI usually exact), amplified by discrete near-tied-record deletion. jl's
    deferred-tripling gives ls_balmod an untripled density basis vs FVS's pre-density triple, but this is
    output-inert (NOTRIPLE proof). Documented + accepted like the SN COMPRESS-eigensolver class.

- **★★★★★ FINAL FIX — the last real LS residual CLOSED (QMDGE5 cumulative cap).** A site-productivity
  sweep (NOTRIPLE, SITECODE 40/50/70/80) exposed a Δ2-5 TPA terminal residual at non-canonical sites.
  Traced end-to-end: terminal mortality residual → BAMAX BA-cap amplification (proven bit-faithful) →
  per-record diameter growth → sp5 (white pine) DGF off 16% → **QMDGE5 jl 14.5 vs live 13.0**. ROOT:
  FVS caps the STAND-WIDE QMDGE5 in place as it walks species in INDEX order (ls/dgf.f:362-390), so a
  species is capped by all lower-indexed present cap-species — white pine (uncapped) sees the 13" cap
  jack pine (sp1, cap-13) applied upstream. jl had applied a LOCAL per-tree cap. FIX: `dgf!` precomputes
  a per-species effective QMDGE5 via the cumulative species-order walk. RESULT: **LS NOTRIPLE now
  BIT-EXACT vs live across ALL site indices AND all cycles** (deterministic model fully faithful). The
  identical latent bug was fixed in CS (`centralstates/diameter_growth.jl`) too. The LS mortality
  (background/VARMRT/BAMAX) was proven bit-faithful throughout. Suite 7658 pass / 2 broken (the 2 =
  accepted SN COMPRESS + NOHTDREG). ⇒ `docs/LS_COMPLETE` flipped; tag `FVSsn+ne+cs+ls-done`.
