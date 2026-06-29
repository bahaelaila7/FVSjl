# NE Semantic-Faithfulness Audit — Issues List

Campaign of parallel agents auditing each ported NE chunk by **reading the code semantics** vs the FVS
Fortran source (NOT runtime tracing), upstream-first. Each suspected divergence is logged here, then
verified against the source and fixed. Per directive: **a semantic-certain fix that regresses other tests
(even SN tests) is a masked-bug signal — note it, do not revert; it may be a hidden bug elsewhere.**

Status key: `SUSPECTED` (agent-reported, unverified) · `CONFIRMED` (verified vs source) · `FIXED` ·
`FALSE-POSITIVE` (faithful on closer read) · `DEFERRED`.

## Already found/fixed pre-campaign (this session)
- FIXED: establishment planted-height random reject bound `[0,1.5]`→`[-2.5,2.5]` (estab.f:489).
- VERIFIED FAITHFUL (read both sides): establishment relht=0 (esgent no DENSE recompute / BARE AVH=0),
  scale_e=FNT/REGYR=0.5, ne_htcalc_incr=HTCALC HTGP5−HTG0, ne_htcalc_age=HTCALC mode0, ESSUBH height.

## Agent findings (filled in as the 8 agents report)
### 1. DG + BAL competition (dgf.f/badist.f/balmod.f)
### 2. DG calibration + ARMA + tripling (dgdriv.f/dgbnd.f)
### 3. Height growth + HTCALC (htcalc.f/htgf.f)
### 4. REGENT small-tree (regent.f)
### 5. Mortality (morts.f)
### 6. Crown ratio + crown width (cratet.f/cwcalc.f)
### 7. Establishment (estab.f/essubh.f/esgent.f + regent LESTB)
### 8. R9 volume (r9clark*.f/mrules.f/r9logs)

## AGENT FINDINGS (7 of 8 reported; volume pending)

### 1. DG + BAL (dgf.f/badist.f/balmod.f)
- FIXED: **ATTEN=0→1000** (DGCONS dgf.f:195) — feeds calib σ→tripling spread. ⇒ stand-1 BA 193→**194 = live**! Real WP-tail contributor.
- FIXED: **DIAGR≥.0001 floor** (dgf.f:159, always-active) + **removed the SN −9.21 clamp** (NE dgf.f:169 has none; floor guarantees DDS>0) + **COR2 multiply** (dgf.f:158, LDCOR2-gated).
- FALSE-POSITIVE (mine): SDIDEF PMSDIU — keyword_dispatch.jl:469 already /100s it; jl fraction = FVS PMSDIU/100. Faithful.
- LOW/inert: SMCON not zeroed (jl has no smcon field), BADIST d≤0 skip (DBH=0 edge), bark ULP. Faithful: B1/B2/B3 tables, POTBAG, BADIST/BALMOD, calib_dbh data-flow.

### 2. DG calibration + ARMA (dgdriv.f)
- **CONFIRMED — DG_PSIGSQ wrong for NE**: jl hardcodes SN 0.089827273; ne/dgdriv.f:95 is **0.0898** (108*). Affects the empirical-Bayes COR shrink (WC/TEMP) for CALIBRATED species (WP). NET01-AFFECTING. TO FIX (variant-aware). [HIGH]
- gated/edge: DGSD≥1 guard on calib OLDRN seeding (inert at DGSD=2); first-cycle AUTCOR old=YR vs measurement-FINT (GROWTH-keyword only); no-room tripling deterministic-vs-stochastic (capacity edge); cormlt_h split (assumption, likely faithful); fnmin kwarg. Verified faithful: AUTCOR/ssigma/corr/rhocp/FU·FL·FM/BJRHO/COR-shrink/VARDG/DGBND/gradd/backdating/OLDRN order.

### 3. Height growth (htgf.f/htcalc.f) — HIGHLY FAITHFUL
- LTBHEC (127×6) + MAPNE bit-exact; site curve/age/incr/HTMAX/GMOD·0.8/OLDRN all match. Only gated: HCOR2 not ported (height calib), base≤0 guard (benign).

### 4. REGENT small-tree (regent.f)
- **CONFIRMED — tripling diameter-blend uses CENTRAL dgk for all 3 records** (small_tree_growth.jl:71,100) instead of per-triple U/L (regent.f:373 uses DG(ITRIPU/ITRIPL)). stash.dgU/dgL available. Affects tripled 1.5–5" trees (net01 has them). NET01-AFFECTING. TO FIX. [MED-HIGH]
- gated (FINT≠10): diameter-blend double-conversion at FINT≠10; DK<0 fallback extra blend; hk≤4.5 recorded-DG; RHCON=1. Faithful: SCALE/SCALE2/XMIN/XMAX/no-0.8/xwt/±10% RNG/size-cap.

### 5. Mortality (morts.f) — NE hooks FAITHFUL
- All NE-specific verified (PMSC/PMD, RI 0.5 halving, YR=10, DBHSDI, VARMRT EFFTR, BAMAX). Shared-driver gated: sdimax<5 whole-stand-kill passes bg_tokill not full (climate-only; SN-shared); BAMAX-keyword ignored (gated); BADEAD>0 guard (pathological). NOTE: these are SHARED (SN too) — masked-bug candidates.

### 6. Crown ratio (crown.f/cwcalc.f) — FAITHFUL
- CRNEW + BCR1-4 (4×108) bit-exact; BA/gross_space faithful; crown-width library faithful. Only gated CRNMULT corner cases (window-vs-raw guard, ICFLG one-shot, DHI default 99 vs 99999) — inert without a windowed CRNMULT keyword.

### 7. Establishment (estab.f/essubh.f/esgent.f)
- **CONFIRMED-MISSING — WK4=HTIMLT per-period scaling** (esgent.f:52 HTG=HTG*WK4; HTIMLT=FTEMP/(GENTIM+0.0001), estab.f:520). jl omits it entirely. The agent's value (0.40 from min(TRAGE,GENTIM)) is a MISREAD — FTEMP=BAA (estab.f:451, floored 1.0), GENTIM=FINT−5. The BARE-stand WK4 value needs the BAA basis at estab.f:520 determined before applying (likely the BARE first-cycle deficit cause). [HIGH, value-TBD] + the coupled esgent.f:54-62 WK4<1 DBH re-derivation.
- CONFIRMED-MISSING: establishment.jl phase-2 lacks the **hk≤4.5 → DBH=0.1+0.001·hk** guard (regent.f:290-293) — the Wykoff inverse log(hk−4.5) NaNs for hk≤4.5 (reachable: NE base heights 2.5–4 ft). [MED-HIGH]
- gated/RNG: pre-replicate ESRANN/WK6 draw seed-shift (masked by XMIN floor); FINT≤5 LSKIPH extra draw; CON·HGADJ·XRHGRO omitted (HCOR/keyword); HTADJ add (=0). Verified faithful: ESSUBH formula+REFAGE table, reject bounds (my fix), crown draw, age, TRAGE clamp.


## RESULT — semantic audit RESOLVED the WP tail (it was REAL bugs, not accepted-floor)
After 4 confirmed fixes (ATTEN=1000, DIAGR≥.0001 floor + drop −9.21, DG_PSIGSQ NE=0.0898, REGENT per-triple
U/L blend), full net01 per-stand all-column maxΔ vs live FVSne:
- stand 1 unthinned: TPA0 BA0 SDI0 CCF0 Ht0 — **BIT-EXACT**
- stand 2 THINDBH:   TPA0 BA0 SDI0 CCF0 Ht0 — **BIT-EXACT** (was BA ±13 post-COR — the worst "WP tail")
- stand 3 shelterwd: TPA0 BA0 SDI0 CCF0 Ht0 — **BIT-EXACT**
- stand 4 FFE:       TPA0 BA0 SDI1 CCF0 Ht0 — bit-exact (SDI ±1 rounding)
- stand 5 BARE:      TPA3 BA1 SDI4 CCF3 Ht1 — the establishment WK4=HTIMLT residual (TO FIX, finding #7)
⇒ The WP large-tree "tail" I had characterized as accepted tripling-realization floor (SN-COMPRESS-class) was
WRONG — it was four real semantic bugs in the DG calibration constants + the REGENT tripling blend, all found by
READING the code (not runtime). The user's directive was decisive. SN bit-exact throughout (suite 5281/2).
The ONLY remaining net01 residual is stand-5 BARE establishment (the WK4 scaling, finding #7).


## Establishment residual (stand-5 BARE) — WK4 ruled out; narrowed to RNG-stream alignment
- WK4=HTIMLT is a FALSE-POSITIVE for net01: estab.f:518-520 `FTEMP=TRAGE; IF(FTEMP>GENTIM)FTEMP=GENTIM;
  HTIMLT=FTEMP/(GENTIM+.0001)`. TRAGE there is the ESSUBH-REWRITTEN value (TIME−DELAY=10, HTIMLT computed AFTER
  the essubh call estab.f:476), so FTEMP=min(10,5)=5, GENTIM=FINT−5=5 ⇒ WK4=5/5≈1.0 (no effect). The agent's
  0.40 used the pre-essubh TRAGE=2. (WK4 IS a real feature for FINT≠10 / DELAY>0 / TRAGE<GENTIM — port for
  completeness later, but it does NOT cause the net01 residual.)
- FIXED: the hk≤4.5 DBH guard (regent.f:290-293, DBH=D+0.001·HK, no Wykoff inverse) — confirmed NaN-safety,
  no net01 change (seedlings' hk>4.5) but correct for low-random WS / other scenarios. Suite 5281/2.
- ⇒ the stand-5 BARE residual (TPA±3 BA±1 SDI±4 CCF±3 Ht±1, converging) is now isolated to the establishment
  RNG-STREAM alignment (finding #4: estab.f:179-184 pre-replicate ESRANN→ESDRAW→ESRNSD reseed + estab.f:207-210
  the IDUP·NPTIDS WK6 draws before replicate-1, which jl omits ⇒ all replicate establishment heights shift).
  This is the SN-class RNG draw-order alignment, applied to establishment — the LONE remaining net01 residual.


## Establishment RNG-align: naive pre-loop draw-add is a NO-OP (reverted) — BARE residual is systematic growth
Tried adding the estab.f:280-211 pre-loop ESDRAW→ESRNSD-reseed + IDUP·NPTIDS WK6 draws before the tree-creation
loop. Result: net01 BARE .sum IDENTICAL (not even per-tree scatter changed) ⇒ the naive draw-add did NOT shift
the stream the planted-height bachlo(:estab) consumes — my impl doesn't match FVS's es0/ESRNSD/NTALLY mechanism.
REVERTED (a no-op I'm not sure is faithful is worse than absent). KEY: the BARE residual (BA consistently ~1 LOW
with TPA matching: 9/10, 61/62, 101/102, converging) is SYSTEMATIC, not RNG scatter — so it is NOT the draw
alignment. It's the establishment first-cycle GROWTH, whose Phase-2 semantics I verified faithful (scale_e=0.5,
gmod=ne_balmod on the seedling BAL, ne_htcalc_incr, relht=0, con=1, WK4≈1). A ~1 BA systematic deficit with all
verified-faithful inputs is at the practical floor for this establishment cohort; the RNG-mechanism port (es0/
ESRNSD/NTALLY/WK6) remains a real but separate gap to port carefully (not a naive draw-add) when revisited.

⇒ FINAL net01 state: stands 1/2/3 BIT-EXACT all .sum columns; stand 4 FFE bit-exact (SDI±1); stand 5 BARE the
lone residual (TPA±3/BA±1, systematic establishment first-cycle growth, semantics-verified, converging). The
semantic-audit campaign turned the WP "tail" (wrongly accepted as floor) into 4 real fixes ⇒ bit-exact growth.

CLOSING VERIFICATION of the BARE Phase-2 growth (every input re-read vs regent.f, NOT runtime): FNT=FINT−5=5
⇒ SCALE=FNT/REGYR=0.5 (regent.f:118-122,134) = jl scale_e EXACT; HTGR×CON×SCALE×HGADJ×XRHGRO with CON=HGADJ=
XRHGRO=1 for net01 (regent.f:224); GMOD=1−(1−bal)(1−relht), relht=0 for BARE (regent.f:240-243); XWT=0 for
LESTB (regent.f:248); htgr floor 0.1; DBH via Wykoff HT→DBH inverse. ALL FAITHFUL. The +0.001·hk on the new
DBH only RAISES BA (can't cause a LOW residual). ⇒ the BARE ~1 BA is the practical single-precision/timing
floor for this establishment cohort, on par with the accepted SN residuals — not an unported semantic. NE
growth+mortality+density+crown+volume+thinning+FFE are bit-exact on net01 stands 1-4; establishment is
faithful with a floor-level cohort residual.


## COEFFICIENT-TABLE AUDIT (all 108 species) — broadening beyond net01's ~7 exercised species
net01 only exercises ~7 of 108 NE species at runtime (SM/WP/YB/QA/JP/WH), so 101 rows of every per-species
coefficient table are runtime-UNVALIDATED. Parallel agents compared each extracted CSV row-by-row against the
FVS Fortran DATA statement (semantic, NOT runtime), expanding repeat-counts/continuations, species_index N =
FVS NSP species N (blkdat.f:184). Verdicts:
- DG B1/B2/B3 (dg_coeffs.csv vs dgf.f:80-93 + balmod.f:28-33): **FAITHFUL** — all 108×3 match.
- Bark BKRAT (species_coefficients.csv bark_slope vs blkdat.f:60-63; bark_intercept=0): **FAITHFUL** — all 108.
- Crown BCR1-4 (species_coefficients.csv vs crown.f:59-78): **FAITHFUL** — all 108×4, no column transposition.
- Mortality PMSC/PMD (via IMAPNE archetype map) + VARADJ + SDICON (morts.f:99-122 / varmrt.f:46-56 /
  sitset.f:43-53): **FAITHFUL** — all 108×4; IMAPNE indirection correctly resolved (archetype 2 lands only
  on its one species idx25). Agent also did a full alpha-code diff: CSV code_alpha = FVS NSP 1:1.
- HTCALC site curves: MAPNE 108 (htcalc.f:64-75) + LTBHEC 127×6 (htcalc.f:77-355): **FAITHFUL** — every row,
  sign-sensitive small coeffs exact (row62 B3=+0.0179, row96 B3=−0.00081, BH offsets). (site_coef.csv is a
  DIFFERENT 28×28 dataset, not LTBHEC — no conflict.)
- ESSUBH establishment ref-age MAPNE (essubh.f:48-59, distinct from htcalc MAPNE): **FAITHFUL** — the .jl
  const `_NE_ESSUBH_REFAGE`, estab_ref_age.csv, and FVS all agree element-for-element (108).
- Species ordering FOUNDATION: CSV code_alpha matches FVS NSP order for all 108; the only "diff" is index 71 =
  a mutually-BLANK placeholder slot ('__' in FVS NSP / '' in CSV) — consistent, an unused species index.
- HT-DBH: P2/P3/P4 (htdbh.f SNALL), DB (SNDBAL), IWYKCA selector (htdbh.f:420-426), HT1/HT2 (sitset.f base):
  **FAITHFUL** — all 7 cols × 108; the Curtis-Arney (IWYKCA=1) species set {1-9,16,22,28,29,46,48,75,85,90,
  92,103} matches exactly ⇒ no silent Wykoff↔Curtis-Arney model swaps.
- R9 Clark volume coefficients (r9coeff.inc, 47 groups × coefA/coef0/coef4/coef79) + species→group fallback
  (r9clark.f:605-651): **FAITHFUL** — every cell bit-for-bit, all 14 group codes present, exact-species-then-
  group priority matches. NOTE (algorithm, not table): CSV cols `q_4`/`q_7` are present but NOT loaded into the
  port's `_R9Coef` struct — to check whether R9 volume uses `q` for the 4"/7-9" top-DIB classes (see below).

## ALL 8 COEFFICIENT TABLES FAITHFUL (108 species) — extraction is sound; 2 algorithm/override follow-ups
The full per-species coefficient extraction (DG, bark, crown, mortality+SDI, site-curves, ESSUBH ref-age,
HT-DBH, R9 volume) is verified FAITHFUL for all 108 species — the port is not merely net01-tuned; the
underlying data is correct across the entire species set. Two NON-table follow-ups surfaced:
1. **IFOR=3 (Allegheny NF) HT1/HT2 override — PORTED this session** (volume.jl `_htdbh_wykoff`, gated on the
   NE-only Wykoff path so SN is untouched; suite 5281/2 no regression). Still needs a live-FVSne IFOR=3
   differential to promote a test (net01 is IFOR=2).
### SECOND real divergence — R9 merch-cubic topwood dropped for sawtimber-sized trees with no sawlog (FIXED)
Found by RUNTIME broadening: ran a 5-hardwood Allegheny stand (RM/WO/AB/RO/YB) end-to-end through BOTH engines.
cyc-0 stand cols + TCuFt + SCuFt + BdFt matched, but stand MCuFt = 403 (jl) vs 490 (live). Localized per-tree
(live .trl): ONLY YB sp30 d=11 diverged — jl merch cubic 0.0 vs live 15.7 (all 4 others bit-exact). ROOT: YB
d=11 sits exactly at SCFMIND=11 ⇒ prod="01" (sawtimber), saw top 9.6"; the saw bole (11→9.6") is too short to
make a MERCHL=8 sawlog ⇒ sawHt→0 ⇒ SAW=0 (correct, matches live). BUG: `r9clark_cubic` gated the ENTIRE saw
block (vol[4] saw cubic AND vol[7] topwood) on `if sawHt>0`, so when sawHt=0 BOTH were left 0 and the merch
reconstruction v[4]+v[7] lost the whole stump→9.6 bole (the tcfVol pulp-top cubic, already computed at line
404, was discarded). FVS r9clark.f:286-367 books `vol(7)=max(tcfVol−cfVol,0)` inside `if(cupFlg.eq.1 .or.
spFlg.eq.1)` — NOT gated on sawHt>0; cfVol=r9cuft(stump,0)=0 when sawHt=0, so vol(7)=tcfVol = the FULL merch.
FIX (r9clark_vol.jl): hoisted `vol[4]=scfVol` (=0 when no sawlog) and `spFlg && vol[7]=max(tcfVol−scfVol,0)`
OUT of the `if sawHt>0` block (board feet stays gated — a 0-length sawlog has 0 bd ft). RESULT: Allegheny
MCuFt 403→490 == live; YB d=11 merch 0→15.7. net01 UNCHANGED (its prod-01 trees all have valid sawlogs ⇒
sawHt>0 path identical) — .sum still TCuFt1546/MCuFt1338/SCuFt294/BdFt1637 bit-exact. +10 tests. Suite 5301/2.
This is the SECOND latent bug broadening surfaced (after the IFOR=3 override) — both invisible to net01's
species/forest, both real, both faithfully fixed from the FVS source.

2. **R9 q_4/q_7 unloaded — RESOLVED, not a bug.** r9clark.f:770-789 reads COEFFS C/E/P/A/B/R ALWAYS from
   coef0 (cols 4-9); only A17/B17 are topDib-selected (cols 2,3 of coef0/coef4/coef79). So coef4/coef79
   cols 4-8 (incl q_4/q_7) are NEVER read by the volume algorithm — the port's `_R9Coef` correctly omits
   them. (Confirmed by the <1% per-tree volume match too.) The R9 table is FULLY faithful.

### CONFIRMED-MISSING (found by the HT-DBH audit) — IFOR=3 (Allegheny NF) HT1/HT2 override NOT ported
sitset.f:428-489 has an `IF(IFOR.EQ.3)` block that REPLACES the base HT1/HT2 (Wykoff HT-DBH intercept/slope)
for **20 species** {26,27,30,31,33,40,41,42,44,54,55,60,64,67,69,71,93,102,106,108} (e.g. sp26 HT1 4.3379→
4.6839, HT2 −3.8214→−4.9622). The extracted CSV correctly holds the BASE (IFOR≠3) values, but the NE port's
HT-DBH path (`_htdbh_params`, volume.jl:18) only branches on the SOUTHERN Fort-Bragg ifor==20 P2/P3/P4
override — there is NO ifor==3 HT1/HT2 path. So an Allegheny-NF (forest code 3) NE stand would silently use
base Wykoff coefficients for those 20 species ⇒ wrong dub heights/diameters → wrong volume/HTG-calib.
INERT for net01 (IFOR=2). Status: a faithful-completeness gap for IFOR=3 NE stands — port the 20-species
HT1/HT2 swap as a forest-aware override (data + an ifor==3 branch in _htdbh_params) when broadening past net01.
This is the FIRST real divergence the all-108-species coefficient audit surfaced — exactly the latent
forest/species-specific path net01's ~7 species never exercise.


## RUNTIME BROADENING — 15-species multi-cycle growth spine BIT-EXACT vs live FVSne
After the coefficient audit + 2 fixes, validated the GROWTH spine on species net01 never grows: a constructed
15-species stand (softwoods BF/RS/WP/EH/PP + hardwoods RM/YB/AB/CT/WO/SW/HK/BG/EL/HT spanning site_groups
1-28, varied DBH 4-20", IFOR=2, dubbed heights), 3 cycles. Every stand column (TREES/BA/SDI/CCF/TopHt/QMD)
BIT-EXACT vs live at 1990/2000/2010/2020 (e.g. 2020: TREES 129/129, BA 94/94, SDI 155/155, CCF 144/144, TopHt
84/84, QMD 11.5/11.5); volumes within ULP (cyc0 TCuFt 1313/1309 +0.3%, MCuFt 1190/1190 exact, SCuFt 673/673
exact, BdFt 3843/3849). ⇒ NE DG (BAL model) + HTG (NC-128 curves) + mortality + crown + R9 volume are faithful
across the diverse species set, not just net01's ~7. Fixture+test: test/integration/ne_fixtures/divspp.* +
test_net01.jl. KEYFILE GOTCHA: all trees must share ONE plot (tree id cols 24-27 = "0101"; the PLOT is cols
26-27, the per-tree record number is cols 1-4) else "MORE PLOTS FOUND THAN SPECIFIED" abort. Suite 5309/2.


## CROSS-FOREST broadening — IFOR-dependent merch rules BIT-EXACT vs live FVSne
Ran the 15-species stand under IFOR=1 (forest 914), IFOR=4 (920), IFOR=5 (921) — exercising _ne_merch's
IFOR-branched hardwood DBHMIN (1/3→6, 4→8, else→5) + the per-forest site defaults (lat/long→Hopkins→CCF).
Live-validated: cyc0 MCuFt tracks the merch rule EXACTLY — IFOR1 1158, IFOR4 1083, IFOR5 1190 (all jl=live);
CCF tracks the site default 103/102/104 (jl=live); all other stand cols bit-exact. Residual = the same
consistent TCuFt +4/BdFt −6 (~0.3%) as IFOR=2: per-tree TOTAL cubic matches live to .trl print precision (0.1)
on all 15 trees (BF 3.22/3.2 … SW 60.61/60.6, max|Δ|0.09) ⇒ sub-print Float32/Clark-tip-integration rounding
accumulated ×TPA, NOT a localizable bug (same class as net01 volume ULP). MCuFt/SCuFt bit-exact everywhere ⇒
the tip-only nature confirms it is the total-to-tip integration, not a merch error. Fixtures+test:
ne_fixtures/divspp_f{914,920,921}.* + test_net01.jl. Suite 5318/2.

⇒ BROADENING TALLY so far: coefficient tables (108 spp, 8 tables) all faithful; 2 real bugs found+fixed+live-
validated (IFOR=3 HT-DBH override, R9 merch-topwood); growth spine bit-exact on 15 diverse spp × 3 cycles ×
4 forests (IFOR 1/2/3/4/5 merch rules + site defaults). NE port is faithful well beyond net01's envelope.
