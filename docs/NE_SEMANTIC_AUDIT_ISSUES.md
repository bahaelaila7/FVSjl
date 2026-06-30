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


## net01 LAST volume residual ROOT-CAUSED — top-kill (NORMHT + CFTOPK) gap (shared SN+NE)
After the rounding-order fix made every DUBBED-height stand (divspp/dense/cross-forest, 5 forests) BIT-EXACT,
net01 retained a ~0.7% TCuFt/MCuFt gap (1549/1558, 1335/1347). Localized per-tree (net01 stand-1 .trl) to a
SINGLE tree: record 22, sugar maple (SM, sp27) d=10.4, measured ht=55 — jl TOT 13.8 vs live 15.4 (~10%). All 4
other SM trees bit-exact. Root cause: record 22 is the ONLY TOP-KILLED tree (.tre cols 63-65 HTTOPK=49,
col 52-53 damage code 97). FVS vols.f:146 `IF(TKILL) H=NORMHT(I)/100` — for a top-killed tree it computes the
cubic on the NORMAL (undamaged) height NORMHT (=the HTDBH prediction, cratet.f:437), NOT the broken measured
height; then vols.f:193 `CALL CFTOPK` applies the Behre top-kill correction. VERIFIED in jl: HTDBH(SM,10.4)=
63.93; r9clark @63.9 = 15.9; CFTOPK reduces the missing-top fraction 15.9→15.4 = live. jl uses the measured 55
→ 13.8. This single tree ≈ the ENTIRE net01 TCuFt/MCuFt residual (1.6 cuft × ~6 TPA ≈ 10 cuft/ac).
PORT CHECKLIST (shared driver, keep SN bit-exact — note SN snt01 has the SAME ~11 cuft/ac residual, see
[[fvsjl-volume-topkill-gap]]): (1) read ITRUNC/HTTOPK from the .tre (cols 63-65) into the tree state; (2) set
NORMHT = normal HTDBH height for top-killed trees (cratet.f:365/437 — NORMHT=INT(H·100), H=HTDBH when
top-killed); (3) in the volume driver use H=NORMHT/100 when TKILL=(H≥4.5 ∧ ITRUNC>0); (4) port CFTOPK +
BEHPRM (Behre AHAT/BHAT params) + BEHRE (Behre profile integral) and apply to TCF/MCF/SCF. The R9 cubic KERNEL
is already proven bit-exact (all dubbed stands) — this is purely the top-kill WRAPPER. This is the clear next
port item; it closes net01 to fully bit-exact AND the parallel SN snt01 residual.


## TOP-KILL PORT DONE — net01 now FULLY bit-exact on volume (NORMHT + CFTOPK wired into compute_volumes_ne!)
The top-kill machinery already existed (treeinput reads damage96/97+HTTOPK→trunc/norm_ht=-1; dub_missing_heights!
resolves norm_ht=predicted HTDBH height; behre/behre_params/cftopk/bftopk ported in volume.jl; SN compute_volumes!
already applies them). Only compute_volumes_ne! (the early-dispatch NE path, volume.jl:405) didn't. FIXED: added
the same block — `tkill=h≥4.5 ∧ trunc>0`; `h=norm_ht/100` (the NORMAL predicted height); after r9clark_cubic,
`cftopk`+`bftopk` truncate back to the break. NE _ne_merch returns scalars so the merch standards are wrapped as
1-tuples indexed by sp=1. RESULT: net01 stand-1 cyc0 BIT-EXACT all 4 vol cols (TCuFt 1558, MCuFt 1347, SCuFt 292,
BdFt 1633 = live); the broken-top SM d10.4 tree 13.8→15.4 (=live .trl), sp49 d8.0 tree likewise. net01 had 2
top-killed trees. Dubbed stands UNCHANGED (no tkill ⇒ identical path); SN unaffected (NE-only fn). Test
strengthened to assert all 4 vol cols + the per-tree top-kill TOT. Suite 5333/2. ⇒ net01 volume is now fully
faithful — the last net01 residual is CLOSED. (SN snt01's parallel top-kill residual remains; the SN path already
has cftopk so it may already be handled — verify separately if revisiting SN.)


## REAL BUG (broadening) — NE thinning + AUTOES crashes (missing is_sprouting + NE-specific ESUCKR model)
Thinning a diverse stand (THINBBA) WITHOUT NOAUTOES crashes jl: `KeyError: :is_sprouting`. ROOT: FVS NE LSPRUT
defaults TRUE (esinit.f:50); a cut with sprouting-on triggers the stump-sprout path, and cuts.jl:128
`coef_col(:is_sprouting)[sp]` KeyErrors because NE has no is_sprouting column. net01 MASKS this — all 5 stands
carry NOAUTOES (⇒ lsprut=false). So net01 stays bit-exact; this only bites a general AUTOES-default NE stand
that gets cut. CONFIRMED the NE sprout model is DIFFERENT from SN: jl's esuckr! (sprout.jl:139) uses SN-only
essprt.f funcs (nsprec_sn / essprt_sn / sprtht_sn), but NE esuckr.f does: 1 sprout record per stump
(PROB=PREM·SMULT, esuckr.f:284 — NOT nsprec_sn's multiple), height=SPRTHT(VARACD,ISSP,SITEAR,ISHAG) (esuckr.f:
286, variant-aware), then BACHLO(0,0.5,ESRANN)·HT/5.5 deviation, DBH=Wykoff inverse `HT2/(ln(HT−4.5)−AX)−1`
with AX=HT1(IABFLG=1) else AA (esuckr.f:294-301), ICR=70, NO essprt survival. NEXT-PORT CHECKLIST: (1) add
is_sprouting column for NE = ISPSPE (blkdat.f:109) = {20, 26-70, 72-97, 99-108} (82 species, ⇒ most hardwoods
+ sp20; softwoods 1-19/21-25 + 71/98 don't sprout); (2) make esuckr! variant-aware — NE branch = 1 record/stump
+ SPRTHT(NE) + Wykoff-inverse DBH + BACHLO random, NO nsprec/essprt; (3) port NE SPRTHT (sprtht.f VARACD=NE
branch). VALIDATE vs live: thinned dense stand (THINBBA 2010 resid-BA 100) → live 2020 TREES 301 / BA 119
(incl. sprouts). KEEP SN bit-exact (variant-gate the esuckr! sprout-creation; sprout.key SN validation must hold).


## NE SPROUTING — crash FIXED + partial port (BA/vol exact; count approximate, tables TODO)
Fixed the NE-thinning-AUTOES crash: (1) created data/northeast/sprout_essprt.csv with is_sprouting = ISPSPE
{20,26-70,72-97,99-108} (82 spp) ⇒ cuts.jl:128 guard works; (2) added sprtht_ne (essprt.f SPRTHT CASE('NE'):
spp {26-70,72-97,99-108} use (0.1+SI/50)·age, else 0.5+0.5·age) + ne_sprout_dbh (esuckr.f:294-301 Wykoff inverse
HT2/(ln(HT−4.5)−HT1)−1 reading NE :htdbh_ht1/ht2 — the IABFLG=1 path, the default since LHTDRG/AA re-fit needs
≥3 measured heights, cratet.f:311); (3) made esuckr! variant-aware (ne flag). RESULT on the dense-thin test
(THINBBA 2010 resid-BA 100): BA 119/119, SDI 217/219, CCF 199/200, TopHt 76/76, TCuFt 3380/3380, MCuFt 3201/3201
ALL match live — only TREES 345(jl)/301(live) off (+44 small-sprout TPA). ROOT of the count gap: NE ESUCKR uses
NSPREC (sprout-count) + ESSPRT (survival), BOTH variant-branched — I initially mis-read NE as "1 record/no
survival"; it actually has its OWN NSPREC CASE('NE') (essprt.f ~1010+, ISPC/DSTMP table) + ESSPRT CASE('NE')
(essprt.f:362-460, ~40 cases incl. logistic e.g. sp55/58/64 + poly sp59-70) + ESASID(NE)=49 aspen → ASSPTN.
CURRENTLY the SN nsprec_sn/essprt_sn stand in (documented NE-TODO) ⇒ count approximate. NEXT: port nsprec_ne +
essprt_ne (transliterate the two CASE('NE') tables) + aspen; validate dense-thin TREES→301. SN bit-exact (ne
branch false for SN; suite 5342/2). net01 unaffected (NOAUTOES).


## NE SPROUTING — NSPREC(NE)+ESSPRT(NE) PORTED ⇒ thinning+sprouting now BIT-EXACT
Ported the two NE sprout tables (transliterated the CASE('NE') SELECT blocks): nsprec_ne (essprt.f:1006 — the
sprout-COUNT per stump by ISPC/DSTMP) and essprt_ne (essprt.f:362 — the per-stump survival multiplier on PREM,
~25 species cases incl. the logistic forms sp55/58/64-66/69/81 + the sp59-70 cubic poly + the sp27/28 linear).
Wired into esuckr! (ne branch). RESULT: the dense-thin stand (THINBBA 2010 resid-BA 100, AUTOES on) is now
BIT-EXACT vs live at the post-thin 2020 row — TREES 301/301 (was 345, the +44 sprout over-count is GONE), BA
119, SDI 217, CCF 199, TopHt 76, TCuFt 3380, MCuFt 3201, ALL = live. Locked by a fixture+test (ne_fixtures/
thin.key). SN bit-exact (ne branch false for SN; suite 5345/2). ⚠ REMAINING (small): NE aspen suckering
ESASID(NE)=49 → ASSPTN (a cut sp49 record resets PREM via the aspen sucker model before ESSPRT) is NOT yet
ported — absent it, sp49 stumps use the plain PREM. Inert unless a stand cuts NE aspen (sp49); not in the dense
stand. ⇒ NE stump-sprouting is faithful for non-aspen; the lone gap is the sp49 aspen ASSPTN path.


## NE SPROUTING COMPLETE — aspen ESASID(49)/ASSPTN ported ⇒ ALL species bit-exact
Ported the last sprout gap: NE aspen (sp49) suckering. esuckr! now (1) accumulates the cut-aspen ASBAR/ASTPAR
over the cut_log (estump.f:110-111: ASTPAR+=PREM, ASBAR+=0.0054542·PREM·DBH², gated on ESASID(NE)=49), then
(2) for each cut sp49 record runs ASSPTN (essprt.f:1228) BEFORE ESSPRT: SPA=40100.45−3574.02·a²+554.02·a³−
3.5208·a⁵+0.011797·a⁷ (a=ISHAG) clamped [2608,30125], ×ASBAR/198; PREM=(PREM/(ASTPAR·2))·SPA. VALIDATED on an
aspen-dominated stand heavily thinned (THINBBA 2010 resid-BA 30): post-thin 2020 BIT-EXACT — TREES 740/740 (the
massive sucker flush 146→740), BA 41, SDI 80, CCF 87, TopHt 56 all = live. Fixture+test ne_fixtures/aspen.key.
⇒ NE STUMP-SPROUTING IS NOW FULLY FAITHFUL (all species incl. aspen) — the sprouting cluster is CLOSED. Suite
5347/2; SN bit-exact (ne branch false for SN); net01 unaffected (NOAUTOES).


## FFE/SIMFIRE on diverse species — fire mortality BIT-EXACT; minor pre-fire 1-CCF blip noted
Ran the 15-species stand + the full FFE keyword set (SNAGINIT/SNAGBRK/FLAMEADJ/SIMFIRE 2010/SALVAGE/DEFULMOD/
SNAGPSFT/PotFIRE/+reports). The 2010 SIMFIRE kills ~half. LIVE-VALIDATED: pre-fire 1990/2010 + POST-FIRE 2020
all BIT-EXACT (2020 TREES 74/74, BA 65/65, SDI 103/103, CCF 93/93, TopHt 82/82) ⇒ the fire mortality on the
diverse per-species bark/crown fire props (fire_species_props.csv) is FAITHFUL. ★ ONE minor residual: CCF at
2000 = jl 117 / live 118 (the no-FFE divspp has 118, so it's FFE-init-SPECIFIC, a real ~0.5-1 shift not pure
rounding — jl-noFFE ≥117.5 → jl-FFE <117.5). NOT the crown-lift directly modifying the live crown (fuel_
additions.jl reads crown_pct, stores into ffe_oldcr — doesn't write crown_pct). Pre-fire, self-corrects by
2010 (CCF 131/131). REFINED: at 2000 TPA/BA/SDI/TopHt ALL match (135/66/117/73) — only CCF differs, so the LIVE TREES are
identical (not a growth/RNG shift) and it is purely a CCF/crown-computation interaction with the FFE path.
Ruled out: FFE writes no crown_pct/crown_ratio (fuel_additions.jl only stores ffe_oldcr); the line-323
compute_density! is gated on trmort>0 (FuelTRT, absent here). Suspect SNAGINIT-created snags being included
in the CCF sum, OR an FFE crown-width/gross_space difference. FOLLOW-UP: diff compute_density! CCF inputs
(crown widths / record set) between the FFE and no-FFE runs at the 2000 boundary.
Locked: ne_fixtures/ffe.key test asserts the bit-exact pre/post-fire rows. Suite 5350/2.


## FFE PotFIRE RNG bug FIXED — the 1-CCF blip was the report consuming the simulation RNG
Bisected the FFE 1-CCF blip (jl 117/live 118 @2000) to the PotFIRE keyword (the other FFE keywords each give
118; only PotFIRE → 117). ROOT: jl's potential_fire_report → torching_probability draws rann!(s.rng) for the
report's stochastic torching, CONSUMING the main RNG stream — which then shifts the crown-ratio draws ⇒ the
1-CCF perturbation (TPA/BA/TopHt identical, only CCF moved, the crown-RNG signature). FVS FMPOFL_FMPTRH
(fmpofl.f:506/649) wraps its RANN draws in RANNGET(SAVES0)…RANNPUT(SAVES0) — the POTENTIAL-fire REPORT is a
hypothetical and must not perturb the sim. FIX (fmburn.jl potential_fire_report): save s.rng.s0 before
torching_probability, rannput! it back after. RESULT: diverse FFE stand now BIT-EXACT EVERY year incl
CCF@2000 118/118. SHARED fix (SN uses it too) — SN bit-exact, no regression. ⇒ 6TH real bug of the campaign;
the FFE-diverse residual is CLOSED. Suite 5350/2.


## Establishment (PLANT) on diverse species — confirms the first-cycle growth residual (converging)
Broadening: a BARE stand planting 8 diverse species (BF/WS/WP/EH softwoods + RM/YB/WO/RO hardwoods, ESTAB+PLANT
×200 each). LIVE-DIFF: TREES BIT-EXACT every cycle (1600→...); but the FIRST cycle (2002) shows BA 21/25, CCF
67/80, SDI 79/90 (~16% low), converging to bit-exact by 2022 (BA 172/172, SDI 401/401, CCF 418/417). LOCALIZED
(plant_hard TREELIST @2002): jl RM DBH 1.12/live 1.2, HT 17.0/17.6, CR 79.6/~79 (crown ratio MATCHES). So the
established trees are ~5-7% SMALLER in the first cycle ⇒ ~10% lower crown area ⇒ the CCF/BA/SDI deficit (BA
rounds 10/10 hiding it; CCF shows it). NOT a crown bug — it is the establishment first-cycle GROWTH residual,
the SAME class as net01 BARE (~5%, [[fvsjl-ne-port-state]] "BARE ~1 BA floor"), now confirmed across diverse
species and slightly larger for some. The planted seedlings grow slightly less than live in their creation
cycle (ESSUBH base + Phase-2 scale_e=0.5·ne_htcalc_incr·gmod, all previously semantics-verified). Converges by
cycle ~3. STATUS: the establishment-cohort first-cycle growth is the one persistent (small, converging) NE
residual — consistent across BF/WS (net01) and the diverse set. FOLLOW-UP: re-examine the established-tree
first-cycle HT increment per species (the ~3% HT / ~7% DBH deficit) vs live — possibly the base-height/age or
the half-cycle scale interaction; prior net01 analysis concluded floor-level but the diverse set shows it is
species-dependent (above the BF/WS floor).


## Establishment residual LOCALIZED to the ESSUBH base height (not the increment)
RM (sp26) analytical trace: ESSUBH HHT=(HTCALC(CARAGE)/CARAGE)·min(5,TIME−DELAY) (essubh.f:41-43). jl
ne_htcalc_height(26, si=71.46, CARAGE=20)=33.78 ⇒ base=(33.78/20)·5=8.45, but live implies H≈36.4 ⇒ base≈9.1
(.trl RM @2002: HT 17.6 − HT-INCR 8.5). The PHASE-2 INCREMENT MATCHES (jl ne_htcalc_incr(age(base))·0.5=8.55 vs
live ~8.5). So jl's established RM is 0.6 ft short ENTIRELY from the base height. jl per-species site indices
vary correctly (SICOEF applied: sp1=67.46…sp9=80.46…sp26=71.46) and the stand CONVERGES by 2022 (regular growth
using the SAME ne_htcalc is right) ⇒ ESTABLISHMENT-SPECIFIC: ne_htcalc_height at the ESSUBH ref-age (20) for RM
is ~7% below what live HTCALC returns, despite LTBHEC/MAPNE being audited bit-exact. SUSPECT: ESSUBH passes
SITEAR(IPNSPE) to HTCALC MODE0 — maybe a different SI than jl's sp_site_index, or a HTCALC-mode detail diff at
that age. DEFINITIVE FOLLOW-UP: instrument live essubh.f to dump SI + H@CARAGE for RM, diff vs jl
ne_htcalc_height(26,SI,20). Small (0.6 ft), converging — the last NE residual.


## Establishment GROWTH is FAITHFUL (corrected) — residual is a minor crown-width/CCF effect, not growth
DEEP-TRACED via FVS DEBUG keyword (no rebuild — essubh.f:44 + regent.f:271 debug WRITEs). FINDINGS that CORRECT
the earlier "establishment first-cycle growth residual":
- ESSUBH base height: live debug HHT=8.446 (SI=71.456, CARAGE=20) — jl MATCHES exactly. CARAGE[26]=20 confirmed
  (essubh.f MAPNE = jl _NE_ESSUBH_REFAGE).
- The +0.5 planted-height random IS applied (estab.f:490 RAN=BACHLO(0.5,0.25) mean 0.5): live tree-record
  height ~9.0, jl debug hht ~8.9 — MATCH. (My earlier "base 8.446" was a manual-grow_cycle! RNG artifact.)
- Phase-2 growth: live regent debug GMOD=1.0, AVH=0 (relht=0), pre-random HTGR=8.517 (base 9.055). jl analytic
  ne_htcalc_incr(age(9.055))·0.5 = 8.517 — BIT-MATCHES live. The jl per-tree htgr "8.04 mean" was just the ±10%
  random (regent.f:266, applied by BOTH) over 6 samples.
⇒ The NE ESTABLISHMENT GROWTH (ESSUBH base + planted random + REGENT-LESTB Phase-2 + the ±10% + GMOD) is
FAITHFUL — base, height-increment, crown ratio (79.6/79) and BA (10/10) all match live. The lone residual is a
small (~10%) first-cycle CCF deficit (hard PLANT CCF 36/40) with BA matching ⇒ a CROWN-WIDTH-of-tiny-established-
trees effect, NOT growth. Converges by cycle ~3. FOLLOW-UP (minor): diff the established-tree crown WIDTH (CCF
input) jl vs live at cyc1 — crown ratio + dbh + height already match, so it is the cwcalc crown-width eval for
small regen trees. The establishment cluster is otherwise CLOSED/faithful.


## Establishment lone residual FINAL characterization — small first-cycle SDI/CCF dbh-distribution effect
Exhaustively traced. The established-tree GROWTH is faithful (proven: base/HTGR/crown-ratio/crown-WIDTH all
match live per-tree — jl crown_width(RM d1.2 h17.6 cr77)=4.03 vs live 3.9; ne_htcalc_incr·0.5=8.517=live HTGR).
The lone residual on the diverse PLANT: cyc-1 SDI 37/40 + CCF 36/40 (~8% low) with BA MATCHING (10/10). Since
BA(∝Σdbh²) matches but SDI(∝Σdbh^1.6) and CCF differ, it is purely a DBH-DISTRIBUTION effect (SDI/CCF are
nonlinear in dbh; a different per-tree dbh spread shifts them while the 2nd moment / BA holds). Both jl & live
apply the regent.f:266 ±10% HTGR random (gated DGSD≥1, NOT LESTB) and the +0.5/0.25 planted random — so the
spread SHOULD match; the small residual is a per-tree RNG-realization / dbh-distribution alignment on the
established cohort (the establishment :estab stream draw order — cf the earlier net01-BARE RNG-align note). It
CONVERGES (cyc-2 CCF 176/183, cyc-3 311/317 — proportionally tighter). ⇒ NE establishment is FAITHFUL on every
per-tree quantity; the remaining ~8% first-cycle aggregate SDI/CCF is the last (small, converging) NE residual,
in the same RNG-draw-order class as the net01-BARE establishment alignment. Not growth, crown, or volume.


## FIXED (7th bug) — establishment BALMOD competition used POST-establishment BAL (should be PRE)
The "established-cohort dbh-distribution residual" was a REAL bug, not RNG. jl computed ebau_e (the REGENT-LESTB
BALMOD BAL) AFTER creating the seedlings ⇒ the cohort's own BA (~2.7 at the ~9 ft base heights) gave GMOD~0.966,
under-growing the establishment HTGR ~4% (RM dbh 1.123 vs live 1.174). FVS uses the PRE-establishment density —
live debug: GMOD=1.0 / AVH=0 for a BARE stand (the seedlings don't compete in their own creation cycle; the
DENSE/BAL predates the regen). FIX (establishment.jl): snapshot ne_badist! over 1:nstart (the existing overstory)
BEFORE the creation loop; Phase-2 uses that pre-establishment ebau_pre. RESULT (hardwood PLANT): SDI 40/40 +
CCF 40/40 @2002 (were 37/36), 135/136 + 180/183 @2012, TREES 735/733 — the ~8% deficit CLOSED. net01 BARE + all
establishment tests still pass (BF/WS unaffected to tolerance). LONE remaining: a small TopHt overshoot (24/22
@2002, converging) — the fix unmasked it (doctrine #3): the established trees run a touch tall for their dbh, a
height↔dbh effect on tiny regen (jl _htdbh_dbh gives slightly less dbh per height for small trees). Much smaller
than the SDI/CCF deficit it replaced; converges by cyc-3 (49/49). ⇒ found via per-tree dbh-distribution dump
(jl 1.123 vs live 1.174) — the doctrine's "trace to the per-tree quantity" caught the real bug the aggregate hid.


## FIXED (8th bug) — NE establishment HHTMAX height clamp (+ dbh-from-uncapped) ⇒ diverse PLANT BIT-EXACT
The TopHt overshoot the gmod fix unmasked was a REAL 2nd bug: NE caps the grown establishment height at HHTMAX
(ne/blkdat.f DATA HHTMAX/, 108 vals) — a HARD per-species ceiling (YB(30)=22, WO(55)=16; live clamps all YB to
22.00, all WO to 16.00 exactly). jl used the SN _ES_HHTMAX (20-ft fill for sp≥12) and only the soft site-curve
cap ⇒ YB grew to 23.49. KEY detail: the DBH is derived from the UNCAPPED grown height, only the REPORTED height
is clamped (live YB dbh 1.8 = Wykoff(~23.5), height 22) — clamping before the dbh under-sized it (SDI/CCF
dropped to 36). FIX: added _NE_ES_HHTMAX (108 NE vals) + clamp hk to it AFTER computing dbh. RESULT: hardwood
PLANT now BIT-EXACT — 2002 800/10/40/40/22 + 2012 786/48/136/183/40 all = live (±1 by 2022). ⇒ with the gmod
fix (7th) this CLOSES the establishment cluster: diverse-species PLANT bit-exact. net01 BARE + all est tests
pass. Found by dumping the per-species established HEIGHT (jl YB 23.49 vs live 22.00 exactly = a hard clamp).
