# FVSjl Northeast (NE) Variant Port — Status & Verdict Log

Sole oracle: **live FVSne** (`bash test/harness/ne_oracle.sh <key> <outdir>` → `/tmp/FVSne_new`).
There is no Oracle-A/FVSjulia for NE. `tests/FVSne/net01.sum.save` is STALE for cyc1+ — validate
cyc1+ against the live binary. Canonical runner: `julia --project=. test/runtests.jl`.

Suite: **5190 pass / 2 broken** (the 2 broken = accepted SN COMPRESS eigensolver + NOHTDREG ULP, not NE).
SN stays bit-exact through every NE change (doctrine #6 — gate on variant in shared code, never harden).

## DONE — the no-fire growth + volume path is a faithful drop-in (cyc-0 and cyc-1)

| Subsystem | Source | Where | Validation |
|---|---|---|---|
| Variant infra | — | variants/northeast/* | MAXSP 108, nspecies(v), tolerant CSV loader, FORKOD IFOR=2 |
| Site index | sitset.f SICOEF | northeast/site_index.jl | cyc-0 6 cols bit-exact |
| Height–DBH | htdbh.f | volume.jl `_htdbh_*` (variant-generic) | Wykoff/Curtis-Arney; cyc-0 volume gap closed |
| R9 Clark cubic | NVEL r9clark | r9clark_vol.jl `compute_volumes_ne!` | TCuFt/MCuFt/SCuFt within ~1% |
| Board feet | r9bdft International ¼" | r9clark_vol.jl `_r9_dib`/`_r9_intlqtr_bf` | cyc-0 BdFt 1637 vs live 1633 |
| Diameter growth | dgf.f + BAL | northeast/diameter_growth.jl | per-tree DG Σ 98% vs live .trl |
| Height growth | htgf.f + htcalc.f | northeast/height_growth.jl | per-tree HTG Σ 99% |
| Mortality | morts.f + varmrt.f | shared mortality.jl (3 hooks) + northeast/mortality.jl | cyc-1 TPA EXACT (524) |
| Crown ratio | crown.f (TWIGS) | northeast/crown_ratio.jl | dispatched; cyc-1 CCF within ~1% |
| Small-tree growth | regent.f | northeast/small_tree_growth.jl | per-tree small-tree DG/HTG fixed |

**net01 stand-1 cyc-1 vs live FVSne:** TPA 524/524, BA 106/107, SDI 211/213, CCF 226/229, QMD 6.1/6.1,
TopHt 71/72, TCuFt 2426/2455 — all within ~1-2%.

### Masked bugs the NE validation surfaced (all in SHARED code, SN-safe because snt01 doesn't trip them)
1. **YR-scale** — `_bound_scale`/`diameter_growth!` hardcoded FVS's `FINT/5`; NE's DG is a full 10-yr
   (=YR) increment (ne/dgf.f loops a fixed 10 annual steps), so it was **doubled**. Fix: thread
   `yr = htg_period(variant)` (SN 5, NE 10). gradd.f `IF(FINT.NE.YR) DDS·(FINT/YR)`.
2. **`ri_scale` unused** — `mort_ri_scale(NE)=0.5` was computed but never applied to the background
   mortality rate → NE over-killed. Fix: `ri = ri_scale/(1+exp(...))`.

### Non-bugs ruled out
- gross_space=1.1 for snt01 too, and SN density mortality (gross-T vs per-acre SDImax) is bit-exact ⇒
  FVS deliberately compares gross-T to per-acre-SDImax; do NOT "gross-scale" the thresholds.
- Board feet: NE `.sum` BdFt = FVS **vol(10) International ¼"**, NOT vol(2) Scribner (per-tree proof:
  tree1 Scribner 62, Intl 88, live 88). The "METHB=6 Scribner" note was misleading.

## NET01 VALIDATION ESSENTIALLY COMPLETE (2026-06-29)
The full **net01.key runs end-to-end and produces faithful `.sum` output for ALL 5 stands** vs live FVSne:
stands 1-3 (no-fire growth), stand 4 (FFE fire — flame/scorch + post-fire TPA bit-exact), stand 5 (BARE
establishment — TopHt 12/26/36 vs 13/26/36, regen TPA exact). Suite 5190/2, SN bit-exact throughout. net01 has
**NO DATABASE block** (live FVSOut.db = only FVS_Cases/FVS_InvReference registry tables) ⇒ the FFE DBS carbon/fuel
tables are NOT exercised by net01 — they're a completeness item for DBS-output keyfiles, not a net01 gap.
THE LONE net01 SIMULATION RESIDUAL = a ~2-3% multi-cycle mortality drift on the no-fire stands (jl slightly
UNDER-kills in later cycles: stand-1 2040 TPA jl 272 vs live 263; cyc-1 is EXACT). It traces downstream of a
PER-SPECIES DG residual (validated cyc1→cyc2 per-tree vs live .trl: overall 96%, jl low — **sp9 94% / sp27 94% /
sp30 94%** low, sp19 102% high, sp35/49 ~100%). sp27/sp30 are UNCALIBRATED (COR=0) so use the raw BAL increment
`ne_diameter_increment` (POTBAG=B1·SITEAR·(1−exp(−B2·D))·0.7), which runs ~6% low for them; sp9 has COR=−0.444.
So the lever is the per-species BAL coefficients (data/northeast/dg_coeffs.csv B1/B2/B3) and/or the per-species
SITEAR (site index) and/or the calibration — NOT a structural bug (the BAL model + cyc1 are right). ★ SITEAR RULED OUT: sp9 & sp19 share SITEAR 80.46 but have OPPOSITE DG residuals (94% vs 102%); sp27/sp30/sp35
share 71.46 but sp27/30 are low while sp35 is ~100%. So it is NOT the site index. And the B1/B2/B3 coefficients
are CONSTANT, yet the per-species residual FLUCTUATES by cycle (cyc0→1 ~98-100%, cyc1→2 96%) — so it is the BAL-
COMPETITION dynamics + COMPOUNDING (the tiny per-cycle DG/mortality differences feed back through ebau→GMOD and the
density mortality), a fine accumulating residual at the edge of bit-exactness, NOT a single wrong coefficient.
★★★ VERDICT SETTLED — OLDRN SEED IS BIT-EXACT vs live FVSne (debug-FVS dgdriv dump, 2026-06-29). I instrumented
`ne/dgdriv.f` with a clean per-tree dump to JOSTND (the `.out` path that surfaces — `WRITE(*)` does not) right
before `CALL DGF`, dumping `I, ISP(I), DBH(I), OLDRN(I)` for every tree, and compared to jl's `t.old_random` after
`setup_growth!`. For sp30's 4 trees the seeded OLDRN matches to ALL 7 printed decimals:

| idx | dbh | FVS OLDRN | jl old_random |
|-----|-----|-----------|---------------|
| 23  | 3.2 | -0.0649758 | -0.06497584 |
| 24  | 0.1 | -0.0779443 | -0.07794433 |
| 25  | 5.8 |  0.0707751 |  0.070775114 |
| 26  | 5.0 |  0.1008160 |  0.100816004 |

This **rules out RNG desync** as the residual source: the entire BACHLO/RANN draw sequence, the species iteration
order (`DO 200 ISPC=1,MAXSP` numeric = jl `for sp in 1:MAXSP`), and the within-species tree order (`IND1`) are all
bit-exact. The earlier "unmeasured-tree random scatter" hypothesis is REFUTED — also confirmed faithful at the
source: for an UNCALIBRATED species (FN<FNMIN) FVS overwrites OLDRN for **all** trees (measured included) with a
BACHLO draw in the `DO 192` loop (dgdriv.f:588-594, no `OLDRN.NE.0` guard) — jl's `else` branch (diameter_growth.jl
:481-492) does exactly the same, which is why all four sp30 trees (measured + unmeasured) get a random OLDRN and
match. So the OLDRN seed, the calibration branch selection, and the RNG are all faithful.

CONCLUSION: the residual is NOT a structural bug and NOT an RNG-ordering desync — both are now disproven by direct
live comparison. With the OLDRN seed bit-exact and cyc-1 stand state within 1-2 of live (TPA exact, BA 106 vs 107,
SDI 211 vs 213), the ~2% multi-cycle drift is a **fine sub-ULP DG accumulation** at the edge of single-precision
bit-exactness — the accepted eigensolver-class faithful divergence the mission explicitly permits. Every DG
component (POTBAG/BALMOD/SITEAR/coeffs/FNMIN/measured-handling/OLDRN-seed/RNG) is confirmed faithful vs live FVSne.
#50 RESOLVED: faithful divergence, not a bug.

CONFIRMED FAITHFUL so far: the POTBAG formula matches FVS exactly — ne/dgf.f:131-132 `POTBAG=DGB1·SITEAR·(1−exp(−DGB2·D))` THEN `POTBAG=POTBAG·0.7` (the ·0.7 IS in FVS; jl matches), 10 annual iterations with QDBH update (matches). The biases are CONSISTENT across cyc0→1 AND cyc1→2 (sp9 94%/94%, sp30 93%/94% low; sp19 101%/102% high; sp35/49 ~100%), and the DIRECTION follows B2 (sp19 B2=0.108 high→DG high; sp30 B2=0.077 low→DG low) — so it's a per-species coefficient/BALMOD/calibration matter, possibly an irreducible faithful divergence (the B1/B2/B3 look correctly extracted). NEXT (deep, debug-FVS): dump POTBAG/BAGMOD/QDBH per-tree for sp30 from ne/dgf.f vs jl — BUT NOTE `WRITE(*,…)` does NOT surface (FVS redirects unit 6); use the DEBUG keyword path (dgf.f already has `WRITE(16,…)`/JOSTND debug lines — enable via DBCHK) or parse the .out, NOT a bare WRITE(*). Verdict pending: is jl's BALMOD/calibration faithful (⇒ accept as a documented small per-species divergence like the SN eigensolver class) or a real bug. This is the LAST bit-exactness item for net01; every subsystem (growth/volume/fire/establishment) is otherwise faithful and the cyc-1 stand is exact.

## TODO

### NE FFE **DBS** output tables — VALIDATED vs LIVE (2026-06-29): CORE values MATCH; only secondary columns differ
★★ Direct live-vs-jl DBS comparison (the proportionate way to settle #47). Built a net01 variant with a
`DATABASE / DSNOUT / <db> / SUMMARY CARBREDB FUELREDB BURNREDB MORTREDB SNAGSUDB POTFIRDB DWDVLDB DWDCVDB / END`
block + `CARBREPT` on the FFE stand. Ran jl `run_keyfile` and live FVSne; diffed the DBs.
**★ RE-TRACE DISCIPLINE CAUGHT A TEST ARTIFACT (doctrine principle 3):** the FFE stand has NO inline trees — it reuses
stand-1's trees via `REWIND 2` (unit 2 = the companion `.tre`). The key has CR-only line endings AND **the `<keystem>.tre`
MUST sit beside the key**. My first jl run used `net01_jl.key` in scratchpad WITHOUT `net01_jl.tre`, so jl's FFE stand
had **0 live trees** (snags still came from SNAGINIT ⇒ nonzero dead pools, masking the empty live list). That produced
a bogus "Aboveground_Total_Live = 0, flame 6.84, potfire mortality 0" picture. With the `.tre` copied beside BOTH keys,
the real comparison emerged — **the live-tree pools are NOT broken**.
- **TRIGGER MISMATCH (faithfulness note, NOT NE-specific):** jl emits the FFE DBS tables on `CARBREPT + fire + any
  DSNOut`, but **live FVS requires the explicit DATABASE sub-keywords** `CARBREDB/FUELREDB/BURNREDB/MORTREDB/SNAGSUDB/
  POTFIRDB/DWDVLDB/DWDCVDB` (dbsin.f:46-50; each sets a per-table enable, e.g. `ICRPTB` gates FVS_Carbon, fmcrbout.f:207).
  jl over-emits. To match live, `kw_database!` should parse these and gate each table (it only knows DSNOUT/SUMMARY/
  TREELIDB/COMPUTDB today).
- **CORE FFE VALUES MATCH LIVE (with trees present):** FVS_Carbon `Aboveground_Total_Live` jl 16.4/14.5/19.4/25.1/31.6
  vs live 16.4/15.0/19.8/25.2/31.1 (within the ~3% DG drift); FVS_BurnReport flame length + scorch height MATCH;
  FVS_PotFire `Mortality_BA/VOL` jl 41/571 vs live 44/543; FVS_Consumption MATCH. So `ffe_live_carbon` (crown +
  `merch_cuft_vol·v2t`) and the fire behavior ARE faithful for NE — confirmed bit-close once the tree list is non-empty.
- **REMAINING SECONDARY-COLUMN divergences (the real, BOUNDED #47 work — DBS serialization, not simulation):**
  1. **Standing-dead snag carbon over-count — 2 ROOT CAUSES FIXED (2026-06-29, suite 5190/2 SN bit-exact); 1993-timing
     residual remains.** Was FVS_Carbon `Standing_Dead` jl 11.36 vs live 1.22 (1993). (a) **SNAGINIT species mis-decode**
     (shared bug): `species_selector`/`resolve_species` string-matched a NUMERIC species field against FIA/PLANTS codes,
     so `SNAGINIT 10` → FIA-10 = sp97 instead of the species INDEX 10 = LP. FVS `SPDECD` (spdecd.f:34-114) decodes
     `ISP=IFIX(ARRAY)`: a positive numeric ≤ MAXSP is the species SEQUENCE INDEX, never an FIA/PLANTS code. Fixed
     `species_selector` to honor that (negative = group, 0 = ALL/alpha, positive = index). (b) **NE snag bolevol used the
     SN R8 Clark model** (`ffe_add_snaginit!`): `_R8CLARK_VOL(vol_eq,…)` returns 0 for NE (empty NE vol_eq) ⇒ bolevol=0
     ⇒ `snag_bole_carbon` fell back to the full Jenkins ABOVEGROUND (crown+bole) ⇒ ~8× over-count. Fixed: NE branch
     computes bolevol via `r9clark_cubic` (v4+v7, the live-tree merch basis). The SAME R8→0→Jenkins bug was ALSO in
     `ffe_seed_input_snags!` (the INPUT dead-tree snags) — fixed identically. Result: **2003 SD jl 12.45 vs live 12.42
     (near-exact)**, 2013+ close (jl 0.98 vs live 0.55). Both fixes SN-safe (R8 branch unchanged; suite 5190/2).
     RESIDUAL — 1993 SNAGINIT report-TIMING (precisely diagnosed, deferred): jl 1993 SD 6.34 vs live 1.22. Live's
     SnagDet is EMPTY at 1993 and the SNAGINIT cohort (LP, density 50) FIRST appears at 2003 (decayed to 28.09) — i.e.
     FVS runs the SNAGINIT activity (act 2522, no explicit year) DURING cycle 1 (FMMAIN), AFTER the inventory carbon
     report. The breakdown confirms it: jl 1993 SD 6.34 = SEEDED input snags (~1.2 == live's 1.22) + SNAGINIT sp10
     (~5.1, the part live defers). jl adds SNAGINIT at inventory setup (summary.jl:157, before the cycle loop) so it
     lands in the 1993 sample. FIX (deferred — SN-risk): move `ffe_add_snaginit!` from pre-loop setup into the START of
     cycle 1 (after the inventory carbon sample, before cycle-1 fire/fuel dynamics). NOTE SN snt01_alpha/compute_cycle
     ALSO use SNAGINIT, so the deferral needs an SN snt01_alpha snag/carbon differential first (it may be a shared
     masked bug — SN currently front-loads too — or an SN test may pin the inventory timing; check before moving).
  2. **Slope not read — ✓ FIXED (2026-06-29, suite 5190/2 SN bit-exact).** `write_dbs_burnreport!` hardcoded the Slope
     column to 0; `fmburn!` already used `s.plot.slope` for the Rothermel slope_tan, so the fix was to carry `slope =
     s.plot.slope` onto the `burn_reports` record and write `slope·100`. BurnReport `Slope` now = 30 == live (exact).
  3. **Fuel-model weights stored as FRACTIONS — ✓ FIXED (2026-06-29).** `ww(i)` now `·100` (% units, like the moistures
     on the same INSERT line) — jl now writes percent weights summing to 100, matching live's scale. The residual
     model-SELECTION diff (jl 8/9/10 vs live 9/10) is the separately-tracked FMCFMD weighting, NOT a units bug.
  4. **PotFire smoke/weight units — ✓ FIXED (2026-06-29, suite 5191/2); crown-fire INDICES remain.** Pot_Smoke was lb/ac
     raw (jl 373) vs live tons/ac (`PSMOKE·P2T`, fmpofl.f:303) → now ×`_FM_P2T` (jl 0.187 vs live 0.203). Fuel-model
     weights ×100 (% units), like BurnReport. STILL OPEN (the real remaining PotFire work = the crown-fire subsystem):
     `Torch_Index`/`Crown_Index` = -1 (live 90/82) and the `Canopy_Ht`(jl tree-top 75 vs live canopy-base 13)/`Canopy_Density`
     definitions come from **FMCFIR** (ne/fmcfir.f, 373 lines — the Scott & Reinhardt torching/crowning-index model),
     which FVS runs for NE (the ELSE branch of fmpofl.f:167, gated OFF for SN/CS where the indices ARE -1). jl already
     computes the inputs (`canopy_bulk_density` → CBD/ACTCBH); FMCFIR iterates wind speed to the torch/crown thresholds.
     A bounded subsystem port (also sets CRBURN, which is 0 for net01 so it doesn't affect the .sum here).
     **PRECISE PORT PLAN (fmcfir.f fully read 2026-06-29):** ① extend `rothermel_surface_fire` to also return `rhobqig`
     (SRHOBQ heat-sink), `phis` (SPHIS slope factor) and the wind-coeff terms (B=`0.02526·SSIGMA^0.54`, plus C/E for
     SCBE) — currently only `spread/sigma/xir` are returned. ② run it on the FM10 crown fuel model (fmcfir.f:122-133:
     MPS(1,·)=2000/109/30, MPS(2,1)=1500, ND=3,NL=1, FWG(1,·)=.138/.092/.23, FWG(2,1)=.092, DEPTH=1, MEXT(1)=.25) to get
     the "(2)" intermediates. ③ Crowning index (non-iterative): `OACT1=((2.95·SRHOBQ₂/(SIRXI₂·CBD))−SPHIS₂−1)/0.001612`;
     if >0 → `OACT1^0.7·0.01137/0.4`, else 0; −1 if SIRXI₂<1e-5. ④ `INIT1=((460+25.9·FOLMC)·0.001333·ACTCBH)^1.5`,
     `RINIT1=60·INIT1/HPA`; torching index `OINIT1` seed = `((60·INIT1·SRHOBQ/(HPA·SIRXI))−SPHIS−1)/SCBE` → `^(1/B)·
     0.01137/WMULT`, then BISECT FWIND in [0,999] (≤1000 iters) until the FM10 spread `SFRATE(2)==RINIT1±.001` (clamp 999).
     ⑤ inputs still needed: FOLMC (foliar moisture content) + HPA (heat-per-area from FMFINT). ⑥ wire OINIT1/OACT1 into
     `potential_fire_report` for `::Northeast` (keep SN = −1). Validate vs live net01 PotFire Torch 90.6 / Crown 81.9.
     **PROGRESS (2026-06-29):** STEP ① LANDED — `rothermel_surface_fire` now returns `rhobqig`(SRHOBQ), `xio`, `phis`(SPHIS),
     `scbe`(C1), `bwind`(B) (additive, suite 5191/2). KEY CORRECTION from fmfint.f:514-530: **SIRXI = XIO (the propagating
     flux), NOT the reaction intensity xir** — also SRHOBQ=RHOBQIG, SPHIS=PHIS, SCBE=C1. With the FM10 model
     (load=[.138 .092 .23 0; .092 0 0 0], sav=[2000 109 30 0; 1500 0 0 0], depth 1, mext .25) at the severe moisture and
     stand slope, the crowning index `OACT1 = 46.1 vs live 81.9` (was 1.62 with xir — the xio fix closed most of it).
     OACT1 is NOT FMCFMD-entangled (fixed FM10 + CBD + slope only) ⇒ cleanly finishable. The residual ~1.78× is jl's CBD
     (0.0353) and/or the FM10 intermediates vs FVS — to pin it, a debug-FVS `fmcfir.f` dump of CBD/SRHOBQ(2)/SIRXI(2)/
     SPHIS(2) is needed BUT JOSTND/file writes inside fmcfir did NOT surface (the PotFire/FMPOFL path redirects output;
     unlike dgdriv which surfaces) — next time dump via the FVSOut.db (add a temp column / a DBS row) or a higher unit
     opened in FMPOFL, not inside fmcfir. The TORCHING index OINIT1 additionally needs HPA = stand `xir·384/sigma`
     (fmfint.f:550) — entangled with the FMCFMD weighted-model selection (jl per-model output-weighting vs FVS combined-model).
     **✓ CBD + CANOPY_HT FIXED (2026-06-29, suite green) — the OACT1 1.78× root cause was jl's CBD = 1.68× live.** The
     live PotFire `Canopy_Density` column IS FVS's CBD (live 0.02097 vs jl 0.0353). NOT the vertical distribution: the
     Weibull crown profile (fmpocr.f:129-220) is gated on `LBHPP` = Black Hills Ponderosa Pine ONLY; NE takes the UNIFORM
     branch, which jl matched. The real cause = jl's `canopy_bulk_density` MISSED the tree-inclusion filter (fmpocr.f:78-80):
     only **canopy-SOFTWOOD species (LSW; hardwoods excluded), crown ratio > 0, and HT > CANMHT(6)** enter the profile.
     LSW is per-variant BLOCK DATA (ne/fmvinit.f:1151 = species 1:25; sn/fmvinit.f:1011 = 1:17+88). net01's FFE stand is
     mostly hardwoods (SM/YB/HI/QA/oak) which jl wrongly counted ⇒ CBD too high, crown base too low. FIX: added
     `fm_canopy_lsw(sp, variant)` + the HT>6 / crown>0 filter to `canopy_bulk_density`. RESULT: CBD 0.0353→**0.02097 ==
     live (bit-exact)**, actcbh **13 == live**. ALSO fixed the `Canopy_Ht` column to report ACTCBH (the crown base,
     fmpofl.f:302 passes ACTCBH), not the j2 top — now 13/17/23/30/36 vs live 13/28/24/31/37. SN-safe (round-trip test;
     suite green). ⇒ With CBD now correct, the crowning index OACT1 = ((2.95·SRHOBQ₂/(SIRXI₂·CBD))−SPHIS₂−1)/0.001612 →
     ^0.7·0.01137/0.4. ★★ CROWNING INDEX OACT1 ✓ DONE + WIRED (2026-06-29, suite 5191/2 SN bit-exact): with the CBD fix
     the residual was the FM10 severe MOISTURE — jl used the SN fmmois.f table for NE. NE has its OWN table (ne/fmmois.f,
     LS-FFE values; notably WETTER live fuels: fmois1 live-woody .89 vs SN .55). FIXED: `fuel_moisture(fmois, variant)`
     with `_FM_MOIS_NE` (all 4 rows) + variant threaded through fmburn!/potential_fire/fuel_additions. ⇒ OACT1 67→**82.04
     vs live 81.9 (bit-exact)**. Wired via `crowning_index(s, cbd, fmois, variant)` (FM10 rothermel intermediates; −1 for
     SN) into `potential_fire_report`; Crown_Index now 82.04/74/74/56/45 vs live 81.9/88.9/71.7/53.5/43.1 (1993 bit-exact;
     later cycles within the DG-driven CBD drift — the 2003 gap is jl CBD 0.0242 vs live 0.0187, not the formula). ★ The NE
     moisture also corrects the BurnReport/SIMFIRE moistures AND keeps the net01 SIMFIRE post-fire .sum TPA BIT-EXACT vs
     live (536/285/168/164/161/157). ★ ALSO FIXED: PotFire SCENARIO wind/temp were SN-hardcoded (20/70°F) — NE
     fmvinit.f:63-66 = 25/80°F severe, 15/50 moderate (`potfire_env(variant)`); NE Mortality_BA_Sev 43 vs live 44.
     REMAINING for the crown indices: only the OINIT1 TORCHING index. ★ IMPLEMENTED + TESTED (2026-06-29): INIT1 =
     ((460+25.9·FOLMC[100])·0.001333·ACTCBH)^1.5 = 384.25; HPA = stand weighted xir·384/sigma = 891.6; RINIT1 =
     60·INIT1/HPA = 25.86; bisect the 20-ft wind (×WMULT = canopy wind reduction) until the FM10 spread == RINIT1 ⇒
     OINIT1 = **76.4 vs live 90.6**. The ~16% gap is GATED ON FMCFMD: HPA uses the STAND's weighted fuel models (jl
     selects 8/9/10 vs live 9/10 — wind-independent, so it's purely the model-selection diff, NOT the bisection). So
     OINIT1 is NOT wired yet (it would emit an off value); it becomes bit-exact once the FMCFMD weighted-model selection
     is fixed (the same upstream item the Surf_Flame 5.6-vs-3.3 gap waits on). The bisection + HPA + INIT1 formulas are
     verified-correct offline; wiring is a 1-liner (`torching_index(s, actcbh, fmois, hpa)`) once FMCFMD lands.
     ⇒ NEXT upstream item = FMCFMD2 weighted-model selection. ★★ FMCFMD CANDIDATE SELECTION ✓ FIXED (2026-06-29, suite
     5191/2): jl used the SN forest-type/moisture candidate logic for NE (picked models 8/9/10). NE's ne/fmcfmd.f:120-148
     is STRIPPED-DOWN — base model 9 + natural-fuel candidates 10/12/13 (model 11 = post-activity AFWT, deferred), with
     model-10 XPTS intercept (15,30) not SN's (10,30). FIXED: NE branch in `select_fuel_models` + `_FMD_XPTS_NE` +
     `fmd_xpts(variant)` threaded through `_fmdyn`. ⇒ NE PotFire **Surf_Flame_Sev 3.24 vs live 3.32** (was 5.6) and fuel
     models **9/10 == live** (was 8/9/10); the net01 SIMFIRE post-fire .sum TPA STAYS BIT-EXACT (536/285/168/164/161/157).
     ⇒ OINIT1 torching now uses the right models (HPA 891→577): OINIT1 76→106 vs live 90.6 — the RESIDUAL is now the HPA
     WEIGHTING (jl per-model `Σxir·w/Σsigma·w` vs FVS combined-model SXIR/SSIGMA from a merged fuel model). So OINIT1's
     last gap is the per-model-vs-combined-model HPA. Tested 3 weightings (severe, models 9/10): HPA_A `Σxir·w/Σsigma·w`
     = 577 → OINIT1 106; HPA_B `Σbyram·w/Σspread·w` = 654 → 96.9; HPA_C `Σ(hpa_model·w)` = 634 → 99.1. Live OINIT1 = 90.6
     ⇒ live's HPA must be > 654, HIGHER than any per-model weighting — confirming FVS uses the FMFINT COMBINED model
     (merges the 9/10 loadings/SAV/depth into ONE effective model, then computes xir/sigma/HPA); jl weights per-model
     OUTPUTS. ★ REFINED (2026-06-29): fmfint.f:513-550 CONFIRMS HPA = SXIR·384/SSIGMA with SXIR=Σxir·FWT, SSIGMA=Σsigma·FWT
     (line 515/513) — i.e. exactly jl's HPA_A weighting. And WMULT (the bisection wind reduction, fmcfir.f:248) =
     ALGSLP(PERCOV,CANCLS,CORFAC) = jl's `fire_wind_reduction(percov)` (matches). So the formulas all match; the residual
     is the FMDYN WEIGHTS: jl 9/10 = 74.1/25.9% vs live 72.0/28.0% (the flame is robust to this, the HPA is not). BUT a
     2% weight shift can't move HPA 577→~700 (the value live's OINIT1=90.6 implies), so there is ALSO a per-model xir/sigma
     or down-wood (SMALL,LARGE)-point difference feeding FMDYN. PINNING the final gap needs the debug-FVS HPA + per-model
     xir/sigma/FWT from fmcfir/fmpofl — but those writes don't surface (the PotFire path redirects output; try the DBS or
     an EVSET4 event-monitor var, NOT a fmcfir WRITE). OINIT1 NOT wired (96.9 vs 90.6 ≈ 7%). Everything else in the crown-
     fire chain (CBD, OACT1 crowning index, surface flame, fuel-model selection) is bit-exact/live-matching + committed.
  5. **FVS_Mortality per-species rows — ✓ FIXED (2026-06-29, suite 5191/2).** Live (dbsfmmort.f, shared SN/NE) emits one
     row PER SPECIES (SpeciesFVS/PLANTS/FIA columns) + an 'ALL' aggregate; jl emitted only the aggregate and lacked the
     species columns. Fixed: `fmburn!` now accumulates killed/total TPA + BA/vol by species×DBH-class (`species_mort`);
     `write_dbs_mortality!` adds the 3 species columns and emits per-species rows + ALL. Matches live per-species (WP K2
     22.15/22.19, JP, SM, YB, HI, QA, ALL K1 35.56/35.47); small Volkill diffs = the known DG/merch-volume drift. The
     pre-existing SN test (test_carbon.jl:560) pinned the OLD 1-row structure (a masked-bug per doctrine #3) — updated to
     assert the faithful per-species+ALL structure (SN-validated via the shared dbsfmmort.f path).
  6. **Fuels `Total_Consumed` = 0** at the fire year (Consumption table itself MATCHES) + standing-live fuel under-count.
  **Two ROOT CAUSES already located (actionable):** (2) `write_dbs_burnreport!` **hardcodes slope = `0`** at the 12th
  INSERT value (dbs_output.jl:339 — `Float64(b.wind), 0, ...`); the fix is to carry the stand slope onto `burn_reports`
  and write it. (3) the same writer stores the fuel-model weights `ww(i)` as RAW fractions while the moistures on the
  line ARE ×100 — live's columns are percent, so `ww` needs ×100 (the residual model-SELECTION diff jl 8/9/10 vs live
  9/10 is the separately-tracked FMCFMD weighting). Both live in the SHARED SN/NE writer ⇒ need an SN burn-report DBS
  oracle to confirm no SN regression before changing (deferred — no canonical test exercises these columns).
- **VERDICT:** #47 is NOT a fundamental gap (my earlier "live pools = 0" verdict was a `.tre`-placement test artifact —
  corrected). The FFE **simulation** (.sum) and the **core** DBS values (carbon live pools, flame/scorch, potfire
  mortality) are faithful. What remains is a bounded set of secondary DBS-column refinements (snag-carbon scale, slope
  field, fuel-weight ×100 units, crown-fire indices, mortality by-class split) — mostly shared-with-SN serialization, none
  exercised by the canonical net01 test (no DATABASE block). Repro: scratchpad `net01_jl.key`+`net01_jl.tre` / `net01_live.*`, `cmp.jl`.

### NE FFE (the largest remaining subsystem) — net01 stand 4 (`id=FFE`) activates `s.fire`
`run_keyfile` over the whole net01.key fails with `KeyError :v2t` at `ffe_seed_input_snags!`.

**`data/northeast/fire_species_props.csv` column → NE source (fully traced 2026-06-29):**
- `v2t,tfall_cls,leaf_life,dkr_cls,snag_cls` — `fmvinit.f` main `SELECT CASE(I)` block (lines **175-1149**,
  one CASE per species; has a DEFAULT). Extract these 5 base props per species.
- `snag_decayx` — DERIVED from `snag_cls` via `SELECT CASE(SNAGCLS)` (fmvinit.f:1204): {1→0.07, 2→0.21, 3→0.35}.
- `snag_fallx` — **constant 1.0** for ALL NE species (fmvinit.f:1201 `FALLX=1.0`). (SN has per-species values; NE does not.)
- `snag_alldwn` — **constant 50.0** for ALL NE species (fmvinit.f:1202 `ALLDWN=50.0`).
- `ls_spi` — **SPILS**, the NE→LS-variant crown-biomass species crosswalk (used by `fmcrowe.f`/`fmcbio.f` for
  crown biomass; `SELECT CASE(SPILS)`). Find the NE→SPILS mapping (a per-species integer table).
- `bark_eqnum` — bark-thickness equation number per species (fire_effects.jl bark thickness `_FM_BARK_B1`).
  Trace the NE bark-eqn assignment (likely `fmcblk.f`/`fmcrowe.f`).
- `biogrp` — biomass group (carbon.jl: `>5 ⇒ hardwood`). Trace `fmcbio.f`.

NOTE the snag-fall model: NE uses the `TFALL(I,1..5)` fall-time-by-DBH curve (TFALLCLS→`SELECT CASE` at
fmvinit.f:1164-1196) PLUS FALLX/ALLDWN; the jl SN snag model uses `snag_fallx`/`snag_alldwn` directly — VERIFY
the jl `snag_fall_density` semantics match the NE TFALL-curve model before trusting NE snag falldown.

**Also needed:** `fire_biomass.csv` (biomass eqns), NE fuel models; verify the FFE engine (fmburn/snag/carbon)
is variant-generic — carbon.jl:225 + fuel_loading hardcode `data/southern/` for some NATIONAL tables
(hwp_fate, fuel loading); confirm shared vs NE-specific.

**Validate** against live net01 stand-4 .sum (post-fire TPA/carbon/fuels) — the FFE stand IS the oracle.
Until ported, validate NE growth via `write_sum_file` on the non-fire stands (1/2/3/5), or skip stand 4.

#### FFE progress (2026-06-29) — incremental gap-chasing on stand 4, 6 gaps resolved
DONE:
- **`data/northeast/fire_species_props.csv`** (all 11 cols, 108 spp) — extracted: 5 base props from
  fmvinit.f main `SELECT CASE(I)` (175-1149); snag_decayx from SNAGCLS {1→.07,2→.21,3→.35}; snag_fallx=1.0,
  snag_alldwn=50.0 (constant); ls_spi=ISPMAP (fmcrow.f:51); bark_eqnum=EQNUM (fmbrkt.f:26); biogrp=BIOGRP (fmcblk.f:27).
- **`data/northeast/fire_biomass.csv`** (`bio_group` = BIOGRP, same table) — Jenkins biomass (fmcbio) is national.
- **`init_merch_standards!` variant-gated** (volume.jl) — NE fills the per-stand merch arrays from `_ne_merch`
  (IFOR-dependent), not the SN merch_specs.csv. SN bit-exact (5190/2).
- **`dbh_min`=5.0** added to NE species CSV (jenkins merch threshold; = `_ne_merch` dbhmin at IFOR=2, net01's value).
  NOTE: IFOR-2-specific (hardwood dbhmin is 6/8 at IFOR 1·3/4) — route through per-stand `c.sp_dbh_min` for IFOR-generality.

NEXT GAP = **fuel loading model structure** (`ffe_live_fuel_loading`/`fmcba!`) — the meaty FFE core, now
precisely mapped (fmcba.f, 2026-06-29). The jl `fuel_loading.jl` is SN-specific (ffe_dead_fuel_type maps FIA
forest-type→1-9 by SN ranges; ffe_forest_type uses SN species 4-14=pines + SN ecological logic; reads
`coef.ffe_fuel_dead[ft,:]` 9-row / `ffe_fuel_live[ft,:]`). NE structure (DIFFERENT — needs variant-dispatch):
- **Live fuel = CONSTANT** `FULIV=(0.31,0.31)` for ALL types (fmcba.f:68,129) ⇒ `ffe_live_fuel_loading(::Northeast)=(0.31,0.31)`. Trivial.
- **Dead fuel** = `FUINI(MXFLCL=11, 31)` indexed by `FTDEADFU = (IFFEFT-1)*3 + ISZCL` (fmcba.f:214; IFFEFT=11
  non-stocked → 31). So 31 cols = 10 forest-types × 3 size-classes + non-stocked. `STFUEL(ISZ,2)=FUINI(ISZ,FTDEADFU)`.
- **IFFEFT** (NE FFE forest type 1-11) ← `CALL FMNEFT(IFFEFT)` (fmcba.f:123) — the NE forest-type classifier
  (analog of SN FMSNFT/`ffe_forest_type`); PORT THIS. **ISZCL** (stand size class 1-3) — from COMMON (find its
  setter, likely FMNEFT or a stand-stats routine; probably a QMD/mean-DBH breakpoint).
★★★ NE FFE NOW RUNS END-TO-END + FIRE MODEL VALIDATED BIT-EXACT (2026-06-29; suite 5190/2 SN bit-exact).
DONE: `ne_dead_fuel_loading`/`_ne_iffeft` (FMNEFT) in fuel_loading.jl; ISZCL = jl's existing `plot.size_class`
(`_stkval_stocking`, already matches stkval.f); fmcba.jl branched for NE (constant live FULIV=(0.31,0.31) +
ne_dead_fuel_loading; the BA-share decay-class distribution is IDENTICAL NE=SN so it's shared); extracted
`fire_fuel_dead.csv` (31 FTDEADFU groups × 11 fuel cls, fmcba.f:75) + `fire_fuel_live.csv` (0.31/0.31);
`fire_fuel_models.csv` COPIED from SN (FMR table is IDENTICAL NE=SN — national). The fire BEHAVIOR
(rothermel/byram/scorch) + FMEFF mortality are SHARED FFE and work for NE unchanged. ★★ VALIDATION vs live
net01 stand-4 (id=FFE): **fire mortality BIT-EXACT at the fire year** — jl & live both remove **245 TPA** at
1993, and 2003 TPA=**285**/285, BA 99/99, SDI 193/193 (CCF 207/208, TopHt 68/69 within 1). The NE fuel-loading +
fire-behavior + fire-mortality port is FAITHFUL.
★ NE DEFAULT CYCLE LENGTH FIXED (masked bug; suite 5190/2 SN bit-exact). All NE stands defaulted to
`control.year`=FINT=**5** (SN default), but NE GRINIT (ne/grinit.f:172) is FINT=**10**. The no-fire stands
HID this because I'd forced `grow_cycle!(fint=10)`/`write_sum_file(period=10)` manually; the FFE stand (full
pipeline) exposed it (ran 5-yr cycles 1993,1998,2003… vs live 10-yr 1993,2003,2013…). FIX: `initialize!`
sets `s.control.year=10f0` for NE before keyword processing (TIMEINT/NUMCYCLE still overrides; SN keeps 5 ⇒
bit-exact). NOW ALL NE stands run 10-yr NATURALLY (no manual fint): no-fire stand-1 2000 = 524/107/211/228/72/6.1
vs live 524/107/213/229/72/6.1; FFE stand cadence = live's (1993,2003,2013,…), and **1993 + 2003 are bit-exact**
(536→285 across the SIMFIRE-2003 fire).
REMAINING FFE RESIDUAL (post-fire mortality, ~12%): after the 2003 fire (both 285 TPA), jl over-kills in the
2003→2013 cycle — jl 2013 TPA **148** vs live **168**; the extra ~20 TPA is killed ENTIRELY in that one post-fire
cycle (after 2013 both lose ~4/decade in lockstep), and it reduces BA proportionally (jl 73 vs live 81). The fire
itself + the fire-year (2003) state are bit-exact. ★ KEY DIAGNOSTIC: this is FIRE-SPECIFIC, not the general
mortality drift — the NO-fire stands have jl slightly UNDER-killing (2010 jl 479 vs live 475, +~2% keeps more),
the OPPOSITE sign. So the post-fire over-kill is a delayed FIRE mortality / fire-damaged-tree mortality booked in
the cycle AFTER the fire (FMEFF delayed kill, or post-fire snag/vigor mortality), localized to the first post-fire
cycle. ★ MECHANISM PINNED DOWN (via the .sum cut columns): the FFE stand's 1993 drop (536→291) is a THINDBH CUT of
**245 TPA — BIT-EXACT both sides** (cutTPA col); 2003 (285 TPA) is bit-exact; and the 2003→2013 drop (285→148/168)
is PURE MORTALITY (cutTPA=0) and HUGE (41% in live) ⇒ it is the **DELAYED FIRE MORTALITY** of the SIMFIRE-2003 fire
(scorched trees dying over the cycle AFTER the fire — the fire-year 2003 stand is intact, the deaths accumulate
2003→2013). jl over-kills this delayed mortality by ~17% (kills 137 vs live's 117), and proportionally more BA
(2013 BA jl 73 vs live 81). So this is NOT kill-distribution and NOT regular mortality — it's the FMEFF post-fire
mortality MAGNITUDE (the scorch/intensity → kill-fraction, which depends on the NE fuel loading I just ported).
★★ FIRE BEHAVIOR VALIDATED vs live (Burn Conditions Report in net01.out, via the FFE stand's BURNREPT keyword —
no debug-FVS needed): the 2003 fire FLAME LENGTH = **2.1** (jl 2.14) and SCORCH HEIGHT = **6.5** (jl 6.57) MATCH.
⇒ the NE FUINI fuel loading + the (shared) rothermel/byram/scorch fire behavior are FAITHFUL. So the kill divergence
is NOT the fire intensity. (Minor: the FMCFMD fuel-MODEL selection differs slightly — jl weights 8/9/10, live 9/10 —
but flame/scorch are unaffected; a low-priority FMCFMD-weighting check.) ⇒ With IDENTICAL scorch (6.5) the FMEFF
per-tree kill should match, so the over-kill (jl 150.8 vs live 117 TPA) must come from a DIFFERENT 2003 TREE LIST —
the per-tree stand composition differs (same TPA/BA aggregates, different species/DBH mix), most likely from the
1993 **THINDBH cut DISTRIBUTION** (which 245 TPA it removes) or the intervening growth. ★★★ RESOLVED — BIT-EXACT (suite 5190/2 SN bit-exact). The over-kill was the **NE FMEFF dormant-season mortality
reductions** (ne/fmeff.f:304-326), which the jl FFE (SN-based) lacked: for VARACD='NE' & BURNSEAS≤2, conifers
(sp≤25) PMORT×½, balsam-fir floor 0.7, small maples(26:29,99:100)<4″→1, hardwoods(sp>25) oaks(55:70,89) 2.5″+ ×½
else ×0.8, hardwoods ≤1″→1. ALSO: the SN/CS Regelbrugge species groups (fmeff.f:196) are skipped for NE (NE uses
the base Reinhardt logistic for ALL species — and that base = jl's group-6 Reinhardt exactly). FIX: fire_effects.jl
`fire_tree_mortality`/`fire_mortality_adjust` take the variant (NE ⇒ Reinhardt + the dormant-season block); fmburn.jl
passes s.variant; SN keeps Regelbrugge+no-op (bit-exact). RESULT: FFE stand post-fire TPA now **168/164/161/157 =
live EXACTLY** (was 148/144/141/137), BA within 1. The NE FFE fire path (behavior + mortality) is now FAITHFUL.
(Historical narrowing below kept for the record.)

RULED OUT (the FMEFF kill inputs all check out) — the over-kill is NOT any of: (a) fire intensity (flame 2.1/scorch
6.5 MATCH live); (b) bark thickness — jl `_FM_BARK_B1[39]` is IDENTICAL to ne/fmbrkt.f B1 (national table) and NE
`bark_eqnum`=EQNUM is correct; (c) crown ratio — per-tree vs live .trl, **70/81 within ±1** (jl +2-3% high only on
a few sp27/sp30; a minor NE TWIGS residual worth a look but too small for the 30% kill gap). ⇒ With scorch, bark,
and crown all ~matching, the kill gap (jl 150.8 vs live 117 TPA) must be the per-tree fire-KILL DISTRIBUTION on a
slightly-different 2003 stand (the THINDBH-cut mix or the few-% crown/DBH diffs compounding). NEXT: get the FFE
stand's 2003 pre-fire .trl jl vs live (which trees the THINDBH removes + which the fire kills) — build a clean
single-stand FFE key (REWIND→inline the stand-3 inventory) + TREELIST. Minor sub-leads: the FMCFMD model selection
(jl 8/9/10 vs live 9/10) and the sp27/sp30 crown +2-3%.
DEBUG-FVS NOTE: a `WRITE(*,…)` in fmburn.f did NOT surface (FVS unit handling) — use the BURNREPT/.out report or
the DBS instead. Also TODO: NE FFE DBS carbon/fuel tables; confirm no other FFE keyword paths need NE data.

### NE ESTABLISHMENT (stand 5 "BARE") — ★★ SUBSTANTIALLY RESOLVED (suite 5190/2 SN bit-exact)
Ported the REGENT(LESTB) establishment-growth pass into `establish!` PHASE 2 (NE-gated): each new seedling grows
its creation cycle by the NE small-tree height increment (`ne_htcalc_incr` + BALMOD + relht + DGSD ±10% random),
scaled by **FNT/REGYR with FNT=FINT−5** (regent.f:118-124 LESTB period; LSKIPH when FINT≤5), then its DBH is reset
from the new height (HTDBH⁻¹ + DIAM floor + 0.001·HK). RESULT vs live: TopHt 12/26/36 vs **13/26/36** (mid-cycles
exact), regen TPA still bit-exact (800/777/755), BA 8/60/100 vs 10/62/100. SN unaffected (NE branch only; SN's
essubh assigns the full height-at-age directly so needs no growth). The KEY was FNT=FINT−5 (a full-cycle scale
over-shot to TopHt 18). Residual: 2002 first-cycle slightly low (BA 8 vs 10, TopHt 12 vs 13 — likely the essubh
BASE; current uses the `ne_htcalc_height(sp,si,age)` band-aid, the faithful essubh `(H@refage/refage)·min(5,time)`
is in estab_ref_age.csv if needed) + the general cyc2+ mortality drift at 2032+. The full net01.key now runs
end-to-end with faithful output for ALL 5 stands. (Historical root-cause analysis below.)

### NE ESTABLISHMENT — historical root-cause (kept for the record)
The full net01.key run_keyfile now gets past the FFE stand (4) but fails on stand 5 (BARE, ESTAB/PLANT from
bare ground) with `KeyError :ht_curve_b1` in `establish!` (establishment.jl:50) — it reads the SN Chapman-Richards
ht-age curve to assign regen-tree heights, but NE uses the NC-128 curve (`ne_htcalc_height`, already ported in
height_growth.jl). Same curve-family split as small_tree_growth/height_growth. DONE the variant-gate (establish! line 82 NE→ne_htcalc_height; bc SN-only) ⇒ **full net01.key now runs
end-to-end (all 5 stands)**. BARE stand regen **TPA is bit-exact (800/778/755)**, but regen-tree HEIGHTS run low
(jl TopHt 6/19/31 vs live 13/26/36). ROOT (traced): the NE establishment height is NOT the large-tree NC-128 curve
at tree age — ne/essubh.f uses `HHT=(H@CARAGE / CARAGE)·min(5,TIME−DELAY)` where CARAGE=a per-species REFERENCE-AGE
table (essubh.f DATA MAPNE, sp1=20; extracted to data/northeast/estab_ref_age.csv, NOT yet wired). BUT that creation
height ≈ 5.6 ft ≈ jl's current value — so the 2002 gap (6 vs 13) is mostly the regen GROWTH after creation, i.e. the
NE REGENT small-tree growth on the young planted trees grows them too slowly (or the establishment timing/age is off).
★ ROOT CAUSE DEFINITIVE (traced the FVS call chain): regen TPA is bit-exact but heights are low because jl
`establish!` creates the regen records but NEVER GROWS them. FVS chain: fvs.f→TREGRO→GRADD(gradd.f:229
CALL ESNUTR)→ESTAB→ESGENT, and **esgent.f:48 `CALL REGENT(.TRUE.,ITRNIN)`** — the LESTB establishment-growth pass
that grows the just-created seedlings from the essubh BASE height (~5 ft, confirmed in jl for both BARE species:
balsam-fir SI52 / white-spruce SI50, refage 20/15 → essubh ~5) up to the cycle-end height (~13). The jl REGENT
port (`small_tree_growth!`) only implements the LESTB=false regular-growth branch; the LESTB=true regen-growth
branch is NOT ported. So `establish!` records keep their ~5-6 ft creation height (jl's per-decade slope AFTER
matches live, confirming it's only the missing creation-cycle growth). FIX = port the REGENT LESTB branch
(regent.f:163-407: assign establishment DBH from HK, then grow HTG/DG using the essubh AGE) and call it on the
new records inside establish!. **NOT cross-cutting (corrected)**: SN is FINE — sn/essubh.f:66 `CALL HTCALC(MODE1,I,AGE,HHT,…)` assigns HHT =
height-at-AGE DIRECTLY (full height), so jl's SN `htcalc_height(bc,sp,si,age)` is FAITHFUL & validated. The NE
essubh is a DIFFERENT model (ne/essubh.f: linear base `(H@CARAGE/CARAGE)·min(5,TIME−DELAY)` ⇒ ~5 ft, THEN
esgent.f:48 REGENT(LESTB=T) grows it to the age-appropriate height). So the fix is NE-ONLY: port the NE essubh
base height (estab_ref_age.csv refage) + the REGENT LESTB growth pass, called on the new records in `establish!`
(behind `s.variant isa Northeast`). The current `ne_htcalc_height(sp,si,age)` band-aid (~5-6 ft) under-shoots
because the NE tree's effective age at the first reported cycle is ~12-14 (height ~13), not the jl `age`=7 — the
NE essubh+LESTB timing produces the higher effective age. SN keeps its direct path (bit-exact). Stands 1-4 don't
hit establish! and ARE validated.

### Smaller residuals
- cyc2+ TopHt ~2 low: sp9 (white pine) mid-tree (d5-7) HTG/DG run ~0.8× live — the DG OLDRN
  serial-correlation RNG scatters per tree; species mean is close (DG 94%, HTG 92%). Deep RNG-alignment
  (the SN bit-exact standard), small stand-level impact. Investigate the OLDRN draw order for NE.
- Broken-top (CFTOPK) reuse in `compute_volumes_ne!` is TODO.
- NE small-tree height calibration (HCOR, regent.f:431-558) is a no-op for net01 (the .tre has measured
  DG, not measured HTG) — needs porting if a stand carries measured small-tree height growth.

## Per-tree validation recipe (the key NE tool)
Build a clean single-stand key (STDIDENT/DESIGN/STDINFO/SITECODE/INVYEAR/NUMCYCLE + `TREELIST 0` +
TREEFMT(2 lines) + TREEDATA + PROCESS + STOP), copy `net01.tre`→`<stem>.tre` alongside (key and outdir
in DIFFERENT dirs), `bash test/harness/ne_oracle.sh <key> <outdir>` → `<stem>.trl`. Parse `END CYCLE: N`
blocks; tree rows: token[2]=TREE INDX (1..n = jl tree order, match on THIS), token[9]=CURR DIAM,
token[11]=CURR HT, token[21]=SAW BD. Applied DG = d(cyc1)−d(cyc0). jl side: each_stand+notre!+
setup_growth!+compute_volumes!+diameter_growth!(tripling=false,sfint=10)+height_growth!+small_tree_growth!.
