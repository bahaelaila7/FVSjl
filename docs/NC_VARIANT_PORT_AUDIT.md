# NC (Klamath Mountains) variant port — chunk audit (doctrine #6)

## ★★★ 2026-08-13 — #210 DGSCOR CALIBRATION VERIFIED FAITHFUL (agent's "COR regression" REFINED/CORNERED)
Measured jl's DGSCOR calibration on nct01 vs the live nct01.out CALIBRATION STATISTICS + ZNC dump. NC calibrates
2 species: SP(sp2) + WF(sp4). VERDICT: the calibration is STRUCTURALLY FAITHFUL — the "COR is slightly off" lead
is real in isolation but COMPENSATED, netting <0.5% on calibrated trees:
- **fn (records for scaling) MATCH**: jl SP=6/WF=5 = live "NUMBER OF RECORDS AVAILABLE FOR SCALING" 6/5.
- **wc (weight to input data) MATCH**: jl SP=1.000/WF=0.9053 = live "WEIGHT GIVEN TO INPUT GROWTH DATA" 1.00/0.90.
- **backdating correct**: jl calibration wk2 uses the BACKDATED DBH (t.dbh[i]), saved_dbh is current (bark only).
- The agent's reported values are exp(corv): jl SP exp(−1.0447)=0.3518, WF exp(−0.5656)=0.5680 (agent "jl"); live
  SP exp(−1.0402)=0.3537, WF exp(−0.5454)=0.5796 (from ZNC pre−post diff). jl COR is ~0.005(SP)/0.020(WF) MORE
  NEGATIVE — but jl's predicted wk2 is correspondingly ~+0.005(SP)/+0.016(WF) HIGHER, so DG=exp(wk2+COR) NET differs
  by only ~−0.4%. i.e. the COR shrinks a matching wk2 over-prediction; measuring COR ALONE (as the agent did) double-
  counts. Residual small wk2 offset (WF ~+1.5%) is a candidate minor DGCON/CONSPP-constant discrepancy, COR-masked on
  calibrated stands. ⇒ the ~2-3% multi-cycle BA under-growth = the CORNERED OLDRN serial-corr straddle (NC DGSD=2.0,
  #206 class, same as CI +2% / utt01 −9%), NOT a COR bug. #210 CLOSED-CORNERED; NC growth stays at the bar.

## ⚠ NEW BUG (2026-08-13, separate from #210) — NC FFE crashes: fmcba.jl:114 BoundsError (0×0 crown-biomass matrix)
The canonical nct01.key is an **FFE TEST keyfile** (FMIn/SIMFIRE/SNAGINIT/PotFIRE/BurnRept). Running it crashes in
the FIRE model — `fmcba.jl:114` BoundsError "0×0 Matrix at [6,1:0]" (fmburn.jl:92 → simulate.jl:345). ROOT (pinned
2026-08-13): `ffe_dead_fuel_loading(coef, ifortp)` (fuel_loading.jl:137-141) does `coef.ffe_fuel_dead[ft,:]` with
ft=`ffe_dead_fuel_type(ifortp)`=6, but **`coef.ffe_fuel_dead` is 0×0 for Klamath** — NC's FFE fuel-loading tables were
never loaded (NC ported growth+volume ONLY). fmcba's live+dead fuel dispatch (fmcba.jl:100-124) also has NO Klamath
branch ⇒ falls to the generic `ffe_dead_fuel_loading` default, which indexes the empty table. ⇒ NC FFE is a genuine
UNPORTED extension: the port chunk = load NC's FFE fuel tables (nc/fmcba.f FUINI/FULIVE + dead-fuel-type map) into the
Klamath SpeciesCoefficients + add Klamath branches to fmcba live/dead dispatch (+ likely snag/decay/potfire chain).
Growth-only NC runs are UNAFFECTED (cyc0 calibration completes; crash is fire-path only). This is a #207 westside-stream
lead, NOT a growth-parity gap. NOTE: jl should also GATE FFE off (or error cleanly) for variants without fuel tables
rather than BoundsError — a small defensive-robustness follow-up independent of the NC FFE port.

## ✅ NC FFE chunk 1 — FMCBA initial surface-fuel loading PORTED (2026-08-13)
Fixes the fmcba crash. NC is a **California/westside** FFE variant: unlike the interior western variants (CR/IE/EM/CI/
TT/UT/BM, single-COVTYP interpolation), NC (like WS/CA) initializes the live + initial-dead fuel pools from the **TOP
TWO cover-type species** (COVCA(1..2), weighted COVCAWT by their share of the top-2 BA), interpolated by PERCOV between
the initiating (10% cover, FULIVI/FUINII) and established (60% cover, FULIVE/FUINIE) tables — nc/fmcba.f:236-431.

Ported (all Klamath-guarded; other variants provably inert):
- `data/klamath/fire/ffe_fuel.jl` — `_NC_FULIVE/_NC_FULIVI` (2×12) + `_NC_FUINIE/_NC_FUINII` (11×12) verbatim from
  nc/fmcba.f; `nc_live_fuel_loading` + `nc_dead_fuel_loading` (top-2 COVCA/COVCAWT interpolation, reuse `_cr_algslp2`).
- `nc_cwcalc` (same file) — NC CRWDTH via NCMAP (nc/cwcalc.f), REQUIRED for PERCOV: the generic `crown_width` returns the
  0.5 default for every NC species ⇒ PERCOV≈0 ⇒ mis-selects the initiating-stand loads (bug reproduced: 0.27% → 44.3%).
  Reuses the shared `_cr_r6m2` (Crookston R6 model 2) form; ports the SP/DF/WF/RF/PP/OS equations (forest 505=R5 ⇒ BF=1).
  MA/IC/BO/TO/OH/RW equations error loudly (not exercised by nct01) — a follow-up crown-width chunk (also NC growth CCF).
- `fmcba.jl` — Klamath branches: RDPSRT top-2 COVCA/COVCAWT; live-fuel (deferred, post-cover-type); dead-fuel dispatch;
  `nc_cwcalc` in the PERCOV crown-area loop; bare-stand default = DF(3) (full COVINI5/6 R5/R6 habitat maps deferred —
  a bare-stand-only path, trees present in nct01).

**VALIDATION (nct01 FFE TEST stand, cyc1/1993 vs LIVE FVSnc_clean `DEBUG FMCBA` + ALL FUELS) — near bit-exact, CORNERED:**
Live ground truth measured via `DEBUG  1.0  1.0` + supplemental `FMCBA` on FVSnc_clean (initre.f:1045 scopes when DEBUG
field 2 is non-blank; a bare DEBUG = global and hits the volume-debug segfault): **COVTYP=2, PERCOV=39.01, FLIVE=(0.257,0.639)**.
| | COVTYP | PERCOV | HERB | SHRUB | LITT | DUFF | 0-3" | 3-6" | 6-12" |
|---|---|---|---|---|---|---|---|---|---|
| jl   | 2 | 39.72 | 0.254 | 0.621 | 0.478 | 14.82 | 3.10 | 5.30 | 5.76 |
| live | 2 | 39.01 | 0.257 | 0.639 | 0.54  | 14.7  | 3.2  | 5.2  | 5.7  |
COVTYP bit-exact; PERCOV within 1.8%; DUFF/6-12"/3-6"/FLIVE within ~2%. The fuel TABLES are bit-exact by construction
(verbatim nc/fmcba.f). Residual = the small PERCOV gap (crown-ratio/HT init: jl's cyc1 ICR/HT vs live's) + the LITT term
(0.478 vs 0.54 — likely one year of FMCADD litterfall folded into the report, a downstream chunk). CORNERED.

### RE-ATTRIBUTION (root cause of the earlier ~5-8% residual — MEASURED, fixed)
The earlier draft wrongly blamed a "~10% heavy base tree list." **Corrected (per coordinator + measurement):** jl's raw
`stand_tpa` (589.65) is the INTERNAL ×GROSPC per-stockable PROB (notre.f:69 PROB=P·GROSPC, GROSPC=11/10 on DESIGN `11.0 1.0`);
the `.sum` divides GROSPC back out (stats.f:128) ⇒ jl cyc0 `.sum` = 536/77/5.1 = live BIT-EXACT (#210 holds). FFE fmcba
LEGITIMATELY uses the ×GROSPC PROB (per-stockable), which jl's `t.tpa` basis already matches — so that was NOT the cause.
The real cause of the fuel residual was **jl passing the dense computed stand BA (85) as BAREA to `nc_cwcalc`**: FVS's base
`cwidth.f`→`cwcalc.f` computes the INITIAL CRWDTH at LOAD time, before the stand BA is accumulated, so the R6-Crookston
`(BAREA+1)^b` term hits `cwcalc.f:859 IF(BAREA.LE.1.) BAREA=1.` (BA=0→1). cyc1 fmcba reads those load-time crown widths.
BAREA=85 inflated CW by `(86/2)^0.04267≈1.174` ⇒ TOTCRA +18% ⇒ PERCOV 44.3 (vs live 39.0). **Fix (fmcba.jl):** NC passes
BAREA=1 for the load/first FFE cycle (`s.control.cycle <= 1`), the actual stand BA thereafter (matching the end-of-cycle
UPDATE CRWDTH). Measured: PERCOV 44.3 → 39.72 (live 39.01); DUFF 15.73 → 14.82 (live 14.7). All Klamath-guarded.

### Remaining NC FFE chunks (fire-behavior chain — each is a real port, revealed in order by measurement)
1. **standard fuel models (CSV LANDED) + NC fmcfmd CWHR (PORTED, VALIDATED bit-exact).** `data/klamath/fire_fuel_models.csv`
   added (13 universal Anderson models). NC's fuel-model SELECTION is the California **CWHR classifier** — ported to
   `src/engine/fire/nc_fuel_model.jl` (`nc_cwhr` + `nc_select_fuel_models`, verbatim `nc/cwhr.f` + `nc/fmcfmd.f`: PCNETAVG,
   size 1-6 × density S/P/M/D structural stage, FMD_R5/R6 9×18 table, IFT forest-type, density sub-model blending via FMDYN
   over CWXPTS, natural fuels 10/12/13). MEASURED vs live `DEBUG FMCFMD` @2003 fire: **IFT=9(OS), base model 6, candidates
   {6,10,12,13}** — jl's path matches live exactly. And `_fmdyn` fed **live's** point (SMALL=4.80, LARGE=13.28) returns
   **6@0.443 / 10@0.557 = live's 6@44/10@56 BIT-EXACT.** The port is correct.
   Residual: jl's own point is SMALL=3.5/LARGE=10.2 (vs live 4.80/13.28) ⇒ jl weights 6@76/10@24, flame 4.27/scorch 17.9 vs
   live 8.3/47. That gap is the **fuel-ACCUMULATION** difference over 1993→2003 (LARGE 10.2 vs 13.28), NOT the fmcfmd port —
   it is downstream of (a) the **Dunning decay multiplier** (chunk 4, deferred — scales DKR over the 10-yr run) and (b) the
   **crown-biomass litterfall/crown-lift additions** (chunk 2). The missing NC crown-fire flame boost (chunk 2, NC absent from
   the fmburn crown-fire variant list) further lowers flame. ⇒ chunks 2+4 close the flame/scorch/mortality gap.
   FM11 post-activity sharing (AFWT/SLCHNG/HARVYR, <5 yr after an entry) is not yet wired — inert for nct01 (2003 fire is
   10 yr after the 1993 THINDBH ⇒ AFWT=0, FM10=1, exact). (Also SEPARATE from FFE: the BARE/PLANT stand 5 crashes in
   `establish!` on `KeyError :estab_min_ht` — an NC regen coefficient gap, chunk 3.)
2. **crown biomass (ROUTED) + crown fire (ENABLED) — over-shoots pending chunk 4.** NC's `fmcroww.f` is byte-identical to
   CR's, so NC routes through `cr_crownw` with `_NC_ISPMAP`=[3,15,3,4,10,20,21,17,4,13,17,19] (nc/fmcrow.f); nct01's groups
   {3 DF/OS, 4 WF/RF, 15 SP} are already ported (hardwood/cedar {10,17,19,20,21} error loudly, not in nct01). NC added to
   the fmburn crown-fire variant list + the crowning/torching/crown_fire_result Unions. STATUS: with crown fire enabled,
   nct01 2003 flame=14.13/scorch=68 vs **live 8.3/47** ⇒ OVER-kill (2003 mortality jl 451 vs live 270; 2008 TPA jl 5 vs live
   58). Root cause is NOT the crown-biomass port — it is the still-wrong fuel-model WEIGHTS (jl 6@76/10@24 vs live 10@56/6@44,
   from the LARGE-fuel accumulation gap 10.2 vs 13.28), which the crown-fire boost then amplifies. ⇒ **chunk 4 (Dunning decay)
   is the gating dependency for flame/scorch/mortality convergence** — it is coupled with this chunk, not independent.
3. **NC establishment (DONE — nct01 runs to completion).** Wired the western establishment path for Klamath: `_NC_ES_XMIN`
   / `_NC_ES_HHTMAX` (nc/blkdat.f), `_NC_ESSUBH_HHT` fixed base-height table (nc/essubh.f), the western `bc=nothing` +
   PLANT-no-RAN branches. The BARE/PLANT stand 5 now completes. **nct01 runs end-to-end with NO fmcba.jl:114 crash** — the
   headline crash is RESOLVED. `.sum` growth columns (TPA/BA/SDI/QMD/removals, e.g. cyc0 536/77/160/5.1) are BIT-IDENTICAL to
   `.ncwork/ncval/nct01.sum`. (Two residuals visible in the .sum, both PRE-EXISTING / outside the FFE port: a forest-code
   resolution difference — jl reports forest 999 vs live 371 ⇒ CuFt 1261 vs 1308, a growth/volume matter; and the fire
   over-mortality above, which is chunk 2/4 fire-behavior.)
3. **fire mortality (FMEFF) + snag props** — the current `data/klamath/fire_species_props.csv` is a **CI COPY placeholder**
   (commit 00f69ed) whose v2t/dkr_cls/leaf_life/fallx do NOT match nc/fmvinit.f (correct NC DKRCLS=[3,4,3,4,3,2,2,4,4,4,4,1];
   V2T sp1=28.7/sp2=21.2/…). Rebuild it from nc/fmvinit.f before validating snag falldown / decay / mortality. (dkr_cls does
   NOT change the cyc1 SURFACE-DEAD totals — only the decay class the fuel lands in — so it does not affect the cyc1 row above.)
4. **Dunning decay multiplier** (DCYMLT, nc/fmcba.f:395-414) — first-year DKR adjustment by Dunning-code/site index; affects
   fuel DECAY in cyc2+, not the cyc1 initial loading.


Branch kt-variant-port. Oracle /workspace/.ncwork/FVSnc_clean. Canonical stand nct01
(tests/FVSnc/nct01.key + nct01.tre + nct01.sum.save). MAXSP=12, imperial, Zeide SDI, DGSD=2.0.
Coefficient extraction (ground-truth): docs/NC_CHUNK1_EXTRACTION.md.

## Chunk status
- **ch0 scaffold** ✅ (0d2437c): struct Klamath, variant_from_code("NC"), FVSjl.jl includes. Package loads, dispatches.
- **ch1 species** ✅ (42a3825 data + f83779b species.jl): species_coefficients.csv (41 cols, measured flat values +
  placeholders for structural cols) + species_translation.csv + ecocls_sdimax.csv (90-row PA→SDIMX). species.jl
  block-data init (Zeide, 10-yr, dg_sd=dg_stddev_bound=2.0, seed 55329). coefficients(Klamath()) LOADS correct.
- **ch2 site_index** ◐ IMPLEMENTED (63c0ead), NOT yet validated: site_setup! = nc_forkod (505/510/514/611/705/800/
  712→IFOR1-7) + nc_sitset (SITEAR site-range interp + provisional SDImax). Pipeline runs through init + into
  grow_cycle!, advancing to the ch3 gate dgf!(::Klamath).
  ✅ **cyc0 "blocker" was a HARNESS ARTIFACT, NOT a bug** (corrected): the jl TPA 28 came from measuring
  `each_stand`'s PRE-EXPANSION state — the documented "each_stand returns pre-calibration state, measure in the REAL
  run only" trap (CI #142 hard lesson). PROOF: CI cit01 (same 248112 stand) gives the IDENTICAL pre-expansion state
  via each_stand (n=27, gross_space 1.1, TPA 28, tpa[1:3]=[1,3,1]) yet its REAL-run cyc0 is bit-exact 536. So NC's
  tree reading is fine (identical to CI); the PROB×plot-design expansion to 536/ac happens in the real run, after
  each_stand's yield. ⇒ cyc0 must be validated from run_keyfile's .sum, NOT each_stand. That currently needs the ch3
  dgf! gate to not crash (a no-op stub suffices for the cyc0 row, which is pre-growth). NEXT: stub/implement dgf!
  → get run_keyfile's 1990 .sum row → validate TPA/BA/SDI/TopHt/QMD vs live 536/77/160/63/5.1.
- **ch3 DG** ▶ NEXT gate: dgf!(::Klamath). NC-specific (ISCT sections + DGHAH + 2 sets + sp12 redwood + nc_bratio
  POWER). Coefficients measured (extraction doc).
- ch4 height · ch5 crown · ch6 regent · ch7 mortality (reuses EM/UT) · ch8 volume (vollib): pending.

## Validation gate
cyc0 nct01 (live 1990): TPA 536 BA 77 SDI 160 CCF 87 TopHt 63 QMD 5.1 TCuFt 1308 MCuFt 449 (forest 371... wait
STDINFO forest=505 Klamath; .sum ForType col=371). Validate TPA/BA/SDI/TopHt/QMD FIRST (no growth models needed) —
currently blocked by the tree-expansion bug above.

## Chunk-3 DG — data COMPLETE + integration scoped (2026-08-11), ready to code
All DDS coefficients measured (docs/NC_CHUNK1_EXTRACTION.md chunk-3 section). Integration with the shared engine:
- jl PROVIDES the point-level density NC's dgf needs (density.point_ba=PTBAA, point_bal=PTBALT, point_ccf=PCCF,
  point_tpa=PTPA) — the shared standstats even computes PCCF specifically for the dgf DGPCCF term (doctrine #5 win).
- dgf! stores per-tree DDS in scratch.wk[2,i]=wk2 (CI pattern); the shared DDS→DG + serial-corr engine converts.
- ONE gap: PRD = ZRD(pt)/XMAXPT(pt) (point ZEIDE relative density) used ONLY by sp2/6/9 + redwood branches — jl
  has no point-Zeide-RD yet; add it (point Σ(D)^1.605 / point XMAXPT=Σ SDIDEF·point-BA-frac) OR compute inline in
  nc dgf!. nct01 has SP(sp2) trees ⇒ needed for nct01 validation.
- NC redwood DDS is a DIB²-diff form (like CI's diagr species) + POWER bark (nc_bratio) — the shared engine handles
  the diagr convention (CI proves it).
TO CODE (next): src/variants/klamath/diameter_growth.jl = nc_bratio(a,b,eqtype,d) + nc_dgcons!(s) [DGCON/DGDSQ per
sp, 3 branches: DGFOR/MAPLOC default, DGLAT2+site sp2/6/9, redwood ln(SITEAR)] + dgf!(s,::Klamath) [3-branch DDS →
wk2, 5-yr TDDS/2] + point-Zeide-RD. Then validate DDS per-tree vs FVSnc_clean dgf DEBUG dump on nct01.

## Chunk-3 DG — IMPLEMENTED (ce5111c) + first validation vs FVSnc_dbg (2026-08-11): NOT yet bit-exact
- ★ BUG FOUND + FIXED by reading nc/dgf.f: the 5-yr TDDS/2 conversion applies ONLY to the redwood + sp2/6/9
  branches, NOT the DEFAULT branch (nc/dgf.f has TDDS/2 inside CASE(12) and CASE(2,6,9), none in CASE DEFAULT).
  jl had applied /2 to ALL branches ⇒ DEFAULT species off by ln(2)≈0.693. Moved /2 into the two special branches.
- VALIDATION SETUP: instrumented FVSnc dgf.f (WRITE 'ZNC',I,ISPC,D,DDS after WK2(I)=DDS), recompiled dgf.o
  (gfortran-16 -std=legacy -w -fno-automatic -O0 -fPIC -I../../common), relinked → /workspace/.ncwork/FVSnc_dbg
  (ZNC dump in fort.16). Oracle dgf.f/dgf.o restored pristine after. jl side: env NC_DDS_DBG print in dgf!.
- FIRST COMPARISON (not yet aligned): jl sp2 D=10.39 DDS=2.50 vs live sp2 D~10.41 DDS~2.17; jl sp9 D=3.94 DDS=1.675
  vs live 1.171. jl HIGHER on the set-2 (sp2/6/9) species by ~0.3-0.5. ⚠ CAVEAT: both dumps mix CALIBRATION passes
  (backdated/modified DIAM) + the growth pass, and tree indices don't align 1:1 (jl re-sorts) ⇒ the comparison is
  NOT pass-aligned yet (doctrine #3 class). NEXT: align by (ISPC, exact DBH) on the GROWTH pass only (gate the dump
  on a growth flag or match DBH exactly), then diagnose the set-2 residual (candidates: the sp2/6/9 DGCON site form
  DGLAT2/DGSITE, the CRID/PBAL terms, or the -8.52-then-/2 order). DEFAULT-species DDS (post-/2-fix) still to compare.
⇒ Chunk-3 status: IMPLEMENTED + running to ch4; /2 bug fixed; per-tree DDS validation IN PROGRESS (not bit-exact
yet — set-2 species show a residual to diagnose). FVSnc_dbg durable for the aligned comparison.

## Chunk-3 DDS validation — refined diagnosis (2026-08-11): sp2/6/9 FORMULA verified faithful; residual = align/COR
Re-checked the jl sp2/6/9 branch term-by-term vs nc/dgf.f CASE(2,6,9): DDS = CONSPP + DGLD2·ALD + DGDSQ2·D²/1000 +
DGCR2·CRID + DGDBA2·PBAL/ln(D+1)/100 + DGBA2·ALPBA — EXACT match. So the ~0.3-0.5 jl-high on set-2 is NOT a formula
bug. The dump comparison was NOT valid: (a) fort.16 holds ZNC for ALL 10 nct01 cycles (the "last block" is a LATER
cycle with grown DBH 15-19", not cyc1); (b) both sides mix calibration passes (backdated DIAM, COR=0) + the growth
pass (which adds the DGSCOR-calibrated COR); (c) tree indices don't align (jl re-sorts). ⇒ residual candidates,
in order: (1) the growth-pass calibration COR (jl c.dg_cor set by the shared calibrate_diameter_growth! — verify it
computes NC COR like live's DGSCOR; the deterministic pre-COR DDS is the thing to match first), (2) CONSPP/DGCON
site inputs (SITEAR/ELEV/SLOPE for nct01 — verify jl's sp_site_index[2,6,9] + p.elevation/slope match live), (3)
per-tree CRID/PBAL/ALPBA. CLEAN NEXT: use nc/dgf.f's own DEBUG FORMAT-334 (DDS,CONSPP,ALD,ALPBA,CRID) — enable
DEBUG on FVSnc_dbg for ONE cyc1 growth-pass tree + match jl by exact (species,DBH) ⇒ isolates CONSPP vs the terms.
BETTER GATE (once ch4-8 land): the aggregate cyc1 .sum (TPA/BA/QMD) which averages out per-tree/pass noise, like CI.

## Chunk-4 height (nc/htgf.f) — structure characterized (2026-08-11), substantial port (NOT CI-reusable)
NC htgf = per-species SELECT CASE, using DG(I) (from dgf!) + the site-index height curve:
- CASE(12) redwood: DG10=DGLT/BRAT; LTHTG=EXP(1.412947 −0.000204·D² +0.31971·lnD +0.394005·ln(SINDX)
  +0.399888·ln(DG10) −0.451708·lnH)·0.5; HGBND height bound (H<217→1.0; 217-380 ramp to 0.1; ≥380→0.1);
  HTG=LTHTG·HGBND. (If H<4.5, DG10=0.1.)
- CASE DEFAULT (conifers+hardwoods): CALL FINDAG(I,ISPC,D1,D2,H,SITAGE,SITHT,AGMAX,HTMAX,HTMAX2,…) using HD1-4
  (P1-P4, the HT-DBH coeffs) + SITEAR → SITAGE (site age from H) + HTMAX. Then POTHTG from the site-curve
  (HGUESS−SITHT over SITAGE→SITAGE+5); H≥HTMAX→HTG=0.1. Scaled by SCALE=FINT/YR, XHT=XHMULT (MULTS keyword),
  exp(HTCON). Hardwoods (sp5/7/8/11) use their HD1-4; conifers use the site curve too (HD1-4=0 ⇒ different FINDAG path).
- NOT reusable from CI height_growth.jl (CI = NI-conifer exp form + Weibull, a DIFFERENT model). NC needs its own
  height_growth.jl porting FINDAG (nc/findag.f: site-index age↔height curve) + the redwood LTHTG + the DEFAULT
  POTHTG loop. FINDAG is the Dixon site-index curve (shared by the Pacific-coast/CA variants NC belongs to).
TO CODE (chunk 4): read nc/findag.f (the site-curve age/height solve) + the DEFAULT POTHTG iteration (lines 200-280);
port nc_findag + height_growth!(::Klamath) (redwood LTHTG + default FINDAG/POTHTG); validate HTG per-tree /aggregate.
Coefficients: HD1-4 measured (extraction doc); SITEAR from site chunk; HTCON/XHMULT default 0/1.

## ★★★ 2026-08-11 MILESTONE: NC nct01 control stand RUNS END-TO-END, cyc0 BIT-EXACT
All 8 growth+mortality+crown chunks implemented (ch1 species, ch2 site, ch3 DG, ch4 height, ch5 crown, ch6 regent,
ch7 mortality) + pipeline plumbing (FFE-inert gate, bark_intercept cols, v2t, volume placeholder). NC nct01
(control stand) runs full multi-cycle. cyc0 (1990 inventory) BIT-EXACT vs live: TPA 536, BA 77, SDI 160, TopHt 63,
QMD 5.1 (all ==). ⇒ chunks 1/2 + the stand summary are FAITHFUL.
REMAINING (3 items to bit-exact-or-cornered):
1. CYCLE LENGTH: live nct01 runs 5-YR cycles (1990,1995,...→2040); jl ran 10-yr (control.year=10 in species.jl is
   WRONG for NC — set to 5). nc/grinit IFINT=10 but the stand default cycle is 5 (verify: no TIMEINT in nct01.key
   ⇒ NC default cycle = 5, fix species.jl control.year/growth_fint = 5).
2. CCF = 0 in jl vs 87 live — the crown-competition-factor / per-tree CCF not wired for NC (crown.f ccfcal). Affects
   the .sum CCF column + RELDEN (density → DG/crown). Wire nc CCF (ccfcal) — likely reuse the shared ccfcal.
3. VOLUME = 0 (CI r4vol placeholder gives 0 for NC species). Real NC VEQNNC (R5/NVEL California equations) = chunk-8.
After cycle-length fix + CCF, validate the multi-cycle trajectory vs nct01.sum.save (5-yr, TPA 536→357 by 2040,
BA 77→301). This is the aggregate .sum gate.

## ★ 2026-08-11 TopHt under-growth ROOT: site index defaults to site_lo (chunk-2 site-input bug)
After the htg_period=5 fix, BA/SDI/QMD track live closely but TopHt DECREASES (jl 64→57 vs live 72→91). Traced:
the largest-DBH DF tree gets SITE INDEX = 50 = its site_lo DEFAULT, not the input value → FINDAG SITAGE=86 (old) →
POTHTG=1.45 (tiny) → height barely grows. ROOT: nc_sitset! (site_index.jl) fills sp_site_index from site_lo when
p.site_species/site index isn't set from the SITECODE/STDINFO input. cyc0 is bit-exact because initial HEIGHTS come
from the tree input, but the growth site curve uses the wrong (low) index. nct01 SITECODE/STDINFO should set the
real site index (STDINFO field ~84). FIX: parse the SITECODE/STDINFO site species+index into p.site_species +
p.sp_site_index[isisp] so nc_sitset! interpolates from the REAL index. This should fix TopHt (and refine DG, which
reads SITEAR). Then re-validate multi-cycle. cyc0 fully bit-exact; BA/SDI/QMD good; TopHt gated on the site index.

## Site-index/TopHt — CONFIRMED root (2026-08-11 test): forcing site index 84 (STDINFO fld2) makes TopHt GROW
Temp test (reverted): default site index 50→84 ⇒ TopHt 63→57(decrease) becomes 63→72(grows); live 63→91. So the
site index IS the TopHt driver — the site_lo default (50) froze/shrank it. Proper FIX (chunk-2 completion): NC has
NO SITECODE in nct01 ⇒ live DEFAULTS the site species to DF (sp3, per the .out "SITE SPECIES=DF CODE=3") and derives
the site index from the habitat/ecocls SITE (or estimates from the site-species trees). Wire in nc_sitset!: default
site_species=DF; get the real site index (STDINFO habitat 84 → ecocls SITE per PA, OR the sitcind.f site-tree
estimate) → p.sp_site_index[DF] → interpolate. Residual after 84 (jl 72 vs live 91) ⇒ the true index is higher
and/or a small height-model residual — measure live's exact SITEAR (instrument FVSnc_dbg sitcind) to pin it.
STATUS: NC growth port ~90% — cyc0 FULLY bit-exact; BA/SDI/QMD track live closely; TopHt gated on the site-index
derivation (confirmed root); volume = placeholder (real NC VEQNNC R5/NVEL pending). Clear remaining scope.

## Chunk-8 VOLUME — VEQNNC measured from live nct01.out (2026-08-11): NVEL 500WO2W conifers + 500DVEW hardwoods
NC volume = National Volume Estimator Library, region/forest 500 (California). Per-species equation numbers:
- OS 500WO2W108 · SP 500WO2W117 · DF 500WO2W202 · WF 500WO2W015 · IC 500WO2W081 · RF 500WO2W020 · PP 500WO2W122 ·
  RW 500WO2W211  (conifers + IC → 500WO2W* profile — the R5 West-side 2-point taper).
- MA 500DVEW361 · BO 500DVEW818 · TO 500DVEW631 · OH 500DVEW981  (hardwoods → 500DVEW* woodland D²H — the SAME
  DVEW/r4d2h family CI/UT/TT already implement; reuse r4d2h_vol with the 500DVEW coefficients).
⇒ chunk-8 port: (a) hardwoods = REUSE the shared DVEW/r4d2h (like CI's compute_volumes_ci! DVEW branch) with the
NC 500DVEW eqnums; (b) conifers = port the 500WO2W NVEL profile (West-side 2-point — check if it's the Flewelling
FW2 family jl already has, or a distinct NVEL call). Replace the CI-r4vol placeholder (compute_volumes!(::Klamath)
-> compute_volumes_ci!) with compute_volumes!(::Klamath) using these. Validate TCuFt/MCuFt vs live nct01.sum
(1990 TCuFt 1308/MCuFt 449). All VEQNNC measured; DVEW reusable, WO2W is the new piece.

## Chunk-8 volume — implementation scoping (2026-08-11):
- DVEW hardwoods (MA361/BO818/TO631/OH981): jl r4d2h_vol (Chojnacky INT-339 D²H) exists (TT/UT/CI/CR) but is keyed
  by species code — the NC R5 hardwood codes 361/818/631/981 are NOT in the existing table (those variants use
  066/475/998/133/065 etc.). ⇒ add the 500DVEW coefficients for the 4 NC hardwoods to r4d2h (or a NC table).
- WO2W conifers (OS/SP/DF/WF/IC/RF/PP/RW): "WO2W" is a DISTINCT NVEL eqtype vs jl's "FW2W" (Flewelling, cr_fw2_vol)
  and "MATW" (Matney, r4vol). VERIFY whether 500WO2W maps to an existing jl profile or is a new NVEL taper to port
  (check ForestVegetationSimulator/volume/NVEL for the WO2W equation). This is the substantial chunk-8 piece.
⇒ compute_volumes!(::Klamath) = DVEW branch (r4d2h + NC hardwood coeffs) + WO2W branch (verify/port), validate
TCuFt/MCuFt per-tree + .sum vs FVSnc_dbg (1990 1308/449). Replaces the CI-r4vol placeholder. All eqnums measured.

## Chunk-8 volume — REFINED (2026-08-11): NC 500* = Region-5 California NVEL, NEW to jl (not R4 reuse)
The NC DVEW codes 361/818/631/981 are NOT in jl's R4 r4d2h_vol1 branches (which cover 064/066/106/475/... R4
species). ⇒ NC's 500DVEW (hardwoods) + 500WO2W (conifers) are the REGION-5 (California) NVEL equations, distinct
from the R4 Matney/Chojnacky and R1 Flewelling jl already has. chunk-8 = a genuine NVEL R5 port: locate the WO2W +
R5-DVE coefficients in ForestVegetationSimulator/volume/NVEL (voleqdef/the R5 taper source), port the taper + the
D²H woodland eq for the 12 NC species, wire compute_volumes!(::Klamath), validate TCuFt/MCuFt vs FVSnc_dbg. All 12
eqnums measured (500WO2W108/117/202/015/081/020/122/211 conifers+IC; 500DVEW361/818/631/981 hardwoods). Substantial
but self-contained; the growth port is unaffected (volume is a reporting/.sum column, output-only).

## Chunk-8 volume — NVEL sources located (2026-08-11): WO2W=Flewelling(r6vol/fwinit), DVE=dvest.f
voleqdef.f confirms NC uses 500WO2W* (conifers) + 500DVEW* (hardwoods). The NVEL profile routines:
- WO2W → fwinit.f + r6vol.f/r6vol1.f (Region-6 FLEWELLING West-side 2-point taper) — the SAME Flewelling family
  jl's cr_fw2_vol already implements (jl uses it for CI/EM/BM/UT "FW2W"). ⇒ the WO2W conifer volume likely REUSES
  cr_fw2_vol with the WO2W (Flewelling-West) profile coefficients — verify the fwinit coeff set for WO2W vs FW2W.
- DVEW → dvest.f (the DVE woodland/D²H estimator) — region-5 species coefficients (361/818/631/981) in dvest.f.
⇒ chunk-8 is MORE reusable than feared: WO2W ≈ jl cr_fw2_vol (Flewelling) + WO2W coeffs; DVEW = port dvest.f R5
coeffs (small). compute_volumes!(::Klamath) = cr_fw2_vol(conifers, WO2W) + dvest(hardwoods). Measure the WO2W/DVE
coefficients from fwinit.f/dvest.f, validate TCuFt/MCuFt per-tree vs FVSnc_dbg (1990 1308/449). Doctrine-#5 reuse.

## Chunk-8 volume — FINAL scope (2026-08-11): genuine NVEL R5/R6 port (cr_fw2_vol INGY-only, needs West-side coeffs)
Tested: cr_fw2_vol("500WO2W202",...) returns 0 — jl's Flewelling is the INGY (Inland NW) + R2/R3 subregion set only
(_fw2_is_ingy || 22≤jsp≤29 gate). NC's 500WO2W is the Region-6 Flewelling WEST-SIDE profile (fwinit.f SHP_C2 west
coeffs), a DIFFERENT coefficient set + species map than jl has. ⇒ chunk-8 = a real NVEL port:
- WO2W conifers: extend cr_fw2_vol (or a nc_fw2) with the R6 West-side SHP/taper coefficients from fwinit.f + the
  NC species→Flewelling-jsp map (framework reused; coefficients new).
- DVE hardwoods (MA/BO/TO/OH): port dvest.f's DVE D²H estimator + the R5 coefficients for 361/818/631/981.
Substantial but SELF-CONTAINED + OUTPUT-ONLY (volume is a .sum reporting column; the growth port is complete and
unaffected). All eqnums measured; sources = volume/NVEL/{fwinit.f, r6vol.f, dvest.f}. Validate TCuFt/MCuFt per-tree
+ .sum vs FVSnc_dbg (1990 1308/449). This is the last NC chunk; growth (ch1-7) tracks live end-to-end.

## NC growth verdict (2026-08-11): cyc0 BIT-EXACT; multi-cycle ~2-5% = plausibly the cornered DGSCOR class
Aggregate .sum nct01 control: cyc0 all 6 density cols BIT-EXACT; multi-cycle BA ~2% high, QMD ~1-2%, TopHt ~5% low,
CCF ~12% low (CCF/TopHt downstream of the DG-distribution + height realization). The BA/QMD ~2% is the SAME
magnitude as CI #142's accepted cornered DGSCOR/ZZRAN straddle (~2% over-kill) — so NC growth is plausibly
bit-exact-or-cornered, NOT a gross bug. TO CONFIRM cornered-vs-small-bug (doctrine #1): per-tree deterministic DDS
jl-vs-FVSnc_dbg on ONE cyc1 growth-pass tree matched by (species, exact DBH) — gate the dgf DEBUG dump on the
growth pass (not calibration) to avoid the earlier mixed-pass alignment error; if the deterministic DDS matches,
the residual is the DGSCOR realization (cornered). The TopHt-5%/CCF-12% may need a 2nd look (height POTHTG or the
site-index HTCALC-conversion approximation — I hardcoded the DF=90-default SI() rather than porting the full
sitset.f DO-30 HTCALC conversion; exact for nct01, approximate if a SITECODE sets a non-default site species/index).
⇒ NC growth port: cyc0 bit-exact + multi-cycle tracks live in the cornered range. Remaining to fully close:
(1) per-tree DDS confirm-cornered, (2) full sitset HTCALC site-conversion (for SITECODE stands), (3) chunk-8 volume.

## Chunk-8 volume — implementing NVEL routines traced (2026-08-11), ready for a focused port:
- DVEW hardwoods (500DVEW361/818/631/981): dvest.f VOLEQ(1:1)='5' → CALL R5HARV (volume/NVEL/r5harv.f, Region-5
  D²H harvest-volume routine) → VOL(4) cubic, VOL(6)=VOL(4)/90 cordwood. Port r5harv.f + its R5 coefficients.
- WO2W conifers (500WO2W*): fwinit.f + r6vol.f/r6vol1.f (R6 Flewelling West-side taper). Port the West-side SHP
  coeffs (extend jl cr_fw2_vol's INGY set) + species map.
⇒ compute_volumes!(::Klamath) = r5harv-port (hardwoods) + cr_fw2_vol-extended (conifers). Both routines located;
this is a focused multi-file NVEL port (r5harv.f + fwinit.f + coeffs). OUTPUT-ONLY; growth port (ch1-7) is complete
+ tracks live. All eqnums measured, sources traced end-to-end (dvest→R5HARV, voleqdef→fwinit/r6vol). Validate vs
FVSnc_dbg .sum (1990 TCuFt 1308/MCuFt 449).

## Chunk 8 — Volume (2026-08-11): R5TAP total-cubic BIT-EXACT; merch-driver measured

**VEQNNC (12 species) CONFIRMED bit-exact** vs live oracle (nct01.out:68-70, from voleqdef.f R5_EQN):
OS 500WO2W108 · SP 500WO2W117 · DF 500WO2W202 · WF 500WO2W015 · MA 500DVEW361 · IC 500WO2W081 ·
BO 500DVEW818 · TO 500DVEW631 · RF 500WO2W020 · PP 500WO2W122 · OH 500DVEW981 · RW 500WO2W211.
⇒ 8 conifers+redwood = WO2W (Wensel-Krumland R5 profile, r5tap.f — NOT Flewelling); 4 hardwoods
MA/BO/TO/OH = DVEW (R5 CA-hardwood D²H, r5harv.f). nct01 = 100% conifers ⇒ only WO2W exercised.

**Dispatch (MEASURED, agent trace):** volinit.f MDL='WO2' → PROFILE → TAPERMODEL → R5TAP (profile.f:1368).
FWINIT is SKIPPED (VOLEQ(4:4)='W' not 'F') ⇒ JSP=0, no Flewelling coeffs. R5TAP gives inside-bark DIB;
DOB stays 0 (no bark in the taper). TCUBIC (stump 1-ft cyl + 4-ft Smalian + tip) = VOL(1). Merch = GETDIB
1-inch-class log bucking to MTOPP.

**R5TAP taper — coefficients (r5tap.f:16-55), VOLEQ(8:10)→SP 1-9:** DF202=1,PP122=2,SP117=3,WF015=4,
RF020=5,IC081=6,JP116=7,LP108=8,RW211=9. R5WKC(5,sp)=C1..C5, R5WKB(2,sp)=B1,B2. White fir (sp4) has the
TERM2≥-1 clamp. Formula: upper stem (htup≥4.499): DIB=DBH·(C1 - TERM2·log(1-TERM3·(1-exp(C1/TERM2)))),
TERM2=C3+C4·DBH+C5·TOTHT, TERM3=((htup-1)/(TOTHT-1))^C2; stump (htup<4.499): DIB=(1-B1)·DBH·exp(B2·(4.5-htup)).
All coeffs + the port in src/variants/klamath/volume.jl (NC_R5WKC/NC_R5WKB/nc_r5tap_dib).

**VALIDATION (per-tree vs standalone live R5TAP driver /workspace/.ncwork/r5drv + FVSnc treelist):**
- ★ TOTAL CUBIC VOL(1): **BIT-EXACT** every tree — SP D11.5 H73 =17.6, DF D12.7 H67 =20.9, SP D9.5 =10.0,
  DF D10 =12.8, WF D10.9 =15.9, RF D6.5 =2.7, big-SP D34.6 =277.4≈277.3. R5TAP DIBs match live to 4dp
  (4.5→9.9332, 17.5→8.4859, 28→7.2509, 38.5→5.9022 — jl identical).
- MERCH: jl VOL(4) to MTOPP=6" IB gives SP#1=12.6; live treelist 'MCH CU FT'=13.9 because the FVS
  reported merch = **MCF = VOL(4)+VOL(7)** (primary-to-6" + topwood 6"→4"), fvsvol.f:512, gated D≥DBHMIN(ISPC).
  MERLEN merch length == jl (_fw2_hs) to 2dp — the gap is the missing VOL(7) topwood, not the length.
- Board VOL(2): jl Scribner to 6"; live uses BFTOPD·BARK top (fvsvol.f:382) — bark-adjusted.

**.sum (control stand):** growth BIT-EXACT (TPA 536/BA 77/SDI 160/CCF 87/TopHt 63/QMD 5.1 == oracle all
cycles). Volume 1990: jl TCuFt 1261 / MCuFt 522 / BdFt 2249 vs oracle 1308 / 449 / 2000. (Old CI-r4vol
placeholder was 1540/833 — R5TAP is far closer.) TCuFt aggregate −3.6% despite per-tree bit-exact =
each_stand pre-expansion/mortality-tree summation nuance (CI #142 class; the big D34.6 mortality tree #5
carries 277 ft³). MCuFt +16% = jl uses VOL(4)-only without the DBHMIN(ISPC) gate (small trees D<9 that
live's fvsvol zeros still get merch in jl) — the DBHMIN gate DROPS jl toward 449; VOL(7) topwood adds back.

**REMAINING chunk-8 driver (fully scoped, MEASURED — port these to reach bit-exact merch/board):**
1. VOL(7) topwood (profile.f secondary-product loop, 6"→4" = MTOPP→MTOPS bucking) → MCF=VOL4+VOL7.
2. DBHMIN(ISPC) gate on merch (fvsvol.f:337,512) — find NC's per-species value (grinit=0 ⇒ a default
   applies; treelist shows the merch cutoff at D≈9).
3. Board recompute with BFTOPD·BARK top (fvsvol.f:362-383) — needs the NC volume BARK ratio.
4. Resolve the TCuFt aggregate via a REAL-run per-tree dump (not each_stand) — pre-expansion summation.
Oracle: /workspace/.ncwork/FVSnc_clean; standalone taper driver: /workspace/.ncwork/r5drv.f (r5tap.f linked).

### Chunk 8 UPDATE (2026-08-11 cont.): DBHMIN=9 gate LANDED; per-tree cubic bit-exact confirmed

**★ DBHMIN(ISPC)=9.0 FOUND + APPLIED** (nc/sitset.f:196-224, forest-default merch specs by IFOR):
  IFOR 4 (Siskiyou/R6): DBHMIN=9, TOPD=4.5. IFOR 5,7 (Simpson/BLM): DBHMIN=9, TOPD=5.0.
  DEFAULT (IFOR 1=Klamath 505, our nct01): **DBHMIN=9.0, TOPD=BFTOPD=SCFTOPD=6.0**, BFMIND=SCFMIND=9.0.
  Applied as the merch/board gate (fvsvol.f:337,512: MCF/board=0 for D<DBHMIN) in compute_volumes_nc!.
  ⇒ .sum MCuFt 1990: 522 → **427** (oracle 449; +16% → −4.9%). Matches the .out MERCH report (smallest
  merch quantile 9.4"). No effect on TCuFt (gate is merch-only).

**★ PER-TREE TOTAL CUBIC BIT-EXACT for ALL 27 trees** (not just large): measured jl nc_wo2w_vol vs FVSnc
  treelist across D=1.2..34.6 incl. every small tree — TPA-weighted deficit = 0.0. The live treelist LIVE
  trees sum to exactly 1308 = oracle TCuFt.

**TCuFt aggregate 1261 vs 1308 (−3.6%) EXPLAINED — NOT a volume bug:** per-tree cubic is bit-exact + TPA-
  weighted per-tree deficit is 0, so the gap is the each_stand-vs-real-run TPA normalization (each_stand
  returns BA=85 PRE-expansion; the real run correctly normalizes to BA=77 bit-exact — the CI #142 trap).
  Pinning the real-run per-tree TPA needs a real-run cycle-0 hook (each_stand is INVALID for this). The
  volume port itself is faithful. NOTE the pre-existing growth residual (TopHt ~5% low / BA drift by 1995:
  jl BA 101 vs oracle 96) also perturbs multi-cycle volume — a GROWTH refinement (#158 class), not volume.

**REMAINING (task #162):** VOL(7) topwood 6"→4" (profile.f:684-770 secondary-product bucking: MERLEN to
  MTOPS=4, LMERCH_topwood = LMERCH_to4 − LENMS, NUMLOG/SEGMNT/GETDIB with LOGST=NUMSEG_primary) → MCF=VOL4+VOL7
  (would raise 427 toward 449). Board recompute BFTOPD·BARK (needs NC volume BARK). Real-run TPA hook for TCuFt.

### Chunk 8 VERDICT (2026-08-11 final): NC volume FAITHFUL — all 3 cols track one TPA residual

Added the exact MERLEN (profile.f:1011-1052, tenth-inch-truncated 0.1-ft binary search) + nc_wo2w_merch
(primary VOL4 to MTOPP=6; optional topwood VOL7 to MTOPS=4 for SPFLG=1 harvest reports). The .sum uses
VOL4-only (SPFLG=0): VOL4+VOL7 overshoots to MCuFt 524 (vs oracle 449), VOL4-only = 431 — so the summary
merch is SPFLG=0 primary product (fvsvol.f:512 MCF=VOL4+VOL7 only when SPFLG=1).

**nct01 control .sum 1990 (growth BIT-EXACT — TPA 536/BA 77/SDI 160/CCF 87/TopHt 63/QMD 5.1 all == oracle):**
  col     jl    oracle   Δ
  TCuFt  1261   1308   −3.6%
  MCuFt   431    449   −4.0%
  BdFt   1939   2000   −3.1%
All three within ~3.7%, the SAME ratio ⇒ a single upstream cause (the each_stand-vs-real-run TPA/expansion
normalization: each_stand returns BA=85 pre-expansion, the real run normalizes to BA=77 bit-exact). The
volume EQUATIONS are faithful — VOL(1) total cubic is bit-exact per-tree (all 27 trees, deficit 0.0), and
the merch/board columns scale with the same TPA factor. ⇒ NC volume is BIT-EXACT-OR-CORNERED: per-tree
bit-exact; aggregate residual = the CI #142-class real-run TPA normalization (measure in the real run, not
each_stand) compounded downstream by the pre-existing TopHt-low growth residual. NOT a volume-port bug.

Remaining (low priority): board BFTOPD·BARK top (−3.1% already ≈ TPA residual ⇒ bark effect is minor here);
DVEW hardwood path unexercised on nct01 (ported faithfully, needs a hardwood stand to validate); the
real-run TPA hook to close the aggregate is a GROWTH/expansion item, not volume.

## GROWTH BUG FOUND (2026-08-11): NC large-tree DG over-grows EXACTLY 10/9 uniformly

While validating volume I found NC growth OVER-predicts on nct01: cyc0 1990 bit-exact, but 1995 BA jl 101
vs oracle 96 (+26% BA-growth: jl 77→101 vs live 77→96). MEASURED (real-run per-tree DG dump at the Klamath
branch, applied increment = diam_growth/bark, vs the live cycle-0 treelist DIAM INCR column — pre-tripling,
valid per-tree):

  tree        live DG   jl DG   ratio
  SP D11.5     1.00     1.111   1.111
  SP D9.5      1.10     1.222   1.111
  SP D9.6      0.50     0.556   1.111
  DF D10.0     1.00     1.111   1.111
  DF D12.7     1.60     1.778   1.111
  DF D9.4      1.80     2.00    1.111
  WF D10.9     1.00     1.111   1.111
  RF D6.5      2.30     2.556   1.111

jl DG = live × 1.111 (= 10/9) for EVERY tree — UNIFORM across the DEFAULT branch (DF sp3, WF sp4) AND the
special sp2/6/9 branch (SP sp2, RF sp9). A perfectly uniform multiplicative factor ⇒ NOT a per-tree/species
DG-equation coefficient error, NOT the bark (nc_bratio verified CORRECT vs live nc/bratio.f: SP 0.874, DF
0.824 — matches; and a bark error would vary by species), NOT the 5-vs-10-yr period (that would be 2×, and
would differ between the TDDS/2 specials and the no-/2 default). ⇒ isolated to a GLOBAL 10/9 multiplier in
the shared DDS→increment path (calibrate_diameter_growth! / the DGDRIV stochastic realization: XDGROW
variance/Jensen correction exp(σ²/2), or a dgscale). Candidate: the E[ln DDS]→E[DDS] variance correction or
a DGSD/period constant applied for NC. NEXT: instrument calibrate_diameter_growth! for NC — dump wk2 (ln DDS)
+ XDGROW + WK4 + the final sqrt, compare each factor to live DGDRIV; the 10/9 (ln=0.10536) is a constant
additive offset in ln(DDS)-space. Repro: scratchpad/nct01_ctl.key + the FVSJL_NC_DG_DEBUG dump (reinsert at
the Klamath dgcons branch, simulate.jl:118). This is DISTINCT from (and larger than) the TopHt-low residual;
it is the dominant NC multi-cycle divergence. Volume port is unaffected (VOL(1) bit-exact per-tree at cyc0).

## CORRECTION (2026-08-11): DG bug is PER-SPECIES (calibration), NOT uniform 10/9

The prior "uniform 10/9" finding was a MEASUREMENT ARTIFACT: the FVSJL_NC_DG_DEBUG probe read
t.diam_growth at the calibrate branch (simulate.jl:118) — but the REAL DG realization happens later
(diameter_growth! at simulate.jl:474). t.diam_growth at line 118 held STALE data. Re-measured IN the
realization loop (variants/southern/diameter_growth.jl:1138, computing dgc/bark), the ACTUAL per-tree DG
is WILDLY per-species (NOT uniform):
  SP (sp2) D11.5: jl ~1.08× live  (base DDS ~right)
  DF (sp3) D10.0: jl ~0.08× live  (base DDS=0.224/exp1.25, ~12× TOO SMALL — massive UNDER-grow)
  WF (sp4) D10.9: jl ~1.96× live  (DDS=3.386/exp29.5, ~2× TOO BIG — OVER-grow)
The aggregate BA (+5% by 1995) MASKED this by mixed-sign cancellation (DF under vs WF/SP over) — the
same "measure per-species, aggregate hides it" doctrine lesson as EM/CI.

ROOT (isolated by per-term dump of the DEFAULT branch): the base (uncalibrated) DDS is too small for the
default species (DF terms sum to 0.224; DGHAH/DGPCCF verified CORRECTLY 0 vs nc/dgf.f DATA:108-114 — NOT a
missing term). **Live CALIBRATES nct01** (nct01.out:155 "DBH GROWTH MODEL SCALE FACTORS WERE COMPUTED" —
the trees carry a measured past-DG F2.1 field in nct01.tre, e.g. SP D11.5 DG=1.0). Live's per-species COR
boosts the small base DDS to match the measured growth (~1.0"/5yr). jl's dg_cor[sp] comes out 0 in the
growth pass (DF) ⇒ the uncalibrated (too-small) base DDS is used ⇒ DF under-grows; WF's COR path over-shoots.
⇒ the bug is jl's NC DGSCOR CALIBRATION (per-species COR wrong/not applied), OR a base DGCON coefficient —
distinguish by comparing jl's pre-cal base DDS + dg_cor[sp] per species to live's ZNC calibration dump
(nct01.out:157-185) and the computed scale factors. NEXT: verify jl reads nct01.tre's measured past-DG into
the NC calibration (WK1), and that calibrate_diameter_growth! sets dg_cor[sp] for NC's default species.
Repro: reinsert the per-term probe at diameter_growth.jl:163 (Klamath default branch). SUPERSEDES the 10/9
framing in fcfdba3. NC growth NOT complete; volume unaffected (VOL1 bit-exact per-tree).

## DGCON DEFAULT-branch missing terms + calibration interaction (2026-08-11)

The nc_dgcons! DEFAULT branch (diameter_growth.jl) is `dgcon = DGFOR[sp,MAPLOC]` ONLY — but live
nc/dgf.f:478-485 DGCON(default) = DGFOR + DGEL2·ELEV² + (DGSASP·sinAsp + DGCASP·cosAsp + DGSLOP)·SLOPE
+ DGSLSQ·SLOPE² + **DGSITE·ln(SITEAR(3))**. The site term alone is DGSITE(3)·ln(90)=0.56356·4.4998=**+2.536**
for DF (SITEAR(3)=DF site index=90, confirmed sitset.f:120,152; used for ALL default species). So jl's
default-species base DDS is ~2.5 too LOW in ln-space (DF base exp 0.224/1.25 vs faithful ~2.87/17.7).
Coefficients captured (ready to apply): DGSASP/DGCASP/DGSLOP/DGSLSQ = nc/dgf.f DATA (per species):
  DGSASP[3,4]=-0.040708,-0.01560 · DGCASP[3,4]=-0.16836,-0.15630 · DGSLOP[3,4]=0.46468,0.58937 ·
  DGSLSQ[3,4]=-0.87145,-1.05045 · DGSITE (jl NC_DGSITE, already present) · DGEL2 (present, 0 for defaults).
  Formula: dgcon += DGEL2·elev² + (DGSASP·sin(asp)+DGCASP·cos(asp)+DGSLOP)·slope + DGSLSQ·slope² + DGSITE·log(SITEAR(3)).

⚠ APPLYING IT ALONE REGRESSES the .sum: 1995 UNCHANGED (101, the calibration pins cycle-1 growth to the
measured past-DG in nct01.tre regardless of base), but 2010 BA 171→191 (oracle 167) — WORSE. ⇒ a masked
interaction: jl's DGSCOR calibration cor was tuned against the too-low base; the faithful base + the same
calibration double-counts in later cycles. Since live HAS the site term and matches (167) while jl WITH it
overshoots (191), jl's calibration cor (or DGFOR value, or cor decay) differs from live's. REVERTED to
preserve the validated ~5%-over state (was assessed plausibly-cornered) rather than ship an unresolved
14%-over regression.

⚠ DOCTRINE LESSON (twice this session): per-species DG measured via a probe at the calibrate branch
(simulate.jl:118) OR during calibration passes is CONFOUNDED — calibrate_diameter_growth! calls dgf! MANY
times (backdated calibration iterations) with different ba/relden/cor. Measure the REALIZED DG only in the
growth-pass realization (diameter_growth! @ simulate.jl:474), and even there the cor differs from the base.
NEXT (needs a working live dgf debug — the DEBUG keyword SEGFAULTs FVSnc_clean/dbg; use a g16-instrumented
nc/dgf.f writing WK2 unconditionally, OR a standalone dgf driver): compare jl base DDS (pre-cor) + dg_cor[sp]
per species to live's, WITH the DGCON site term applied, to find the calibration/DGFOR discrepancy. THEN
apply the DGCON term + the calibration fix together. NC growth ~5% multi-cycle over is likely DGSCOR-class
(cornered like EM/CI) once the base is faithful; volume remains done (VOL1 bit-exact).

## DGCON fix VALIDATED BIT-EXACT via live dgf debug (2026-08-11) — APPLIED

UNBLOCKED the live dgf debug: the FVS `DEBUG <cycle> <nonblank>` keyword + a `DGF` supplemental record
debugs ONLY dgf (the bare `DEBUG` → DBALL segfaults). Output → fort.16. dgf.f:435 prints per-tree
I/ISPC/DBH/BAL/CR/RELDEN/BA/**LN(DDS)**. The 3rd pass (RELDEN 95.7/BA 85.1) is the growth pass (current DBH).

Applied the DGCON default-branch fix (nc_dgcons!: + DGEL2·elev² + slope/aspect + DGSITE·ln(SITEAR(3))).
**DF now BIT-EXACT vs live** (growth-pass LN(DDS), ba=85.1):
  DF D10.0  jl 2.794  == live 2.7940 ✓   DF D12.7  jl 3.0123 == live 3.0123 ✓
  DF D10.4  jl 3.0428 == live 3.0428 ✓   DF D9.4   jl 2.7204 == live 2.7204 ✓
The +2.570 the port was missing = exactly DGSITE(3)·ln(90)+slope = 2.536+0.034. Fix is CORRECT & KEPT.

**Remaining: WF (sp4) calibration COR** — jl WF D10.9 LN(DDS)=3.3406 (cor=+0.0591) vs live 2.7362. The WF
BASE matches (computed from live coeffs = 3.282 == jl 3.2815; all sp4 coeffs DGLD/DGCR/DGCRSQ/DGDBAL/DGBA/
DGSITE/DGSASP.. verified == nc/dgf.f DATA). So the 0.60 gap is PURELY the cor: live WF cor = 2.7362−3.282 =
**−0.546**, jl computes **+0.059**. DF cor matches (both ~0). ⇒ jl's DGSCOR calibration produces a wrong
per-species cor for WF (reads WF measured-DG or the residual/shrinkage differently). The DGCON fix EXPOSED
this pre-existing WF-cor bug (previously masked: the too-small base + a compensating cor netted ~5%; now the
faithful base + the wrong cor nets ~14% over). .sum regressed to 14% over TEMPORARILY — will resolve when the
WF cor is fixed. NEXT: instrument calibrate_diameter_growth! for NC WF — compare jl's measured-DG (WK1) +
residual + dg_cor[4] to live's calibration (nct01.out ZNC dump); DF works so it's WF-species-specific.
Live dgf debug recipe (durable): /workspace/.ncwork/nctree.key (DEBUG 1./1. + DGF record) → fort.16.

## NC PSIGSQ missing-branch FIXED (2026-08-11) — CI/IE-class

NC was MISSING from the shared PSIGSQ dispatch (southern/diameter_growth.jl:538) → fell through to the SN
default DG_PSIGSQ=0.0898. Real gap (same class as the CI 3b9aa35 / IE 96cde22 fixes). Added NC_PSIGSQ from
nc/dgdriv.f:95 DATA: [0.0408,0.0586,0.1556,0.0970,0.0858,0.1433,0.0636,0.0970,0.0970,0.0636,0.0858,0.0898].
WF(sp4)=0.0970 vs 0.0898. Wired `s.variant isa Klamath ? NC_PSIGSQ[sp]`. NC-only, no cross-variant regression.

⚠ But it did NOT resolve the WF cor gap (.sum still ~14% over): the shrunk cor `corv = wc·cornew` follows the
SIGN of the RAW cornew, and PSIGSQ only scales the shrinkage weight `wc∈[0,1]` — it cannot flip +0.059→−0.546.
So the WF divergence is the RAW cornew (the calibration regression residual): jl WF cornew is ~+ (measured≈base)
while live's is ~− (measured<base). Since the growth-pass BASE now matches live (DF bit-exact, WF base 3.282==),
the raw-residual difference is in the CALIBRATION pass — either jl reads WF's measured past-DG (the .tre F2.1
field) differently, or the backdated-DBH base at the calibration state differs. DF cornew≈0 matches (so it's
WF-specific, likely a per-species measured-DG or backdated-state handling). NEXT: instrument
calibrate_diameter_growth! for NC — dump per-WF-tree measured-DG (WK1) + backdated base + the accumulated
cornew, compare to live (nct01.out ZNC dump / an instrumented nc/dgdriv.f). The DGCON + PSIGSQ fixes are
faithful and KEPT; the .sum will converge once the WF calibration raw-cornew is corrected.

## NC bark cache = real nc_bratio (was 0/0.9 placeholder) — faithful; WF cornew still the blocker (2026-08-11)

nc_dgcons! set c.bark_a/bark_b to a 0/0.9 PLACEHOLDER (constant 0.9) with a comment "NC uses nc_bratio
directly in dgf!". But the SHARED code uses this cache: the DDS→increment realization (dbh += dg/bark,
simulate.jl:556) AND the calibration term (2·bark·wk3). Live (dgdriv.f:205) uses BRATIO(ISPC,D,HT) = the
real varying nc bark. FIXED: encode nc/bratio.f into the linear cache (eqtype 1: bark_a=−a,bark_b=1−b;
eqtype 2: bark_a=a,bark_b=b; sp12 RW eqtype 3 POWER left as 0.9 — RW-only, no nct01). Verified: (bark_a+
bark_b·D)/D == nc_bratio (DF@10 0.8235, WF@10.9 0.8765). Faithful — matches live BRATIO.

⚠ This barely moved the WF cornew (0.059→0.039) — confirming the WF divergence is NOT bark, it's the
calibration RAW RESIDUAL. And because the correct (lower) bark increases the realized increment (OB_incr =
sqrt(d²+dds/bark²)−d), it AMPLIFIES the still-unfixed WF over-growth ⇒ .sum 2010 191→201 over (oracle 167).
⇒ THE .sum WILL NOT CONVERGE until the WF calibration cornew is fixed; each faithful fix (DGCON/PSIGSQ/bark)
UNMASKS/amplifies it. This is the doctrine's "faithful fix regresses ⇒ masked bug" — the masked bug is the
WF (sp4) calibration raw cornew: jl +0.039 vs live −0.546, WF base matches live (3.282), DF cornew≈0 matches.
NEXT (the ONE remaining NC-growth bug): instrument calibrate_diameter_growth!'s residual accumulation for NC
WF — the measured-DG (WK1/`dg`) + backdated DIB (`wk3`) + reslog=log(dg·(2·bark·wk3+dg))−wk2 per WF tree vs
live (an instrumented nc/dgdriv.f writing the per-tree residual, same DEBUG-DGF technique). DF works ⇒ it's a
WF-species-specific measured-DG or backdated-state handling diff. All of DGCON+PSIGSQ+bark are FAITHFUL & KEPT.

## ★★ ROOT CAUSE — calibration SCALE=0.5 (measured-DG 10yr → model 5yr) FIXED (2026-08-11)

Instrumented live DGDRIV (DEBUG + `DGDRIV` supplemental record → fort.16, dgdriv.f:449 prints per-tree
OBS.DG/TERM/RESLOG + the SNX/SNY sums). Found: DF has FN=4 measured-DG trees (<FNMIN=5 ⇒ NOT calibrated,
cor=0 — matches jl) but WF has FN=5 (calibrated, cor=−0.5454). jl's WF FN=5 too, BUT jl's calibration TERM
was EXACTLY 2× live's for every WF tree (jl 17.462 vs live 8.731, etc.) ⇒ RESLOG wrong sign ⇒ SNY +6.0 vs
live −53.5 ⇒ WF cor +0.04 vs −0.55.

ROOT: NC's DG MODEL basis is 5-yr (blkdat DATA YR/5.0/) but the measured past-DG in nct01.tre is a 10-YEAR
measurement (nct01.out:12879 "TALLY 2 AT 10 YEARS"). Live scales the measured DDS by SCALE=YR/FINT=5/10=0.5
(dgdriv.f:419,328); jl used scale=1. FIX: NC calibrate scale = 0.5·dgscale (simulate.jl Klamath branch).
Other western variants have YR=10 + 10-yr measurement ⇒ scale 1; NC is the unique YR=5-with-10yr-measurement.

RESULT (nct01 control, all 4 NC-DG fixes DGCON+PSIGSQ+bark+scale): WF SNY −54.92 == live −53.51; .sum
1990 BIT-EXACT, then 1995 BA 94/96, 2000 115/119, 2005 138/142, 2010 162/167 — all within ~2-3% (was 14%
OVER before scale). Early cycles now bit-exact-or-cornered. Later cycles (2015+) drift to ~+5-9% BA / −SDI
(mixed-sign compounding, the accepted DGSCOR/tie-break tail like EM/CI) — residual, not the dominant bug.
The DEBUG-DGDRIV technique + the 4 fixes together resolve the NC large-tree DG. NC growth now tracks live.

## MULTI-STAND FIA SWEEP (2026-08-11) — the real validation gate; scale refinement fixes 4/7 divergences

Ran a 12-stand NC FIA sweep (scratchpad/nc_sweep.jl, live oracle FVSnc_clean via DATABASE keyword). Initial
(with the hardcoded 0.5·dgscale): bit-exact=0, cornered(≤3%)=5, DIVERGED(>3%)=7, 0 crashes. ⇒ nct01
single-stand was INSUFFICIENT (doctrine: multi-stand FIA is the real gate).

SCALE REFINEMENT (growth_dg_set ? dgscale : 0.5): the hardcoded 0.5·dgscale DOUBLE-SCALED any FIA stand that
carries a GROWTH card (which supplies the remeasurement FINT ⇒ dgscale=yr/dfint IS already YR/FINT_meas). The
no-GROWTH default keeps the 10-yr-measurement 0.5. nct01-inert (no GROWTH ⇒ still 0.5, bit-exact unchanged).
Re-ran the 7 divergent CNs — FIXED 4: 850447806 6.2%→0.4%, 248613816 18%→2.8%, 248615564 4.4%→0.8%,
1123874415 3.2%→0.0%. ⇒ NC now 9/12 cornered.

REMAINING 3 divergences (task #164):
- cn 850447807 jl BA=1 vs live 59 (98%): DEGENERATE near-total growth failure (TPA 2390 present but BA=1 ⇒
  the trees exist but don't grow — likely a species/site edge case, a NaN, or all-large-trees-excluded). REAL bug.
- cn 1288130126 +8.4% DG-over (TPA 1598 vs 1474): DG over-prediction / minor mortality; species mix ≠ nct01 conifers.
- cn 1123874220 jl TPA 1168 vs live 625 (+27.6% BA): MORTALITY UNDER-KILL ~2× — NC Zeide self-thin under-kills,
  the BM #140 class. Distinct from DG.
NEXT: instrument the degenerate BA=1 stand (per-tree DG dump — why zero growth) + the mortality-2× stand
(DEBUG-DGDRIV/self-thin vs live). NC species beyond nct01's DF/WF/SP/RF need per-tree DG validation.

## RESOLUTION (2026-08-11) — species-crosswalk completion + BAMAX residual-BA cap

Two SYSTEMIC fixes drove the 3 remaining divergences down; NC FIA now ~10/12 cornered (0 crashes).

### Fix 1 — species_translation.csv completeness (commit d0d40a9)
jl's NC crosswalk (data/klamath/species_translation.csv) carried ONLY the 12 PRIMARY species' FIA codes;
all other FIA species (150) DEFAULTED to OS (softwood) instead of live's spctrn.f NC mapping. Parsed
spctrn.f (NC = COLUMN 13: CASE('NC') SPCOUT=ASPT(I,13); two DATA blocks, FIA@col2 in J=1-10, NC@col13 in
J=11-21, aligned by I) and appended the 150 missing FIA→NC rows. VALIDATED:
- cn 1288130126 (801 coast-live-oak → BO, was OS): BA 181→172 vs live 167 (8.4%→3.0%, CORNERED).
- cn 850447807 (768 → OH, was OS): BA 1→36 vs 59 (98%→39%; the mis-mapped SEEDLING now grows).
META: MEASURED MA(sp5) DG bit-exact FIRST (ruled out the assumed madrone-DG culprit) → the real bug was the
crosswalk. Systemic — every hardwood/mixed-species FIA stand had been mis-grown as softwood.

### Fix 2 — BAMAX residual-BA mortality cap (klamath/mortality.jl) — the redwood self-thin root
nc/morts.f header: "SDI-BASED MORTALITY IS USED AS LONG AS QMD < 10 INCHES, AT WHICH TIME BAMAX-BASED
MORTALITY TAKES OVER. IF NOT SET BY THE USER, BAMAX IS DETERMINED FROM MAX SDI AT 10 INCH DBH." Mechanism:
- morts.f:308 CALL SDICAL(0,SDIMAX) → vbase/sdical.f:203-204 (.NOT.LBAMAX branch): BAMAX = SDIMAX·0.5454154·PMSDIU.
  (PMSDIU=85 in grinit is converted 85→0.85 at morts.f:167 BEFORE the SDICAL call, so BAMAX uses 0.85.)
- morts.f:685-754: after the SDI/Zeide self-thin AND size-cap loops, if residual BANEW > BAMAX+1, scale ALL
  per-tree mortality up by ADJFAC=(BANEW-BAMAX)/BADEAD, iterate ≤100× until residual BA ≤ BAMAX.
jl had NO such cap → dense stands (redwood, SDIMAX~1000) under-killed ~2×. Ported the loop into
mortality!(::Klamath) after the per-tree kill loop, before apply_fixmort!. Inert (immediate break) when BA is
already ≤ BAMAX, so it CANNOT touch below-cap stands. VALIDATED cn 1123874220 (99% RW1): BA 597→464 vs live
468 (27.6%→0.9%, CORNERED); the 2 cornered control stands unchanged (no regression). Residual: TPA 909 vs 625
(QMD 9.7 vs 11.7) — finer per-tree kill realization; BA (the cap variable + classified metric) meets the bar.

★ CROSS-VARIANT LEAD: the identical BANEW-BAMAX cap block exists in bm/em/ut/tt morts.f (all CALL vbase/sdical.f,
all convert PMSDIU 85→0.85 before SDICAL) but is ABSENT from the shared jl mortality (southern/mortality.jl).
Latent under-kill on any >BAMAX dense stand — matches the BM #140 "under-thin on actively-self-thinning stands"
signature exactly. NEXT CHUNK: port the cap into the shared mortality apply-path (after MSB + size-cap), validate
per-variant vs FVSbm/em/ut/tt_clean on a dense self-thinning stand each.

### Remaining NC divergences
- cn 850447807: 39% residual (seedling regent small-tree tail beyond the crosswalk fix).
- cn 504618389: 3.0% borderline (accepted DGSCOR/tie-break tail).
