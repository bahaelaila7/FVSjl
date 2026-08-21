# WWPB synthetic-PPE-harness port plan (USER-approved 2026-08-21)

Goal: an actual end-to-end WWPB (Westwide Pine Beetle) outbreak in FVSjl. Beetle
kernels port bit-exact vs pristine wwpb/*.f drivers; the thin absent outer layer
(PPMAIN/ALSTD1/ALSTD2/SPLAEX/SPLAAR/GPGET/GPNEW) is faithfully reconstructed from
the archive/PPEcommons/*.F77 common-block interfaces. See [[fvsjl-wwpb-ppe-oracle-absent]].

## COMPLETE ARCHITECTURE (mapped 2026-08-21 from pristine wwpb/*.f)

### The PPE call chain (what the absent harness does)
```
PPMAIN  (absent — reconstruct)  ── per master cycle:
  ├─ per stand:  BMSDIT              ← FVS→BM tree BRIDGE (exists, bmsdit.f)
  ├─ once:       BMDRV               ← per-YEAR outbreak loop (exists, bmdrv.f)
  └─ per stand:  BMKILL              ← BM→FVS mortality HANDBACK (exists, bmkill.f)
BMSETP (exists, bmsetp.f) once per cycle BEFORE PPMAIN's stand loop:
  ├─ GPGET(301,…)  → IBMYR1 outbreak start year (from DISPERSE kw)   [scheduler — reconstruct]
  ├─ build stand list BMSDIX/BMSTDS, BMSTND count
  ├─ SPLAEX(…)     → per-stand spatial location+area                 [spatial — reconstruct]
  └─ LBMSPR = (IRC==0)   ← gates the whole outbreak (bmdrv.f:35 `IF LBMSPR .AND. IBMYR1>0`)
```

### BMSDIT — FVS→BM tree bridge (bmsdit.f, exists, ground-truth-portable)
Bins the FVS treelist into NSCL=10 DBH size classes × {host=1, nonhost=2}:
- per FVS tree I (via species-index ISCT(ISPC,1..2) + IND1(II)): K=BMDBHC(DBH(I));
  LX = host? via HSPEC(PBSPEC,ISPC). Accumulate BA=DBH²·(π/576)·PROB, TREE+=PROB,
  HTS+=HT·PROB, CRS+=ICR·PROB, HGS+=HTG·PROB, TVOL+=CFV·PROB into [K, host/nonhost].
- ISPH(,1/2) = species with max host/nonhost BA; IQPTYP = ISPFLL(ISPH).
- then HTS/CRS/HGS/TVOL /= TREE (size-class averages). OTPA = initial TREE snapshot.
- MICYC==2: seed PBKILL/ALLKLL from inventory damage LBMDAM(I).

### BMDBHC — leaf (bmdbhc.f, TRIVIAL): INDEX=1; DO I=1,NSCL-1: if DBH<UPSIZ(I) break; INDEX++.

### BMDRV — per-year loop (bmdrv.f) years IBMYR1..IBMYR2, gated LBMSPR & IBMYR1>0:
  mgmt: BMPHER,BMAPH,BMPSTC,BMSALV,BMSANI,BMSMGT → BMDRGT(drought,stoch 2 draws)
  stand-loop-1 {STOCK only}: BMCWIN(windthrow),BMLITE(stoch),BMFIRE→BMFMRT,BMOBB,
     BMDFOL,BMQMRT,BMMORT(.FALSE. fast),BMCGRF(→OLDGRF),BMCBKP(OLDGRF),BMCNUM
  BMATCT (landscape attractiveness + BKP redistribution; calls SPLAAR)
  stand-loop-2 {STOCK only}: BMIPS(stoch, if PBSPEC==3||IPSON) | BMISTD(stoch, fills
     PBKILL), BMOUT, BMMORT(.TRUE. beetle kills), BMAGDW(age dead-wood pools)

### BMKILL — BM→FVS handback (bmkill.f): TPBK ledger → per-record WK2(I) mortality,
  bounded PROB(I)−WK2(I) ≥ 1e-6; SDWP/DDWP dead-wood; SVMORT. (feeds FVS growth.)

## RECONSTRUCT (thin outer layer — no pristine source, use archive/PPEcommons/*.F77)
- **GPGET(act,iyr,…)/GPGET2(act,iyr,maxprm,nprms,prms,mxstnd,scnt,mylst,lok)** — PPE
  activity scheduler *get*: return scheduled params for activity code (301=DISPERSE,
  303..323 mgmt, 305 windthrow, 307 qmrt, 308 salv, 310/311 lite/fire…). FVSjl has an
  Event Monitor; the GP-scheduler is its PPE analogue. Reconstruct a small activity
  table keyed (activity, year) → PRMS, populated by the BMPPIN keyword block.
- **GPADD/GPNEW(kode,idt,act,nparms,prms[,j,mylst])** — schedule/reschedule an activity.
- **SPLAEX(bmstds,bmstnd,mylist,irc)** — load per-stand spatial location; **SPLAAR(istd,
  area,irc)** — return stand area. Single-stand degenerate: area = stand's EXPAND/acres,
  location = a single point; IRC=0 ⇒ LBMSPR=T. (spatial dispersal between stands is moot
  at MXSTND=1 — BMATCT's redistribution self-loops.)
- **PPMAIN top loop + BMSETP wiring** — reconstruct as an FVSjl landscape seam.

## PORT ORDER (dependency-bottom-up; each kernel driver-validated vs pristine hex)
0. ✓ RNG wwpb_rand!/seed (done) + BMIN block (done) + defaults (done).
   ✓ CHUNK 1-2 DONE (215d2b22): wwpb_dbh_class (BMDBHC) + WwpbStand state + bmsdit! (FVS→BM bridge). 57 tests.
   ✓ CHUNK 3a DONE (3691658c): bmmort! (fast/slow tree decrement + FASTK/TPBK ledgers).
   ✓ CHUNK 3b DONE (e569de35): bmcgrf! (GRF/GRFSTD/RVDNST) BIT-EXACT vs gfortran-16 driver (glibc expf/powf).
     Driver-golden recipe PROVEN for transcendental beetle biology. driver_bmcgrf.f reusable template.
  ✓ CHUNK 3c DONE (68726505): wwpb_init_coeffs (BMINIT MSBA/UPBA/INC) — 31 values BIT-EXACT (driver_bminit.f);
    also fixed WWPB_PI24 to Float32 (was 1-ULP off) ⇒ bmsdit BA now bit-exact too.
  ✓ CHUNK 3d DONE (bmcbkp): BKP brood core BIT-EXACT (driver_bmcbkp.f, GPGET2/GPADD stubbed → normal path).
    ⇒ the GRF→BKP reproduction chain is now bit-exact end-to-end.
  ✓ CHUNK 3e DONE (bmcspt): special-tree proportion BIT-EXACT (measured; overlap correction is DEAD code — JK loop
    never fires ⇒ SPCLT=clamp(ΣSP,1)).
  ✓ CHUNK 3f DONE (bmcnum): attractiveness numerator NUMER/TFOOD BIT-EXACT (driver_bmcnum.f). GRF→BKP→NUMER bit-exact.
  ✓ CHUNK 3g DONE (bmatct_single): single-stand BKP saturation BIT-EXACT (driver_bmatct.f). GRF→BKP→NUMER→BKP(sat) exact.
    NEXT: BMISTD (stochastic PBKILL from BKP, BMRANN stream) → BMKILL handback → harness (BMSETP/PPMAIN) + simulate.jl seam.
   → NEXT: the transcendental kernels need gfortran-16 driver-goldens (Float32 exp/logistic bit-exactness):
     BMFMRT (fire mort logistic), BMCGRF→BMCBKP→BMCNUM (susceptibility/BKP/attractiveness core).
1. State design: WwpbLandscape Julia struct mirroring BMCOM/BMFCOM/BMPCOM (MXSTND-dim
   arrays; single-stand first). BMDBHC (trivial).
2. BMSDIT bridge + HSPEC host-species matrix + ISPFLL quality-pool map (from bmblkd*.f).
3. Deterministic leaf kernels (driver-validate each): BMFMRT (fire mort logistic),
   BMQMRT, BMCGRF→BMCBKP→BMCNUM (GRF/BKP/attractiveness), BMMORT (tree decrement),
   BMCWIN (windthrow), BMOBB, BMDFOL, BMAGDW.
4. BMATCT (landscape attractiveness; SPLAAR).
5. Stochastic kernels (exact BMRANN stream + call order): BMDRGT(2), BMLITE(2),
   BMIPS(1), BMISTD(2). ← the RNG order is load-bearing; MEASURE both sides.
6. Management: BMPHER/BMAPH/BMPSTC/BMSALV/BMSANI/BMSMGT (scheduler-driven).
7. BMKILL handback + BMSETP + PPMAIN loop + BMPPIN kw (DISPERSE/HOST/PBSPEC/RANNSEED).
8. simulate.jl landscape seam; wire mortality into FVS growth; end-to-end single-stand
   pine-host DISPERSE run. Corner the composition (self-referential harness).

## DRIVER-GOLDEN RECIPE (per kernel, mirrors driver_bmrann.f)
gfortran-16 -std=legacy -w -fno-automatic -I/workspace/ForestVegetationSimulator/wwpb \
  -I<variant>/common (PRGPRM) -Iarchive/PPEcommons (PPEPRM) driver_<k>.f wwpb/<k>.f … \
  set the common-block inputs, CALL <K>, print outputs as Float32 hex (Z8.8). Compare
  to the Julia port bit-for-bit. Commons all present in wwpb/*.F77 + archive/PPEcommons.

Gate every chunk: multicycle 339/11. Additive/inert until the seam lands.

## FOUNDATIONAL DATA (extracted 2026-08-21)
- **HSPEC host designations** (bmblkd*.f): shipped DEFAULT = all 0 (no host) — set by the
  HOST keyword in BMPPIN. The commented "natural" defaults (bmblkdni.f) confirm the model
  doc: MPB(PBSPEC=1) host = **sp7 LP** (lodgepole); WPB(PBSPEC=2) host = **sp10 PP**
  (ponderosa); Ips(PBSPEC=3) host = sp7 LP. IE 11-sp order: WP L DF GF WH C LP S AF PP OTH.
  ⇒ FVSjl: wire HOST kw → HSPEC(pbsp, sp); default the natural LP/PP hosts when DISPERSE
  is given without HOST.
- **ISPFLL** (falldown-rate 1fast/2med/3slow, IE 11-sp): /2,3,2,1,2,3,2,1,2,1,2/.
- **UPSIZ**=/3,6,9,12,15,18,21,25,30,50/, **WPSIZ**=/10,20,60/, **ISCMIN**=/3,3,2/,
  seed 55329, NBGEN=1, NIBGEN=2, PFSLSH=0.9, IPSON=F/IPSMIN=2/IPSMAX=5. (all in wwpb.jl
  defaults already except HSPEC/ISPFLL — add those next.)
- Per-variant bmblkd: bm(generic)/ca/cr/nc/ni/so/wc. NI=IE-family. BA in bmsdit uses
  π/576 (=PIE/(24·24), PI24) — the FVS BA constant, per-tree DBH²·PI24·PROB.

## CONTINUATION POINT (next session)
Start step 1-2: add WWPB_HSPEC + WWPB_ISPFLL data + wwpb_dbh_class to wwpb.jl; design the
WwpbLandscape struct (single-stand MXSTND=1 first); port BMSDIT (the FVS→BM bridge) reading
FVSjl's treelist (dbh/ht/icr/htg/cfv/prob/species) into the size-class×host/nonhost table;
unit-test BMDBHC + the binning. THEN step 3 (deterministic kernels, driver-validated).
All additive/inert (no simulate.jl seam) until step 8 — gate 339/11 must hold each chunk.

## CONTINUATION (next session) — the BKP outbreak-dynamics core
CHUNK 3b (bmcgrf) done + BIT-EXACT. NEXT dependency chain for the outbreak:
  • BMINIT (bminit.f) — the coefficient SETUP: MSBA(isiz)=MID²·(π/576) [size-class-midpoint BA];
    INC(1,i)=RSLOPE·DBHMID+B clamped to REPMAX, INC(2,i)=INC(1,i), INC(3,i)=INC(1,1)·0.1 [Ips].
    RSLOPE/B/REPMAX are keyword/default params (bminit.f). Port this FIRST — bmcbkp/bmcnum/bmistd all read INC/MSBA.
  • BMCBKP (bmcbkp.f) — BKP "brood" from last year's PBKILL: BKP += MSBA·PBKILL·INC(pbspec,isiz) + STRIP + FINAL;
    ×REPRD (gen mult from NBGEN). Bad-year branch = GPGET2(317) scheduler stub (LBAD=false default ⇒ normal path).
    Needs INC/MSBA (BMINIT) + FINAL/STRIP fields + BKP/BKPIPS/OLDBKP + a bad-year table stub. Driver-golden-able.
  • BMCNUM (bmcnum.f) — attractiveness numerator; BMATCT — landscape BKP redistribution (SPLAAR); BMISTD — fills
    PBKILL from BKP (stochastic, BMRANN). Then the harness (BMSETP/PPMAIN/GP-scheduler/SPLA) + BMKILL handback + seam.
Driver-golden recipe: scratchpad/wwpb/driver_bmcgrf.f is the reusable template (INCLUDE the wwpb/*.F77 + ie/common
PRGPRM + archive/PPEcommons PPEPRM; set the BMCOM inputs; CALL; print Z8.8; route Float32 exp/pow via glibc).

## CONTINUATION (2026-08-21) — the BMISTD stochastic kill kernel + BMCBET beta dist
CHUNK 3g (bmatct_single) done. NEXT chain = the actual beetle kills:
  • BMCBET (bmcbet.f) — beta-distribution weights BETA(NSCL) over size classes, keyed on ABETA (a function of
    BKP/acre). **DEPENDS on ALNGAM (log-gamma, AS245) + BETAIN (incomplete beta, AS63) which are ABSENT from the
    tree** (only DECLARED `REAL ALNGAM, BETAIN` in bmcbet.f — like the PPE harness, reconstruct-category). FVSjl's
    LPMPB mpb_betin (betin.f AS-package, Float64) is a DIFFERENT algorithm ⇒ not bit-identical. APPROACH: port the
    canonical AS245 alngam + AS63 betain in BOTH Julia and as Fortran stubs, then driver-validate BMCBET's OWN
    weight-computation logic bit-exact (jl vs Fortran-with-same-stubs) — the AS functions are the shared
    reconstruction (documented; no pristine oracle exists for them). ABETA rules (bmistd.f:199-204):
    BKP>6→15; ≤6&>3.6→2.5+4.46·(BKP−3.6); ≤3.6&≥1.6→1.2+0.65·(BKP−1.6); <1.6→1. B=2.0 fixed. BETA(i)=BETAIN(x_i)−
    BETAIN(x_{i-1}) over [MINSIZE−0.5, MAXSIZE+0.5]; BETA init 1e-5.
  • BMISTD (bmistd.f) — STOCHASTIC (BMRANN=wwpb_rand!) kill allocation: SPRAY reduce → find MXISIZ/ISIZ1/ABETA →
    BMCBET → special-tree kills (BMRANN loop) → "group-kill" (deterministic bulk, DO 555 while TOTKL>2·MXISIZ) →
    "individual-kill" (BMRANN loop, DO 888) → PBKILL proportion→TPA + dead-wood. Fills PBKILL/PITCH/STRIP/FINAL.
    The exact BMRANN call order is load-bearing. SAREA from SPLAAR (single-stand = the FVS stand area). Driver-
    validatable (seed the RNG, set TREE/BKP/GRF/MSBA/SPCLT, compare PBKILL/PITCH/STRIP/FINAL) once BMCBET lands.
  • Then BMKILL (bmkill.f) TPBK→WK2 handback, BMSETP/PPMAIN single-stand harness, simulate.jl seam.
Doctrine note: BMCBET/AS-functions are reconstruct-category (source absent, like the harness — USER-approved);
every other kernel so far is bit-exact-vs-pristine. BMISTD's logic IS bit-exact-vs-pristine given the same BETA.


## CHUNK 3h DONE (3001ea02) — CORRECTION: BMCBET is BIT-EXACT, not reconstruct
The prior note (ALNGAM/BETAIN absent) was WRONG — both AS functions are DEFINED inside bmcbet.f (lines 136 AS245
log-gamma, 246 AS63 incomplete beta; I'd read only 120 of 323 lines). Ported _wwpb_alngam (AS245) + _wwpb_betain
(AS63, ACU=0.1e-14) + bmcbet! — BIT-EXACT vs the self-contained pristine bmcbet.f (ABETA 15/1..5 + 2.5/1..8, all
10 BETA match). ⇒ NINE kernels + AS funcs, ALL bit-exact-vs-pristine (nothing is reconstruct-category after all).
NEXT: BMISTD — now fully unblocked (has BETA via bmcbet!, BMRANN via wwpb_rand!, MSBA/GRF/TREE/SPCLT). The stochastic
kill kernel: SPRAY reduce → MXISIZ/ISIZ1/ABETA → bmcbet! → special-tree kills (BMRANN) → group-kill (deterministic
bulk) → individual-kill (BMRANN) → PBKILL/PITCH/STRIP/FINAL + dead-wood. Driver-validate (seed RNG, set TREE/BKP/
GRF/MSBA/SPCLT, compare PBKILL/PITCH/STRIP/FINAL). Then BMKILL handback + BMSETP/PPMAIN single-stand + simulate.jl seam.


## CHUNK 3i DONE (e2b64e46) — BMISTD stochastic kill kernel BIT-EXACT ⇒ BEETLE BIOLOGY COMPLETE
bmistd! validated bit-exact (driver_bmistd.f) on BOTH the individual-kill (BKP=3) AND group-kill (BKP=40) paths —
PBKILL/PITCH/STRIP/BKP/FINAL all match, the BMRANN stream order is correct. ⇒ ALL 10 beetle kernels + AS245/AS63
are bit-exact-vs-pristine. NEXT (the only remaining biology→FVS + orchestration):
  • BMKILL (bmkill.f) — the BM→FVS handback: TPBK ledger → per-record WK2(I) mortality, bounded PROB−WK2≥1e-6;
    SDWP/DDWP dead-wood; SVMORT. Reads the per-FVS-tree treelist (needs the FVS→BM tree-index map from bmsdit).
  • BMSETP + PPMAIN single-stand harness (reconstruct thin outer layer): BMSETP gets DISPERSE→IBMYR1, SPLAEX
    (single point/area, LBMSPR=T), BMSTND=1. PPMAIN loop: per cycle → bmsdit! (bridge) → bmdrv per-year loop
    (management stubs → bmdrgt drought → stand-loop-1 {bmcwin/bmfire/bmmort-fast/bmcgrf/bmcbkp/bmcnum} → bmatct_single
    → stand-loop-2 {bmistd → bmmort-slow → bmagdw}) → bmkill handback.
  • simulate.jl landscape seam (single-stand MXSTND=1); DISPERSE/HOST/PBSPEC keyword block (BMPPIN).
Note: many bmdrv sub-calls (bmpher/bmaph/bmpstc/bmsalv/bmsani/bmsmgt mgmt, bmobb/bmdfol other-agents, bmfire/bmfmrt
fire, bmcwin windthrow, bmlite lightning, bmagdw aging) are OPTIONAL for a minimal outbreak — their inputs are 0/off
by default (no fire/wind/mgmt keywords) ⇒ the minimal chain is bmsdit→bmcgrf→bmcbkp→bmcnum→bmatct→bmistd→bmmort→bmkill.


## CHUNK 3j DONE — BMKILL handback BIT-EXACT ⇒ full biology↔FVS bridge complete
bmkill! (WK2 handback) validated bit-exact (driver_bmkill.f) on both the kill case + the PROB-bound case. ⇒ bmsdit!
(FVS→BM) + 10 beetle kernels + bmkill! (BM→FVS) all bit-exact. REMAINING = ONLY the orchestration (no more per-tree
biology):
  • BMSETP + PPMAIN single-stand harness (reconstruct — the absent outer layer): per FVS cycle, bmsdit!(bridge) →
    bmdrv per-year loop → bmkill!(handback). For a MINIMAL outbreak (no fire/wind/mgmt keywords) the per-year loop
    is: bmcgrf! → bmcbkp! → bmcnum! → bmatct_single! → bmistd! → bmmort!(slow). GRF stressors default neutral,
    OAKILL=0 (no fast agents), so bmmort-fast/bmcwin/bmfire/bmobb/bmdfol/bmqmrt are no-ops.
  • BMPPIN keyword block: DISPERSE (→IBMYR1 outbreak start year + duration IBMYR2), HOST (→HSPEC, default LP/PP),
    PBSPEC, RANNSEED. Wire into keyword_dispatch.jl.
  • simulate.jl landscape seam (single-stand MXSTND=1): call the PPMAIN-equivalent per cycle if a DISPERSE outbreak
    is active + the year is in [IBMYR1,IBMYR2]; apply the returned WK2 to FVS mortality (pre-tripling, like the other
    insect seams). ⇒ end-to-end outbreak. Corner the composition (reconstructed harness, per the USER decision).
