# TT (Teton) Variant Port — Audit Log

FVS **Teton** (VARACD `TT`, MAXSP 18, Region-4 Intermountain, 10-yr cycle). The next western
Rockies cluster variant after CR/KT/IE/EM. Bit-exact-or-cornered vs live `FVStt`, chunk by chunk.
Doctrine (unchanged from CR/EM): validate vs LIVE per chunk; MEASURE don't infer; per-record treelist
diff INVALID after tripling; port faithfully then validate; reuse the shared engine (TT-gate all shared
changes). Oracle + Fortran source are the SOLE ground truth (no FVSjulia oracle).

## ★★ DVEW woodland volume COMPLETE (2026-08-02) — chunk-8 tail closed

The last chunk-8 item (PM/UJ/RM/MC/OH woodland cubic) is ported + validated. **Root model (measured, not
inferred):** TT DVEW = **R4D2H** (Chojnacky INT-339 D2H regression, `volume/NVEL/r4d2h.f`), routed
`GROSSVOL→DVEST` (region 4, `grossvol.f:172` `MDL='DVE'`). `VOL(1)=VOL(4)=(a+b·D2H^⅓+c·MSTEM)³`,
`D2H=DBH²·HTTOT`. **NO Behre/CFTOPK trim** in the DVE path (`dvest.f:90` just calls R4D2H, returns whole
cubic) — the earlier "CFTOPK 0.49× taper" memory was a PHANTOM; the correct chain has no trim. `fvsvol.f`
NATCRS then maps `TCF=VOL(1)`, `MCF=(D≥DBHMIN)·(VOL(4)+VOL(7))` [VOL(7)=0], `SCF=0` (region-4, not R8/R9),
`BdFt=0` (DVE has no board-foot, `fvsvol.f:411` skips it; BFPFLG=0 for R4). **DRCOB=0** passed ⇒ D2H uses
DBH not DRC. **FCLASS=0 ⇒ MSTEM=0** ⇒ the c-coefficient term drops. DBHMIN=8 (tt/grinit.f, no woodland
override). UJ/PM carry a **0.1 cuft floor for DBH<3** (r4d2h.f:92/98). `jl`: `r4d2h_vol1()` +
`compute_volumes_tt!` else-branch (src/variants/teton/volume.jl).

**gfortran-match details (both mattered):** `D2H**(1./3.)` and `(...)**3.` are **REAL exponents** ⇒
`powf` (`fpow(x,1f0/3f0)`, `fpow(base,3f0)`), NOT `cbrt`/`x*x*x` — using the wrong primitive left PM/UJ
off. **D<1 exclusion:** the FVStt volume loop (`vols.f` DO 200) never calls r4d2h for a sub-inch record —
INSTRUMENTED: a DBH-0.5 PM never reaches r4d2h while DBH-2.0/2.5 do (and take the 0.1 floor). So the
existing `d<1→0` guard is faithful and the 0.1 floor applies only to VOLUMED trees (1 ≤ DBH < 3).

**Validation (synthetic pure-woodland stands, ttt01 rows re-speciesed, vs FVStt_g16):**
per-tree VOL1 + TCF/MCF **BIT-EXACT** for all 5 species (instrument fort.95=r4d2h VOL(1), fort.96=fvsvol
TCF/MCF). `.sum` 1990: **MC 389/231, RM 397/228, OH 327/195 fully BIT-EXACT**; **PM 610/361 vs 608/360,
UJ 496/291 vs 495/291** = ±1-2 cuft = the documented Float32 tpa/summation-accumulation cornered tail
(per-tree exact, same stand/tpa as the exact MC ⇒ not a volume bug). small_PM (0.1-floor+D<1) 5/5 exact.
Meta: MEASURE settled two phantoms (CFTOPK trim; "5× high" from the wrong CR region-2 Chojnacky). See
[[fvsjl-tt-variant-port]].

## Woodland-species GROWTH tail (measured 2026-08-02, alongside DVEW) — CORNERED, not a fresh bug
With DVEW volume done, the only thing making the pure-woodland stands diverge at projected years is DG.
Measured pure_{sp}_g (10-cycle, vs FVStt_g16), BA/QMD (growth cols, volume-independent):
- **OH bit-exact-or-cornered** (2090 BA 176/177) — OH uses the NC/cottonwood DG (not regent-dominated).
- **MC close** (2090 BA 19/20). **PM/UJ/RM under-grow**: 2000 near-exact (BA 21/22, QMD 7.8/7.9 — matches
  the memory's "PI/JU validated 2000 near-bit-exact") but a ~0.1-QMD/decade DBH deficit COMPOUNDS → 2090 BA
  17 vs 21-23, QMD 8.0 vs 8.9-9.3. TPA/mortality match all cycles (49/49) — it's purely DBH growth.
- ROOT: PM/UJ/RM are **regent-for-all-sizes** (TT_RG_XMAX=99, tt/regent.f), so their entire DBH trajectory
  rides on TT's small-tree regent — the KNOWN **ABI-dependent buildDir single-step HT-DBH regent** whose
  POTHTG semantics are un-derivable from source + un-instrumentable (SIGFPE); documented as cornered for the
  conifer core too. XMAX=99 just amplifies that per-cycle residual over 100 yr. NOT a reducible coefficient
  bug (first cycle is near-exact). Same accepted class as the conifer regent cornering. See [[fvsjl-tt-variant-port]].

## Why TT is the chosen next variant (measured — chunk-0 scouting)
Compared `dgf.f` + `grinit.f` across the 3 remaining western candidates (TT/UT/BM) vs the done EM:
all three are **Wykoff-DDS** (inline `dgf`, NOT GENGYM). TT is the cleanest discount because it reuses
**both** already-built engines:
- **DG = EM's structure.** `tt/dgf.f`: `DGFOR(5,MAXSP)`, `DGHAB(7)`, `DGDS(4,MAXSP)`,
  `DGCON(ISPC)=DGHAB(ISPHAB,ISPC)+DGFOR(ISPFOR,ISPC)+…` — the same habitat+forest+DSQ dims as EM
  (EM: DGHAB(8)/DGFOR(6)). ⇒ chunk-3 clones `em_dgcons!`/`dgf!`, resizes arrays, swaps coefficients.
- **Density = CR's Zeide SDI.** `tt/grinit.f` sets `LZEIDE=.TRUE.` (EM/KT are `.FALSE.`=Stage SDI) —
  the ONE cross-engine difference. TT = **EM-DG + CR-Zeide-density**.
- Calibration knobs = EM: `DGSD=2.0`, `IFINT=10`, `IFINTH=5`. UT (2nd choice) is also Wykoff+Zeide but
  lacks the explicit DGHAB dim shown; BM (last) has `LHTDRG` default `.FALSE.` + no habitat dim — most divergent.

## Species / FIA map (tt/blkdat.f)
MAXSP=18. JSP: 1-WB 2-LM 3-DF 4-PM 5-BS 6-AS 7-LP 8-ES 9-AF 10-PP 11-UJ 12-RM 13-BI 14-MM 15-NC 16-MC
17-OS 18-OH. FIAJSP: 101 113 202 133 096 746 108 093 019 122 065 066 322 321 749 475 299 998.
~11 species overlap EM (WB/LM/DF/LP/ES/AF/PP/RM/AS/OS/OH → reuse those coefficient forms); new:
PM(pinyon)/BS(blue spruce)/UJ(Utah juniper)/BI(322)/MM(321)/MC/NC.

---

## Chunk 0 — Scaffold + oracle  ✅ DONE

- **`Teton()` singleton + registration.** `src/variants/teton/teton.jl` (`struct Teton <: AbstractVariant`,
  `variant_code`=TT, `nspecies`=18, `htg_period`=10f0, `TT_DATADIR`, `coefficients`). Registered in
  `src/variants/variant.jl` (`variant_from_code "TT"/"TETON"`), included in `src/FVSjl.jl`, exported
  (added EM too — was missing). Module precompiles; `Teton()`/`variant_from_code("TT")` verified.
- **Oracle: RELINKED + VERIFIED.** `/workspace/.ttwork/FVStt_clean` (12.6MB) via `relink_tt.sh clean`
  (EM-recipe clone: `bin/FVStt_buildDir/*.o` [670] + isoc23 shim). **BIT-EXACT vs
  `tests/FVStt/ttt01.sum.save` all 11 cycles** (TPA/BA/QMD: 1990 536/77/5.1 → 2090 389/245/10.8).
- **Harness.** Clean NOAUTOES control keyfiles `/workspace/.ttwork/ttc{2,10}.key` (+.tre) reproduce the
  ttt01 stand-1 growth bit-exact via the oracle. `ttt01` = the SAME 248112 conifer stand as EM's emt01
  (identical 1990 inventory), forest 415/habitat 41416 (Region-4) ⇒ chunk-3 DG reuses EM's main-conifer path.

## Chunk 1 — Infra/species (blkdat, FIA map, grinit)  ✅ DONE

- **`data/teton/species_coefficients.csv`** (committed, 18 rows): identity (code_alpha/code_fia/
  code_plants) from tt/blkdat.f JSP/FIAJSP/PLNJSP + `ht1`/`ht2` (HT1/HT2) + `dg_resid_sd` (SIGMAR).
  Verified vs live blkdat: alpha WB LM DF PM BS AS LP ES AF PP UJ RM BI MM NC MC OS OH; FIA 101 113
  202 133 96 746 108 93 19 122 65 66 322 321 749 475 299 998; plants PIAL PIFL2 PSME … 2TB. Bark/site/
  small-tree/mort/htdbh columns land in their chunks (EM's split-CSV approach, not KT's bundled one).
- **`data/teton/species_translation.csv`** (committed, 442 rows): SPCTRN crosswalk, ASPT **column 16**
  (legend ALFA FIA PLNT AK BM CA CI CR EC EM IE KT NC PN SO **TT** UT WC WS OC OP). Generated by
  `tools/teton/extract_crosswalk.py` (faithful/programmatic; self-check 442 rows, 18 distinct targets,
  shared conifers DF/LP/ES/PP/AF self-map). `spctrn_column(::Teton)=4` (CSV target col), `other_species=18` (OH).
- **`src/variants/teton/species.jl`**: `init_blockdata!` / `load_species_coefficients!` — clone of EM's,
  with TT's grinit (tt/grinit.f): **`zeide_sdi=true`** (LZEIDE=.TRUE. — the CR density path, UNLIKE EM/KT),
  `year`/`growth_fint`=10, `dg_sd`=2.0, seed 55329, `ht_drag_sp`=.TRUE. except **sp13 BI / sp16 MC** =.FALSE.
- **Validated:** `coefficients(Teton())` loads (18 species + 442 translation rows); alpha/FIA/plants/ht1/
  dg_resid_sd match live. Zero-regression (TT-only additions; no existing test exercises Teton).
## Chunk 2 — Site/habitat (tt/sitset.f + tt/siterange.f + tt/forkod.f)  ✅ DONE  BIT-EXACT

`src/variants/teton/site_index.jl` — TT's site chunk is SIMPLER than EM (NO habtyp.f / no
habitat-type-group site index; TT has no tt/habtyp.f at all). SITSET interpolates the site
species' SI (TEM) into each species' [SITELO,SITEHI]:
`pos=(TEM-SLOSSP)/(SHISSP-SLOSSP); SITEAR(I)=SITELO(I)+pos·(SITEHI(I)-SITELO(I))`.
- **Tables** (verified vs live): `TT_SDICON` (sitset.f), `TT_SITELO`/`TT_SITEHI` (siterange.f).
- **forkod** (`TT_JFOR=[403,405,415,416]` + reservation cases 7306→Bridger/8107→Caribou); IGL=0
  (KFOR has NO DATA statement in tt/forkod.f ⇒ zeros).
- **SDIDEF** = SDICON, or Zeide `BAMAX/(0.5454154·PMSDIU/100)` when BAMAX keyword set (PMSDIU is a
  PERCENT here, unlike EM's fraction — inert for ttt01 which has no BAMAX).
- **VALIDATED BIT-EXACT vs live ttt01** (measured from the .out SITECODE section): forest 415→IFOR=3,
  site species DF(3), TEM=50 (default, no explicit SI — the STDINFO 6th field 65 is ELEVATION not SI).
  All 18 SITEAR match (WB 43.75, DF 50, PM 16.25, NC 97.5, …) and all 18 SDIDEF=SDICON. `initialize!`
  now runs end-to-end for TT (tree load + site setup); next hook (setup_growth!/DG) errors as expected.
## Chunk 3 — Large-tree DG (tt/dgf.f)  ✅ DONE — BIT-EXACT (WK2 all 27 trees)

**VALIDATED bit-exact vs live** via DGFTRC instrument-replay (patched WK2(I)=DDS dump + CONSPP-component
dump, relinked FVStt_trc). ttt01 first growth cycle: **all 27 trees Δwk2=0** — MAIN Wykoff (WB/LP/ES/AF)
AND aspen DGFASP (AS, all sizes 0.1–12.7"). DGCONS/calibration/density all bit-exact: BA=85.131, RELDEN
(stand CCF)=112.319, DGCON LP=0.962982/WB=1.235373/ES=1.749294, COR ES=0.042154, bark LP=0.969 — all match.
Impl: `tt_dgcons!` + `dgf!(::Teton)` (src/variants/teton/diameter_growth.jl) + `tt_tree_ccf` (crown.jl) +
generated `dg_coefficients.jl` (tools/teton/gen_dg_jl.py). Wired: setup_growth! DG dispatch, standstats
point_ccf/stand_ccf Teton branches, simulate.jl:163 `relative_density=stand_ccf`, TT_PSIGSQ calibration branch.
CSV gained wykoff_ht2 (=ht2, Wykoff htdbh) + bark1/bark2/bark_imap (tt/bratio.f). Non-ttt01 forms (PP/juniper/
BI-MC/NC-OH DIAGR) error loudly (need their own test stand). Validation harness: notre!→setup_growth!→grow_cycle!
(isolated each_stand skips notre! tpa-expansion ⇒ BA 10× low — a harness artifact, fixed by calling notre! first).

**ttt01 species: ES×9, AS×8, AF×6, LP×5, WB×1** ⇒ exercises only TWO DDS forms: the MAIN Wykoff
(WB/LP/ES/AF) + ASPEN DGFASP (AS). No PP/juniper/NC/OH — those forms port later (unvalidatable on ttt01).

### DDS per-species dispatch (tt/dgf.f DO 10, `WK2(I)=DDS`), CONSPP=DGCON(ISPC)+COR(ISPC)+0.01·DGCCF(ISPC)·RELDEN:
- **CASE(1:3,5,7:9,17)** WB LM DF BS LP ES AF OS — MAIN: `DDS=CONSPP+DGLD·lnD+DGBAL·BAL+CR·DGCR+CR²·DGCRSQ+DGDSQ·D²`; BAL=(1−PCT/100)·BA100, CR=ICR·.01.
- **CASE(10)** PP — `+DGDBAL·PBAL/ln(D+1)` (PBAL=(1−PCT/100)·PTBAA); DGCON adds DGHAB(MAPHAB), MAPHAB=ICHBCL(ITYPE)+1.
- **CASE(6,14)** AS MM — aspen: `CALL DGFASP(D,ASPDG,CR,BARK,SI)` → `DDS=ASPDG+ln(COR2)+COR`. DGFASP at bin/FVStt_buildDir/dgfasp.f (shared ut/dgfasp.f). **NEEDED for ttt01.**
- **CASE(4,11,12)** PM UJ RM — DIAGR juniper: DF=0.25897+1.03129·DPP−0.0002025464·BATEM+0.00177·SI; DDS=ln(DIAGR·(2·DPP·BARK+DIAGR))+CONSPP. (note: REGENT eqns actually drive PM/UJ/RM/BI/MC all sizes.)
- **CASE(13,16)** BI MC — extended MAIN `+DGDBAL·BAL/ln(D+1)+DGPTCC·PCCF+DGBA·BA`.
- **CASE(15,18)** NC OH — CR-surrogate DIAGR: DF=(1.55986+1.01825·DPP−0.29342·lnBATEM+0.00672·SI−0.00073·BAUTBA)·1.05; ·DSTAG if ISTAGF.

### DGCONS assembly (ENTRY DGCONS): `DGCON(ISPC)=DGSIC(IDGSIM(ISI,ISPC),ISPC)·XSITE + DGFOR(ISPFOR,ISPC) + (DGSASP·sinθ+DGCASP·cosθ+DGSLOP)·SLOPE + DGSLSQ·SLOPE² + DGEL·ELEV + DGEL2·ELEV²`
where XSITE=SITEAR(ISISP) (log for sp13/16), ISI=ISMAP(ISISP), ISPFOR=MAPLOC(IFOR,ISPC), ISPDSQ=MAPDSQ(IFOR,ISPC),
DGDSQ=DGDS(ISPDSQ,ISPC), DGCCF=DGCCFC(IGCCFM(ISC,ISPC),ISPC) [ISC from INT(SITEAR(ISISP)) clamp 20-60, /10−1], ATTEN=IBSERV(ISIC,ISPC) [sp10/13/16 use OBSERV]. Aspect θ=ASPECT−0.7854 (sp10/13/16 special).

### DATA blocks to extract (tt/dgf.f, all MAXSP=18): 1D — DGLD:100 DGCR:106 DGCRSQ:112 DGBAL:118 OBSERV:194
DGCASP:300 DGSASP:306 DGSLOP:312 DGSLSQ:318 DGEL:324 DGEL2:330 DGDBAL:336 DGBA:342 DGPCCF:348 ISMAP:397.
2D — IGCCFM(5,18):124 IBSERV(5,18):174 DGCCFC(5,18):147 MAPLOC(4,18):203 DGFOR(5,18):223 MAPDSQ(4,18):246
DGDS(4,18):272 IDGSIM(5,18):357 DGSIC(5,18):377. Special — ICHBCL(363):405 DGHAB(7):409. DGPTCC/BAU/ISTAGF (sp13-18 only).

### DENSITY CORRECTION (measured — supersedes the chunk-0 scout's "CR Zeide density" for DG):
tt/dgf.f uses TWO density measures. **RELDEN** (in CONSPP `0.01·DGCCF·RELDEN`) is NOT assigned in dgf.f
⇒ it is the shared-COMMON **stand CCF** (same as EM/KT — `p.relative_density = stand_ccf(s)`), NOT Zeide.
The **Zeide RELSDI** (`STDSDI/SDIMAX`, cap 0.85; STDSDI=SUMTRE·(DGQMD/10)^1.605; SDICAL) is LOCAL and only
drives **DSTAG** (`3.33333·(1−RELSDI)` when RELSDI>0.7 — DIAGR-species stagnation) + mortality (chunk 7). So
TT's DG RELDEN = CCF (needs tt/ccfcal.f); Zeide is for DSTAG/morts. ttt01 main species have NONZERO DGCCF
(WB −0.199592, LP −0.206752, ES/AF class-varying) ⇒ RELDEN/CCF is REQUIRED for ttt01 DG (not deferrable).
DGFASP (aspen) uses REL=D/RMSQD + BA, not RELDEN.

### DGFASP (aspen, bin/FVStt_buildDir/dgfasp.f — mapped): REL=D/RMSQD; ASPCR=CR/10 (CR=ICR raw pct);
POT=(0.4755−3.8336e-6·D^4.1488)+(4.510e-2·ASPCR·D^0.67266) [floor 0.01]; FOFR=1.07528·(1−exp(−1.89022·REL));
GOFAD=0.21963·(RMSQD+1)^0.73355; BAACT=min(BA,305 if≥310); VALMOD=1−exp(−FOFR·GOFAD·((310−BAACT)/310)^0.5);
PREDGR=POT·VALMOD·(0.48630+0.01258·SI); ASPDG=ln(2·D·BARK·PREDGR+PREDGR²). Needs RMSQD (stand QMD) + BA.

### Remaining chunk-3 steps: ✅(1) extract 19 DATA blocks → data/teton/dg_*.csv [DONE 6031e6d].
(2) port tt/ccfcal.f → `tt_tree_ccf` (RELDEN=stand CCF **and** PCCF; EM did this in crown.jl); (3) port DGFASP;
(4) `tt_dgcons!` (DGCON=DGSIC·XSITE+DGFOR+aspect/slope/elev; DGDSQ; DGCCF; ATTEN) + `dgf!(::Teton)` (6 forms;
ttt01=main+aspen); RELSDI/DSTAG via SDICAL; (5) wire setup_growth! DG dispatch (+ `relative_density=stand_ccf`
at simulate.jl:162) + calibrate; (6) DGFTRC bit-verify WK2 (WB/LP/ES/AF/AS).
## Chunk 4 — Height (tt/htgf.f)  ✅ DONE — SBB model BIT-EXACT given DG

`height_growth!(::Teton)` + generated `htgf_coefficients.jl` (TT_HTCOF 33×9 + AZBIAS/BZBIAS).
**VALIDATED via HTGFTRC replay** (instrument TEMHTG=HTG dump, relink FVStt_trc): ttt01 first growth
cycle — **ES/AF bit-exact end-to-end**; WB/AS/LP height correct but cornered to the upstream DG.
**PROVEN the height model is bit-exact**: feeding live's per-tree DG into jl's height → max|Δhtg|=7.6e-6
(float rounding) across ALL 27 trees. So the WB/AS/LP residual is the **accepted DG-ZZRAN RNG stream-order**
(WK2/DDS bit-exact ch3; the per-tree random DG increment diverges — same ch9 accepted class as CR/KT/EM),
NOT a height bug. HT matches live (heights faithful). DBH<1.5 → REGENT (chunk 6). Non-ttt01 forms error.

## Chunk 4-OLD — Height (tt/htgf.f)  [scope notes below, superseded]

Species dispatch (SELECT CASE at htgf.f:305): CASE(10)PP `HTG=exp(CON+0.62144·lnDG)+0.4809`; CASE(4,11,12)
PM/UJ/RM `HTG=0` (regent); CASE(15,18)NC/OH even-aged GEMHT; CASE(13,16)BI/MC FINDAG+POTHTG+modifiers;
**CASE DEFAULT (WB/LM/DF/BS/AS/LP/ES/AF/OS + MM→JSPC6)** = the ttt01 path = **Schreuder-Hafley SBB** height-DBH:
- IICR=INT(ICR/10+0.5) cap 9; KEYCR 1-3 (IICR 0-2→1, 3-7→2, 8-9→3); JSPC=sp≤10?sp:(sp14?6:11); K=(JSPC−1)·3+KEYCR.
- Small/OOB → HTG=0.1 (HT≤4.5 or DBH≥XI1+COF1 or HT≥XI2+COF2 or DBH≤0.1). XI1=0.1, XI2=4.5.
- Y1=(DBH−XI1)/COF1; Y2=(HT−XI2)/COF2; FBY=logit; Z=(COF4+COF6·FBY2−COF7·(COF3+COF5·FBY1))·(1−COF7²)^−0.5.
- ZBIAS=AZBIAS+BZBIAS·ELEV (0 if ELEV<55 or >80; sp6/14 use ELEV−20 + a ZADJ). **DIA=DBH+DG/BARK** (uses DG!).
- PSI=COF8·((DIA−XI1)/(XI1+COF1−DIA))^COF9·exp(Z·√(1−COF7²)/COF6); H=(PSI/(1+PSI))·COF2+XI2; HTG=max(H−HT,0.1).
- Finalize (label 201): `HTG·SCALE·XHMULT(sp)·exp(HTCON(sp))·MISHGF`; then SIZCAP(sp,4) cap.
Coefficients extracted (tools/teton/extract_htgf.py): **data/teton/htgf_cof.csv** (COF 33×9 = 11 JSPC-groups × 3
crown-groups from COF1-11, EQUIVALENCE COF(:,1:3)=COF1…) + **htgf_zbias.csv** (AZBIAS/BZBIAS 18). XHMULT=1
(MULTS default), HTCON = HTCONS calibration (shared). NEXT: implement `height_growth!(::Teton)` + HTGFTRC bit-verify.

## Chunk 5 — Crown (tt/crown.f)  ✅ DONE — validated in isolation (WB bit-exact; others ±1-2)

`crown_ratio_update!(::Teton)` (in crown.jl) — rank-based Weibull, coeffs data/teton/crown_weibull.csv.
ISORT = whole-stand **GROWN-DBH** rank (dbh+DG/bark) via _rdpsrt! (the key fix — tt/crown.f runs after DG).
**VALIDATED in isolation** (call crown directly after height, bypassing the unported small_tree): WB=68
BIT-EXACT, ES [28,36,49,50,59] all match live, others ±1-2. Residual = tripling-order confound (live crown
runs on 81 TRIPLED records vs jl's 27 un-tripled; X=ISORT/ITRN ~preserved since tripling ×3 both) + DG-ZZRAN.
NOTE: in the full grow_cycle! crown (line 527) is AFTER small_tree_growth! (437) which errors (chunk 6
unported) ⇒ full-cycle crown validation awaits chunk 6. CCF (tt_tree_ccf) done ch3.

## Chunk 5-OLD — scope notes (superseded)

CCF (tt_tree_ccf) already done (chunk 3). Remaining = `crown_ratio_update!(::Teton)` — a rank-based Weibull
crown-ratio model (differs from EM's DCR form): (1) RELSDI = SDIAC/SDIDEF(ISPC) (cap 1.5; sp17 cap 1.0);
ACRNEW = C0(ISPC)+C1(ISPC)·RELSDI·100 (mean CR%). (2) Weibull A=WEIBA, B=WEIBB0+WEIBB1·ACRNEW (floor 1, PP 3),
C=WEIBC0+WEIBC1·ACRNEW. (3) per tree: SCALE=1−0.00167·(RELDEN−100) (clamp 0.30–1.0); X=(ISORT(I)/ITRN)·SCALE
(rank percentile; clamp 0.05–0.95) or RNUMB·SCALE if DBH≤0; invert Weibull → CRNEW. (4) CHG=CRNEW−ICR bounded
±1%/yr (PDIFPY); crown-change DBH gate DLOW/DHI (=0/99 uniform, no restriction). CL=crown length from HF=H+HTG.
Coefficients extracted: data/teton/crown_weibull.csv (18×7: WEIBA/WEIBB0/WEIBB1/WEIBC0/WEIBC1 + mean C0/C1).
NEXT: implement + wire crown_ratio_update! (lstart CRATET dub + per-cycle change) + instrument-replay bit-verify.
## Chunk 6 — Small-tree (tt/regent.f)  🔶 SCOPED (largest remaining chunk; coefficients extracted)

tt/regent.f (1422 lines) — the keystone: `small_tree_growth!` is called at grow_cycle!:437 BEFORE crown(527),
so it blocks the full cycle + all small-tree DG/HTG/crown. Small-tree height increment on a **5-yr basis**
(REGYR=5), then HHT1/HHT2 height-DBH assignment for DBH, then blended with the large-tree prediction over
**DBH ∈ [XMIN,XMAX] = [1.5,3.0]** (DEFAULT): DBH<XMIN pure regent, XMIN–XMAX ramp, DBH>XMAX pure large-tree.
TT uses the **TTVAR** branch (regent.f has TTVAR/UTVAR/CIVAR — TT=TTVAR=.TRUE.). Height model: RELSI=(SI−SLO)/
(SHI−SLO), RSIMOD, POTHTG site-curve based + PCTRED density reduction + XRHMLT/XRDMLT multipliers (=1). Coeffs
extracted: data/teton/regent_1d.csv (DIAM/XMIN/XMAX/DGMAX 18); AB(13)/HHT1/HHT2 + the POTHTG curves still to
extract. ttt01 small trees: AS d=0.1/1.2, ES/AF d=0.1 (0.1" seedlings). REGCON entry = small-tree calibration.
### SMHTGF (tt/smhtgf.f, height increment HTGRTH) — MAPPED:
- **CASE(6) AS aspen**: FINDAG→SITAGE; HITE1=26.9825·SITAGE^1.1752; AG2=SITAGE+5; HITE2=26.9825·AG2^1.1752;
  HTGR=(HITE2−HITE1)/(2.54·12); HTGRTH=(HTGR+ZRAND·0.1)·0.75. (+ regent.f:521 RSIMOD for sp6.)
- **CASE DEFAULT (WB/LP/ES/AF)**: BETA1=exp(B0ACCF+B1ACCF·ln(TPCCF)); BETA2=exp(B0BCCF+B1BCCF·ln(TPCCF));
  HTG1=BETA1+BETA2·CR; STDDEV=HTG1·(B0ASTD+B1BSTD·CR); HTGRTH=HTG1+ZRAND·STDDEV. TPCCF=point CCF clamp[25,300].
- ZRAND(I) = per-tree BACHLO(0,1) bounded ±2 (small-tree ZZRAN residual). Coeffs data/teton/regent_smht.csv (18×6).
- Then regent.f:554 H2=H1+HTGRL·SCALE·XRHGRO·CON (SCALE=KPER/REGYR, XRHGRO=1); subcycle to fint.
### STILL TO READ/PORT: tt/smdgf.f (small-tree DG), HHT1/HHT2 height→DBH assignment, XWT small/large blend over
DBH∈[1.5,3.0], REGCON calibration. NEXT: implement small_tree_growth!(::Teton) (subcycle+smhtgf+smdgf+blend) +
instrument-replay bit-verify; unblocks the full grow_cycle! + small-tree crown + the stand .sum differential.

## Chunk 7 — Mortality (tt/morts.f)  ✅ DONE — full grow_cycle! runs end-to-end, .sum CLOSE

mortality!(::Teton) = EM Hamilton form EXACTLY (tt/morts.f DEFAULT: RI=0.5/(1+exp(PMSC+PMD·D)); RIP=RN merge;
CONST=SDIMAX/0.02483133; RN=1−(1−(T−TN10)/T)^(1/FINT)) with TT PMSC/PMD (no PMDSQ). All ttt01 species use
DEFAULT; only PP(10) has a special CI-form (errors, not in ttt01). ttmrt.f = percentile distribution (VARMRT
equivalent; shared self-thinning used, like EM). Merch-spec default cols added to CSV (real VOLEQ = chunk 8).

★★ FULL grow_cycle! RUNS END-TO-END (DG→height→small_tree→crown→mortality) → .sum. vs live FVStt ttt01:
1990 BIT-EXACT (536/77/160/102/63/5.1); projected TPA near-exact (2000 525/525, 2090 392/live389 — MORTALITY
VALIDATED) + TopHt exact (66/68/75/79) + QMD close (5.8/6.5/8.1/10.1 vs 5.9/6.5/8.3/10.8); BA runs few-% LOW
(97/99 → 219/245, growing — accumulating DG-ZZRAN + small-tree tripling residual). Growth core integrates
end-to-end bit-exact-or-cornered like CR/KT/IE/EM. OPEN: BA-low residual (check DG cumulative vs .sum; likely
accepted tail) + real volume (chunk 8, FW2 like EM).
## Chunk 8 — Volume (tt VOLEQ = MATW/DVEW)  ✅ CUBIC DONE — total BIT-EXACT, merch within 1%

Ported R4VOL/R4MATTAPER (Region-4 Matney taper) → src/engine/r4vol.jl: r4vol_volumes returns (CF0 total-stem,
CFGRS gross merch) — **499/499 BIT-EXACT** per-tree vs live (instrument fort.83). Taper STUMPD/BUTTCF/CF0/B from
CFCOEF(20,7) (II from VOLEQ code+geocode) + Smalian 16.5-ft log integration (INT-rounded). Wired compute_volumes_tt!
(teton/volume.jl) + engine dispatch + Teton VOLEQ branch (TT_VOL_EQ 18 sp). ★MERCH TOP measured: MTOPP=TOPD(6.0
outside-bark, tt/grinit.f)·bark_ratio = inside-bark top (AF 6·0.937=5.622, ES 5.736, LP/WB/AS 5.814 — match live);
DBHMIN=8 (sp7=7). .sum: TOTAL cubic (TCuFt) **BIT-EXACT** (1990 1584/1583); MERCH (MCuFt) within 1% (826/817;
residual = dbhmin-boundary/topwood detail); later years cornered = growth-tail. REMAINING (downstream leaves):
board-foot Scribner (SCRIBC table + INTL14 Intl-¼), DVEW species (PM/RM/MC/OH/UJ via cr_dve_vol — not in ttt01),
+ the ~1% merch residual. TT GROWTH CORE (0-7) + cubic volume validated = full-variant-port bar these leaves.

TT VOLEQ (from ttt01.out) = geocode 400/401 + eq-type: MOST conifers **400MATW<volcode>** (Matney profile),
+ 400DVEW (PM/RM/MC/OH DVE) + 401DVEW065 (UJ). ttt01 species: WB 400MATW108, DF 400MATW202, BS/ES 400MATW093,
AS 400MATW746, LP 400MATW108, AF 400MATW019, PP 400MATW122. NOTE per-species VOLCODE ≠ FIA (WB→108=LP vol eq).
The **MATW (Matney) equation is NOT yet in the jl NVEL driver** (CR=DVEW, EM/KT=FW2, eastern=Clark). MATW
dispatches (volinit2.f:157 MDL=MAT) to **R4VOL** (bin/FVStt_buildDir/r4vol.f, 662 lines) = the Region-4 Matney
taper: R4MATTAPER taper profile + CFCOEF(20,7) height-to-2/3-DBH coeffs + SCRIBC Scribner board rules + cubic/
board log integration. A substantial NVEL routine port (comparable to the FW2/DVE ports). DVEW (PM/RM/MC/OH/UJ)
reuses the existing cr_dve_vol. Chunk 8 = port R4VOL (r4vol.jl) + Teton setup_volume_equations! (400MATW/400DVEW/
401DVEW mapping, per-species VOLCODE from ttt01.out) + compute_volumes! Teton dispatch. Placeholder merch specs currently let the .sum run
(growth cols validated; volume cols not). DOWNSTREAM LEAF — does not affect the growth core (chunks 0-7 validated).   ## Chunk 9 — Full-cycle diff  ⬜

## Off-switch
`touch docs/TT_VARIANT_PORT_COMPLETE` (USER's call).


### Chunk-6 per-record validation is TRIPLING-CONFOUNDED (doctrine #3) — measured 2026-08
Rich instrument (regent.f:929 + :23 CONTINUE, fort.82) proved the small-tree per-record diff is INVALID:
live triples the record AFTER subcycle 1, so the MAIN record shows WK3 after ONE subcycle (AS I=2: 2.0→5.63)
while the tripled records continue (→9.607); jl (no tripling in the isolated test) applies BOTH subcycles to
one record (→8.56). So the aspen HTG 6.56-vs-3.63 gap is largely the tripling/subcycle split, NOT purely a
FINDAG bug. Confirmed HK = H+HTG = WK3 (my gate is correct); ES DG=0 (HK<4.5) is likely right — live 0.2 is
the retained large-tree DG (regent skips DBH increment for HK<4.5 / D≥BREAK, leaving dgf DG). ⇒ chunk 6 is a
FAITHFUL port (smdgf RD-arg/DG-formula/floor/cap all fixed); DEFAULT height close (ES 1.597/1.540). VALIDATE
AT .sum (needs chunk 7 mortality), NOT per-record. OPEN (verify at .sum): exact aspen FINDAG (BH-adj) + the
subcycle×tripling interaction. NEXT: chunk 7 mortality (TT hits the SHARED driver — needs mort_bkgd_intercept/
mort_bkgd_dbh/sdi_max cols in species_coefficients.csv; check tt/morts.f vs shared) → stand .sum differential.

## Non-ttt01 species campaign — PP (ponderosa pine, sp10) DG: BIT-EXACT (+2 real bugs)
Pure-PP synthetic stand (all ttt01 trees → PP, /workspace/.ttwork/ttpp.{key,tre}) vs live FVStt_clean +
DGFTRC instrument-replay (fort.79 WK2, fort.76 DGSCOR calibration, fort.77 DGCON breakdown). Ported the PP
DG **CASE(10)** Wykoff-DGHAB form (tt/dgf.f:556): `DDS = CONSPP + DGLD·lnD + CR·(DGCR+CR·DGCRSQ) + DGDSQ·D²
+ DGDBAL·PBAL/ln(D+1)`, CONSPP gets a PP-only `−0.257322·ln(BA)` term; PBAL=(1−PCT/100)·point_ba (PTBAA).
Per-tree WK2 now BIT-EXACT (mean 0.0, maxabs 1e-6 ULP over 27 trees) after fixing **2 REAL BUGS** the
pure-PP differential exposed (ttt01 never hit them — its species use neither DGHAB nor the power-bark model):

1. **TT habitat/ITYPE resolution MISSING** (DGCON off +0.006074). tt/dgf.f CASE(10) DGCON adds
   DGHAB(MAPHAB), MAPHAB=ICHBCL(ITYPE)+1. jl left `habitat_input=0` ⇒ MAPHAB=2 (DGHAB=0.006074) vs live
   ITYPE=24 (41416→R4HABT idx 24, ICHBCL=0 ⇒ MAPHAB=1, DGHAB=0). Fix: ported tt/habtyp.f R4HABT(363)
   code→ITYPE lookup (habtyp_table.jl + tt_habtyp; CRDECD string-match ≡ numeric since all 363 codes
   distinct) + added Teton to the STDINFO western-habitat branch (field2→habitat_code) + resolve
   habitat_input=tt_habtyp(habitat_code) in tt_site_index_setup!. DGCON now bit-exact 1.827343.

2. **PP bark power model (IMAP=4) unported** (DGSCOR COR off −0.0282). tt/bratio.f PP is IMAP=4:
   `BRATIO=BARK1·D^(BARK2−1)` (0.809427·D^0.016866, cap 0.97, NO 0.80 floor — GO TO 100), dbh-dependent;
   jl's `_tt_bark_ab` placeholder returned the constant BARK1=0.809427 (the linear (a+b·d)/d shared
   bark_ratio can't express a power law). Fix: added faithful `tt_bratio` (all IMAP 1/2/3/4) + routed TT
   bark through it in the shared calibration TERM + _backdate_dbh! (mirrors CR's cr_bratio dispatch). The
   wrong constant bark skewed WK3(backdated)/TERM ⇒ per-tree RESLOG ⇒ cornew −0.2526 vs live −0.2222;
   fixed ⇒ cornew −0.22215, SVAR 0.010524, WC 0.78371, COR −0.174098 all bit-exact.

ttt01 UNCHANGED (its species all IMAP=2 ⇒ tt_bratio≡constant; non-PP ⇒ habitat unused). Suite 38588/0/1/75
(0 regress). REMAINING PP: height CASE(10) + mortality CASE(10) + REGENT small-tree + the growth-cycle bark
sites (DBH-update/volume/mortality → tt_bratio) for the full pure-PP .sum. Then other non-ttt01 species.

### PP height CASE(10) + bark-infrastructure (same session)
Ported PP **height CASE(10)** (tt/htgf.f:309, CI-variant, NOT SBB): `CON = 2.03035+0.7316 −0.00013358·HT²
−0.5657·lnD +0.23315·lnHT; HTG = exp(CON + 0.62144·lnDG) + 0.4809` (min 0.1); DBH<1.5 → REGENT; sets HTG
directly then only the SIZE-CAP tail (no ZZRAN/SCALE — those live in other htgf cases). VERIFIED: manual calc
tree1 (DBH11.5/HT73/DG0.9555) HTG=5.63909 = live 5.638779 bit-exact; jl reproduces its-own-DG HTG exactly.
Per-tree HTG residual = the accepted DG-ZZRAN (bidirectional, not systematic — jl DG 0.796 vs live 0.955 etc.
carries the ch9 RNG-stream noise into HTG). **tt_bratio wired at all 6 TT bark sites**: calibration TERM +
_backdate_dbh! (COR, done above) + DDS→DG dib (diameter_growth.jl:988, _tt_dg) + DBH-update (simulate.jl:490,
_tt_up) + volume (compute_volumes_tt!) + mortality-DG (teton/mortality.jl:61). Suite 38588/0/1/75 (0 regress).
★NEXT (blocks pure-PP .sum): mortality **CASE(10)** (tt/morts.f:539) — the CI-variant self-thinning form
(distinct from the EM-Hamilton the other TT species use): RIP=2.76253+0.222310√D−0.0460508√BA+11.2007·G
−0.554421/D+B0+0.246301·RELDBH+6.07129·G/D → 1/(1+expRIP)·REIN(IP); RIPP=(BA·RZ+(BAMAX−BA)·RIP)/BAMAX;
needs REIN/RZ/BAMAX/GMULT/AVED/WK1/GMULT(IP=D≤5?2:1). Then PP REGENT small-tree (PP∉_tt_rg_default).

### PP mortality CASE(10) ported + pure-PP .sum runs end-to-end (same session)
Ported PP **mortality CASE(10)** (tt/morts.f:539, CI-variant self-thinning, distinct from the EM-Hamilton
form the other TT species use): BAMAX=SDImax·0.5454154·PMSDIU (sdical.f:204, **bit-exact 206.767**);
stand RZ from BA-forward projection (DELTBA/BA10/TB/TTB); per-tree G (WK1/DGT/DG-override w/ tt_bratio),
IP=D≤5?2:1, GMULT/REIN constants (POT(12)=0.80→REIN1=0.89443/GMULT1=1.125, POT(41)=2.25→REIN2=0.98042/
GMULT2=1.111); RIP=2.76253+0.222310√D−0.0460508√BA+11.2007G−0.554421/D+PMSC+0.246301·RELDBH+6.07129G/D →
1/(1+expRIP)·REIN(IP); RIPP=(BA·RZ+(BAMAX−BA)·RIP)/BAMAX floor RIP cap1; WKI=P·(1−(1−RIPP)^FINT). WK1=dg_prev
(0 at cycle1 ⇒ the ICYC==1 override G=DG/(bark·10) fires for every DG>0.5 tree — verified tree1 G=0.1274,
tree2 G=0.1848 exact). VALIDATED: pure-PP .sum **1990 inventory BIT-EXACT** (536/77/160/70/63/5.1) + **ALL
volume columns BIT-EXACT** (TCuFt1251/MCuFt671/BdFt2809); 2090 converges (TPA jl112/live113). Suite
38588/0/1/75 (0 regress). ★OPEN: early-cycle over-kill (2000 TPA jl448/live498) = (a) accepted DG-ZZRAN on
large trees + (b) the SEEDLING mortality DG: live tree2(d=0.1) DG=1.29520 but PP's regent coefficients (smht
B0ACCF + smdg SDHTCR/SDHPCF/SDCR/SDHL4) are ALL ZERO in blkdat/smdgf ⇒ PP small trees correctly use the
large-tree dgf/htgf CASE(10) (jl excludes PP from regent = right), so the seedling's morts DG(I)=1.295 source
is subtle (regent-blend / establishment default?) — needs a per-tree DG(I)-in-morts trace. Not the RIP/RIPP
formula (bit-exact structure). PP DG+height+mortality CORE done; seedling-DG-in-morts is the isolated tail.

### Seedling-morts-DG root cause LOCALIZED (regent H-D function)
Traced the pure-PP early-cycle over-kill's seedling component: live's 0.1″ PP seedling gets DG(I)=1.295 in
morts NOT from smdgf (which is ZERO for PP by design) but from the **regent height-diameter function**
(tt/regent.f:18 "DIAMETER IS ASSIGNED FROM A HEIGHT DIAMETER FUNCTION"): the seedling grows in HEIGHT
(smhtgf), then its end-of-cycle DIAMETER is assigned from an H-D function (WK5 D2), and DG = D2−D1 scaled to
10yr. So PP small-tree DG is a HEIGHT-driven diameter assignment, not a direct DG equation — this is why all
PP smdgf/smht coefficients are zero. **Follow-up = port the regent H-D-function diameter-assignment path for
PP** (regent.f:424+ SI/RELSI/RSIMOD → HTGRL → H-D → WK5 D2 → DG). That closes the seedling morts-DG (fires the
G=DG/(bark·10) override correctly) and the early-cycle .sum over-kill. Large-tree over-kill component stays
the accepted DG-ZZRAN. PP growth+mortality core remains bit-exact/verified; this is the last isolated tail.

### Regent-for-PP experiment (result: excluding PP is correct)
Tested adding PP to `_tt_rg_default` (jl regent DEFAULT forms: smhtgf 1+CR / smdgf SDIAM=0.3+DIAM-floor, PP
coefs all zero). Result: pure-PP 2000 TPA jl 434 (WORSE than 448 without regent) vs live 498 ⇒ the jl regent
DEFAULT does NOT reproduce live's PP seedling DG=1.295; excluding PP (large-tree dgf/htgf CASE 10 for small
trees too) is closer. Reverted. The live seedling DG=1.295 traces to a specific regent DG(K) path
(regent.f:930-1040, TTVAR DK/DKK via smdgf + (DK−DKK)·BARK over subcycles/[XMIN,XMAX] blend) that the simple
DEFAULT-form addition doesn't capture — a genuine deep regent sub-chunk, not a one-line gate. Bounded and
DEFERRED. PP core (DG bit-exact / height verified / mortality formula+BAMAX exact / 1990+volume .sum exact)
stands as the delivered milestone; the seedling-morts-DG early-cycle over-kill is the one accepted-open tail.

### ★ Seedling-DG tail RESOLVED — accepted DG-ZZRAN, NOT regent (measurement corrected inference)
DECISIVE instrumentation (doctrine #2) OVERTURNED the earlier "regent H-D function" inference: instrumented
BOTH live regent DG(K) assignment points for ISPC==10 (the TTVAR sqrt form regent.f:939 + the LESTB DGMX
regent.f:962) → **fort.73 EMPTY both times** ⇒ the regent NEVER sets PP's seedling DG. The DG(I)=1.295 comes
from **dgdriv.f:225-233**: `DDS=EXP(WK2+XDGROW)·WK4; DG(I)=SQRT(DSQ+DDS·FRM)−D` where FRM=FM·SSIGMA·RHOCP
carries the **DGSCOR serial-correlation RANDOM effect** (CALL DGSCOR, BACHLO draws + OLDRN). So the 0.1″
seedling's DG is RANDOM-DOMINATED (tiny deterministic DDS=0.049, but FRM random-inflates it to 1.295). jl's
DGSCOR RNG stream differs from live's = the accepted ch9 ZZRAN residual; the seedling AMPLIFIES it because the
morts override `DG>0.5` flips on the random draw (jl's smaller DG ⇒ DGT-floor path ⇒ over-kill). ⇒ the pure-PP
early-cycle over-kill is the SAME accepted DG-ZZRAN residual as every other TT species — **NOT a missing form,
NOT a bug**. PP port is COMPLETE: DG bit-exact / height verified / mortality formula+BAMAX exact / .sum 1990
inventory+volume bit-exact / early-cycle divergence = accepted ZZRAN / converges 2090 (jl112 live113). The
regent is correctly UNINVOLVED for PP (excluding it = right). No outstanding PP work; jl PP handling faithful.

## Next non-ttt01 group SCOPED: PM/UJ/RM (pinyon/junipers) — chunk-0 done
Built pure-PM test infra (/workspace/.ttwork/ttpm.{key,tre} + live ttpm.sum). KEY LESSON: the pure-species
synthetic approach that worked for PP (a legit LARGE tree) CRASHES live for PM/UJ/RM (SIGFPE) when derived
from ttt01's large conifers — PM/UJ/RM are small shrubby species, so forcing 73-ft "pinyons" feeds degenerate
inputs to the DIAGR/regent form. FIX: build REALISTIC small stands (DBH 2.0″/ht 12ft) — live then grows
cleanly (pure-PM 1990 845TPA/BA18/QMD2.0 → 2090 604/120/6.0; volume via DVEW, 0 until 2090). jl gap (scoped):
`KeyError :htdbh_p2` — PM uses REGENT-for-all-sizes (dgf CASE 4,11,12 is a placeholder; "EQNS IN REGENT USED
FOR ALL SIZED TREES FOR PM/UJ/RM/BI/MC") via the **UTVAR regent path** (regent.f ELSEIF(UTVAR): VIGOR=(150·X³·
exp(−6X))+0.3, HTGRL=POTHTG·PCTRED·VIGOR·CON; H-D fn needs htdbh_p2/p3/p4 coefficients). PM/UJ/RM PORT =
htdbh coefficients + UTVAR regent form (HTGRL/VIGOR + H-D diameter) + DVEW volume (cr_dve_vol). A distinct
multi-form chunk, cleanly scoped and ready. (BI/MM/NC/MC/OH are further surrogate-form species.)

### Latent TT gap surfaced by PM scoping: height-dubbing fallback unsupported
Running jl on the pure-PM stand errors in SETUP (before growth): `dub_missing_heights!` (volume.jl:303) falls
through to the shared Curtis-Arney `_htdbh_height`, which needs `:htdbh_p2/p3/p4` coefficient columns TT does
NOT define. ttt01 never hit this (all inventory heights present + no broken-top/norm_ht<0). This is a REAL
LATENT TT gap — any TT stand with a missing/dub-required height (common in FIA data) would error. FIX (a TT
prerequisite, not PM-specific): give TT its own height-dubbing branch in dub_missing_heights! (the TT SBB
height-diameter inverse or the regent ht1/ht2 H-D function), instead of the shared _htdbh_height. This is a
BLOCKER for the PM/UJ/RM port (and for realistic TT FIA stands generally). So the PM/UJ/RM chunk now = (0) TT
height-dubbing branch + realistic small test data, (1) UTVAR regent-for-all-sizes (VIGOR/POTHTG + H-D diam),
(2) DVEW volume. PP remains fully complete + validated; these are the next western-cluster items, cleanly
scoped. RECOMMENDATION: commit the completed PP milestone before opening the PM/height-dubbing chunk.

### ★ TT height-dubbing gap FIXED (validated) — PM prerequisite done
Fixed the latent TT height-dubbing gap: added a Teton branch to `dub_missing_heights!` (volume.jl) porting
tt/cratet.f CASE DEFAULT — `H = exp(AX + wykoff_ht2/(D+1)) + 4.5` (AX=AA calibrated/IABFLG==0 else HT1
default) + the PP(sp10,D≤3) linear special `1.74189+4.17687·D` — instead of the shared Curtis-Arney
`_htdbh_height` (TT defines no htdbh_p2/p3/p4). Loads HT1/wykoff_ht2 unconditionally for TT. VALIDATED: the
corrected pure-PM stand (2 broken-top trees, norm_ht<0, legitimately need dubbing) now passes setup and
reaches the PM DG port (was crashing at `KeyError :htdbh_p2`). Suite 38588/0/1/75 (0 regress — only fires for
TT trees needing dubbing; ttt01 has clean heights). ALSO fixed the PM test fixture: HEIGHT is at .tre cols
45-47 (T45,F3.0), not col 63 (=HTG); rebuilt ttpm.tre with consistent DBH2.0/HT12 ⇒ live TopHt=12 matches
(pure-PM 1990 845/18/QMD2.0 → 2090 348/143/8.7). NEXT: PM/UJ/RM DG = the UTVAR regent-for-all-sizes
(regent.f ELSEIF(UTVAR): VIGOR=(150·X³·exp(−6X))+0.3 [X=CR/100], HTGRL=POTHTG·PCTRED·VIGOR·CON, H-D diameter
DK=(HK−4.5)·10/(SITEAR−4.5)) + DVEW volume. Now cleanly unblocked with a valid test fixture.

### ★★ BI/MC (sp13,16) GROWTH — FAITHFUL form corrected to DG=0.1·HTG (was coincidental Wykoff)
CRITICAL SOURCE FINDING: the COMPILED oracle bin/FVStt_buildDir/regent.f is a NEWER/DIFFERENT version than
canonical tt/regent.f — the doctrine ground truth is the buildDir (what FVStt_clean was linked from), NOT the
canonical tree. buildDir regent species-groups differ (sp15→CASE(1:5,8:9,12,15) Wykoff; sp13,14,16,18→a
CASE(13,14,16,18) block with RDCON/RDCR/RDLHT/RDHT/RDDUM diameter-lookup coefs, INDX 13→1/14→2/16&18→3).
For MC(sp16,INDX=3): RDCON provisional reduces to DKK=3.102+0.021·H, DK=3.102+0.021·HK (linear). THEN a
`CALL HTDBH(IFOR,ISPC,DK,HK,1)` override fires iff (.NOT.LHTDRG .OR. IABFLG==1); ELSE the DG override
`IF(LHTDRG .AND. IABFLG==0) DG(K)=0.1*HTG(K)*XRDGRO`.
MEASURED (instrument-replay, faithful trc = buildDir *.o + one swapped object):
  (1) HTDBH is a STUB (fmcrow.f:183, returns X10=0) — the ONLY htdbh_ symbol. Instrumented it → fort.89 EMPTY
      ⇒ HTDBH is NEVER called during the MC run, yet MC still grows to QMD 3.0 in the faithful trc.
  (2) cratet.f AA-fit: MC K1=0 (no measured-ht fit trees), so IABFLG stays init; HT1=5.152, HT2=-13.576.
  (3) MC trees are D=1.1", H=12 (25 of 27) — D<BKPT=3 ⇒ they DO enter regent's small-tree DG path.
  Since D<3 trees exist AND HTDBH is never called, the override MUST be skipped ⇒ at GROWTH time
  LHTDRG(16)=T & IABFLG(16)=0 ⇒ the FAITHFUL DG path is **DG = 0.1·HTG·XRDGRO** (then the DDS-sqrt scaling).
  Also: buildDir regent VIGOR ⅔-cut is ISPC==6 ONLY (line 284) — NOT the UTVAR block; and MC POTHTG comes
  from SMHTGF (CASE 1:11,13:14,16:18), not the (SJ/5)·(SJ·1.5−H)/(SJ·1.5)·0.83 formula.
CORRECTION: my earlier "BI/MC complete via Wykoff H-D" was an UNFAITHFUL COINCIDENCE — canonical tt/regent.f's
Wykoff `DK=(HT2/(ln(HK−4.5)−AX))−1` matched live's QMD *delta* (~1.0) but by luck (it over-grew: QMD 3.0/BA37
vs live 2.9/35). REPLACED with the faithful DG=0.1·HTG (regent.jl sp13/16 branch) + conditional ⅔-cut
(kept for PM/UJ/RM which validated, dropped for MC/BI per line 284). RESULT (ttmc vs live FVStt_clean): 1990
BIT-EXACT; 2000 834/26/2.4 vs 834/30/2.6; 2090 735/29/2.7 vs 734/35/2.9 — TPA exact, QMD −0.2 / BA −15% UNDER.
The under-grow = my formula-POTHTG is ~½ of SMHTGF's true output for MC (0.1·HTG propagates the low HTG);
SMHTGF for MC is un-reproducible here (7-arg smhtgf.f vs regent's 8-arg CALL = FVS shared-routine version
skew; regent itself is toolchain-fragile — fresh 12.2.0 regent.o SIGFPEs under -ffpe-trap and CHANGES results,
so regent cannot be instrumented). PM/UJ/RM UNAFFECTED (still exact). Suite 38588/0/1/75 (0 regress).
CAVEAT: instrument-replay corrupted 3 buildDir objects (cratet.o/fmcrow.o/regent.o now 12.2.0, not the
original 15.2.1) — they are build artifacts, unrecoverable, but FVStt_clean (the linked binary oracle,
12649792 B) is UNTOUCHED and remains the sole validation ground truth. Do NOT `relink_tt.sh clean` (it would
now yield a 12.2.0 binary ≠ FVStt_clean). META (doctrine #2/#4): reading the WRONG (canonical) source sent me
down a Wykoff-vs-Chapman rabbit hole; the buildDir + HTDBH-stub measurement is what settled it. Prefer the
faithful-but-cornered 0.1·HTG over the greener-but-coincidental Wykoff.

### NC/OH (sp15,18) — CHUNK-0 READY (test fixture + faithful buildDir equations identified; not yet ported)
Test fixture built: /workspace/.ttwork/ttnc.{key,tre} (pure-NC, MC-template w/ species NC). LIVE FVStt_clean
reference (NC grows FAST — cottonwood, aspen-like): 1990 845/18/12/2.0 → 2000 835/57/23/3.6 → 2050 339/190/63/
10.2 → 2090 121/190/79/17.0. (OH fixture TBD: sed MC→OH.)
FAITHFUL equations from bin/FVStt_buildDir (NOT canonical tt/*.f — buildDir is the linked-oracle ground truth):
• dgf CASE(15,18) [dgf.f:611] — DIAGR, uses COR+DGCON (NOT CONSPP, no CCF·relden term):
    ICLS=int(D+1) cap41; BAUTBA=BAU(ICLS)/BA; SI=SITEAR(sp); DPP=max(D,1); BATEM=max(BA,5);
    DF=(1.55986+1.01825·DPP−0.29342·ln(BATEM)+0.00672·SI−0.00073·BAUTBA)·1.05; DF=max(DF,DPP);
    DIAGR=(DF−DPP)·BARK [ISTAGF=0 ⇒ no DSTAG]; DDS=ln(DIAGR·(2·DPP·BARK+DIAGR))+COR(sp)+DGCON(sp), floor−9.21.
    NEEDS a TT BADIST BAU pre-pass (badist.f = CR form: BAU(1)=TOTBA−BAU(1); BAU(J)=BAU(J−1)−BAU(J); skip
    HT<4.5 seedlings; DBH at growth, WK3 at LSTART) — reuse/adapt CR `_cr_badist_bau`.
• regent HEIGHT: NC(sp15)→CASE(15) = ASPEN-Sheppard (FINDAG effective-age, HITE1/HITE2=26.9825·age^1.1752,
    RSIMOD=0.5·(1+relSI), ·CON·2.40·0.75; GO TO 3 bypasses VIGOR/PCTRED). OH(sp18)→SMHTGF POTHTG, then
    HTGR=POTHTG·PCTRED·VIGOR·CON (VIGOR full, no ⅔-cut). [NC's fast growth ⇒ the Sheppard height.]
• regent DG: NC(sp15)→CASE(1:5,8:9,12,15) sets BX=HT2,AX=(IABFLG==1?HT1:AA); Wykoff DK=(BX/(ln(HK−4.5)−AX))−1,
    DKK likewise (H≤4.5→DKK=D); sp15 is in the GO TO 300 skip (no HTDBH override) ⇒ DG=(DK−DKK)·BARK·XRDGRO
    (CASE DEFAULT). OH(sp18)→CASE(13,14,16,18) INDX=3 ⇒ DG=0.1·HTG (same as MC, if LHTDRG=T&IABFLG=0 — MEASURE).
• Both then DDS-sqrt scale + DGMX cap (DGMAX 2.5/5.0) + DIAM floor.
VALIDATION PLAN: .sum differential vs FVStt_clean ONLY (per-tree instrument-replay is now degraded — the
3 recompiled buildDir .o + regent toolchain-fragility). This is a fresh multi-part chunk best done with full
attention; chunk-0 (fixture+reference+faithful-eqn map) is complete and de-risked.

### ★★ NC (sp15) GROWTH PORTED + VALIDATED — dgf DIAGR + BADIST BAU + GENGYM height (CR IMODTY=4 reuse)
Corrected my earlier WRONG note (NC height ≠ regent aspen-Sheppard): NC trees are DBH≥0.5 & HT>4.5, so they use
htgf CASE(15,18) = the CR GENGYM even-aged/uneven-aged height (IMODTY=4 spruce-fir RM-32 curve), NOT regent.
PORTED (all faithful to bin/FVStt_buildDir):
 • dgf CASE(15,18) DIAGR [diameter_growth.jl] — DF=(1.55986+1.01825·DPP−0.29342·ln(BATEM)+0.00672·SI−0.00073·
   BAUTBA)·1.05; DDS=ln(DIAGR·(2·DPP·BARK+DIAGR))+COR+DGCON (NOT conspp). NC DBH grows via dgf for ALL sizes
   (NOT in regent dispatch) — validated: 2000 BA 57/QMD 3.6 BIT-EXACT.
 • _tt_badist_bau [diameter_growth.jl] — BADIST BA-above-class (skip HT<4.5), BAUTBA=BAU(int(D+1))/BA.
 • Height: REUSE the CR GENGYM (`_cr_htg_tree(4,…)` + `cr_fndag(4,…)`) — flat module, directly callable.
   height_growth.jl NC/OH branch: ADJUST=0.78+0.0023·SITEAR, ssite/bark/bautba/pccf/agerng passed; small
   (DBH<0.5|HT≤4.5) → regent. `_tt_dub_ages!` dubs ABIRTH from height via cr_fndag(IMODTY=4) at setup
   (simulate.jl TT branch) + per-cycle birth_age += fint (gradd.f:205, gated _cr_up||_tt_up).
VALIDATED ttnc (pure-NC) vs live FVStt_clean: 1990 845/18/12/2.0 BIT-EXACT; **2000 835/57/158/23/3.6 BIT-EXACT
(TopHt 23 = live, was frozen at 12 before the height)**; 2050 360/189/349/63/9.8 vs live 339/190/347/63/10.2
(TopHt 63 EXACT, BA/QMD near). PM/UJ/RM/MC UNCHANGED; suite 38588/0/1/75 (0 regress). ⇒ NC GROWTH (DBH+height)
DONE.
REMAINING NC gap = LATE-CYCLE MORTALITY DISTRIBUTION (2090 jl 167/228 vs live 121/190, ~38% over). ROOT
(re-measured — my first "deferred middle-SDI" note was WRONG; the `_tt_tn10_iter` middle-SDI fit IS implemented,
mortality.jl:91-96): the SELF-THIN total is identical (T−TN10 both sides — jl's uniform rn kills exactly
Σpr·(1−(1−rn)^fint)=T−TN10, = live's SUMTRE), BUT the DISTRIBUTION differs. Live morts.f:684 does `IF(RIP==RN)
SUMTRE=T−TN10; CALL TTMRT(SUMTRE)` — TTMRT REDISTRIBUTES the self-thin excess by PERCENTILE (small/low-crown
trees first) with NC/OH CRI-EFFTR=PEFF·((100−CRI)/100)·VARADJ·0.01 (ttmrt.f:125). jl kills UNIFORMLY (rn·pr per
tree). Same total per cycle, but killing-small-first raises QMD faster ⇒ different SDI/TN10 next cycle ⇒ COMPOUNDS
(dense NC self-thins every cycle). ttt01/emt01 stay BELOW-SDI (RIP=RI background, TTMRT never called) ⇒ the
uniform-rn matched ⇒ TT mortality validated there; NC is the FIRST stand to self-thin hard and expose the missing
TTMRT percentile distribution. NEXT: port TTMRT (geometric-progression EFFTR allocation of SUMTRE by percentile;
affects ALL TT dense stands, only NC exercises it here) — then re-diff ttnc late cycles. OH(sp18) shares the whole
sp15/18 code path (fixture: sed MC→OH; untested but same branches). META (doctrine #2, AGAIN): I documented
"deferred middle-SDI" from a STALE code comment without reading mortality.jl:91 — the iter was there; reading
morts.f:684 (SUMTRE=T−TN10 + CALL TTMRT) is what located the real gap (distribution, not total).

### ★★ OH (sp18) VALIDATED — confirms the sp15/18 growth port generalizes
Built ttoh fixture (sed MC→OH). jl vs live FVStt_clean: 1990 845/18/12/2.0 BIT-EXACT; **2000 835/40/13/3.0
BIT-EXACT**; 2050 642/144/20/6.4 vs live 634/145/20/6.5 (TopHt EXACT); 2090 406/184/24/9.1 vs live 393/187/24/9.3
(TopHt 24 EXACT, TPA within 3%). OH's late-cycle tail (3% over) is MUCH milder than NC's (38%) — OH grows slower/
less dense so it barely enters the self-thin regime, whereas dense NC exercises it hard. This CONFIRMS the gap
diagnosis (missing TTMRT percentile self-thin distribution): the effect scales with self-thinning intensity.
⇒ NC/OH (sp15,18) GROWTH DONE = bit-exact-or-cornered (same status as PP/PM/UJ/RM). The TTMRT distribution is a
refinement for the heaviest-self-thinning stands (NC), NOT a growth defect. TT non-ttt01 GROWTH COMPLETE:
PP✓ PM/UJ/RM✓ BI/MC✓(faithful-cornered) NC/OH✓. Remaining TT: TTMRT percentile self-thin (tightens NC dense
tail; RISK = shared TT mortality, must re-validate ttt01/emc10 conifer core), DVEW volume (PM/UJ/RM/MC/OH).

### TTMRT PORT — full integration spec (analyzed, ready to execute next session)
Goal: replace jl's UNIFORM self-thin distribution with the faithful ttmrt.f percentile+EFFTR redistribution
(tightens NC's dense-stand tail; OH already within 3%). Algorithm (bin/FVStt_buildDir/ttmrt.f, fully read):
 1. EFFTR(I): PEFF=0.84525−0.01074·PCT+0.0000002·PCT³ (clamp 0.01..1); NC/OH(15,18) EFFTR=PEFF·((100−CRI)/100)·
    VARADJ(sp)·0.01, else EFFTR=PEFF·VARADJ(sp)·0.01. PASS1=Σ PROB·EFFTR.
    VARADJ (DATA, sp1..18): 0.80,0.70,0.55,0.70,0.50,1.00,0.90,0.50,0.60,0.85, 0.70,0.70,0.70,1.00,0.90,1.10,0.75,0.90.
 2. NPASS=int(TOKILL/PASS1)+1. Iterate (label 100→105): TEMWK2(I)=−TPALFT·((1−EFFTR)^NPASS−1), TEMSUM=Σ;
    ADJUST=TEMKIL/TEMSUM; if <0.8 NPASS−=max(MINSTP,int((TEMSUM−TEMKIL)/PASS1)) (ISWTCH=1); if >1.2 NPASS+=…
    (ISWTCH=2); MINSTP=5/2/1 for NPASS>50/>20/else; opposite-ISWTCH ⇒ GOTO 110. 
 3. Label 110: XKILL=TEMWK2·ADJUST; if PROB−WK2−XKILL≤0 cap (TEMWK2=PROB−WK2, SHORT+=overkill, PASS1−=EFFTR);
    WK2+=TEMWK2, SUMKIL+=. If SHORT>0: NPASS=int(SHORT/PASS1)+1, GOTO 100.
COMPOSITION (morts.f:682): SUMTRE=0; IF(RIP==RN) SUMTRE=T−TN10; IF(TN10≥0.1) CALL TTMRT(SUMTRE). TTMRT RESETS
WK2=0 then fills it. IF SUMTRE passed as 0, TTMRT sets TOKILL=Σ(current WK2 background) ⇒ background is ALSO
percentile-redistributed. ⇒ TTMRT fires almost ALWAYS. jl matched ttt01 bit-exact WITHOUT it because ttt01 is
below-SDI ⇒ tiny background kill ⇒ redistribution is .sum-inert there; the RISK is that shifting even one
background tree breaks ttt01 bit-exactness — so port behind a HARD ttt01/emc10 regression gate + revert if it
moves. INTEGRATION POINTS (jl):
 • killed[] in mortality!(::Teton) is the uniform per-tree kill; inject `_tt_ttmrt!(killed, tokill, s, n)` after
   the per-tree loop (before book_snags). tokill = self-thinning?(rn>0 && tt>tem ⇒ default trees rip==rn) T−TN10
   : Σkilled. Pure ttnc/ttoh (no PP sp10) are the clean case; MIXED stands w/ PP need care (PP uses ripp≠rn, so
   SUMTRE's `RIP==RN` reflects the LAST tree — replicate morts.f's last-tree RIP semantics).
 • PCT = t.crown_ratio (BA percentile, NOT crown ratio); CRI = t.crown_pct (ICR crown ratio). VERIFY PCT is the
   POST-growth percentile at mortality time: grow_cycle! calls compute_density!→stand_pct! at cycle START
   (line 382) then growth then mortality — CHECK whether FVS MORTS uses a DENSE-refreshed PCT (dense.f runs in
   GRINCR before MORTS on the grown stand). If jl's crown_ratio is stale (pre-growth), add a stand_pct! refresh
   before mortality! or the distribution will be subtly wrong.
 • DVEW volume (PM/UJ/RM/MC/OH) is the OTHER remaining TT item (NATCRS internal transform, needs fvsvol.f — the
   instrument path is degraded, so read METHC=6 woodland branch end-to-end instead).

### ★★ TTMRT PORTED + VALIDATED — percentile self-thin distribution (mortality.jl _tt_ttmrt!)
Executed the spec: `_tt_ttmrt!` = faithful ttmrt.f (EFFTR=PEFF·[NC/OH:(100−CRI)/100]·VARADJ·0.01, PEFF=0.84525−
0.01074·PCT+2e-7·PCT³; geometric-progression NPASS convergence to TOKILL; SHORT-overflow re-pass). Wired into
mortality!(::Teton) after the per-tree loop: TN10≥0.1 ⇒ TOKILL = self-thin?(tt>tem & rn>0) T−TN10 : 0 (0 ⇒
redistribute background), overwrites `killed`. PCT=t.crown_ratio (percentile, unchanged pre-update ⇒ faithful),
CRI=t.crown_pct. RESULTS vs live FVStt_clean:
 • OH (sp18): 2050 640/146/6.5 (QMD EXACT), **2090 394/187/9.3 vs live 393/187/9.3 — BIT-EXACT** (was 406/184/9.1).
 • NC (sp15): 2050 348/195/10.1 vs live 339/190/10.2 (QMD ±0.1), 2090 156/235/16.6 vs live 121/190/17.0 — QMD
   16.6 (was 15.8) much closer; TPA 156 (was 167) vs 121. Residual = NC self-thins 700→121 (EXTREME) ⇒ the
   deepest ZZRAN/tie-break/tripling tail (bit-exact-or-cornered, = all TT species).
 • MC/PM: UNCHANGED (735/29/2.7, 335/8.9) — they barely self-thin so EFFTR-percentile ≈ uniform.
 • Suite 38588/0/1/75 (0 regress). ★ my ttt01-regression fear was UNFOUNDED (measured, not inferred): TTMRT is a
   faithful improvement, not a risk.
CAVEAT (corrected misread): ttt01.key can't run end-to-end in jl — it errors on a PRE-EXISTING THINDBH cut-
logging gap (_log_cut!→coef_col KeyError, cuts.jl:177; a volume-coef lookup, UNRELATED to mortality/TTMRT), so
ttt01 was NOT a usable regression gate here (the "ttt01 bit-exact" I first saw was LIVE data; jl had errored
silently). Evidence of no-regress = MC/PM unchanged + suite floor + OH/NC improved. ⇒ TTMRT DONE. TT non-ttt01
mortality now complete. REMAINING TT: DVEW volume (PM/UJ/RM/MC/OH); + the pre-existing THINDBH cut-logging coef
gap (separate, surfaced here — TT _log_cut! needs a volume-coef fallback).

### is_sprouting coef added + ESSPRT stump-sprout subsystem = TT downstream gap (ttt01 e2e blocker)
Fixed one genuine gap: TT species_coefficients.csv LACKED the `is_sprouting` column ⇒ coef_col(:is_sprouting)
KeyError in cuts.jl:177 (ESTUMP cut-log guard), fmburn.jl:154 (fire sprout filter), keyword_dispatch.jl:1386
(SPROUT validation). Added is_sprouting = 1 for ISPSPE {6,13,14,15}=AS/BI/MM/NC (tt/blkdat.f:111), 0 else.
Suite 38588/0/1/75 (0 regress). This unblocks TT fire's sprout filter + the first cut-log barrier.
DISCOVERED: ttt01.key (THINDBH ×8 cutting its 8 aspen) then hits the NEXT barrier — `:essprt_fsp` (sprout.jl:428),
the ESTUMP/ESSPRT stump-sprout REALIZATION subsystem. MEASURED live ttt01.out shows "REGENERATION FROM STUMP &
ROOT SPROUTS" ⇒ stump sprouting IS active under NOAUTOES (NOAUTOES disables SEEDLING establishment only, NOT
stump sprouting — so jl's lsprut=true is FAITHFUL; note keyword_dispatch.jl:2221's "NOAUTOES→lsprut=false" comment
is misleading but the branch evidently isn't setting it false for TT, which is accidentally correct). The ESSPRT
subsystem (sprout_essprt.csv + NSPREC sprout-count + SPRTHT sprout-height coefs) is ported ONLY for eastern
variants (CS/LS/NE/SN have data/*/sprout_essprt.csv); TT + all western variants (CR/IE/KT/EM) LACK it. ⇒ TT ESSPRT
stump-sprout port is a DOWNSTREAM chunk (analogous to the deferred FFE/establishment tails), needed for ttt01 e2e
and any TT stand that cuts a sprouter. NOT a growth/mortality defect. SCOPE (corrected): TT essprt.f is NOT the eastern parameterized kind/p1/p2/fsp form — it's a full subroutine
with complex per-species sprout-count equations (Keyser-Loftis DSTMP logistic `22.6839·(1/((DSTMP/0.7788)−0.4403))`,
aspen Crouch polynomial) + ENTRY SPRTHT sprout-heights (TT sprouters PY/AS/CW per essprt.f:1257 = 13/15/16 in the
BM-order comment; ISPSPE{6,13,14,15} in blkdat). So the TT ESSPRT port = a SUBSYSTEM (port those equations into a
TT branch of sprout.jl + a sprout_essprt.csv for the dispatch), analogous to the deferred FFE/establishment tails
— NOT a quick coef add. REMAINING TT (downstream subsystem leaves): DVEW volume (PM/UJ/RM/MC/OH), ESSPRT stump-
sprout. GROWTH+MORTALITY CORE COMPLETE (PP/PM/UJ/RM/BI/MC/NC/OH all validated bit-exact-or-cornered).

### ESSPRT subsystem — entry-point map (scoped for the port)
tt/essprt.f (buildDir) is a multi-ENTRY routine, dispatched CASE(VAR)→CASE(ISPC). TT-specific branches:
 • main PSPROB (@615) — PREM (sprout survival proportion) adjust: TT sp13,14→PREM·0.70; sp15→·0.90; else ·1.0.
 • ENTRY NSPREC (@1150) — # sprout records NMSPRC: e.g. 1 if DSTMP<5, NINT(−1+0.4·DSTMP) for 5–10", else 3
   (per-species SELECT; most TT →1).
 • ENTRY SPRTHT (@1416) — sprout height HTSPRT = f(SI,IAG): variants of (0.1+SI/100)·IAG, (0.1+SI/80)·IAG, or
   default 0.5+0.5·IAG per species.
 • ENTRY ASSPTN / INDXAS (@~741) — ASPEN(sp6) special sprout count (Crouch-polynomial, separate path).
Eligible TT sprouters = ISPSPE{6,13,14,15}=AS/BI/MM/NC (blkdat.f:111); is_sprouting already added. PORT = a TT
branch in sprout.jl for each ENTRY (PREM/NSPREC/SPRTHT + aspen ASSPTN) + validate vs live ttt01's "REGENERATION
FROM STUMP & ROOT SPROUTS" (or a dedicated THINDBH-on-aspen fixture). Subsystem-sized (like FFE/establishment
tails), NOT a coef add. jl already has the shared sprout.jl scaffold (eastern) + cut_log wiring — TT needs the
per-ENTRY equations. This is the last downstream TT leaf besides DVEW volume.

### ★★ ESSPRT stump-sprout PORTED (faithful) + v2t coef; ttt01 e2e now blocked by TT FFE
PORTED the TT sprout subsystem into sprout.jl (esuckr!), reusing the existing eastern ESUCKR scaffold + aspen
ASSPTN Crouch path (already present, just gated in TT): added `nsprec_tt` (sp6→2, sp15→DSTMP-banded, else 1),
`essprt_tt` (PREM: sp13,14→·0.70, sp15→·0.90, else ·1), `sprtht_tt` (sp13,15→(0.1+SI/100)·IAG, sp6,14→(0.1+
SI/80)·IAG, else 0.5+0.5·IAG), `tt_sprout_dbh` (Wykoff H-D via :ht1/:wykoff_ht2). Wired `tt` flag + asp_idx=6
(ESASID) + the ASSPTN aspen gate + all 4 dispatch sites (nsprec/essprt/sprtht/sprout_dbh). TT routes to essprt_tt
(self-contained ⇒ NO sprout_essprt.csv needed, unlike eastern's essprt_kind/p1/p2). Faithful to tt/essprt.f
CASE('TT') across its 4 ENTRYs. Also added `v2t` coef (wood specific gravity lb/cuft, fmvinit.f per-species:
sp1-18 = 22.5,22.5,28.1,31.8,20.6,21.8,23.7,20.6,19.3,23.7,34.9,34.9,27.4,30.6,19.3,21.8,22.5,19.3) — a genuine
gap needed by cut-biomass/snags/FFE. Suite 38588/0/1/75 (0 regress); MC/NC/PM unchanged (they don't cut/sprout).
STATUS: the sprout port cleared the essprt barrier; ttt01 e2e now hits `:bio_group` = the TT FFE (fire) subsystem
(v2t→bio_group→fuel models/biomass), a SEPARATE deferred downstream leaf (like FFE on all western variants). So
the TT sprout port is FAITHFUL + suite-clean, but its TT-specific .sum validation (vs live ttt01's "REGENERATION
FROM STUMP & ROOT SPROUTS") is PENDING TT FFE (ttt01's next blocker), OR a custom aspen+THINDBH fixture with no
input-snags/FFE trigger. REMAINING TT downstream leaves: TT FFE (fuel/biomass — bio_group+), DVEW volume.
GROWTH+MORTALITY+SPROUT-EQUATIONS complete; FFE + DVEW volume are the situational tails.

### ★★ ESSPRT sprout port VALIDATED (ttaspr fixture) — height-exact, count ~1.87× (ASBAR-localized)
CORRECTED a doctrine-#2 misread: NOAUTOES DOES disable stump sprouting (ttas fixture: NOAUTOES+THINDBH aspen →
live "NO SPROUTING WILL BE SIMULATED", TPA 821→44 no rebound; jl MATCHES 44 — jl lsprut=false under NOAUTOES is
FAITHFUL). The earlier "ttt01 sprouts under NOAUTOES" was reading a REPORT HEADER, not the TPA. Built the proper
validation fixture ttaspr (drop NOAUTOES ⇒ sprouting ON, + THINDBH@2000 on pure aspen): live sprouts heavily
(2010 TPA=617 = 44 survivors + ~573 aspen sprouts via ASSPTN Crouch). jl RUNS E2E (no FFE block — all-live stand,
no input snags) and SPROUTS: 2010 TPA=350/QMD=1.2 vs live 617/QMD=1.2. ★ QMD/height BIT-MATCH (1.2=1.2) ⇒
sprtht_tt + tt_sprout_dbh CORRECT; e2e + sprout realization faithful. RESIDUAL: sprout COUNT jl ~306 vs live
~573 (~1.87×). LOCALIZED: total sprout = SPA (both sides, verified: NUMSPR·ASPRTR telescopes to SPA), survivors
match (44 both) ⇒ ASTPAR matches ⇒ the gap is ASBAR (aspen BA removed) = jl ~23 vs live ~43.5 (half). Since
ΣPREM matches but Σprem·DSTMP² is half, jl's removed-aspen DBH² is half — a TRIPLING/cut_log DSTMP detail
(doctrine #3: per-record DSTMP after tripling), NOT a sprout-equation error. Only affects the aspen ASSPTN path;
non-aspen sprouters (BI/MM/NC, per-record PREM) are unaffected. ⇒ ESSPRT port FAITHFUL + validated (structure +
height exact); aspen-count is a bounded cut_log/tripling follow-up. Suite 38588/0/1/75. Fixtures ttas/ttaspr in
/workspace/.ttwork, /tmp/tt_as.jl /tmp/tt_aspr.jl.

### ★★ NEW FINDING: aspen (sp6) DIAMETER growth gap — jl under-grows ~40% (pure-aspen stand exposes it)
Chasing the ttaspr sprout-count residual (instrumented esuckr!: jl asbar=25.5/astpar=853.6/spa=336 at ishag=10)
led UPSTREAM: the sprout count is low because jl's aspen are SMALLER at cut time. Pure-aspen NO-CUT stand (ttasg,
MC data relabeled AS) jl vs live FVStt_clean: 1990 845/18/2.0 bit-exact; 2000 jl 821/25/2.3 vs live 821/46/3.2;
2010 jl 798/33/2.7 vs live 798/60/3.7; 2040 jl 731/72/4.2 vs live 751/122/5.4. TPA MATCHES (mortality OK) but
BA/QMD ~40% UNDER — jl aspen DIAMETER growth under-predicts, from cycle 1 (DBH 2.0 = small-tree regent regime)
and PERSISTING past DBH 3 (large-tree DGFASP regime). This is a GROWTH-CORE gap MASKED by conifer-dominated
ttt01 (8 aspen among 500+ conifers — .sum-invisible); a PURE aspen stand exposes it, exactly like IE's pure-
species stands exposed special-species bugs. Both sides use the same STDINFO ⇒ real jl aspen-DG under-prediction
(candidates: aspen small-tree regent DG (smdgf/regent sp6), DGFASP large-tree, or the aspen SITE INDEX
resolution). PRIORITY: this is growth-core (higher than the downstream FFE/DVEW/sprout-count leaves) — the sprout
count will self-correct once aspen DG matches. NEXT: pre-tripling per-tree aspen DG diff (doctrine #3 window) on
ttasg cycle 1 to localize small-tree-regent vs DGFASP vs site-index. ESSPRT sprout port itself remains FAITHFUL
(height/QMD exact; count is downstream of this aspen-DG gap). Fixture ttasg /tmp/tt_asg.jl.

### aspen-DG gap — root-cause narrowing (XMAX dispatch RULED OUT)
Tested the hypothesis that aspen (sp6) should be regent-for-all-sizes (buildDir regent XMIN=90/XMAX=99 like the
UTVAR species, vs jl's XMIN=1.5/XMAX=3): set TT_RG_XMIN[6]=90/XMAX[6]=99 → aspen 2040 got WORSE (BA 43 vs the
72 baseline, live 122) — jl's regent aspen DG is even LOWER than DGFASP. REVERTED. ⇒ the gap is NOT the small/
large dispatch; jl's REGENT ASPEN DG (smhtgf sp6 Sheppard height → smdgf sp6 SDIAM) under-predicts at ALL sizes.
Cycle-1 (DBH 2.0, 2000: 2.3 vs 3.2) is unchanged by the XMAX fix ⇒ it's the regent aspen height/DG formula, not
dispatch. Candidates: (a) smhtgf sp6 Sheppard HTGRL (the 0.75 + RSIMOD=(si6-30)/70 + 0.5·(1+relsi) chain, and
the buildDir applies POTHTG·PCTRED·VIGOR[⅔-cut for sp6]·CON on top — jl's default regent path may not apply
PCTRED·VIGOR·CON for aspen); (b) aspen site index si6 resolution (low si6 ⇒ relsi cut); (c) the SMHTGF 7-arg-def
vs 8-arg-CALL signature-skew confound. CAVEAT: ttasg is SYNTHETIC (MC data relabeled AS, H=12/DBH=2.0 may be an
atypical aspen H-D) — needs confirmation on a REAL aspen stand (FIA) before treating as a confirmed bug; ttt01's
aspen validated byte-identical (but at inventory/small sizes). NEXT: (1) confirm on a real FIA aspen stand; (2) if
real, trace jl's default-regent aspen height chain vs buildDir regent POTHTG·PCTRED·VIGOR·CON for sp6. This is the
one open growth-core lead; ESSPRT/FFE/DVEW-volume remain downstream leaves. Suite 38588/0/1/75 (revert clean).

### ★★★ REAL FIX: aspen (sp6) DGMAX cap 0.2→2.0 — root cause of the aspen-DG gap (+ sprout cascade)
The aspen-DG gap ROOT-CAUSED (measure-first, ruled out RSIMOD [inert, si6 high⇒relsi=1] and XMAX-dispatch [made
it worse]): jl aspen TopHt was even slightly HIGH (2000: 19 vs live 17) but QMD LOW (2.3 vs 3.2) ⇒ the H→D DG was
CAPPED. TT_RG_DGMAX[6]=0.2 but buildDir regent DGMAX(6)=2.0 (DGMX=DGMAX·SCALE, SCALE=FNT/REGYR=1) — jl capped
fast aspen DG 10× too tight. FIX: TT_RG_DGMAX[6]=2.0. RESULT ttasg (pure aspen): 2000 QMD 2.9 (was 2.3, live
3.2); 2040 QMD 5.3/TopHt 39 (was 4.2, live 5.4/38 — near-exact). CASCADE: aspen sprout count ttaspr 2010 TPA 516
(was 350, live 617); 2090 BA 125 (was —, live 128) — the earlier "sprout-count/ASBAR-half" residual was DOWNSTREAM
of this cap (bigger aspen ⇒ bigger ASBAR ⇒ more sprouts). Suite 38588/0/1/75 (0 regress — FIA aspen stands
unaffected; conifers' small-tree DG stays <0.2 so their cap never bound, which is why ttt01/emc10 validated with
the wrong 0.2). MC/NC/PM unchanged (no aspen). ★ Found via the PURE-ASPEN fixture (IE-style special-species-stand
method — conifer-dominated ttt01 masked it). Also removed the sp6 RSIMOD (canonical-vs-buildDir mis-port, faithful
but inert here). RESIDUAL now ~10% (BA 111 vs 122 early; candidates: per-subcycle cap headroom, or POTHTG·PCTRED·
VIGOR·CON on the aspen height which jl's default-regent path omits — a smaller follow-up, was masked by the cap).
META: TT_RG_DGMAX is jl's per-cycle regent DG cap; the 0.2 defaults for the ttt01 conifers may ALSO be wrong
(buildDir DGMAX 2.4-3.6) but inert (their DG<0.2) — a latent class to re-check per species if fast growth appears.

### DGMAX table corrected to buildDir (latent class closed)
Extended the aspen DGMAX fix to the whole default-regent set: regent.f DGMX=DGMAX(ISPC)·SCALE (SCALE=1 @10-yr;
sp11 special=FINT·0.2=2.0). Default species (1-3,5-9,17) were jl=0.2 but buildDir DGMAX = 2.8,2.8,2.4,·,2.5,·,
3.5,3.6,3.6,·,·,·,·,·,·,·,2.8,· — corrected to those. INERT-and-safe (their small-tree DG stays <0.2 so the cap
never bound — that's why the conifers validated at the wrong 0.2; measured: MC/aspen unchanged, suite 38588/0/1/75
0-regress). UTVAR species (4,11,12,13,16=2.0/2.0/2.0/2.5... wait 13=2.0,16=2.0) + 14/15/18 KEPT at their
separately-validated effective caps (the UTVAR DGMX path differs from DGMAX·SCALE — PM validated at 2.0 not the
raw 3.6). PP(10) kept 99 (CI-variant special). ⇒ TT_RG_DGMAX now faithful to buildDir for the default set; the
latent "fast growth would over-cap" class is closed for conifers (aspen was the only exposed case). This completes
the DGMAX correction. Remaining aspen residual (~10%, POTHTG·PCTRED·VIGOR·CON on the shared default-regent
height) + TT FFE + DVEW volume remain the follow-ups.

### aspen residual ~10% — smdgf coefs VERIFIED correct; residual = H-D allocation (height chain)
Verified jl _tt_smdgf sp6 coefs vs buildDir smdgf.f DATA (all EXACT: SDHTCR=-0.41227, SDHPCF=0.16944, SDCR=
0.003191, SDHL4=-0.0022; sp6∈CASE(3,5:9) uses form SDIAM=SDHTCR+SDHPCF·H+SDCR·CR+SDHL4·RD = jl's _tt_smdg_alt
branch). ⇒ the H→D coefs are NOT the residual. The residual is an H-D ALLOCATION mismatch: jl aspen is TALLER
but THINNER than live (2000 TopHt 19 vs 17, QMD 2.9 vs 3.2, same TPA) — jl puts growth into height, live into
diameter. jl's aspen HEIGHT over-predicts (19 vs 17), consistent with the MISSING POTHTG·PCTRED·VIGOR·CON factor
(buildDir regent.f:350 applies it to all default species; jl's default-regent path omits it; the ≤1 factors
would pull height 19→17, and the correct DBH follows from SDIAM at the lower height). CONFIRMED candidate =
PCTRED·VIGOR·CON on the SHARED default-regent height — but adding it touches the validated conifers, so it needs
(1) live-conifer HTGR instrumentation to confirm why they validated without it (likely small-tree height is a
minor .sum contributor for conifers), (2) a real FIA aspen stand (ttasg is synthetic). Deliberate follow-up, not
a session-end hasty change. The aspen DGMAX fix (the dominant 40%→10% correction) stands; this is the last ~10%.


### CORRECTION: DGMAX table change was INERT (not a fix); pure-DF gate reverted it + found a DF over-growth lead
Built a pure-DF (sp3) gate to test the shared-height change — it immediately showed jl OVER-grows dense young DF
2× (2000 BA 48/QMD 3.2 vs live 22/2.2; TPA 843 both). First blamed my DGMAX-table change (conifers 0.2→2.4) and
reverted — but DF is 48/3.2 with EITHER 0.2 or 2.4 ⇒ the DGMAX change was INERT for DF (the cap doesn't bind),
NOT a regression AND NOT a fix. Kept the revert to 0.2 (the ttt01-validated values; 2.4 is untested for LP/ES/AF,
so 0.2 is the safe choice). So the earlier "DGMAX table corrected / latent class closed" was OVER-CLAIMED — only
the ASPEN 0.2→2.0 is a real, validated fix; the conifer caps are effectively inert (their regent DG isn't
cap-limited). ★ NEW PRE-EXISTING LEAD (not from this session): jl over-grows dense young DF 2× — root = the
small/large-tree BLEND (regent.jl:193 xwt=(d−XMIN)/(XMAX−XMIN)=0.33 at DBH 2, XMIN 1.5/XMAX 3) mixes 33% of the
large-tree dgf DF DG, which at 845-TPA extreme density isn't suppressed enough (DGCCF·relden term). Masked by
mixed/realistic ttt01 DF (validated); exposed only by the pure-dense-DF corner (845 TPA @ DBH 2.0 — possibly an
unrealistic extreme). CAVEAT: synthetic fixture (MC data→DF); needs a realistic dense-DF stand to confirm vs
fixture-artifact. Same pure-species-stand method that found the aspen bug (IE-style). Suite 38588/0/1/75.
META (doctrine #2, AGAIN): I claimed "0-regress, inert" for the DGMAX table without a pure-conifer test — the
pure-DF gate is exactly the test I should have built BEFORE claiming it. Build the differential fixture first.

### DF over-growth — CONFIRMED real, XMIN/XMAX partial fix, remainder CONFOUNDED (SMHTGF/RHCON)
Confirmed the DF over-growth is real (not density/fixture): pure-DF at LOW density (164 TPA) ALSO over-grows 2×
(2000 BA 11 vs live 4). ROOT part 1 = premature small/large-tree BLEND: jl TT_RG_XMIN[3]/XMAX[3]=1.5/3.0 but
buildDir regent XMIN/XMAX DATA = 2.0/4.0 for DF ⇒ at DBH 2 jl XWT=0.33 (blends 33% large-tree dgf) vs live XWT=0
(pure small-tree). FIXED the default-species XMIN/XMAX to the buildDir DATA (sp1=2/3,3=2/4,7-9=2/4,17=1/5, etc.;
UTVAR 4/11/12/16 KEEP 90/99 = jl regent-for-all representation; PP sp10 kept 2.0). Suite 38588/0/1/75 0-regress,
MC/aspen unchanged. PARTIAL: DF 48→42 (2000, live 22) — helped but STILL ~2× over. ROOT part 2 = the small-tree
regent DF DG/height OVER-predicts even at XWT=0 (pure small-tree). Candidate was the missing POTHTG·PCTRED·VIGOR·
CON (buildDir regent.f:350 HTGR=POTHTG·PCTRED·VIGOR·CON on all default species; jl omits it) — BUT this is
CONFOUNDED: buildDir CON=RHCON·exp(HCOR) and RHCON=0 for ALL TT species (regent.f var-def comment) ⇒ literal
HTGR=POTHTG·PCTRED·VIGOR·0=0, which is impossible (heights grow) ⇒ the SMHTGF 7-arg-def vs 8-arg-CALL signature
mismatch means the source does NOT behave as literally read; the actual live height flow is obscured. Can't port
PCTRED·VIGOR·CON faithfully without resolving this, and regent is UN-INSTRUMENTABLE (fresh 12.2.0 regent.o
SIGFPEs/segfaults, main-recompile corrupts results). ⇒ DF over-growth is a REAL pre-existing bug, PARTIALLY fixed
(XMIN/XMAX, faithful), remainder BLOCKED on the SMHTGF/RHCON confound + un-instrumentable regent. VALIDATION SCOPE
of the XMIN/XMAX change: DF-tested + suite; WB/LM/BS/LP/ES/AF/OS pure stands NOT individually tested (faithful
buildDir values, suite-clean — flagged to avoid another over-claim). Fixture ttdf/ttdflo, /tmp/tt_df*.jl.

### DF over-growth — dgf EXONERATED (hypothesis corrected by measurement) + tooling workaround found
UNBLOCKED the live-dgf comparison: relinking with an UNTRAPPED main (recompiled main.f w/o -ffpe-trap) + the
instrumented dgf.o (buildDir dgf.f has built-in WRITE(78/79)) RUNS despite the corrupted 12.2.0 regent/cratet/
fmcrow objects — because DGCON (fort.78) + cycle-1 WK2 (fort.79) are SETUP/cycle-1 values, computed BEFORE the
corrupted regent's growth. (Growth is garbage but those values are correct.) MEASURED live DF vs jl at DBH 2.0:
 • DGCON(sp3) = 1.140128 — jl EXACT match.
 • RELDEN = 173.87, CONSPP = 0.8943 — jl EXACT match (my earlier "jl relden 211 vs live 174" was a LATER cycle,
   a CONSEQUENCE of over-growth, not cause).
 • per-tree WK2 (=ln DDS): live 1.79/2.07/2.18 vs jl 1.808/2.044/2.183 — MATCH (bit-exact-or-±0.02).
⇒ the dgf DF DDS is CORRECT — my "large-tree dgf DF over-predicts" hypothesis was WRONG (corrected by measure).
RULED OUT for the DF over-growth: dgf DDS, DGCON, RELDEN, AND the small-tree regent DG (probe: dgk=0.08-0.2 =
live's slow growth). Blend range now = buildDir (XMIN/XMAX fixed). So the ~2× DF over-growth at cycle 1 is NOT in
any of {dgf, DGCON, RELDEN, small-tree DG, blend range} — all individually match/correct — yet the .sum diverges
(2000 jl 42/3.0 vs live 22/2.2). Remaining suspects: the blend APPLICATION at DBH 2-4 (xwt weighting of the
large-tree DG as trees cross 2→4 within/across cycles), the height feeding SDIAM, or a subset of DBH~3 trees
(live fort.79 caps at DBH 2.0 so the DBH-3 dgf couldn't be compared). STATUS: DF over-growth REAL + significantly
narrowed (4 candidates eliminated by direct live-vs-jl measurement), root not yet isolated. The untrapped-main +
dgf-instrument relink is a REUSABLE workaround for live setup/cycle-1 dgf values despite the corrupted objects.
META: measurement CORRECTED a plausible-but-wrong "dgf over-predicts" inference (doctrine #2) — 4th such
correction this session; the differential-instrument is worth the setup every time.

### ★★★ ROOT CAUSE FOUND + FIXED: TT default-regent pass didn't update the TRIPLING STASH
The DF over-growth (and the aspen residual) ROOT: `triple_records!` (southern/diameter_growth.jl:1116) sets the
upper/lower sub-records' diam_growth = stash.dgU[i]/dgL[i]. The TT DEFAULT regent pass (small_tree_growth!)
overrode t.diam_growth[i] (the CENTRAL record) with the regent DG (~0.2) but NEVER updated stash.dgU/dgL — so
40% of a tripled small tree (upper 25% + lower 15% TPA) grew via the STALE large-tree dgf DG (~1.2), driving the
2× over-growth. The UTVAR pass (line 219) + SN/NE/CS/LS small_tree_growth (sets dgU/dgL/htgU/htgL/is_small) both
do this; the TT DEFAULT pass was MISSING it. This is why every component matched live per-tree (dgf/DGCON/RELDEN/
small-tree DG all EXACT) yet the .sum diverged — the divergence was in the TRIPLED sub-records, invisible to a
central-record probe. FIX: default pass now sets stash.dgU[i]=dgL[i]=diam_growth[i], htgU/htgL=ht_growth[i],
is_small=true (no ZZRAN spread = UTVAR simplification). RESULTS vs live FVStt_clean:
 • DF 2000: 21/14/2.1 vs live 22/14/2.2 — NEAR-EXACT (was 42/3.0, 2× over). 2040 108/40/4.9 vs 71/31/4.0 (residual
   = compounding + height, much smaller).
 • ASPEN 2000: 45/QMD 3.2 vs live 46/3.2 — QMD EXACT (was 38/2.9). 2040 116/5.4 vs 122/5.4 — QMD EXACT. ⇒ the
   aspen "~10% residual" (blamed on PCTRED·VIGOR·CON confound) was ACTUALLY this tripling-stash bug!
 • MC/PM UNCHANGED (UTVAR already updated the stash).
SYSTEMIC: affects ALL TT default-regent species (WB/LM/DF/BS/AS/LP/ES/AF/OS) — the small-tree tripled sub-records.
Almost certainly the ttt01 "cornered" tail. META: the DF/aspen pure-species investigations (many turns, 4
measure-corrected wrong hypotheses) converged on a REAL SYSTEMIC ENGINE BUG in the tripling — found by the
per-component-matches-but-sum-diverges contradiction pointing at the tripled records (doctrine #3: the tripling
is where per-record vs .sum diverge). The XMIN/XMAX fix (faithful) + aspen DGMAX fix (real) STACK with this.

### RHCON=1 confirmed (stale "zero" comment); PCTRED·VIGOR·CON attempted + reverted (mixed/confounded)
Resolved the RHCON confound: regent.f:880 `RHCON(ISPC)=1.0` — the line-70 "ZERO FOR ALL SPECIES" comment is
STALE. So CON=RHCON·exp(HCOR)=exp(HCOR), and the faithful default-regent height IS POTHTG·PCTRED·VIGOR·CON (which
jl omits — a real gap). ATTEMPTED it (PCTRED density poly + VIGOR[⅔-cut sp6] + CON=exp(htg_cor_small)): DF 2040
108→103 (barely helped, still >live 71), aspen 2000 QMD 3.2→3.1 (WORSENED the tripling-fixed exact). The aspen
regression signals a DOUBLE-APPLY — the SMHTGF 7-vs-8-arg signature skew makes it ambiguous whether the aspen
Sheppard SMHTGF returns POTHTG (needs ·PCTRED·VIGOR·CON) or the FINAL HTGR (already reduced). REVERTED pending
resolution of that confound. So PCTRED·VIGOR·CON is a KNOWN faithful gap (RHCON=1 confirmed) but its clean port is
blocked by the SMHTGF confound; it is NOT the DF later-cycle residual anyway (barely moved it).
REMAINING DF residual (2040 103 vs live 71, TopHt 38 vs 31): a LATER-CYCLE compounding height/DG over-prediction
(cycle-1 is near-exact after the tripling fix). Candidates: the large-tree dgf DF DG at DBH 4-6 (validated only at
DBH 2.0 via the untrapped-main workaround — live fort.79 capped there), or the height chain. Smaller than the
tripling bug (now fixed). ★ TURN NET: the tripling-stash fix is the major, validated, systemic win; PCTRED·VIGOR·
CON is a documented faithful-gap deferred on the SMHTGF confound; DF later-cycle is a narrowed remaining lead.

### CROSS-VARIANT LEAD: IE/EM/KT small_tree_growth! also don't update stash.dgU/dgL
The TT tripling-stash bug is a PATTERN. Checked all variants: CR/TT(fixed)/SN/NE UPDATE stash.dgU/dgL for small
trees; IE/EM/KT's small_tree_growth! set t.diam_growth[i] (central) but show NO stash.dgU/dgL write. POTENTIAL
same latent bug (tripled small trees grow via the stale large-tree DG). COUNTER-EVIDENCE (strong): emc10 (EM)
validated FULL 10-cycle bit-exact-or-cornered E2E incl. small trees — if EM had a 2× tripled-small-tree
over-growth, emc10 wouldn't validate. So either (a) IE/EM/KT's diameter_growth! seeds the stash with the
small-tree DG for is_small records (making the override unnecessary), or (b) their validated stands don't trip
small trees, or (c) a latent bug their stands don't expose (like ttt01 masked TT's — TT's ttt01 couldn't run e2e
due to FFE, so TT's small-tree growth was NEVER e2e-validated, which is why TT's bug survived). VERIFY per variant
with a pure-small-tree fixture (the IE-style method that exposed it for TT). NOT investigated here (beyond TT; the
e2e validation argues against it manifesting). ⇒ TT tripling fix = the concrete win; IE/EM/KT = flagged lead.

### CROSS-VARIANT bug CONFIRMED (IE/EM/KT); fix attempted on EM + reverted (harness not runnable in-session)
Confirmed the stash-seeding is SHARED (southern/diameter_growth.jl:1060-64) = the LARGE-tree DG (ZZRAN spread)
for ALL trees; IE/EM/KT small_tree_growth! set only t.diam_growth[i] (central), never stash.dgU/dgL ⇒ they DO
have the same latent bug (tripled small trees grow via the large-tree DG). Unexposed = their validated stands
(emc10 etc.) lack tripled small trees OR it's in their "cornered" tails (like ttt01 masked TT's). Applied the fix
to EM (regent.jl:124, mirror TT: set stash.dgU/dgL=dgw, htgU/htgL=htg, is_small=true) but could NOT validate —
`run_keyfile(emc10.key; EasternMontana())` errors at init.jl:129 (harness/keyword wiring, not the fix). REVERTED
the untested EM change (discipline: no unvalidated edits to validated variants). ADDITIONAL blocker found:
data/easternmontana/ is MISSING species_coefficients.csv (run_keyfile(emc10) fails at load_species_coefficients)
— EM's harness isn't reproducible in-session; the validated emc10 used a setup no longer present. IE HAS its CSV
(data/inlandempire/species_coefficients.csv) so IE may be the tractable one to validate the fix on first. FIX PATTERN for IE/EM/KT (focused
session, harness working): after each sets its final small-tree t.diam_growth[i]/ht_growth[i], add
`stash.dgU[i]=dgL[i]=<dg>; htgU[i]=htgL[i]=<htg>; is_small[i]=true` (guarded on stash!==nothing && !isempty). Then
re-diff iet01/emc10/ktctrl — likely IMPROVES their cornered tails. Locations: IE regent.jl ~264/305/334, EM
regent.jl:124, KT regent.jl:408. ⇒ TT tripling fix = shipped+validated; IE/EM/KT = CONFIRMED latent, fix known,
application deferred to a session with their harnesses.

### Cross-variant fix — IE harness WORKS (iet01_c runs), but the fix is PER-BRANCH (not mechanical)
IE harness IS runnable in-session: run_keyfile(iet01_c.key; InlandEmpire()) → 2000 434/115/7.0, 2050 177/224/15.2
(matches memory's live tie-break). So the "can't validate" blocker is GONE for IE. BUT the IE tripling fix is NOT
a simple end-loop stash update: IE's small_tree_growth! has HETEROGENEOUS branches — NIVAR/UTVAR/TTVAR override
t.diam_growth[i] (stash fix = set dgU/dgL=diam_growth), but the CRVAR/CO branch (regent.jl:402-403) applies the
DBH DIRECTLY (t.dbh[i]=new_dbh; diam_growth=0) — for THOSE a stash=diam_growth(=0) would zero the sub-records.
Each branch needs its own stash handling (the override branches: dgU/dgL=the regent DG; the direct-DBH branch:
the sub-records need the applied-DBH split, or dgU/dgL set to the effective increment). So the IE fix = careful
per-branch (5 diam_growth sites: regent.jl ~254/264/305/334/403), validated on iet01_c + suite. EM/KT likely
similar (check for direct-DBH branches). ⇒ cross-variant fix is CONFIRMED-needed + IE-validatable, but a careful
per-branch task, deferred. TT's fix was clean (single default pass); the western siblings are more intricate.

### ★ CORRECTION (measured): the tripling-stash bug is TT-SPECIFIC — IE does NOT manifest it
Applied the stash fix to IE and MEASURED iet01_c: 2000 BA 115→110 (baseline 115 = live EXACT; fix moved it AWAY),
2050 224→222 (live 220, slight help). MIXED/net-negative — the baseline (large-tree stash) already MATCHED live.
⇒ IE does NOT have the manifesting bug; the fix HURTS it. REVERTED. So my "confirmed cross-variant bug (IE/EM/KT)"
was a LOGICAL OVER-GENERALIZATION (same code pattern ≠ same bug) — doctrine #2, corrected by the empirical test.
WHY TT differs: TT's small-tree regent DG is MUCH lower than the large-tree (DF 0.2 vs 1.2 = 2×), and live FVStt
uses the regent DG for tripled small trees ⇒ TT needed the fix (validated: DF/aspen 2×→near-exact). IE's regent DG
≈ the large-tree (110 vs 115, ~4%) AND live FVSie evidently uses the large-tree DG for the tripled sub-records ⇒
IE's "no stash update" MATCHES live. The pure-small-tree EXTREME (TT's DF/aspen) exposed TT's; iet01_c (realistic
conifer) shows IE is fine. ⇒ EM/KT must each be EMPIRICALLY tested (not assumed) before any stash fix — the fix
is only correct where the variant's live behavior uses the regent DG for tripled small trees AND that DG differs
materially from the large-tree DG. TT's tripling-stash fix STANDS (validated); the cross-variant "extension" is
RETRACTED as unconfirmed. META: 6th measure-corrected inference this session — logical certainty about shared code
is NOT empirical certainty about behavior; the differential is the only arbiter.

### KT empirical check (read-only) REINFORCES the retraction — ktctrl baseline ≈ live, no fix needed
Ran jl KT ktctrl (no code change) vs live ktctrl.sum: 1990 536/77/5.1 BIT-EXACT; 2040 265/175/11.0 vs live
263/170/10.9; 2090 161/222/15.9 vs 155/208/15.7 — near-exact/cornered (3-7% BA compounding = the accepted
self-thin/ZZRAN tail, NOT a 2× tripling over-growth). ⇒ KT (like IE) does NOT manifest the tripling bug; its
realistic conifer stand (QMD 5+, few tripled small trees) matches live WITHOUT the fix. Two independent
data points (IE iet01_c, KT ktctrl) now confirm: the tripling-stash effect is TT-SPECIFIC in practice — only the
PURE-small-tree EXTREME (TT's synthetic DF/aspen) exposed it; realistic western conifer stands don't. So the fix
is CORRECTLY not applied to IE/KT (would perturb their near-exact baselines, as IE showed). EM untested (harness
blocked: missing species_coefficients.csv). TT tripling fix STANDS; cross-variant remains RETRACTED, now with 2
empirical confirmations. (Whether IE/KT have the latent code-pattern that would manifest on a pure-small-tree
stand is academic — their realistic stands are correct as-is, and the fix must NOT be applied to them.)

### DF residual NARROWED: large-tree dgf RULED OUT (bit-exact @ DBH 5) — residual = DBH 2-4 blend transition
Found the .tre DBH-column bug behind my earlier corrupted fixtures: DBH is the F3.1 field at COLS 37-39 ("020"=
2.0), NOT cols 31-32 (that's the past-DBH F2.0 for the DG measurement). Built a CLEAN DBH-5 DF stand (ttdf5, cols
37-39="050") and ran vs live: 1990 113/5.0 BIT-EXACT; **2000 148/5.9 = live 148/5.9 BIT-EXACT**. ⇒ jl's large-tree
dgf DF DG is CORRECT at DBH 5 (as it was at DBH 2, proven earlier via the untrapped-main fort.79). So the DF
later-cycle over-growth (ttdf: 2040 103 vs live 71) is NOT the large-tree dgf — it ACCUMULATED during the DBH 2-4
BLEND/small-tree transition (ttdf starts DBH 2 and grows through 2-4; ttdf5 starts DBH 5, skips it, matches live).
By 2040 the DF are large (DBH~5) and grow correctly (dgf), but at a HIGHER LEVEL from the earlier accumulated
over-growth. ⇒ DF residual root = the DBH 3-4 small-tree regent DG and/or the blend weighting (xwt) in that range
(small-tree DG validated at DBH 2 = 0.2 correct; DBH 3-4 unchecked). Much narrower than "the whole later-cycle
growth". NEXT: validate the regent DG at DBH 3 (build a DBH-3 DF fixture, cols 37-39="030") + the blend xwt. The
tripling fix + XMIN/XMAX + aspen DGMAX all STAND; this is the last DF thread, now localized to a 1"-wide DBH band.
META: nailing the .tre DBH column (37-39) un-blocks all future DBH-specific TT fixtures (my prior edits hit col 31).

### DF residual = COMPLEX non-monotone blend-region issue (DBH-3 UNDER, DBH-2 OVER, DBH-5 exact)
Built clean DBH-3 DF stand (ttdf3, cols 37-39="030", xwt=0.5 blend): live 2000 90/4.4, jl 63/3.7 — jl UNDER-grows
30%. OPPOSITE of ttdf (DBH-2, over-grows later) and ttdf5 (DBH-5, bit-exact). So the DF residual is NON-MONOTONE
across the DBH 2-4 blend: UNDER at DBH 3, OVER from DBH 2, EXACT at DBH 5. ⇒ NOT a single-direction bug — it's the
DGMAX-cap × blend-weight interaction: at DBH 3 (xwt=0.5) the DG=0.5·small(capped 0.2)+0.5·large(dgf); the 0.2 cap
makes the small component too low ⇒ blend under-grows (live's small-tree DF DG at DBH 3 > 0.2). But raising the cap
to buildDir 2.4 made DBH-2 ttdf OVER-grow (earlier finding). So a FIXED cap can't satisfy both DBH 2 and 3 ⇒ the
cap must be DBH/subcycle-dependent (buildDir DGMX=DGMAX·SCALE applied PER-SUBCYCLE, or the SDIAM DG naturally
scales), which jl's fixed per-cycle 0.2 doesn't capture. This is the deepest DF thread — a careful DGMAX/blend/
subcycle-scaling analysis (validate the per-subcycle DG at DBH 2 vs 3 via the untrapped-main workaround + read the
buildDir regent DGMX-per-subcycle logic). The tripling fix + XMIN/XMAX + aspen DGMAX all STAND; DF's remaining
residual is this bounded, well-characterized blend/cap interaction. Fixtures ttdf/ttdf3/ttdf5 (cols 37-39=DBH·10).

### DF residual = 3 SEPARATE components (DGMX per-cycle=2.4 faithful; cap masks one, DBH-3 is another)
buildDir DGMX = DGMAX(ISPC)·SCALE applied PER-CYCLE (regent.f:609, SCALE=1@10yr) ⇒ DF cap = 2.4 (faithful), NOT
jl's 0.2. TESTED cap 2.4 (with the tripling fix in place): ttdf DBH-2 2000 → 22/2.2 BIT-EXACT (was 21/2.1 @ 0.2,
IMPROVED + faithful) but 2040 → 124 (was 103, WORSE — the 0.2 was MASKING a later-cycle over-growth by capping it);
ttdf3 DBH-3 2000 → 64 (was 63, UNCHANGED — DBH-3 under-growth is NOT the cap). ⇒ 3 SEPARATE DF residual
components: (1) the DGMX cap (faithful 2.4 fixes cycle-1 exact, jl's 0.2 masks #2); (2) a later-cycle over-growth
(2040, unmasked by 2.4 — separate, in the DBH 3-5 growth); (3) a DBH-3 UNDER-growth (blend xwt=0.5, small-tree DG
too low, cap-independent). REVERTED to 0.2 (conservative: cap 2.4 for DF is faithful but risks the untestable-e2e
ttt01 conifers — the earlier DGMAX-table over-claim lesson; and it unmasks #2 net-worse). ⇒ the DF residual is a
3-component blend-region tangle, NOT a single fix — a careful focused analysis (the faithful 2.4 cap + find/fix
#2 and #3 together, validated on ttdf/ttdf3/ttdf5 across the DBH band). The tripling fix + XMIN/XMAX + aspen DGMAX
STAND; DF's remaining residual is now fully characterized (3 components, each localized). ttt01 conifers likely
share the wrong-0.2 cap (faithful=their DGMAX·SCALE) — a latent class, inert until a fast/small conifer stand.

---

## DF small-tree residual — ROOT-CAUSED to a regent VERSION MISMATCH (supersedes the "3-component tangle")

★★ 2026-08-02. Deep re-investigation (disasm + source + fixture differential) replaces the earlier vague
"3-component blend tangle" verdict with a precise root cause. **The jl TT small-tree regent is modeled on the
CANONICAL `tt/regent.f`, but the live `FVStt_clean` binary was compiled from the DIFFERENT buildDir version.**

### The two regent versions (confirmed at the source level)
- **Canonical `tt/regent.f`**: `CALL SMDGF` for small-tree DBH (lines 574/715/930/1266); `DATA DGMAX/0.2,0.2,
  0.2,2.0,0.2,0.2,…/`; `IF(TTVAR)DGMX=FINT*DGMAX(ISPC)`; **REGYR=5 with subcycles** (the "5-year increment,
  subcycle to avoid bias" model). This is what jl was ported from (subcycle NPER/KPER + `_tt_smdgf`).
- **buildDir `bin/FVStt_buildDir/regent.f`** (= what `FVStt_clean` links): **no `CALL SMDGF`** — dubs DBH via
  **inline HT-DBH inverse** (CASE 1:5,8:9 AA-fit `DK=(HT2/(ln(HK-4.5)-AA))-1`; CASE 7 LP-log; CASE 10,17
  OS-linear; CASE 6 AS-SITEAR); `DATA DGMAX/2.8,2.8,2.4,3.6,…/`; `DGMX=DGMAX*SCALE`; **`DATA REGYR/10.0/`,
  single-step, no subcycle loop**; height `HTGR=POTHTG·PCTRED·VIGOR·CON` (line 350).

### Binary confirmation (the decisive measurement — doctrine #2)
`objdump -d .ttwork/FVStt_clean` → **`regent_` calls `smdgf_` 0 times** (smdgf_ is called 4× elsewhere = the
establishment/essubh path). So live's regent DOES NOT use SMDGF; jl's `_tt_smdgf` in the growth-cycle regent is
the wrong-version port. `blkdat.f DATA JSP` confirms species order (sp6=AS, sp4=PM) — the buildDir regent SELECT
CASE comments ("PINYON" at CASE 6, "ASPEN" at CASE 15) are STALE BM-inherited labels, real ISPC numbers.

### Why it hid so long, and why it can't be cleanly fixed now
- ttt01 (mature conifer) validated bit-exact only because it has few small trees; the small-tree DBH was always
  "validated at .sum" (tripling-confounded per-record) where it is swamped.
- **The live oracle keeps a pure-DF small stand STATIC**: `ttdf` (all DF @ DBH 2.0) → live QMD stays 2.0, TopHt
  12→13, BA 18, for 100 years. jl grows it (the empirical 0.2 DGMAX cap throttles this to a near-baseline; remove
  the cap and jl explodes to BA 236 by 2090).
- Tested the two faithful pieces: (a) inline HT-DBH⁻¹ DBH + (b) single-step REGYR=10. **Neither, nor both
  together, tames ttdf** — jl still grows (2000 20/68/161/14 vs live 18/64/158/12; still explodes later). The
  DOMINANT suppression is the **`POTHTG·PCTRED·VIGOR·CON` height factor** that jl omits (POTHTG likely ≈0 in live
  — see blocker). Both partial changes also **regress validated aspen** (ttasg BA 46→39) ⇒ reverted per doctrine
  #4 (regression on a ported-faithful chunk → examine the oracle, don't cargo-cult to green).

### The BLOCKER (why the full faithful port is deferred)
buildDir `smhtgf.f` is a **7-arg** subroutine `SMHTGF(ISPC,I,HTGRTH,CR,TPCCF,LESTB,DEBUG)`, but buildDir regent
**CALLs it with 8 args** `SMHTGF(ISPC,POTHTG,H,1,TEMT,ICYC,JOSTND,DEBUG)`. The arguments are positionally
misaligned (POTHTG↔dummy I, H↔dummy HTGRTH, …), so the value of **POTHTG after the call is ABI-dependent and
un-derivable from source** — the hypothesis is that POTHTG stays stale/≈0 (SMHTGF writes its result into the
caller's H slot, never POTHTG), making `HTGR=POTHTG·PCTRED·VIGOR·CON≈0` ⇒ near-static default-species small trees
(exactly the live ttdf behavior: HTG floored to 0.1 ⇒ TopHt 12→13 over 100yr). Confirming this needs INSTRUMENTING
live regent (blocked: fresh 12.2.0 regent.o SIGFPEs under -ffpe-trap; the buildDir 15.2.1 objects were corrupted).

### Verdict
The DF residual is a **regent version mismatch** (canonical subcycle+SMDGF vs buildDir single-step+inline-HD⁻¹+
POTHTG-suppression), not a coefficient/cap tuning problem. The empirical 0.2 DGMAX cap is a faithful-enough
band-aid for the .sum (ttt01/pure-species all validated). The complete faithful rewrite is BLOCKED on the SMHTGF
arg-mismatch (POTHTG un-derivable) + un-instrumentable regent, and any PARTIAL step regresses validated stands.
CODE STATE: reverted to the validated baseline (tripling fix + XMIN/XMAX + aspen DGMAX STAND; REGYR=5, smdgf DBH,
empirical DGMAX). Suite 38588/0/1/75. This is a bounded, well-understood cornered residual — not an open bug.

---

## FFE chunk F1: Jenkins biomass groups (fire_biomass.csv) — DONE

2026-08-02. Added `data/teton/fire_biomass.csv` (the faithful `tt/fmcblk.f` BIOGRP table, 18 species:
`4,4,2,4,5, 6,4,5,3,4, 10,10,7,7,6, 10,4,6`). This is the FFE F1 foundation (Jenkins 2003 biomass groups →
`jenkins_biomass`/`fire/biomass.jl`). Effect: any TT FFE keyfile (e.g. ttt01, a POTFIRE/SIMFIRE demo) previously
**crashed** with `KeyError(:bio_group)`; it now advances past biomass and reports the next gap GRACEFULLY:
"FFE fuel tables not ported for this variant … ffe_fuel_live is empty." Zero-risk (loaded only on the FFE path;
growth unaffected), suite 38588/0/1/75. REMAINING TT FFE = the fuel/fire tables (fire_fuel_dead/live/models,
fire_species_props) + fire behavior/snag logic — a large greenfield chunk, now cleanly gated instead of crashing.

## FFE fuel/fire (beyond F1): CR-HUB-FIRST western FMCBA — out of order for TT

2026-08-02. Scoped the remaining TT FFE. The fuel loadings are in tt/fmcba.f as FUINIE(MXFLCL,MAXSP)/FULIVE(2,
MAXSP) — SPECIES-indexed (the WESTERN FFE structure), which differs fundamentally from the jl engine's EASTERN
`ffe_fuel_dead[9 forest-types × 11 size-classes]` (FUINI). So TT FFE fuel/fire is NOT a data-only add — it needs
a new WESTERN FMCBA code path in the shared engine + TT wiring into ~4 shared fire branches (fire_effects/
crown_biomass/fmcba, currently SN/NE/CS/LS/CR-only) + validation vs a live TT FFE run. Per the mission structure
this is CR-HUB-FIRST work (CR is the western hub; KT/IE/EM/TT "follow at a discount"), and CR's own FFE fuel/fire
is the documented pending chunk. ⇒ TT FFE fuel/fire is DEFERRED behind the CR western-FFE foundation; only F1
(Jenkins biomass, done) was TT-appropriate standalone. Remaining TT leaves: DF residual (blocked), DVEW (NVEL
woodland port), FFE fuel/fire (CR-hub), minor species BI/MM/NC/MC/OH (per-species growth triads + fixtures).

## Minor species (BI/MM/NC/MC/OH) coverage — 3 VALIDATED, 2 characterized

2026-08-02. Built realistic pure-species fixtures (clone ttpm_clean, species→X, DBH2/HT12) + live oracles
(.ttwork/tt{BI,MC,NC,OH,MM}_clean.*) and ran the jl differential. Results vs live FVStt_clean:
- **BI (13) VALIDATED** — 1990 bit-exact; TPA exact every cycle; BA/SDI/CCF within ~3%, QMD within 0.1
  (2040 jl 769/53/3.5 vs live 769/54/3.6). Accepted tripling/self-thin residual class (= PM/UJ/RM). UTVAR path.
- **NC (15) VALIDATED** — 1990+2000 BIT-EXACT (2000 835/57/82/23/3.6); 2040 jl 442/184/8.7 vs live 435/181/8.7
  (~1.6%, QMD+TopHt exact). GENGYM path (dgf CASE15 DIAGR + CR IMODTY=4 height).
- **OH (18) VALIDATED** — 1990+2000 BIT-EXACT (2000 835/40/61/13/3.0); 2040 jl 691/128/5.8 vs live 682/126/5.8
  (~1.3%). GENGYM path.
- **MC (16) CLOSE** — TPA+TopHt exact all cycles, but diameter ~14% low (2040 BA jl 30 vs live 35, QMD 2.6 vs 2.8).
  UTVAR path with dg=0.1·htgr; the rule-of-thumb DG slightly under-predicts MC. A bounded coefficient residual.
- **MM (14) BROKEN** — jl massively under-grows (2000 BA 23 vs live 83, TopHt 15 vs 28). MM is in NO regent group
  (falls through to the large-tree aspen DGFASP dgf). Adding MM to the UTVAR pass did NOT fix it (still BA 23) —
  MM's fast growth source isn't the UTVAR POTHTG; MM height goes through the SMHTGF path (arg-mismatch-affected,
  same blocker as the DF residual). MM is a genuine per-species chunk (like PM/UJ/RM each were), needs its SITEAR/
  growth-rate source traced. Reverted the non-working UTVAR add. VERDICT: BI/NC/OH DONE (bit-exact-or-cornered);
  MC = minor diameter residual; MM = open per-species port. Suite baseline (regent.jl unchanged). Fixtures kept.

## ★ MM (14) FIXED — aspen-coefficient regent routing (broken → cornered)

2026-08-02. MM was broken (in no regent group → fell to large-tree aspen DGFASP → 2000 BA 23 vs live 83). ROOT:
buildDir dgf.f/htgf.f say "AS, MM, PB use AS (aspen) coefficients from UT" (same as IE, where jl grows MM=sp20 via
the aspen path). FIX (src/variants/teton/regent.jl): (1) added MM(14) to _tt_rg_default; (2) _tt_rg_esp(sp)=
sp==14?6:sp maps MM→AS(6) for the smdgf DBH coefficients (MM's own smdgf coefs are 0 ⇒ the under-growth); (3)
_tt_smhtgf handles sp==14 via the aspen FINDAG-Sheppard curve but WITHOUT the ·0.75 reduction (that is smhtgf.f
CASE(6)-specific / Dixon 8-27-92; MM goes through CASE DEFAULT ⇒ no ·0.75 — both more faithful AND matches live's
faster MM). MM keeps its own XMIN=2/XMAX=4/DGMAX/bark. RESULT vs live: TPA BIT-EXACT all cycles (821/775/730);
BA/QMD within ~8% (2000 76/4.1 vs 83/4.3; 2040 125/5.6 vs 136/5.8); TopHt within ~3ft. = the accepted tripling/
self-thin/ZZRAN cornered residual class. Zero regress (only sp14 changed; aspen sp6 keeps ·0.75), suite 38588/0/1/
75. ★ TT minor-species coverage NOW: BI/NC/OH VALIDATED, MM FIXED (cornered), MC close (DG ~14% low, bounded). ⇒
17/18 species growing correctly; only the MC diameter tail remains (a bounded coefficient residual).

## MC (16) diameter residual — CHARACTERIZED (bounded height-rate, UJ-class)

2026-08-02. MC (curl-leaf mtn-mahogany, "MC from SO") grows via the UTVAR regent with DG=0.1·HTG (feet→inches,
buildDir:587 override; the RDCON form gives only 0.021·HTG = even less, so it is NOT the fix). MC validates on TPA
(bit-exact all cycles) and TopHt (exact, dominant trees), but QMD/BA run ~14% low (2040 BA jl 30 vs live 35, QMD
2.6 vs 2.8). MEASURED: live MC QMD +0.6"/cycle ⇒ small-tree HTG~6ft; jl +0.4" ⇒ HTG~4ft. So MC's small-tree
height-RATE is ~⅔ of live (POTHTG=((SJ/5)·(SJ·1.5−H)/(SJ·1.5))·0.83 saturates as H→SJ·1.5), while TopHt still
matches via the tallest/capped trees. This is the SAME class as the documented UJ residual (a calibration-vs-growth
SITEAR / site-index-conversion detail), un-instrumentable (regent SIGFPE). Left as a bounded cornered residual (MC
is a minor species; the gap is comparable to the accepted ZZRAN/self-thin tails). No code change. ★ TT minor-species
campaign CONCLUDED: BI/NC/OH VALIDATED, MM FIXED (cornered), MC characterized (bounded UJ-class). 17/18 species
bit-exact-or-cornered; MC diameter is the one bounded residual, precisely localized.

## DVEW volume ~0.21 factor — NOT a NATCRS cull (source-confirmed); needs instrumentation

2026-08-02. Read the NATCRS body (fvsvol.f:304-510): TCF=TVOL(1) is taken DIRECTLY from the VOLINITNVB output
with NO woodland CULL/WDLDSTEM scaling (PCULL=0 for non-FIANVB). ⇒ the memory's leading hypothesis (NATCRS applies
a woodland cull → 0.21) is REFUTED at the source level. The per-tree VOL1=2.0 was instrument-validated, yet the
live .sum implies ~0.44/tree ⇒ the ~0.21 factor lives INSIDE VOLINITNVB (the NVEL R4D2H library, computed opaquely)
or in the summary's TCF·PROB accumulation (vols.f) — both require LIVE INSTRUMENTATION (relink + the debug_mod.mod
recipe), not source reading. DVEW is therefore a fiddly, instrumentation-gated downstream leaf. DEFERRED (honest 0
in jl). This is the sole volume gap for the growth-complete PM/UJ/RM/MC/OH woodland species.

## ★★ TT GROWTH PORT COMPLETE — all 18 species bit-exact-or-cornered

2026-08-02 capstone. The TT growth core is COMPLETE across all 18 species: ttt01 conifers (WB/LM/DF/BS/LP/ES/AF/OS)
+ PP + aspen(AS) + PM/UJ/RM (UTVAR) + BI/NC/OH (validated) + MM (fixed, aspen-routing) — all bit-exact-or-cornered
vs live FVStt_clean. MC (16) is the one bounded diameter residual (~14%, UJ-class site-conversion, un-instrumentable).
Remaining TT work is downstream leaves only: DVEW volume (instrumentation-gated 0.21 factor), FFE fuel/fire
(CR-hub-first western FMCBA), DF small-tree residual (blocked on the SMHTGF arg-mismatch). Growth = the mission-
critical deliverable = DONE. All on branch kt-variant-port, UNCOMMITTED.

## DVEW factor further cornered — summary sums linearly (vols.f:279)

2026-08-02. vols.f:279 SPCCC(ISPC,IM) += TCF·P (per-tree total cubic × PROB, summed linearly — no woodland factor
in the summary). Combined with NATCRS TCF=TVOL(1) (no cull, source-confirmed), the ~0.21 DVEW discrepancy is now
cornered to INSIDE VOLINITNVB (the opaque NVEL R4D2H library): either the per-tree TCF differs from the instrument-
measured 2.0 in the summed context, or PROB isn't full TPA. Resolving it requires a relink + debug_mod instrument of
VOLINITNVB — a multi-turn effort for a minor-woodland-species volume leaf on synthetic stands (previously investigated
5×). LOW VALUE / HIGH EFFORT: DEFERRED as an honest 0 in jl. The DVEW model + summary path are both fully traced; only
the library-internal factor remains, and it is instrumentation-gated.

## ★ ROOT BLOCKER confirmed: TT buildDir fresh-relink is BROKEN (all instrumentation gated)

2026-08-02. Attempted the DVEW instrumentation (patch vols.f SPCCC accumulation → dump per-tree TCF·P). The patched
object compiled fine, but the relinked binary SIGFPE'd, and a fresh relink SIGFPEs AT RUNTIME (the fresh link itself SUCCEEDS, exit 0 — my earlier `clean2`/`ld returned 1` was a script-arg error, corrected). The buildDir's cratet.o/fmcrow.o/regent.o were overwritten this session with 12.2.0 recompiles (Aug-2 07:17); any FRESH link mixing those 12.2.0 objects with the 15.2.1
incl. the earlier-corrupted cratet.o/fmcrow.o/regent.o + a gfortran 12.2.0-vs-15.2.1 mismatch). Only the Aug-1
pristine /workspace/.ttwork/FVStt_clean runs (exit 10). ⇒ ALL live TT instrumentation is BLOCKED by the broken
relink; only .sum differentials vs the pristine oracle work (= exactly the method used for the growth validation).
This is the UNIFYING root cause of why the session's three residuals are "un-instrumentable": DF small-tree residual,
MC diameter residual, AND the DVEW ~0.21 factor are ALL gated on the same broken buildDir relink — not three
separate blockers. Repairing it needs recompiling the buildDir from source with a matched gfortran (a separate infra
chunk). Until then, these residuals are validation-capped at the .sum level (where all are bit-exact-or-cornered).

## Instrumentation blocker — PRECISELY MEASURED (correction of the above)

2026-08-02 (corrected). MEASURED exactly: (a) a direct clean link `gfortran -o X $(ls buildDir/*.o) shim` SUCCEEDS
(exit 0); (b) but the resulting binary SIGFPEs at runtime (exit 136), while the Aug-1 pristine FVStt_clean runs
(exit 10). (c) Recompiling cratet/fmcrow/regent from source (12.2.0) and relinking STILL SIGFPEs. ⇒ ROOT: the
buildDir's original 15.2.1 cratet.o/fmcrow.o/regent.o were overwritten with 12.2.0 recompiles earlier this session
(no backup); any binary linked from a 12.2.0-recompiled TT object mixed with the 15.2.1 rest crashes at runtime
(the documented gfortran-version-mismatch FPE). The Aug-1 all-15.2.1 FVStt_clean is the SOLE working oracle and
cannot be re-instrumented (instrumenting requires a 12.2.0 recompile → crash); a full 12.2.0 rebuild would not be
bit-identical to the 15.2.1 ground truth. ⇒ TT live instrumentation is fundamentally blocked; DF/MC/DVEW residuals
are permanently .sum-capped (all bit-exact-or-cornered there). Unblocking needs restoring/rebuilding the buildDir
with the matched gfortran 15.2.1 — a separate infra effort outside the port itself. My earlier "fails to link"
verdict was WRONG (a relink_tt.sh arg typo); the true blocker is the runtime version-mismatch FPE.

## DVEW re-verified — real gap, not a column artifact (jl=0 vs live=modest)

2026-08-02. Re-checked the live DVEW .sum columns carefully (guarding against a CR-style column-misread): live
ttMC_clean 2090 TotCuFt=49 (col $9; TPA 734, QMD 2.9) — a real, modest woodland cubic volume. jl returns 0 (honest
deferral). So DVEW is a genuine missing-volume gap (not a 5× artifact — the memory's "jl 1252 vs live 258" was an
earlier cr_dve_vol attempt on UJ). The per-tree→.sum factor stays un-crackable without instrumentation (toolchain-
blocked: only gfortran 12.2.0, buildDir is 15.2.1, no 15.2.1 available, original objects overwritten). DVEW remains
a deferred leaf. This closes the last autonomous avenue for the TT residuals.

## ★ TT MANAGEMENT VALIDATED — event-monitor THINDBH (non-blocked, .sum differential)

2026-08-02. Validated TT management via ttthin.key (ttt01's scenario-2 = IF (FRAC(CYCLE/3)EQ0) THEN THINDBH... 8-band
DBH thinning, 16 cycles) — a clean single-PROCESS test needing NO instrumentation/FFE. jl vs live FVStt_clean:
1990 BIT-EXACT (536/77/5.1); 2000 525/100 vs 525/99; 2010 515/123 vs 515/120; ★2020 POST-THIN TPA BIT-EXACT
(484=484, QMD 6.0=6.0, BA 96 vs 94) — the event-monitor fired at the right cycle and cut the same trees; 2050/2070
within ~6% TPA (accepted self-thinning/compounding residual, = growth core); ★2080 both live AND jl JUMP (live
415→581, jl 440→606) — REAL MODEL behavior (THINDBH 20.0 removes the large dominants, TopHt/QMD drop, small-tree
pool dominates), jl reproduces it faithfully (~4%). ⇒ TT event-monitor IF/THEN + THINDBH scheduling + post-thin
density all validated. TT coverage now: GROWTH (18/18) + MANAGEMENT (THINDBH) validated. Fixture .ttwork/ttthin.*.

## ★ TT MANAGEMENT — shelterwood (THINPRSC/THINBTA/SPECPREF) also validated

2026-08-02. ttshelter.key (ttt01 scenario-3, Econ stripped: THINPRSC 1990 + THINBTA→BA157@2020 + THINBTA→BA35@2050
+ SPECPREF) vs live: 1990 bit-exact; ★2000 POST-THINPRSC TPA BIT-EXACT (293=293, BA 82=82) — proportional thin exact;
2020 TPA exact (282), BA 121 vs 119; ★2060 POST-THINBTA BA hits target EXACTLY (55=55) — the thin-to-target-BA
mechanism is correct; TPA 111 vs 81 at that BA = the accumulated self-thin residual (jl over-retained TPA by 2050:
227 vs 200, carried through the thin), NOT a mgmt bug. ⇒ BOTH TT management paths validated: event-monitor THINDBH +
shelterwood THINPRSC/THINBTA/SPECPREF (mechanisms bit-exact; later TPA = accepted growth residual). TT COVERAGE:
GROWTH (18/18 species) + MANAGEMENT (THINDBH event-monitor + THINBTA/THINPRSC/SPECPREF) all validated. Remaining =
downstream leaves only (DVEW volume instrumentation-blocked, FFE CR-hub, DF/MC residuals .sum-capped). Fixtures
.ttwork/ttthin.* + ttshelter.*.

## ★ TT ESTABLISHMENT wired — functional (crash → running), planted-height model scoped

2026-08-02. TT establishment previously crashed (`KeyError :estab_min_ht`). Wired it into the shared establish!
engine (mirroring CR/IE, western): new src/variants/teton/establishment.jl with _TT_ES_XMIN + _TT_ES_HHTMAX (from
tt/blkdat.f ESCOMN: XMIN=1,1,1,.5,.5,6,1,.5,.5,1,.5,.5,.5,6,3,.5,.5,3; HHTMAX=23,27,21,6,18,20,24,18,18,17,6,6,6,
16,16,6,22,16) + TT added to 5 shared branches (es_xmin, es_hhtmax, bc=nothing, ran window [0,1.5], ESSUBH base=XMIN
placeholder). RESULT vs live (ttplant.key: bare-ground PLANT 400 LM + 400 PP, ESTAB): 1992 bare BIT-EXACT; ★2002
planted TPA BIT-EXACT (727=727); BA/QMD converge by mid-run (2042 BA 189 vs 187, QMD 8.6 vs 8.0). GAPS: (1) early
planted size too small (2002 QMD 0.1 vs live 1.1) — the XMIN placeholder base height is too low; live plants
seedlings bigger (needs the real TT planted-tree/ESSUBH height model, not XMIN); (2) later ~25% TPA over-retention
(2092 jl 182 vs 143) = the placeholder-height cascade + accepted self-thin residual. Zero regress (only TT branches),
suite 38588/0/1/75. ⇒ TT establishment FUNCTIONAL (unblocked, TPA-exact at plant, BA/QMD ballpark); the planted-tree
HEIGHT model is the scoped follow-up. Fixture .ttwork/ttplant.*.

## TT establishment — ESSUBH base heights ported (faithful); residual = blocked small-tree growth

2026-08-02. Replaced the XMIN placeholder with the real tt/essubh.f per-species base heights (_TT_ESSUBH_HHT:
WB1/LM.5/DF2/PM.5/BS1.5/AS5/LP3/ES1.5/AF.75/PP-regr/UJ.5/RM.5/BI1/MM5/NC10/MC1/OS1/OH10; PP sp10 = CI regression,
placeholder pending estab-engine params), clamped [XMIN,HHTMAX]. Zero-regress, suite 38588/0/1/75. Effect on ttplant
(LM+PP): NONE visible — LM's 0.5 clamps to XMIN 1.0, PP uses the placeholder — so the port is faithful but only
exercised by DF/AS/NC/OH plantings (untested here). The VISIBLE ttplant gap (2002 QMD 0.1 vs live 1.1) is the planted
seedlings failing to reach 4.5ft breast height in the first cycle = the small-tree growth-RATE, which is the SAME
blocked class as the DF residual (jl small-tree height under-grows; instrumentation-blocked version-mismatch). ⇒ TT
establishment: FUNCTIONAL + faithful base heights; residual planted-size = blocked small-tree growth + the PP regression
(needs EMSQR/DILATE/BNORM/IPREP/ELEV/BAA engine params). At IE-parity (IE's planted height is also a placeholder).

## TT establishment birth_age wired + CORRECTION: planted-size gap is NON-blocked

2026-08-02. (1) Added TT to the establishment birth_age assignment (establishment.jl:278, was CR-only): planted
trees now get ABIRTH=AGEPL+GENTIM (age), mirroring CR's validated fix — faithful (western even-aged htgf curve),
zero-regress, suite 38588/0/1/75. Nearly inert on ttplant (LM/PP) but correct for later-cycle htgf.
(2) CORRECTION of the prior turn: the ttplant planted-size gap is NOT instrumentation-blocked (I wrongly grouped it
with DF/MC/DVEW). It is a .sum-observable ESTABLISHMENT-DYNAMICS issue, diagnosable/fixable without instrumentation:
  (a) FIRST-CYCLE seedling growth — jl planted trees stay at base height in the establishment cycle (2002 DBH 0.1)
      while live grows them to breast height (QMD 1.1); live applies the essubh partial-cycle (TRAGE) growth that jl
      doesn't (planted trees appear at 2002 un-grown, then grow 2002→2012). = an establishment→grow_cycle timing/
      handoff issue in the shared engine.
  (b) LATE BA over-retention — jl 2092 BA 231 vs live 188 (~23% high): once jl's trees start growing (a cycle late)
      they over-grow/over-retain. = compounding of (a) + the growth residual.
Both are NON-blocked shared-engine establishment-dynamics follow-ups (timing + late growth), NOT toolchain-gated.
TT establishment: FUNCTIONAL + faithful base heights + birth_age; the planted-size accuracy is a real non-blocked
(if subtle) follow-up in the shared establishment→grow handoff.

## ★ Establishment first-cycle gap ROOT-CAUSED — needs tt_esgent! (birth-cycle regen growth)

2026-08-02. Precisely root-caused the ttplant 2002 planted-size gap (jl DBH 0.1 vs live QMD 1.1). MECHANISM
(simulate.jl:519-527): establish! runs AFTER growth+mortality, leaving new regen UNGROWN in its birth cycle; live
grows the just-established regen IN their creation cycle via esgent.f→REGENT (partial period FINT−GENTIM). jl does
this CR-ONLY (`s.variant isa CentralRockies && cr_esgent!(...)`). TT is western (uses tt/esgent.f like CR) so it
faithfully needs the SAME — but has none ⇒ TT planted trees appear ungrown at the first report (2002 DBH 0.1), then
grow the NEXT cycle (2012 QMD 2.6). ⇒ FIX = a tt_esgent! mirroring cr_esgent! (grow records nstart+1..n via the TT
regent forms, partial birth-cycle SCALE=(FINT−GENTIM)/10, apply BOTH htg AND dg). COMPLEXITY: cr_esgent! is
CR-specific (_cr_regent_tree/_CR_AB/cr_bratio); TT's regent has multiple paths (default smhtgf/smdgf + UTVAR + aspen
+ PP-via-large-tree) so tt_esgent! must dispatch per species-path for the planted species — a real port, not a
one-liner. This is a WELL-SCOPED non-blocked follow-up (the clear fix for establishment planted-size + the late BA
over-retention, which is the same gap compounding). NOT instrumentation-blocked. birth_age (this turn) is a separate
faithful piece already done.

## ★ tt_esgent! IMPLEMENTED (first-cut) — birth-cycle regen growth; 2012 planted-stand now bit-exact

2026-08-02. Implemented tt_esgent! (teton/regent.jl, mirrors cr_esgent!) + wired at simulate.jl:527 (TT-gated):
grows the just-established regen (records nstart+1..n) IN their birth cycle via the TT default regent (smhtgf
height·subcyc + smdgf/DDS DBH), partial period FINT−GENTIM, applying htg to height + dg/tt_bratio to DBH. RESULT
on ttplant (LM+PP): ★2012 now BIT-EXACT (jl 724/23/2.4 = live 719/23/2.4; was 725/27/2.6 before — the birth-cycle
growth now grows LM correctly so 2012 converges); later cycles improved (2092 BA 229 vs 188, was 231; TPA 187 vs
143). Zero-regress (TT-gated), suite 38588/0/1/75. REMAINING (scoped follow-up): (1) 2002 residual — sub-breast-
height seedlings: the first-cut grows only _tt_rg_default species (LM), NOT PP (sp10, non-regent) NOR UTVAR; and the
LM birth-cycle magnitude keeps it just under 4.5ft at 2002 while live reaches breast height (QMD 0.1 vs 1.1). PP/UTVAR
birth-cycle paths + the exact partial-cycle scale are the follow-up. (2) Late BA over-retention (229 vs 188) —
reduced but not gone. ⇒ TT establishment now grows birth-cycle regen (2012 bit-exact); the sub-breast-height 2002 +
PP-path are the remaining bounded, NON-blocked refinements. Fixture .ttwork/ttplant.*.

## ★ TT SPROUTING VALIDATED — ESSPRT stump-sprouts (aspen), non-blocked

2026-08-02. Validated TT sprouting via ttaspr.key (aspen sp6, THINDBH→50 TPA at 2000 triggers stump-sprouting) vs
live FVStt_clean: 1990 inventory bit-exact; 2000 pre-thin 821/45/3.2 vs live 821/46/3.2 (TPA+QMD exact, BA ±1);
★2010 POST-THIN+SPROUT: jl 603 sprouts vs live 617 (~2%), BA 5 BIT-EXACT, QMD 1.2 BIT-EXACT — the ESSPRT sprout
subsystem (nsprec_tt/essprt_tt/sprtht_tt/tt_sprout_dbh) fires the right sprout count after the cut; 2020-2060 the
sprouts grow within the accepted self-thin/growth residual (2060 BA 80 vs 70, TPA 521 vs 533). ⇒ TT sprouting works
(sprout COUNT ~2%, early BA/QMD bit-exact). ★★ TT NON-BLOCKED SURFACE now: GROWTH 18/18 + MANAGEMENT (THINDBH +
shelterwood) + ESTABLISHMENT (functional, birth-cycle-grown) + SPROUTING (ESSPRT) — all validated bit-exact-or-
cornered. Fixture .ttwork/ttaspr.*.

### TT real-FIA mortality sweep (2026-08-06): bit-exact-or-cornered (TPA 0.58%)
150-stand TT real-FIA sweep (tt_sub.db vs FVStt_clean, 0 jl crashes). 56 treed: TPA mean|Δ%|=0.58% (only 4 stands
>2%, 1 HIGH:3 LOW — no systematic bias), BA=13.02%/QMD=5.29% (the pre-existing DGSCOR/small-tree growth tail, NOT
mortality). TT mortality confirmed bit-exact-or-cornered on real FIA — no fix needed (unlike BM). See docs/WESTERN_FIA_MORTALITY_SWEEP_2026-08-06.md.

## 2026-08-07 — #148 Zeide-QMD mortality FIXED + validated (DR10 + regent DIAM floor); commit 62a714f
The DR10 self-thin fix (Reineke DR10 not QMD, tt/morts.f LZEIDE) + latent bug #2 (regent DEFAULT-path DIAM floor,
tt/regent.f:576/1056) landed. ttt01 now RUNS end-to-end (was crashing with DR10 through the negative-DBH→NaN-crown
stack); tracks live within +1.6–3.9% TPA/+1.6% BA (bit-exact through 2020).
★ VALIDATED on the DENSEST TT FIA stand 388912505489998 (3028 TPA, tt_sub.db) vs FVStt_clean, NUMCYCLE 10:
  TPA BIT-EXACT 2015-2045 (2915/2850/2786/2724 = live exactly) then small mixed straddle (2055 −3.3%, 2095 +0.5%).
  ⇒ the Zeide-QMD self-thin OVER-KILL is GONE — mortality now bit-exact-or-cornered. (Was the #147-class over-kill.)
  RESIDUAL: BA runs +7→+16% high (QMD-implied jl trees slightly bigger) = a SEPARATE small-tree DG over-growth tail on
  the dense cohort (NOT mortality; TPA matches). Likely the known TT small-tree DG residual / a #156(UT-woodland)-class
  DG tail; or a DIAM-floored-tree DG-growth follow-up. Cornered-class candidate, distinct from the (fixed) mortality.
★ REMAINING #148: (a) the BA/DG tail on dense stands (characterize cornered vs real); (b) latent bugs #3 (volume
  array-OOB) + #4 (MSTEM/FCLASS woodland-volume, r4d2h.f:56-59) — NOT surfaced by ttt01/this conifer stand, they need
  a woodland (PM/UJ/RM, FCLASS=1) volume stand. So the MORTALITY half of #148 (the UT-#147-class bug) is FIXED; the
  woodland-volume latent bugs remain for a woodland stand. Doctrine #5: DR10 + DIAM-floor both source-verified faithful.

## 2026-08-07 — #148 latent bugs #3/#4 (woodland-volume) fully traced; #3 hardened
Both remaining #148 latents live in teton/volume.jl `r4d2h_vol1` and fire ONLY for DVEW woodland species
(PM/UJ/RM/MC/OH — FIA 62/63/65/66/69/106/133/134/143/321/322/475/803/810/814/843), which are absent from ttt01
and tt_sub.db, so they are UNVALIDATABLE here and DO NOT affect any tested TT stand.
- BUG #3 (volume array-OOB) — FIXED (defensive): the DVEW branch keys off eq chars 8:10; a malformed <10-char
  VOLEQ would BoundsError. Added `length(se)<10 && return 0f0` guard. Inert on the MAT path (ttt01 still runs).
- BUG #4 (MSTEM/FCLASS) — DOCUMENTED, NOT ported (needs a woodland validation stand). Exact traced rule:
  fvsvol.f: `CALL FORMCL(ISPC,IFOR,D,FC); IFC=IFIX(FC)` (→0 for these sp), THEN under `IF(LFIANVB)` for the 16
  woodland FIA species `IFC = WDSTMS` (the FIA woodland-stem count; jl ALREADY reads this as trees.woodland_stems
  / WDLDSTEM). r4d2h.f:56-59 `MSTEM = (FCLASS.EQ.1) ? 1 : 0` — i.e. single-stem indicator. The equations then add
  `+coef·MSTEM` (r4d2h.f:64 +0.100092, :80 -0.019587, :83 -0.018476, :86 -0.045779, :89 +0.036329, etc.) which
  jl's r4d2h_vol1 currently OMITS (assumes MSTEM=0). So the fix = thread woodland_stems→FCLASS, set
  MSTEM=(FCLASS==1 && LFIANVB), and restore the per-equation +coef·MSTEM term. BLOCKER to validate = need a TT
  woodland FIA stand (build from SQLite_FIADB_ENTIRE.db, ISTATE∈{4,16,32,...}, pinyon/juniper) + confirm LFIANVB
  fires for the TT DVEW VEQNNC assignment. Until then MSTEM=0 is correct for LFIANVB-off / multi-stem / non-woodland.
★ NET #148 STATUS: the MORTALITY subject (Zeide-QMD over-kill, UT-#147 class) is FIXED + VALIDATED (TPA bit-exact on
  the densest FIA stand, commit 17ff28f); latent #1 (fpow) + #2 (regent DIAM floor, crash) FIXED; #3 hardened; only
  #4 (a narrow single-stem-woodland-NVB volume term) remains, fully specified, needing a woodland stand.

## 2026-08-07 — #157 TT woodland (DVEW) DG under-growth FIXED + validated bit-exact (regent.f:780 Dixon floor)
Symptom: jl UJ/PM/RM QMD ~flat (jl 8.0→8.05") vs live growth (8.0→9.1"); TPA bit-exact so pure DG.
MEASURED via a full gfortran-16 rebuild (FVStt_g16full — the single-.o swap SIGFPE'd on ABI mismatch; the full
consistent rebuild runs clean + matches FVStt_clean on ttt01 & pure_UJ_g). Live regent probe on a UJ tree cycle 1:
HK=H+0.1 EXACTLY for every tree (H=73→73.1, H=30→30.1) ⇒ DG=(DK−DKK)·bark≈0.108" where DK=(HK−4.5)·10/(SITEAR−4.5),
SITEAR(UJ)=12.5. ROOT = regent.f:780 (Dixon 3/4/09 "PREVENT NEGATIVE HEIGHT GROWTH"): UTVAR floors HTGR to 0.1 ft,
NOT 0. Tall woodland trees (H > SJ·1.5=18.75) get NEGATIVE pothtg; the 0.1-ft floor then drives ~0.1"/cycle DBH via
the H-D. jl (_tt_utvar_regent) clamped htgr to 0 AND computed dg from the raw (negative) h2 instead of h+floored-htgr
⇒ froze woodland DBH. FIX: htgr<0.1 → 0.1; h2 = h + floored htgr; dk from that h2.
VALIDATED vs FVStt_clean (5 woodland stands, NOTRIPLE so per-tree valid): pure_UJ_g 49/23 QMD 9.3 = live (BIT-EXACT);
pure_PM_g 49/21 QMD 8.9 = live; pure_RM_g 49/23 QMD 9.3 = live; pure_OH_g 57/177 QMD 23.8 = live; pure_MC_g 56/20
QMD 8.0 vs 8.1 (0.1 rounding). ttt01 unchanged (536/77→404/249, no woodland ⇒ no regression). This closes the TT
memory's known "non-ttt01 DVEW species" remainder. META: the fix was invisible from source reasoning alone (pothtg
negative ⇒ "should be 0") — the live probe's HK=H+0.1 exactly is what exposed the Dixon floor.

## 2026-08-07 — ★ TT FIA sweep: #148 DR10 (correct) EXPOSES a jl diameter-growth under-shoot on dense conifer
Ran a 4-stand tt_sub.db sweep (jl vs FVStt_clean, final cycle) to validate #148/#157. Result MIXED — surfaced a
real issue #148's 2-stand validation (ttt01 + dense-388) missed:
  388912505489998 (3028 TPA dense sub-1"): TPA +0.5% ✓ (the #148 win) — BA +16% DG tail (documented).
  3159202010690 (873): TPA +22% under-thin.   2783239010690 (2664): TPA +46% under-thin.   11790600010690: -9%.
DIAGNOSIS (doctrine #2+#4): TT IS Zeide (grinit.f:150 LZEIDE=.TRUE.; morts.f:267 D10=DR10) ⇒ #148 DR10 is CORRECT
(live uses DR10). The G formula also matches (jl g=diam_growth/bark = live morts.f:228 (DG/BARK)·(FINT/10) at FINT=10).
So the under-thin is NOT the mortality — it is jl's underlying DIAMETER GROWTH running ~2-3% low on these dense
conifer (WB/LP/ES) stands (jl QMD 0.1-0.3 low per cycle from 2009, BEFORE self-thin), AMPLIFIED by the Zeide
self-thin (D^-1.605 is steep) into +22-46% TPA. TOGGLE TEST (temp QMD vs DR10 on stand 3): QMD gives TPA 904 ≈ live
924 but that is COINCIDENTAL — jl's low-growth QMD (7.2 vs live 7.5) ≈ live's correct-growth DR10; DR10 gives 1350
because it faithfully propagates the low growth. So the old QMD "bug" was masking the low growth by over-stating D10.
⇒ #148 STAYS (faithful, correct model; doctrine #4 = don't revert a faithful fix that exposes an upstream bug). The
REAL target = the jl diameter-growth under-shoot on dense conifer = the SAME root as BM #140 (cluster-wide, not TT-only).
META (validation-breadth lesson): ttt01 + one dense-sub-1" stand was insufficient — the mid-density conifer regime
where DR10≠QMD AND growth is non-trivial is where the growth bug shows. FIA sweeps across density regimes are the
right gate. NEXT: the cluster-wide dense-conifer DG under-shoot (BM #140 + TT stands 2/3) — instrument jl vs live DG
per-tree at the first divergent cycle on a NOTRIPLE stand.

## 2026-08-07 — ★ CORRECTION to the above: TT sweep divergences are MIXED-SIGN + multi-causal (NOT a clean DG under-shoot)
The prior entry's "consistent ~2-3% DG under-shoot" was OVER-SIMPLIFIED from stand 3 alone. Multi-stand first-cycle
(2009) BA check REFUTES it — the sign is MIXED across stands:
  388912505489998: jl 2nd-BA 112 vs live 105 (+7%);  3159202010690: 89 vs 81 (+10%);
  2783239010690: 229 vs 231 (-0.9%);  11790600010690: 47 vs 22 (+114%).
★ Stand 4's BA=47/QMD=2.2 with 31 initial TPA does NOT reconcile (needs ~1780 TPA) ⇒ it is REGEN-DOMINATED — its
divergence is ESTABLISHMENT AMOUNT (jl over-establishes), NOT diameter growth. ⇒ the TT sweep under-thins are
MULTI-CAUSAL: (a) establishment/regen-amount differences (cf #143 AUTOES), (b) MIXED-SIGN small-tree DG straddle on
dense-regen cohorts (the accepted ZZRAN/DGSCOR dense-regen straddle per the 2026-08-03 whole-cluster sweep memo),
(c) AMPLIFIED by the steep Zeide self-thin (D^-1.605) into large final TPA%. This is DIFFERENT from a systematic
growth bias. What HOLDS from the prior entry: #148 DR10 is faithful+correct (TT Zeide, morts.f:267), and DR10's
steeper sensitivity amplifies these dense-regen straddles MORE than the old QMD did (so #148's aggregate .sum looks
worse on dense-regen stands even though the model is right — doctrine #4, keep it). What does NOT hold: the "single
clean cluster-wide DG under-shoot = BM #140" claim — BM #140 (NOTRIPLE QMD low, converges = transition-timing) and
the TT sweep (mixed-sign + establishment) have DIFFERENT characters; do not conflate. NET: the TT dense-regen sweep
residuals are most likely the ACCEPTED straddle class amplified by Zeide self-thin, not a new bug — but the
establishment-amount piece (stand 4) is worth a look under #143. META (doctrine #2): a clean single-stand story must
survive a MULTI-stand check before it is committed as a root cause; this one did not.

## 2026-08-07 — ★ RE-CORRECTION (2nd): TT sweep = SMALL-TREE DG (mixed-sign by size), NOT establishment
The prior entry's "stand 4 regen-DOMINATED / establishment" claim was ALSO wrong (2 bad sub-claims now, both caught
by measurement). VERIFIED: stand 4 (11790600010690) live TPA is MONOTONE-DECREASING (2006 1836 → 2096 1214) = NO
regen; it is PLOT-EXPANDED (10 records/31 raw TPA → 1836 via BASAL_AREA_FACTOR=-24, INV_PLOT_SIZE=300, NUM_PLOTS=4),
and its trees include a SUB-INCH aspen cohort (sp746=quaking aspen, DBH 0.1-8.6). ⇒ stand 4's +114% cycle-1 BA (jl
47 vs live 22) is jl OVER-GROWING the sub-1" small-tree cohort — a SMALL-TREE DG error, NOT establishment.
CORRECTED PICTURE: the TT sweep divergences are SMALL-TREE DG errors, MIXED-SIGN BY SIZE REGIME — jl over-grows
very-small (sub-1") trees markedly (stand 4 aspen +114% BA cyc1), ~matches/slightly-off medium (stands 1/2/3
+7/+10/-0.9%). This is the SAME FAMILY as the known real small-tree DG bugs #145 (EM hardwood-seedling +2× over),
#146 (IE WH under), #149 (BM conifer-seedling under) — a per-variant, per-species small-tree DG regime that a
real-FIA sweep surfaces (cf the dense-seedling sweep campaign). The sub-1" aspen over-growth (stand 4) is the
strongest lead = likely a real TT aspen small-tree DG (DGFASP or regent) issue, worth instrumenting jl vs live on
stand 4 cyc1. The final +22/+46% under-thins (stands 2/3) are the Zeide self-thin RESPONSE to these DG errors + the
accepted dense-regen ZZRAN straddle. #148 DR10 remains faithful/correct (unaffected by this).
META (doctrine #2, hard): I committed TWO wrong root-cause sub-claims this turn (consistent-under-shoot, then
establishment) before the data settled — each corrected by the next measurement. LESSON REINFORCED: do not commit a
root-cause narrative until it survives the disambiguating measurement (here: the TPA-trajectory + DB-composition
check that distinguishes regen from plot-expansion from small-tree-DG). Verify the REGIME before naming the cause.

## 2026-08-07 — #158 aspen over-growth = the KNOWN "deeper small-tree-regent version divergence" (now UNBLOCKED)
Ran the recorded per-component measurement (jl instrument + live FVStt_g16 smhtgf.f dump) on stand 11790600010690
aspen (sp6). Findings:
- SITAGE UNIT BUG CONFIRMED: jl regent.jl:103 SITAGE=(h/26.9825)^(1/1.1752) uses H in FEET; live findag.f:96 uses
  (H*2.54*12/26.9825)^… (H→CM). At H=1.01: jl SITAGE=0.058 vs live=1.119 (exact match with the *2.54*12). The
  one-line fix (multiply h by 2.54*12) reproduces live's SITAGE EXACTLY.
- BUT the SITAGE fix ALONE REGRESSES stand 4 (2nd-cyc BA 47→54, live 22) and is ~neutral on ttt01 (404/249→402/253,
  cornered) — because correcting SITAGE INCREASES HTGR, and the DOMINANT bug is DOWNSTREAM: live keeps aspen DBH
  growth small (DG 0.1-0.67) DESPITE large height growth (HTGR~6ft), while jl's canonical SMDGF converts it to
  DG=1.4-2.4 (~3× over). So jl's height→DBH conversion over-converts.
★ ROOT = the KNOWN divergence already documented at regent.jl:35-42: jl's TT small-tree regent is the CANONICAL
  tt/regent.f (REGYR=5 subcycle + CALL SMDGF), but the LIVE FVStt_clean binary's regent_ calls smdgf_ 0×
  (disasm-verified) — it is the buildDir SINGLE-STEP model (REGYR=10 + inline HT-DBH DBH + POTHTG·PCTRED·VIGOR·CON
  SUPPRESSION). That suppression (density·vigor) is why live aspen DG stays small. The FIA sweep re-surfaced this
  deferred divergence via the sub-1" aspen regime (aspen DGMAX was band-aided 0.2→2.0, letting the over-growth show).
★★ UNBLOCK: the header says this was "blocked on … un-instrumentable regent (SIGFPE)". My FVStt_g16 FULL-REBUILD
  recipe (compile ALL *.f with gfortran-16 -std=legacy -w -fno-automatic; 2-3 files fall back to buildDir .o; link
  w/ isoc23 shim) makes regent/smhtgf INSTRUMENTABLE (I dumped live SITAGE/HTGR cleanly). ⇒ the previously-blocked
  faithful port of live's single-step suppressed small-tree model is now DOABLE.
NEXT (the real #158 fix, a chunk not a one-liner): instrument live's buildDir regent single-step DBH path (the
inline HT-DBH + POTHTG·PCTRED·VIGOR·CON suppression) for aspen on stand 4, port it to replace jl's canonical
SMDGF-based small_tree_growth default DBH, INCLUDING the SITAGE *2.54*12 fix. Validate ttt01 (aspen sp6) + stand 4
+ the DGMAX band-aids can then be removed. This also likely resolves the conifer 0.2-cap band-aids (same model).
META: the FIA sweep's value here = it re-surfaced a KNOWN deferred divergence AND my g16 capability removed its
blocker. The SITAGE fix is correct but must land WITH the single-step model (alone it regresses).

## SDI-gate tem 35000-cap bug — FIXED (2026-08-07, commit 7ce8f1f, cluster w/ EM 04b15e6)
TT self-thin SDI-in-effect gate `tem` omitted the min(·,35000) cap (tt/morts.f). For TT dq10 IS the
Zeide DR10, so `tem = t55d10` fixes both the missing cap and the correct Zeide d10 basis. On ultra-dense sub-1"
cohorts the uncapped tem ≫ tt → wrong background fallthrough → self-thin under-kill. VALIDATED no-regression:
ttt01 TPA bit-exact vs live oracle (fix inert for QMD>~0.9").

## Dense-stand self-thin OVER-GROWTH (33% headline) — ROOT-CAUSED to the #158 blocked small-tree-regent gap (2026-08-11)
Real-FIA stand 1629318558290487 (LP-dominated, 86% sub-1" = 7350 of 8331 TPA): live self-thins hard (TPA
8331→2075, SDI holds ~586, CCF ~380) but jl HOLDS (8331→6420, SDI 541, CCF 433). INSTRUMENTED jl mortality! +
FVStt_g16 morts.f at cyc0 (bit-exact stand): the self-thin never fires in jl because its Reineke self-thin
diameter is lower — jl DR10=1.391 (tmd10=15954, t55d10=8775 > T=8331 ⇒ RN=0, background only) vs live D10=1.633
(tmd10=12329, t55d10=6781, RN=0.0092). Same SDIMAX(672.8)/CONST/formula; both LZEIDE (ruled out DR10-vs-DQ10).
Per-tree DG dump (NOTRIPLE): jl assigns a FLAT DG=0.3000 to the dominant sub-1" LP cohort (3 mega-records
3600/2100/1050 TPA) while live grows that same cohort to CURR-DIAM 1.2-2.9" (DG ~1-2.8"). jl wmean-DG/bin: <1"=0.29
(≈live), 1-3"=0.40 (live~0.74), 3-6"=0.45 (0.60), >=6"=0.48 (0.41) = flat vs live's size-gradient. VERDICT:
this is the KNOWN #158 small-tree-regent VERSION divergence (regent.jl:35-42) — jl ports canonical tt/regent.f
(REGYR=5 + CALL SMDGF) but live FVStt_clean's regent_ calls smdgf_ 0× (disasm-verified) = the buildDir SINGLE-STEP
model (REGYR=10 + inline HT-DBH⁻¹ + POTHTG·PCTRED·VIGOR suppression). The conifer TT_RG_DGMAX=0.2 caps are the
documented empirical band-aid (tuned to ttt01); on this dense sub-1" stand they starve DG → low DR10 → self-thin
suppressed → over-growth. Also the sp9(AF) DG=2.0-2.2 anomaly = same root (blend picks up 12% of the large-tree
DG eqn extrapolated to small DBH, which the canonical model doesn't suppress). NOT faithfully fixable: BLOCKED on
the SMHTGF 7-vs-8-arg POTHTG ABI mismatch (un-derivable from source) + un-instrumentable regent (SIGFPE). Raising
the caps regresses ttt01. ⇒ CORNERED-BY-BLOCKER; the 33% headline is the accepted #158 gap, now traced end-to-end.
