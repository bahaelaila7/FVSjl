# EM (Eastern Montana) variant port — audit log

The next western Rockies cluster variant after CR/KT/IE. EM is a **standard western Wykoff DDS**
(`em/dgf.f`, no `CALL GEMDG`) — a **discount of the KT/IE engine**, not a fresh port. Validated
bit-exact-or-cornered vs live `FVSem` per chunk (same doctrine as CR/KT). Oracle:
`/workspace/.emwork/FVSem_clean` (relink_em.sh), verified vs `tests/FVSem/emt01.sum.save`.

## Foundation (established)
- **Oracle** relinked from `bin/FVSem_buildDir/*.o` + isoc23 shim; VERIFIED bit-match to
  `emt01.sum.save` (1990 536/77, 2000 526/96, 2010 517/114).
- **MAXSP=19**, **YR=10** (IFINT=10). grinit IDENTICAL to KT: LZEIDE=.FALSE. (Stage SDI), DGSD=2.0,
  IFINTH=5, LHTDRG=.TRUE., seed 55329.
- **Species (19):** WB WL DF LM LL RM LP ES AF PP GA AS CW BA PW NC PB OS OH.
  FIA 101 073 202 113 072 066 108 093 019 122 544 746 747 741 745 749 375 299 998.
- **Reuse:** `em/dgf.f` uses the identical Wykoff coefficient structure as KT (DGLD/DGBAL/DGCR/
  DGCRSQ/DGDBAL/DGHAB/DGFOR/DGCCF + DGCONS), differing in dims: DGHAB(8,MAXSP), DGFOR(6,MAXSP)
  (KT is 9/7). Distinctive behavior: EM barely self-thins (emt01 536→455 TPA by 2060) vs KT's
  hard thin (536→217) on the same trees — a mortality/density path to validate carefully.

## Chunk 0 — scaffold (DONE)
`EasternMontana <: AbstractVariant` singleton + registration:
- `src/variants/easternmontana/easternmontana.jl`: singleton, `variant_code="EM"`, `nspecies=19`,
  `htg_period=10`, `EM_DATADIR` + `coefficients(::EasternMontana)` (lazy — errors loudly until
  `data/easternmontana/` lands in chunk 1).
- `variant_from_code` accepts "EM"/"EASTERNMONTANA"/"EASTERN MONTANA".
- Included in `FVSjl.jl` after the IE block.
- VERIFIED: package precompiles + loads; `EasternMontana()` constructs; `variant_code`/`nspecies`/
  `htg_period`/`variant_from_code` all correct; unimplemented hooks error loudly (doctrine #5).
- Suite: 0 regressions (scaffold is inert to existing variants — new file + one dispatch branch).

## Chunk 1 prep (findings — ready to execute)
- **Crosswalk = ASPT column 10.** `bin/FVSem_buildDir/spctrn.f` holds the shared `ASPT(442,21)` table;
  its column legend (line ~33) is `ALFA FIA PLNT AK BM CA CI CR EC EM …` → cols 1=alpha, 2=FIA,
  3=PLANTS, 4=AK, 5=BM, 6=CA, 7=CI, 8=CR, 9=EC, **10=EM**. SPTRN matches an input code's ASPT row
  (by alpha/FIA/PLANTS), reads the EM column → SPCOUT 2-char target, then matches SPCOUT against the
  19 `NSP` alpha codes for ISPC1 (fallback ISPC1=MAXSP=19 = OS). Extract the 442-row EM column
  PROGRAMMATICALLY from the DATA blocks (j=1-10 & j=11-21 per 10-row group) → build
  `data/easternmontana/species_translation.csv` in the KT 7-col schema
  (code_alpha,code_fia,code_plants,target×4 — fill all 4 target cols with the EM target; KT uses 4 for
  Kaniksu/TallyLake sub-variants, EM has one). Add `spctrn_column(::EasternMontana)=4`.
- **species_coefficients.csv reader** wants ~40 cols (see data/kootenai/ header). Bootstrap with the 19
  species identity (index,alpha,fia,plants — already known) + coefficients extracted per subsystem as
  its chunk lands. FIRST validatable milestone: EM reads emt01.tre, species-decode correct, 1990
  inventory (TPA/BA/SDI/QMD/TopHt) bit-exact vs FVSem_clean (needs identity + htdbh for missing-height
  trees only — NOT growth coefs). CHECK emt01.tre species vs the 19 NSP (it may use EM-distinct species
  WB/LM/LL/RM/GA that the KT/IE 248112 stand lacked). RESOLVED: emt01 IS the 248112 stand — species DF/WL(L)/LP/PP/ES(S) only, ALL shared conifers (EM idx 3/2/7/10/8), NO EM-distinct species. So the inventory milestone needs only the shared-conifer decode = tractable; the 442-row crosswalk + EM-distinct coefs can follow.

## Chunk 1b progress — coefficient extraction (proven tractable)
Approach confirmed: em/*.f DATA statements are clean + species-order documented, and the engine hooks
REUSE the CR/IE/KT patterns (only the DATA differs). First validated column-group:

**Bark (em/bratio.f → bark1/bark2/bark_imap).** EM's BRATIO is IDENTICAL to CR/IE: IMAP(IS) dispatch
(1/2/3); IMAP=1 with BARK1=BARK2=0 → `0.9002 − 0.3089/min(D,19)`; IMAP=2 → BARK1; IMAP=3 →
BARK1+BARK2/D; clamp [0.80, 0.99]. Reuses `bark_ratio` (southern/bark_and_bounds.jl) +
centralrockies/inlandempire diameter_growth bark-coef load (bark_a=−0.3089,bark_b=0.9002 for the
zero-coef species). Species order 1=WB..19=OH. VALUES (verified vs em/bratio.f DATA):
- BARK1 = 0.934 0.934 0.867 0.969 0.937 0.000 0.969 0.956 0.937 0.890 0.892 0.950 0.892 0.892 0.892 0.892 0.950 0.934 0.892
- BARK2 = 0.000×10, −0.086(11 GA) 0.000(12 AS) −0.086(13 CW) −0.086(14 BA) −0.086(15 PW) −0.086(16 NC) 0.000(17 PB) 0.000(18 OS) −0.086(19 OH)
- IMAP  = 2 2 2 2 2 1 2 2 2 2 3 2 3 3 3 3 2 2 3   (only RM=idx6 uses IMAP=1 w/ BARK1=0 ⇒ the 0.9002−0.3089/D default)
Note the documented species-expansion lineage (em/bratio.f header): LM←IE LM(TT), LL←IE AF(NI),
RM←IE JU(UT, really PP from CR), AS/PB←IE AS(UT), GA/CW/BA/PW/NC/OH←IE CO(CR), OS←EM WB.
Remaining 1b columns (site/sdimax sitset.f, small-tree regent.f, ht/htdbh htcalc.f, mort morts.f,
dg_resid_sd dgf.f DGCONS, volume merch specs) extract the same way — clean DATA + existing hooks.

**species.jl DONE** (mirrors KT blkdat init — EM grinit identical: seed 55329, Stage SDI, year=10,
growth_fint=10, ht_drag_sp=true, dg_sd=2.0; spctrn_column=4, other_species=19=OH). Included in FVSjl.jl.
Package precompiles (lazy coefficients — errors loudly until the CSV lands).

**Dependency order MEASURED** (diagnostic: identity+bark CSV, zeros elsewhere, load emt01): EM gets
past species-load + tree-read; the FIRST missing hook is **`site_setup!`** (MethodError) ⇒ the 1990
inventory REQUIRES chunk 2 (sitset/habtyp) before it can run. So the "first validatable milestone"
needs chunks 1 (species data) + 2 (site_setup!) together. After site_setup! the setup will next need
height-dubbing (htdbh) + CCF (crown) hooks/data for TopHt/CCF columns. Extract order: finish 1b data
→ chunk 2 site_setup! → then emt01 1990 inventory differential vs FVSem_clean.

## Chunk 2 scope (site/habitat — the next required hook `site_setup!`)
`em/sitset.f` (228 lines) + `em/habtyp.f` — mirrors KT's site_index.jl + habitat_tables.jl. EM uses the
**30 NI habitat types** (ITYPE index). Tables to extract (all clean DATA in em/sitset.f):
- `BAMAXA(30)` (line 53) — BA-max per habitat type. `MAPSDI` (58) — habitat→SDICON index (9*2,13*5,…).
- `SDICON(9)` (62) = 467,634,696,768,775,751,707,661,635 — SDImax per SDI group.
- `MAPSIT(30,11)` (64) — site index per habitat type × site-species-group. `MAPSS(30)` (88) — site species
  per habitat type (2*10,7*3,2*8,…).
Logic: ISISP (site species) = MAPSS(ITYPE) if unset; SITEAR(I) = MAPSIT(ITYPE, group) via a per-species
SELECT (DF-default historically, now habitat-mapped); SDIDEF(I) from BAMAX/MAPSDI/SDICON. Implement
`site_setup!(::EasternMontana)` + `em_habtyp` (habitat code → ITYPE, the 30 NI types) + the data files.
Then the diagnostic will advance past site_setup! to the next hook (htdbh/crown) — iterate to the emt01
1990 inventory differential. Also still owed for chunk 1b data: site_lo/hi + dbh_max + small-tree (regent.f)
+ ht1/ht2/wykoff_ht2 (htcalc.f) + mort_bkgd (morts.f) + htdbh (htdbh.f) + volume merch + dg_resid_sd
(dgf.f DGCONS) + varmrt_varadj + is_sprouting.

## Chunk 2 COMPLETE (e7d622b) + measured next-blocker
site_setup!(::EasternMontana) done + validated (see above / commit): ITYPE=4, SDIDEF=696 all sp, DF SI=51,
BAMAX=310 — bit-exact vs live FVSem emt01. Also fixed kw_stdinfo! habitat_code gate to include EM (4th such
gate: SN→eco_unit, KT/IE/EM→habitat_code; found running from-scratch, same class as the KT bug).

**Measured (diagnostic identity+bark CSV, run full .sum):** setup_growth! + site_setup! run clean; the
projection then hits `MethodError diameter_growth!` (dgf!) ⇒ the .sum needs **chunk 3 (DG)** to run all
cycles. The 1990 inventory ROW itself is pre-growth but needs height-dub (htdbh) + crown coefs for TopHt/CCF.
So the READ+SETUP+SITE path is validated; next unblockers are chunk 1b (real species_coefficients.csv so EM
loads with the columns setup/inventory touch) + chunk 3 (diameter_growth!(::EasternMontana) — reuse KT
diameter_growth.jl engine, swap EM DG coefs from em/dgf.f DATA + DGHAB(8)/DGFOR(6) dims). Both substantial.

## Chunk 1b coefficient-source map (measured — all clean DATA)
Per-species coefficients for species_coefficients.csv, by source file (all clean DATA, extract like bark):
- **em/blkdat.f**: HT1/HT2 (Wykoff large-tree height → ht1/ht2/wykoff_ht2), SIGMAR (height resid var),
  XMIN (small-tree min DBH), B0ACCF/B1ACCF/B0BCCF/B1BCCF/B0ASTD/B1BSTD (CCF open/stand crown-width —
  eastern-style, feeds ccfcal), OCURHT(16 habitat-grp × MAXSP)+OCURNF(20 NF × MAXSP) (establishment).
- **em/bratio.f**: bark1/bark2/bark_imap (DONE, validated).
- **em/sitset.f**: site tables (DONE chunk 2, validated).
- **em/dgf.f**: DG coefficients (DGLD/DGBAL/DGCR/DGCRSQ/DGDBAL/DGHAB(8,MAXSP)/DGFOR(6,MAXSP)/DGDS/DGEL/
  DGEL2/DGSASP) + DGCONS (dg_resid_sd) — the chunk-3 bulk.
- **em/regent.f**: small-tree st_dgmax/xmax/xmin/diam/htadj/break + regen coefs (chunk 6).
- **em/morts.f + varmrt.f**: mort_bkgd_intercept/dbh + varmrt_varadj (chunk 7).
- merch/volume specs (stump/top_dib/dbh_min/scf_*/bf_*): shared MRULES defaults + em VOLEQDEF (chunk 8).

**Height group extracted+validated vs em/blkdat.f DATA (species order 1=WB..19=OH):**
- HT1 = 4.1539 4.1539 4.4161 4.192 4.76537 3.2 4.5356 4.7537 4.5788 4.414 4.4421×6(11-16) 4.4421(17) 4.1539(18=OS) 4.4421(19)
- HT2 = -4.212 -4.212 -6.962 -5.1651 -7.61062 -5.0 -5.692 -8.356 -7.138 -8.907 -6.5405×6 -6.5405 -4.212 -6.5405
- SIGMAR = 0.11645 0.11645 0.14465 0.4671 0.4345 0.2 0.14465 0.1585 0.14465 0.1342 0.2 0.375 0.2 0.2 0.2 0.2 0.375 0.11645 0.2
- XMIN = 1 1 1 1 0.5 0.5 1 0.5 0.5 1 3 6 3 3 3 3 6 0.5 3

**Measured inventory note:** emt01's 27 trees ALL have measured heights (ht>0) ⇒ the 1990 inventory TopHt
does NOT need htdbh-dubbing (heights from .tre); it needs correct TPA/BA/SDI/QMD (DBH) + CCF (crown-width
B0*CCF) + volume. So the 1990-inventory validation is gated on: real DBH/expansion path + CCF crown coefs +
volume merch — NOT the full growth set. (A BA=8-vs-77 gap seen with the placeholder CSV must be re-checked
with real coefs, not diagnosed on zeros.)

## Chunk 3 (DG) scope — a HYBRID, not a pure KT reuse (measured from em/dgf.f)
EM's `diameter_growth!` dispatches on species into DISTINCT DDS paths (em/dgf.f DO 20 loop):
- **Standard Wykoff** (main conifers 1-3,5,7-10,18 = WB/WL/DF/LL/LP/ES/AF/PP/OS): the classic form
  DDS = CONSPP + DGLD·lnD + CR·(DGCR+CR·DGCRSQ) + DGBAL·BAL + DGDBAL·BAL/ln(D+1) + DGDSQ·... +
  slope/aspect (DGCASP/DGSASP/DGSLOP/DGSLSQ) + DGEL·elev + DGEL2·elev². CONSPP = DGCON(sp)+COR(sp)+
  0.01·DGCCF·RELDEN. DGCON (ENTRY DGCONS) = DGHAB[MAPHAB[kh,sp]] + DGFOR[MAPLOC[loc,sp]] + DGDS[MAPDSQ]
  + elev. (Fuller than KT's reduced form — KT lacks DGBAL/DGDSQ/slope-aspect.)
- **sp6 (RM juniper)**: GENGYM DIAGR path from UT — DF=0.25897+1.03129·DPP−0.0002025464·BA+0.00177·SI;
  DDS=ln(DIAGR·(2·DPP·BARK+DIAGR))+CONSPP.
- **sp12,17 (AS,PB aspen)**: `CALL DGFASP(D,ASPDG,CR,BARK,SI)` (em/dgfasp.f — separate routine to port);
  DDS=ASPDG+ln(COR2)+COR.
- **sp4 (LM)**: Wykoff CONSPP + extra 0.01·(−0.199592)·RELDEN term + a CO-like DIAGR path.
- **sp11,13-16,19 (CO group GA/CW/BA/PW/NC/OH)**: CR GENGYM DIAGR path — DF=0.24506+1.01291·DPP−
  0.00084659·BA+0.00631·SI (cap 36), DDS=ln(DIAGR·(2·DPP·BARK+DIAGR)).
All DIAGR paths floor DDS at −9.21. ⇒ chunk 3 = KT-Wykoff-main (extend with DGBAL/DGDSQ/slope-aspect) +
RM/CO DIAGR (like CR gemdg) + DGFASP aspen. Data still to extract: 2D DGHAB(8×19)/DGFOR(6×19)/DGDS(4×19)
via MAPHAB(117×19)/MAPLOC(7×19)/MAPDSQ(7×19) + OBSERV(6×19) + DGDSQ(1D) + em/dgfasp.f + dg_resid_sd(DGCONS).
1D DG terms already in data/easternmontana/dg_1d_coeffs.csv.

## Chunk 3 DGCONS derivation (em/dgf.f ENTRY DGCONS, line 615) — the per-stand setup
Per species ISPC (needs IFOR from em_forkod!, ITYPE+JDTYPE from habtyp, ELEV/ASPECT/SLOPE from STDINFO):
- JDTYPE = IDTYPE (habitat idx); if >117 → 30.
- ISPHAB = MAPHAB[JDTYPE,ISPC] for sp 1-3,7-10,18 (Wykoff conifers); else MAPHAB[ITYPE,ISPC].
- ISPFOR = MAPLOC[IFOR,ISPC]; ISPDSQ = MAPDSQ[IFOR,ISPC]; ISPCCF = MAPCCF[ITYPE].
- TMPASP = ASPECT (−0.7854 for sp 4,6,12,17); XSLOPE = SLOPE (÷10 for sp 1-3,7-10,18; =0 for sp18 OS).
- **DGCON(ISPC) = DGHAB[ISPHAB,ISPC] + DGFOR[ISPFOR,ISPC] + DGEL·ELEV + DGEL2·ELEV² +
  (DGSASP·sin(TMPASP) + DGCASP·cos(TMPASP) + DGSLOP)·XSLOPE + DGSLSQ·XSLOPE²**
- **DGDSQ(ISPC) = DGDS[ISPDSQ,ISPC]**; DGCCF(ISPC) species-specific (sp4 special; see line 665+).
- Site class ISIC(1-5) from SITEAR/10 (used downstream).
⇒ NEEDS em_forkod! (IFOR, the JFOR-subscript from em/forkod.f — deferred in chunk 2; MAPLOC/MAPDSQ index it).
Then the main-loop DDS (Wykoff): CONSPP=DGCON+COR+0.01·DGCCF·RELDEN; DDS=CONSPP+DGLD·lnD+CR·(DGCR+CR·DGCRSQ)+
DGBAL·BAL+DGDBAL·BAL/ln(D+1)+DGDSQ·D²(?)+... (verify exact large-tree terms at em/dgf.f:540-585). Then
`DDS=DDS+COR+DGCON` (line 559) for the DIAGR species. All coefficient DATA extracted; forkod + dgfasp remain.

## Chunk 4 (height growth) scope — em/htgf.f (479) + em/pothtg.f (169), measured
Multi-path like DG. MAIN Wykoff conifers (sp1-3,7-10,18 — covers emt01):
- H≤4.5 → small-tree/regen (GO TO 60, chunk 6). Else: BAL=((100−PCT)/100)·BA; CR=ICR/100.
- **CALL POTHTG(I,ISPC,H,SI50,SI100,PHTG)** — potential height growth (em/pothtg.f): per-species site-index
  curves (A=9.72443−0.00091·SI100·CCF−H, B/C forms; TEMSI=SI50−4.5 curves; SI100 base-100 for the conifers,
  SI50 for DF sp3). Returns PHTG.
- **RALPH correction:** PHTG = PHTG·0.706·(1−exp(−10.19·CR))·(1−exp(−0.1·18.158·DG))^0.944 + 0.0265·H.
  (DG = this cycle's diameter growth — so height depends on the DG chunk, already validated.)
- **Modifiers:** RLHTMD = exp(C3MOD·((H/AVH)^C4MOD − 1)); CRMOD=1; HTMOD=min(CRMOD·RLHTMD, 1);
  HTG = max(PHTG·HTMOD, 0.1). C1MOD-C4MOD + HTCON from ENTRY HTCONS (em/htgf.f) — extract.
- Coefficient DATA in em/htgf.f: MAPHAB(30), COFLM(9,3), COFAS(9,3) (LM/aspen), + the mod/HTCON consts.
- SPECIAL forms: sp5 (LL) `CON=HTCON+H2COF·H²−0.1997·lnD+...`; sp4/12/17 (LM/aspen COFLM/COFAS); RM/CO.
Validate via em/htgf.f DEBUG (IN HTGF WRITE line 192/50/901 fire under DEBUG) — same instrument-replay as DG.
emt01 = DF/WL/LP/PP/ES ⇒ only POTHTG main path needed for the first height validation + full-cycle .sum.

## Chunk 6/7 scope — the 2 remaining grow_cycle! hooks for the full .sum (measured)
The full-cycle .sum needs `small_tree_growth!(::EM)` (chunk 6) + `mortality!(::EM)` (chunk 7); grow_cycle!
calls them after height. emt01 has a 0.1" DF seedling so both fire.

**Chunk 7 mortality (em/morts.f, 1266 lines) — HYBRID of KT-density + CS-background:**
- RIP density Hamilton (line 678) = IDENTICAL to KT mortality.jl: `2.76253+0.222310·√D−0.0460508·√BA+
  11.2007·G−0.554421/D+B0+0.246301·RELDBH+6.07129·G/D` (B0=PMSC[sp]). ⇒ reuse KT's mortality! RIP path.
- RI background (line 630) = CS/LS form: `RI = 1/(1+exp(B0+B1·D+B2·D²))` with B0=PMSC/B1=PMD/B2=PMDSQ (DATA
  lines 126/131/136); RI×0.5 halving (line 634). RN = trend-matching rate; RIP=RN unless T≤TEM or RN≤0 → RI.
- Two regimes: SDI-based (default) vs BAMAX-based takeover (header lines 12-14) — like eastern. BAMAX=BAMAXA
  [ITYPE] default. Extract EM PMSC/PMD/PMDSQ + POTEN/BREAK/GMULT/REIN; reuse KT's mortality! structure +
  the shared self-thinning RDPSRT tie-break. Validate via full emt01 .sum (TPA self-thin) vs FVSem_clean.

**Chunk 6 small-tree (em/regent.f):** REGENT small-tree DG/HTG for dbh<threshold (the 0.1" seedling). Reuse
the KT/IE regent pattern (ie_regcons!/small_tree_growth!) + EM coefficients. HTGF's H≤4.5 GO TO 60 path also
lands here. Needed for the seedling's growth in the projection.

## Chunk 7 mortality FULL flow (em/morts.f) — measured, ready to implement
Stand-level (lines 360-486, = KT mortality.jl structure): D10=DQ10 (or DR10); DELTBA→BA10=BA+((BAMAX−BA)/
BAMAX)·DELTBA; TB=BA10/(0.005454154·D10²); TTB=(T−TB)/T (cap .9999); **RZ=1−(1−TTB)^0.1**; CONST=SDIMAX/
0.02483133; TMD10=CONST·D10^−1.605 (cap 35000); T85D10=TMD10·PMSDIU; T55D10=TMD10·PMSDIL. AVED=Σ(D·P)/ΣP.
Per-tree, TWO paths by species:
- **ORIGINAL EM species (1-3,7-10,18 = emt01's WL/DF/LP/ES/PP):** RI=0.5/(1+exp(PMSC+PMD·D+PMDSQ·D²)) [halved];
  RIP=RN (the SDI trend rate from RZ/T85D10 logic, lines 391-600); TEM=CONST·D10^−1.605·PMSDIL;
  IF(T≤TEM OR RN≤0) RIP=RI (background when SDI not yet limiting); WKI=P·(1−(1−RIP)^FINT)·X (X=XMORT in
  [D1,D2] only if RIP==RI). ⇒ a background+SDI-trend model, NOT the KT density Hamilton.
- **ADDED species (4,5,6,11-17,19):** the KT density Hamilton (RIP=2.76253+…, line 678) — reuse KT path.
Implement mortality!(::EM): reuse KT's stand-level (rz/tb/ttb/ba10/CONST/TMD10) + add the RI/RN original-species
branch (PMSC/PMD/PMDSQ in mort_bkgd_coeffs.csv, extracted). Then shared RDPSRT self-thin + full emt01 .sum vs
FVSem_clean (EM barely self-thins — 536→455 by 2060 — so the RN/RI balance is the key validation).
FIND RN's exact assignment (em/morts.f 391-600, the T85D10/T55D10 interpolation) before coding.

## Chunk 6 (regent / small_tree_growth!) scope — the LAST grow_cycle! hook (measured)
em/regent.f (1491) REUSES the KT/IE REGENT pattern: REGYR=5.0 (5-yr subcycles), HSIGMA=0.59,
BACON=0.005454154 (em/regent.f:147) — structurally identical to KT's small_tree_growth! (subcycle count/
lengths from cycle length, per-subcycle stand density from the large trees, small-tree HTG increment model
then DG derived). ⇒ reuse KT's small_tree_growth! structure + EM regent coefficients (extract from
em/regent.f DATA) + em_regcons! (RHCON/RHGL calibration). Small-tree threshold dbh<3" (KT) — emt01 has ONE
0.1" DF seedling, so the stand .sum impact is tiny, but faithful is the goal. Once implemented, grow_cycle!
runs end-to-end → the full-cycle emt01 .sum differential vs FVSem_clean VALIDATES DG+height+CCF+mortality+
regent together (TPA/BA/SDI/CCF/TopHt/QMD across all cycles). This is the FINAL growth-core hook.

## Chunk plan (mirror KT)
1. Species: `data/easternmontana/species_coefficients.csv` (19 species × ~40 cols from em/*.f DATA:
   bark bratio.f, site sitset.f, small-tree/regent, htcalc/htdbh, morts, sdimax, volume specs) +
   `species_translation.csv` (SPCTRN western crosswalk, EM target column) + `species.jl` blkdat init.
2. Site/habitat (sitset + em_habtyp, 8 habitat-type groups). 3. DG (em/dgf.f — reuse KT
   `diameter_growth!`, swap coefs + DGHAB(8)/DGFOR(6) dims). 4. Height (em/htgf.f). 5. Crown/CCF
   (em/crown.f/ccfcal.f). 6. Regent (em/regent.f). 7. Mortality (em/morts.f + varmrt.f — verify vs
   shared driver + KT Hamilton; the low-thinning behavior). 8. Volume (shared NVEL/cr_fw2_vol).
   9. Full-cycle differential vs FVSem_clean on emt01 + native stands.

## Chunk-2/inventory validation: EM 1990 read+expansion+site+height BIT-EXACT
Validated EM's cycle-0 inventory (each_stand → notre! → setup_growth!, diagnostic identity+bark+site CSV,
compare via the .sum normalization stand_X/gross_space): TPA=536 BA=77 SDI=184 TopHt=63 QMD=5.1 — ALL
BIT-EXACT vs live FVSem emt01. Only CCF=0 (needs crown-width coefs B0ACCF/B1ACCF/B0BCCF/B1BCCF/B0ASTD/B1BSTD
from blkdat.f). ⇒ the READ + notre! expansion + site_setup! + measured-height path is VALIDATED end-to-end
at inventory. Design params all match live (BAF=40/fixed=300/break=5/11 plots/1 nonstock/wt=11).

META (2 phantom bugs avoided this session): (1) "BA=8" — isolated setup_growth! skips notre! (the tpa
expansion). (2) "TPA 10% high" — the .sum divides stand_X by gross_space (simulate.jl:456); direct stand_ba
returns the gross_space-inflated internal value. ALWAYS replicate the .sum normalization (notre! + /gross_space)
when validating inventory in isolation — else you chase measurement artifacts (doctrine #2/#3 meta-lesson).

## ★★ CHUNK 3 DG VALIDATED vs live FVSem (em/dgf.f DEBUG instrument-replay)
Enabled em/dgf.f DEBUG (WRITE IN DGF line 598) via the DEBUG keyword on a 1-cycle emt01 keyfile; captured
per-tree DDS (calibration pass). jl em_dgcons!+dgf! (Wykoff-main) DDS vs live, 4 WL trees:
  I=4 DBH7.26 BAL36.61 CR0.25: jl 1.8961 / live 1.8958 (Δ0.0003)
  I=5 DBH7.25 BAL33.32 CR0.25: jl 1.8980 / live 1.8980 (EXACT)
  I=7 DBH6.99 BAL27.71 CR0.35: jl 1.9488 / live 1.9490 (Δ0.0002)
  I=10 DBH6.92 BAL30.48 CR0.45: jl 2.0044 / live 2.0039 (Δ0.0005)
All within F7.4 debug-print rounding (live values are themselves 4-decimal). ⇒ EM DG Wykoff-main DDS +
DGCON(em_dgcons!) + all coefficients VALIDATED bit-exact-or-cornered vs live FVSem. The largest EM chunk (DG)
is correct. RELDEN=90.837/BA=67.097 (calibration-pass backdated stand) also matched (same debug line).
NOTE: live FVSem SEGFAULTS at run-end under DEBUG mode (DGF debug printed fine first) — a DEBUG-only live
crash to trace later (doctrine: fix live crashes; low-pri, non-DEBUG runs are clean).
NEXT: height (ch4 htgf) to enable the full-cycle .sum projection differential; then NI/DIAGR/aspen DG paths
(non-emt01 species) + regent/mortality/volume.

## Chunk 7 RN/TN10 — the SDI self-thinning (measured, more complex than KT)
RN = 1−(1−(T−TN10)/T)^(1/FINT) where TN10 = target tree count after SDI mortality (em/morts.f 486-590):
- T > T85D0 (85% SDI at DIA0): TN10=T85D10 (kill to 85% line).
- T55D0 < T ≤ T85D0: ITERATIVE linear-fn fit between the 55%/85% SDI lines (label 220 loop, ≤100 iters,
  Newton-ish TREEIT+=0.5·DIFF) → TN10=exp(CEPMRT+SLPMRT·ln(D10)), capped at T85D10.
- T ≤ T55D10: TN10=T (HOLD — no SDI mortality). ⇒ ★KEY: a below-55%-SDI stand gets RN=0 → RIP=RI (BACKGROUND
  ONLY). This is WHY emt01 barely self-thins (536→455) — it's below the SDI limit, so PMSC background dominates.
T85D10=TMD10·PMSDIU, T55D10=TMD10·PMSDIL, TMD10=CONST·D10^−1.605, CONST=SDIMAX/0.02483133. DIA0/D10=prev/cur
DQ10. This iterative SDI boundary is NOT in KT (KT used rz=1−(1−ttb)^0.1) — EM mortality is a distinct, larger
port. Implement: stand-level TN10 fit + RN + the RI-background original-species branch (PMSC extracted) +
KT-density-Hamilton added-species branch. For emt01's early cycles RN≈0 so the RI/PMSC path is the first
validation target (full .sum TPA). MSB (mature-stand-boundary SLPMSB/QMDMSB) is an optional refinement.

## Chunk 6 regent small-tree HTG detail (em/regent.f 419-500) — measured
Small-tree height increment dispatches by species-SOURCE variant flag (CRVAR/UTVAR/NIVAR/TTVAR — from the
species-expansion lineage). emt01's conifers (DF/WL/LP/PP/ES) use **NIVAR** (North Idaho, the main EM path):
  RELH=(H1−4.5)/(AH−4.5) [0,1]; DADJ=DELMAX·RELH²−2·DELMAX·RELH+0.65 (bias); TPCCF=clamp(PCCF·PPCCF,25,300);
  **HTGRL = CON + BH·ln(H1) + BCCF·RDJ + BBAL·BAL** (Wykoff small-tree form); then HTG via subcycles (REGYR=5,
  reuse KT). Needs per-species CON/BH/BCCF/BBAL/DELMAX + SLO/SHI (site low/high). DG derived from HTG.
  (CRVAR path: POTHTG=SITEAR/(15−4·RELSI); UT sp6 different — not emt01.)
⇒ chunk 6 = reuse KT small_tree_growth! subcycle scaffold + the NIVAR HTGRL small-tree HTG + EM coefs +
em_regcons! calibration. THE LAST growth-core hook — then full-cycle emt01 .sum differential validates
DG+height+CCF+mortality+regent end-to-end. Growth core (DG validated / height / mortality matching live) done.

## Chunk 6 regent NIVAR coefficients (em/regent.f) — measured, COMPLETE
NIVAR HTGRL = CON + BH·ln(H1) + BCCF·RDJ + BBAL·BAL, then HTG via REGYR=5 subcycles + DADJ bias:
- BH=0.3740, BCCF=−0.00391, BBAL=−0.22957 (FIXED, em/regent.f:391-393).
- CON = RHCON(sp) + HCOR(sp)  (calib mode 40: RHCON·EXP(HCOR); line 383/385). RHCON per-species base +
  HCOR the small-tree height calibration (em_regcons!, like KT/CR/IE htg_cor_small).
- DELMAX = min((AH/36)·(0.01232·R − 1.75), 0); RELH=(H1−4.5)/(AH−4.5) clamp[0,1]; DADJ=DELMAX·RELH²−
  2·DELMAX·RELH+0.65. RDJ = relative density (TPCCF-based point CCF, clamp[25,300]).
- SLO/SHI per-species DATA (em/regent.f:155/161) for the CRVAR POTHTG (not NIVAR); NIVAR uses HTGRL directly.
⇒ EM ENTIRE GROWTH MODEL NOW FULLY MEASURED. Chunk 6 impl = KT small_tree_growth! subcycle scaffold + this
NIVAR HTGRL + RHCON extraction + em_regcons! (RHCON/HCOR calib) + DG-from-HTG. TT/CR/UT/aspen small-tree
branches deferred (not emt01). Then full-cycle emt01 .sum validates the whole growth core end-to-end.
GROWTH CORE STATUS: DG VALIDATED vs live · height DONE · mortality DONE (matches live) · inventory+CCF+site
bit-exact · regent = last hook, fully measured.

## Chunk 6 regent — RHCON + implementation plan (measured, ready)
em_regcons! (em/regent.f:1481): RHCON(sp) = REGCH(sp) + 1.0667 + RHHAB(IRHHAB) [per-species base REGCH +
habitat-indexed RHHAB]. CON = RHCON + HCOR (HCOR = small-tree ht calibration, calib mode-40 RHCON·EXP(HCOR)).
IMPLEMENTATION = mirror KT small_tree_growth! (regent.jl ~250 lines: subcycle count/lengths from FINT,
per-subcycle stand density banext/rdnext from large trees, per-tree HTGRL over subcycles, then HTGR1+ZZRAN
(dgsd≥1, bachlo)+XWT blend with the large-tree HTG + D<3 diameter dub) — but swap KT's HTGRL for EM NIVAR:
HTGRL = CON + 0.3740·ln(H1) + (−0.00391)·RDJ + (−0.22957)·BAL + DADJ-bias-path; DG dub via EM DIAM/HCON/DCON.
Extract EM regent DATA: REGCH/RHHAB (RHCON), XMIN/XMAX (small-tree dbh bounds), DIAM/HCON/DCON (DG dub),
SLO/SHI, HSIGMA=0.59. Reuse KT's ZZRAN/XWT/size-cap logic verbatim. emt01 = 1 DF seedling (NIVAR).
⇒ THIS is the last growth-core implementation; a large-but-scaffolded port (KT reuse). After it: full-cycle
emt01 .sum differential = the end-to-end validation of DG+height+CCF+mortality+regent.

## Chunk 6 regent — REGCH + DG-dub (final measurement)
REGCH = RHGL(IGL) + (RSAB0 + RSAB1·cos(ASPECT) + RSAB2·sin(ASPECT))·SLOPE; RSAB0=−0.10987, RSAB1=0.22157,
RSAB2=−0.12432 (em/regent.f:1440,1446). RHGL(IGL) = geo-location base (IGL from forkod, emt01 IGL=1).
DG-dub (em/regent.f:575): D1 = DIAM(sp)+DADJ if H≤4.5 else AX·(H−4.5)^BX+DADJ — a POWER form (AX/BX per
species), NOT KT's linear HCON·H+DCON ⇒ EM regent DG-dub is EM-SPECIFIC. Extracted DATA: RHHAB(5)=−0.2146,
−0.0941,−0.4916,−0.3582,0; MAPHAB(30)=12*4,3,2,4*5,1,7*4,1,3,3,4; XMAX/XMIN/DIAM(19). Still need: RHGL, AX/BX.
⇒ regent = a LARGE EM-specific port (KT subcycle SCAFFOLD reusable, but HTGRL(NIVAR) + DG-dub(AX/BX power) +
RHCON(REGCH+RHHAB) are EM-specific). EM GROWTH MODEL 100% MEASURED. This is the last growth-core implementation.

## FULL-CYCLE STATUS — all growth chunks implemented; crown_ratio_update! is the LAST hook
grow_cycle! hooks for EM: DG(✓validated) → height(✓) → small_tree_growth!/regent(✓ d9f48ad) → mortality(✓
matches live) → establish!(shared, NOAUTOES=inert) → crown_ratio_update!(✗ — the ONLY missing hook).
crown_ratio_update!(::EM) = em/crown.f Weibull crown-ratio (WEIBA/WEIBB0/WEIBB1/WEIBC0/WEIBC1/C0/C1 +
PARM[IE/NI sp9 form] + CRHAB[habitat] + MAPHAB(30)=13*2,3,4*4,5,6,6,7,6,1,8,1,9,2*10,6 + CRSD=6.35). Reuse
KT/IE crown_ratio_update! structure (Weibull ACR/BCR + habitat CRHAB + ZZRAN) + EM coefficients. Once in,
grow_cycle! runs END-TO-END → full-cycle emt01 .sum differential vs FVSem_clean validates DG+height+CCF+
mortality+regent+crown together (TPA/BA/SDI/CCF/TopHt/QMD × all cycles). THE last growth-core hook.

## Full-cycle divergence LEAD (2000 BA/QMD ~9% high) — under investigation
grow_cycle! runs end-to-end; 1990 BIT-EXACT, 2000 TPA BIT-EXACT (526), but 2000 BA 114 vs live 96 / QMD 6.3
vs 5.8 (~9% high, TopHt 69 vs 68 close). BA divergence is ENTIRELY QMD (DBH growth) — the DG over-applies
~9%/cycle. NOT fixed by the real dg_resid_sd=SIGMAR + EM_PSIGSQ (added — faithful, but 2000 unchanged). The DG
calibration-pass DDS is bit-exact vs live (DEBUG replay 0.0005), so the issue is the GROWTH-pass DG (current-
stand RELDEN/CCF/BAL feeding dgf!, or the DDS→DG conversion/period). NEXT: instrument the live GROWTH-pass DGF
(the DEBUG-keyword run SEGFAULTS at run-end before the growth pass prints — need a non-DEBUG WRITE patch +
relink, OR compare jl growth-pass RELDEN/CCF vs the inventory-validated 97). Candidate causes: (a) growth-pass
RELDEN differs (jl RELDEN during growth vs live); (b) DDS→DG conversion; (c) DGSCOR self-cal firing wrongly.
Deferred: mortality iterative TN10 fit (55-85% SDI, fires ~2020). EM_PSIGSQ wired into the shared psigsq dispatch.

## Residual DG ~7% (post-RELDEN-fix) — DG-precision lead, needs live growth-pass DGF
After the RELDEN fix: emc2 2000 BA 103 vs live 96 (7% high), 2010 BA 133 vs 114 (compounds). jl WL growth-pass:
RELDEN=106.653 (raw CCF; ÷gross_space=97), DDS=1.9146, DG=0.5436, dg_cor[WL]=0. NOT the COR: emt01.tre has
measured DG only on a few records (rec6=5.6") ⇒ jl calibration fires correctly (PP COR=-0.0226; WL/DF/LP/ES
COR=0, matching no-measured-DG). So the residual is DDS-precision, candidates: (a) RELDEN value — jl uses RAW
stand_ccf 106.653; if live's DGF wants per-acre 97 it'd be HIGHER DDS (worse) so probably not; (b) CROWN RATIO
feeding DDS's CR·(DGCR+CR·DGCRSQ) — my new crown_ratio_update! (minimal d<3 dub) may drift CR; (c) the growth-
pass RELDEN timing (pre vs post prev-cycle); (d) DDS→DG bark conversion. TO PIN DOWN: capture the live GROWTH-
pass DGF DDS — the DEBUG-keyword run SEGFAULTS at run-end before it prints, so patch em/dgf.f WK2(I)=DDS site
with an unconditional WRITE(unit) + relink (relink_em.sh <name> patched.o) + compare per-tree DDS at cycle-1
growth. This is the last full-cycle precision gap (+ the deferred mortality iterative TN10 fit at 2020).


## Full 10-cycle runs (iterative TN10 fit done) — compounding residual LEAD
emc10 (control 10cyc) jl vs live: 2000 526/98/5.8 vs 526/96/5.8 (TPA/QMD exact, BA±2); 2020 507/145/7.3 vs
507/132/6.9 (★TPA EXACT — mortality iterative fit CORRECT — but QMD+0.4/BA+13 ~10%); 2040 474/197/8.7 vs
488/168/7.9; 2080 370/279/11.8 vs 423/219/9.7. ⇒ mortality iterative TN10 fit WORKS (2020 TPA exact). Residual
= COMPOUNDING QMD/BA (2000 exact → 2020 +6% → 2040 +10%), driving later TPA divergence (denser stand self-thins
more). This is LARGER than CR/KT/IE DGSCOR tail (~1-3%) ⇒ likely a small FIXABLE feedback, candidates: (a)
CROWN RATIO — my new crown_ratio_update! (minimal d<3 dub + DCR) feeds the DG CR·(DGCR+CR·DGCRSQ) + height RALPH
CR term; if CR drifts over cycles, DG/HT compound; (b) DGSCOR self-cal on PP; (c) a per-species DG coef precision.
growth-pass DDS was BIT-EXACT for WL only — CHECK DF/LP/ES/PP growth-pass DDS via DGFTRC replay across cycles.
EM GROWTH CORE VALIDATED END-TO-END (full 10-cycle runs, early bit-exact-or-cornered, 2020 TPA exact); the
compounding residual is the last precision lead. NEXT: crown-ratio drift check + multi-species DGFTRC; then
volume + non-emt01 DIAGR/aspen/NI species-paths.


## ★★★★ EM GROWTH CORE FULLY VALIDATED END-TO-END — full 10-cycle BIT-EXACT-OR-CORNERED
7th real bug (point_ccf/PCCF): standstats.jl point_ccf dispatch had KT/IE but not EM ⇒ EM fell to crown-width
else (→~0) ⇒ PCCF 0.66 not ~100 ⇒ dgf! DGPCC1/2·PCCF term dropped ⇒ DG high for DGPCC species (LP/ES/PP; WL
immune DGPCC=0). Localized via multi-species DGFTRC replay (WL DDS exact 1.9146=live; ES off 1.5488 vs 1.38423
= the missing PCCF term). Fix: EM uses em_tree_ccf for PCCF. ★ emc10 full 10-cycle: 2000 526/96/5.8 EXACT;
2020 507/133/6.9 (TPA/QMD exact, BA±1); 2040 497/174/8.0 vs 488/168/7.9 (±0.1 QMD); 2080 420/229/10.0 vs
423/219/9.7 (cornered self-thin/DGSCOR tail like CR/KT/IE). EM DONE at growth-core = CR/KT/IE level.
★ 7 REAL BUGS this session: habitat_code-gate, strip_key_ext(shared), IE-ESXCSH, EM-habitat_code, RELDEN-
dispatch, DBH-update-bark, point_ccf-PCCF-dispatch. Instrument-replay recipe (DGFTRC): patch dgf.f WK2=DDS +
WRITE(77,...); compile to /workspace/.emwork/dgf.o (NOT into BD — stray obj → multiple-def dgf_); gfortran
manual link excluding $BD/dgf.o$ + shim. NEXT (downstream leaves, like CR/KT/IE): volume(NVEL) + non-emt01
DIAGR(RM/CO)/aspen(DGFASP)/NI(sp4/5) DG paths + full multi-scenario stands.


## ★ EM VOLUME VALIDATED (= IE/KT standard) — chunk 8
EM = Region-1 (Northern) Flewelling FW2, SAME kernel as KT (cr_fw2_vol). Added VOLEQ I00FW2W<FIA> (em/sitset.f
VAR=EM/IREGN=1) + compute_volumes! EM dispatch → compute_volumes_kt! (variant-agnostic: Region-1, sp7=LP
lodgepole). emc2 1990: Merch=1058 BIT-EXACT, BdFt=5388 BIT-EXACT, TotCuFt 1620 vs 1637 (~1% = INGY geocode-I
taper upper-stem residual, SAME as IE — cornered). ⇒ EM NOW = KT/IE PARITY: growth core validated end-to-end
(full 10-cycle bit-exact-or-cornered) + volume validated (merch/board exact, tot-cubic INGY-cornered).
Remaining (situational, like IE): hardwood DVE volume species + forest-keyed VOLEQDEF; non-emt01 DIAGR(RM/CO)/
aspen(DGFASP)/NI(sp4/5) DG paths; FFE fuel; establishment. EM is a FOURTH western variant to growth+volume
validated this session (after CR done; KT+IE growth-core; EM now growth+volume). 7 real bugs fixed.


## Non-conifer species need FULL paths (DG+height+regent), not just DG — measured
Validated via all-LM synthetic stand (emt01.tre species→LM): live 2000 485/101/6.2 (distinct from conifer
526/96/5.8, confirms LM/NI path fires). jl errors at height_growth! sp4 — so LM/LL/RM/AS/CO need their HEIGHT
path (em/htgf.f COFLM(9,3)/COFAS(9,3) for LM/aspen; the sp5 LL HTCON form) + REGENT small-tree path, in
ADDITION to the DG path (NI done sp4/5; DIAGR RM/CO + aspen DGFASP still todo). ⇒ each non-conifer species-
group is a full DG+height+regent triad. emt01 (conifers) exercises only the main paths (all validated). EM
DONE on conifers (growth end-to-end + volume, = KT/IE); non-conifer species-groups (LM/LL/RM/AS/CO/hardwoods)
are the situational remaining coverage — needs a per-group DG+height+regent port + a test stand. Same tail as IE.

## 2026-08-05 — emt01 CONTROL (NOAUTOES) growth-only: measured ~+7% BA over-growth by 2090 → CROWN-RATIO lead
Re-validated the first emt01 stand (S248112 UNTHINNED CONTROL, NOAUTOES growth-only, 10 cycles) vs live FVSem_clean:
| year | TPA l/j | BA l/j | QMD l/j |
| 1990 | 536/536 | 77/77 | 5.1/5.1 |  (bit-exact)
| 2000 | 526/526 | 96/96 | 5.8/5.8 |  (bit-exact .sum; per-tree dbh Δ~0.002 already present)
| 2010 | 517/517 | 114/115 | 6.4/6.4 |
| 2050 | 474/471 | 183/194 | 8.4/8.7 |
| 2090 | 411/394 | 230/246 | 10.1/10.7 |  (jl BA +7%, QMD +6%, TPA −4%)
jl OVER-grows DBH; the excess compounds ~0.5%/cycle from cyc1. This is the growth-only CONTROL — DISTINCT from
task #137's establishment-stand self-thin over-kill.
MEASURED (instrument-replay, live em/dgf.f WRITE(16) of D/CR/PCT/BAL/RELDEN/BA/DDS at ICYC=2, matched to jl by
sp+dbh):
- DG formula + bark BIT-EXACT at cyc0 (2000 .sum bit-exact; sp2 bark=0.934 jl==live; NI/DIAGR paths not exercised).
- At cyc1 the sp2 per-tree DDS is jl-LOWER than live (jl 1.9279 vs live 1.9553), which RECONCILES with a measured
  CROWN-RATIO divergence: live sp2 crowns are tight 0.27–0.28, jl's are spread 0.25 / 0.25 / 0.32. Crown feeds the
  DDS `cr·(DGCR+cr·DGCRSQ)` term, so a lower jl crown → lower jl DDS (the formula is consistent given the inputs).
- The ~3pp crown difference is FAR too large to be explained by the ~0.002" cyc0 dbh divergence ⇒ the EM crown-ratio
  UPDATE model (crnew/crown.f) is the prime suspect, feeding DG and compounding.
NOT YET ISOLATED (cause-vs-effect open — per the CI each_stand lesson, NOT claiming a root before proving it): per-
tree matching ACROSS the divergence is unreliable (doctrine #3). To confirm the crown model is the CAUSE (not an
effect of dbh drift), the next step is a cyc0-ANCHORED check: pick one input tree (bit-identical at 1990), dump its
crown ratio at the START of cyc1 from BOTH jl and live, and see if the crown-update itself diverges for identical
input. If yes → real crown-model bug (likely cross-variant, shared crnew). If the crowns match for identical input
→ the divergence is a dbh-drift effect and the whole thing is the accepted coupled precision tail (cornered).
Puzzle to resolve alongside: jl sp2 per-tree DDS is LOWER yet aggregate BA is HIGHER — so the aggregate over-growth
is carried by OTHER species / the small-tree cohort / size redistribution, not sp2. Verify species BA contributions.

### CORRECTION (same day, doctrine #2) — retract the crown-specific claim; robust lead is DG-point DENSITY
The crown-ratio pairing above is UNDER-SUPPORTED: the CR-carrying live dump (I3≤6) covered sp2 trees d≈8.31–8.56,
but the d≈8.909 tree I paired with jl's d≈8.907 came from the earlier dump with NO CR captured — so "jl crown 0.32
vs live 0.28" mixed two different trees. Per-tree matching across the divergence is unreliable (doctrine #3: trees
don't correspond 1:1 post-cyc0; CR not fully captured). RETRACTED.
ROBUST (stand-level, not per-tree) fact that stands: at the cyc1 DG point, jl's density inputs are ~0.3% HIGHER than
live's — jl ba=105.53 / relden=123.98 vs live ba=105.2159 / relden=123.2158 — EVEN THOUGH the 2000 .sum stand BA is
bit-exact (96/96). Note the DG-point "BA" (~105) ≠ the reported .sum BA (~96), so it is a different aggregation
(per-point / expansion / small-tree-inclusive). Higher jl density → more competition → lower jl DDS (matches the
observed jl-lower sp2 DDS), yet aggregate BA ends up HIGHER — so the compounding is a density-aggregation +
size-redistribution effect, not a single DG term. PROPER next step (STAND-level, avoids the per-tree trap): compare
how jl vs live aggregate the DG-point BA/RELDEN/PCCF at cyc1 (point_ccf, expansion, small-tree inclusion) — find the
0.3% source. If it's a faithful-but-different rounding in the density aggregation, this is the accepted coupled
precision tail (cornered); if a real aggregation bug, it's cross-variant (shared density code). DG formula + bark
remain bit-exact at cyc0 (unchanged).

## 2026-08-05 — cyc0-DG test (NOTRIPLE, emt01) — the bug-vs-cornered settler, at resolution limit
Ran the definitive test for the compounding over-growth tail: emt01 first stand, NOTRIPLE (trees stay 1:1 with
input ⇒ per-tree matching is doctrine-valid), NUMCYCLE 1 (1990→2000), TREELIST. Compared jl's exact post-growth
2000 DBH (instrumented at simulate.jl:533 `dbh += diam_growth/bark`) to live FVSem's TREELIST 2000 CURR DIAM.
FINDINGS:
- **cyc0 aggregate growth is BIT-EXACT**: NOTRIPLE .sum 2000 = 526/96/5.8 jl == live.
- **Per-tree 2000 DBH matches within the .trl's 1-decimal precision** for most trees (e.g. WL start7.9→8.40 jl ==
  live 8.4; start8.0→8.465 jl == live 8.5), EXCEPT 1-2 boundary cases: WL start8.2 → jl **8.7731** vs live "8.7";
  WL start8.4 → jl **8.9475** vs live "8.9". These sit on the rounding boundary ⇒ AMBIGUOUS between a real
  ≤0.5-0.8% over-growth and a print rounding/truncation artifact (the .trl gives DBH to 1 decimal only).
- ★ TRAP AVOIDED (doctrine #2): the END-CYCLE-0 .trl "DIAM INCR" column (0.60/0.70/1.20…) is NOT the predicted
  cyc0 DG — it is the OBSERVED/calibration increment. Comparing jl's predicted DG (0.47/0.43/0.54) to it fabricated
  a bogus "2× difference". The valid comparison is the grown 2000 DBH.
VERDICT: the compounding tail's cyc0 seed is **at or below the .trl's 1-decimal resolution (≤~0.5-0.8% on a couple
of trees, aggregate bit-exact)** — consistent with the accepted ULP/DGSCOR-accumulation cornered class, but NOT a
definitive "zero cyc0 bias" proof. DEFINITIVE SETTLE (documented next step): instrument live FVSem's gradd.f UPDATE
(or dgf→gradd DG) to dump per-tree 2000 DBH to full Float precision and diff against jl's exact d2000 — if they
match to ~1e-4 the tail is confirmed cornered (downstream ULP accumulation); if jl is systematically higher by
~0.5%, it is a real, cluster-wide, fixable crown/PCT-fed DG bias (the shared crown-ratio/BA-percentile computation
at cyc0 is the prime suspect, per the retracted-but-directionally-suggestive earlier crown/PCT lead).

### DEFINITIVE SETTLE (full-precision, same day) — the tail is CORNERED (tie-break precision), not a fixable bias
Instrumented live FVSem dgf.f to dump D at ICYC=2 (= the full-precision grown 2000 DBH; NOTRIPLE, NUMCYCLE 2) and
diffed against jl's exact d2000. Matched by species+start:
| tree | jl d2000 | live d2000 | Δ |
| WL 7.9 | 8.4038 | 8.3618 | +0.5% |
| WL 8.0 | 8.4650 | 8.5144 | −0.6% |
| WL 8.2 | 8.7731 | 8.7085 | +0.7% |
| WL 8.4 | 8.9475 | 8.8833 | +0.7% |
| DF 10.0 | 10.6688 | 10.7544 | −0.8% |
| DF 10.4 | 11.1038 | 11.1675 | −0.6% |
| DF 1.2 (small) | 1.5379 | 1.3620 | +13% |
| DF 1.9 (small) | 2.1178 | 2.2227 | −4.7% |
⇒ cyc0 per-tree DG differences are **REAL (~0.5-0.8% large-tree, ±5-13% small-tree), NOT ULP** — BUT **mixed-sign
and mostly-cancelling** (aggregate BA bit-exact 96/96). That is precisely the signature of **DDS-input tie-break
precision** — the BA-percentile (PCT→BAL) and crown-ratio fed into the Wykoff DDS differ per-tree by the unstable-
sort/tie ordering, the SAME accepted RDPSRT/AVHT40-tie-break irreducible primitive the largest-FIA-divergence
campaign already verified (see memory fvsjl-largest-div-campaign). The imperfect cancellation leaves a small net
that compounds into the ~7% BA-by-2090 tail. **VERDICT: the EM/IE/cluster growth-only over-growth tail is CORNERED
(accepted tie-break/coupled-precision class), NOT a systematic fixable DG bias.** This closes the bug-vs-cornered
question the earlier entries left open (and supersedes the "needs the definitive measurement" caveat). The small-
tree (D<3) larger deltas are the separate SMHTGF/regent stochastic-ZZRAN class (also cornered). No fix warranted.

### EM real-FIA sweep (2026-08-06, 200-stand) — DOMINANT gap is AUTOES (#143), NOT mortality; "EM validated" was cyc0-only
Ran the BM-style multi-metric real-FIA sweep on EM (200 stratified stands → sub-DB `em_sub.db`, live FVSem_clean,
NUMCYCLE 3, TPA/BA/QMD at final year). jl: 0 crashes. Two measured findings:

1. ★★ **AUTOES establishment gap is the DOMINANT EM real-FIA divergence** — ~180/200 stands: jl produces **0 TPA**
   where live has 100-260. CONFIRMED not mortality: stand 5332701010661 (and the repeated 118-TPA `…010661` cohort)
   START TREELESS (1979: 0 TPA both) and live **auto-establishes to 118** by 1989 via AUTOES natural regen; jl stays
   0 (jl only does explicit PLANT/NATURAL, not AUTOES — the estb/ ~6280-line subsystem, task #143). ⇒ This CORRECTS
   the goal-file "EM real-FIA validated, cyc0 bit-exact" claim: cyc0 IS bit-exact (both treeless), but cycle-1+
   diverges catastrophically because live establishes and jl doesn't. AUTOES (#143) is NOT IE-only — it is the
   dominant real-FIA fidelity gap for EM (and by extension the grassland/regen-heavy western variants). mean|Δ%| on
   the full 200 = 91% purely from these 0-vs-100+ AUTOES stands.

2. **EM mortality UNDER-thins genuinely-treed stands — same signature as BM #140.** Of the ~18 stands with real tree
   data (jl>0), the 7 divergent ones lean under-thin (5 JL-HIGH +4.2%/+5.4%/+6.4%/+6.7%/+7.3% : 2 JL-LOW −1.3%/
   −4.4%). Exemplar 83402334020004: cyc0 TPA 7721 BIT-EXACT (same tree data), then live→7026 vs jl→**7500** (+6.7%),
   BA live 20 vs jl 67, QMD live 0.7 vs jl 1.3 — jl retains the self-thinning cohort, exactly BM's pre-fix pattern.
   EM uses the SAME custom naive single-pass uniform-RN mortality! that BM had (missing BMTMRT distribution + IPASS
   QMD-convergence).

★ COMPLICATION (why the BM shared-driver routing does NOT drop into EM): EM's mortality! is a HYBRID —
easternmontana/mortality.jl:126-165 branches per tree: ORIGINAL species use the SDI self-thin RN (like BM), but
ADDED species (4-6,11-17,19) use a completely separate KT/IE Hamilton potential-mortality regression (GMULT/REIN/RZ/
BAMAX-limited, em/morts.f:661-723) that the shared mortality!(::AbstractVariant) does NOT implement. So fixing EM's
under-thin (#137) is NOT a clean "route through the shared driver like BM/CR" — it needs the BMTMRT distribution +
IPASS convergence added to EM's ORIGINAL-species SDI path while PRESERVING the added-species Hamilton branch (a
hybrid, more involved than BM). Lower priority than #143 for EM real-FIA (mortality affects only the ~18 treed
stands; AUTOES affects ~180). Harness: scratchpad em_sample.txt + em_sub.db + em_live.txt (reusable).

### AUTOES architecture — SCOPE DE-RISKED (2026-08-06): shared estb/ model + variant species-crosswalk, NOT a 6280-line per-variant port
Traced the FVS AUTOES data model to scope the EM extension (#143). KEY FINDING — the establishment model is almost
entirely SHARED estb/ code+DATA, not variant-specific:
- **ESTOCK stocking regression** (estb/estock.f): the SHAB(16)/SSER(5)/SPRE(3,4)/FORDF(20)/FORGF(20) coefficients +
  the PN logistic are hardcoded DATA in the SHARED file — jl's `_IE_ESTOCK_SHAB/SSER/...` consts are actually these
  SHARED coefficients (the `_IE_` prefix is misleading; they apply to all Northern-Rockies variants).
- **Habitat-series maps** (estb/estab.f): MYTYPE(30)/MYHABG(16)/MYHTS(16) — the ITYPE→series/height-series maps —
  are SHARED DATA (`MYTYPE/9*1,2*5,2*2,3*3,4,…/, MYHABG/4*1,4*2,3,4,6*5/`), same across variants.
- **OCURHT occupancy + the species tally** (estb/esnutr.f/estab.f): built on a FIXED Northern-Rockies estb/ species
  set that each variant crosswalks to its own species indices (IE 23 sp, EM 19 sp).
⇒ VARIANT-SPECIFIC pieces are only: (a) habitat_code→ITYPE resolution (already read per-variant), (b) the estb-
species ↔ variant-species crosswalk, (c) the species count. The DATA/model is shared and ALREADY PORTED as jl's
ie_autoes_* (establishment.jl). ⇒ EM AUTOES = GENERALIZE ie_autoes_run/ie_autoes_tally/ie_autoes_establish! to be
variant-agnostic (parameterize the species count + crosswalk, drop the `s.variant isa InlandEmpire` gate at
simulate.jl:561) + wire EM's habitat + species map + validate vs FVSem_clean on the AUTOES-0 stands (em_sub.db). This
is FAR more tractable than "port estb/ ~6280 lines". CAVEAT: jl's IE AUTOES is "v1" (~22% residual, target-
contamination unresolved) so EM would inherit that maturity; the generalization is the enabling step. #143.

### EM AUTOES ground truth captured (2026-08-06, instrument-replay on stand 5332701010661)
Instrumented live estb/estab.f (WRITE at IHAB/ISER derivation + OCURHT dump), relinked FVSem_trc, ran the EM
AUTOES-0 repro stand (starts treeless 1979, live establishes 118 TPA by 1989). GROUND TRUTH for the eventual
EM AUTOES port validation:
```
EMEST plot=1 ITYPE=4 IHAB=3 ISER=1 IHTSER=2 IFO=11 NOFSPE=19 IPREP=1 NTALLY=1 SLO=0 BAA=1 ELEV=55
OCURHT(IHAB=3, sp1..19) = 0 1 1 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0   (occupancy: sp 2,3,7,10)
OCURNF(IFO=11, sp1..19) = 0 0 1 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0   (forest occ: sp 3,7,8,9)
→ 118 TPA established
```
KEY FACTS for the port: (1) NOFSPE=19 = EM's species count ⇒ the tally + OCURHT/OCURNF are over EM's 19 species
(NOT a fixed 23-set) ⇒ OCURHT/OCURNF are VARIANT-SPECIFIC static occupancy tables (corrects the earlier "occupancy
shared" note — the ESTOCK STOCKING regression SHAB/SSER is shared/IHAB-indexed, but the SPECIES-level occupancy is
per-variant). (2) The index chain is ITYPE(4, jl already has via NIHMAP in p.habitat_input) → IHAB(3, via the estb
IPHAB per-plot map) → ISER=MYHABG(IHAB)=1 (MYHABG shared = [4*1,4*2,3,4,6*5], MYHABG[3]=1 ✓). ⇒ the one index jl
still needs is ITYPE→IHAB (the IPHAB derivation) — for ITYPE=4→IHAB=3. REMAINING PORT PIECES for EM AUTOES:
(a) the ITYPE→IHAB map (estb IPHAB), (b) EM's OCURHT(16,19)+OCURNF(6,19) occupancy tables (extract from estb/ —
source location TBD, not a simple DATA in estab.f/esnutr.f/esblkd.f), (c) generalize ie_autoes_* to nspecies +
these EM tables, (d) validate total-TPA vs em_live.txt (118 for this stand). Oracle FVSem_clean; harness em_sub.db.

### EM AUTOES — COMPLETE index derivation (2026-08-06), port fully specified as a near-drop-in
Final trace of the IHAB derivation (the last unknown). For a DB stand (no PLOTINFO), esplt2.f takes the STAND-VALUES
path: `IHTYPE = ICL5` → the SHARED MYGRUP bracket (esplt2.f:46-53, `DO 3: IF(IHTYPE>IEND(I)) skip; IHTYPE=MYGRUP(I)`,
IEND/MYGRUP = jl's _IE_ESTAB_IEND/MYGRUP) → `IPHAB=IHTYPE` (line 266). ICL5 = the FVS habitat code set by em/habtyp
(the NI 3-digit habitat code). For the repro stand ICL5≤269 → MYGRUP[1]=3 → IHAB=3 (matches the ground truth).
⇒ THE EM AUTOES PORT IS FULLY SPECIFIED — a near-drop-in of jl's ie_autoes machinery (shared ESTOCK/CHAB/OCURHT +
shared MYGRUP bracket + fixed 10-species estb set):
  1. jl needs EM's ICL5 (the NI habitat code) — em/habtyp sets it; verify jl computes/stores it (jl has
     habitat_code=IEMTYP + habitat_input=ITYPE; ICL5 is the 3rd value — add it to em site_setup if absent).
  2. em_estab_indices: ICL5 → IHAB via the SHARED bracket (reuse _IE_ESTAB_IEND/MYGRUP); ISER=MYHABG[IHAB]; forest→IFO.
  3. generalize ie_autoes_run/tally/establish! to take nsp (=nspecies, 19 for EM) + drop the InlandEmpire gate.
  4. validate total-TPA vs em_live.txt (118 for stand 5332701010661).
Estimated ~30-50 line coding chunk (the machinery is all shared/ported). Ground truth + oracle in hand
(FVSem_clean, em_sub.db). CAVEAT: inherits IE AUTOES v1 maturity (~22% residual). This closes the AUTOES SCOPING;
the implementation is the next chunk.

### Cluster habitat-type DEFAULT audit (2026-08-06) — no-PV_CODE real-FIA stands defaulted to the wrong ITYPE
While tracing the EM AUTOES habitat chain, found jl mis-defaulted the habitat-type ITYPE for stands with NO PV_CODE
in the DB (grassland/regen-heavy real-FIA — ~182/200 EM stands, many IE stands). Each variant's live default = its
grinit ITYPE:
- **EM** grinit ITYPE=29 → default (IEMTYP=29, ITYPE=NIHMAP[29]=4). jl was (1,1). FIXED 633a114, INSTRUMENT-VERIFIED
  vs live (stand 5332701010661, PV_CODE NULL: live ITYPE=4).
- **IE** grinit ITYPE=4 (ie/grinit.f:201; ie_habtyp keeps-prior on unmatched). jl fell back to 1. FIXED 8c86bc6,
  source-faithful + iet01-validated (byte-unchanged), pending live-FVSie confirmation.
- **UT/TT** grinit ITYPE=0 — jl already defaults to 0. No fix needed.
IMPACT: corrects the ITYPE-keyed DG (DGCONS MAPHAB/MAPCCF) + ITYPE-keyed mortality + the AUTOES indices on every
no-PV_CODE stand. emt01/iet01 (keyword habitats) byte-unchanged. The 4 EM treed no-PV_CODE stands shift +1-4% TPA =
the known #137 under-thin exposed (the wrong itype=1 was a compensating error; doctrine #4 masked-bug signal, keep).
This is essential groundwork for EM AUTOES (#143): the ~182 AUTOES-0 stands now carry the correct IHAB (=3 for the
repro) for when the tally is wired.

### AUTOES tally-AMOUNT — two distinct paths (2026-08-06 measured): disturbance=PLPROB, ingrowth=per-tree-TPA
Instrument-replay on the EM repro (treeless→INGROWTH tally, FVSem_trc estab.f:7009): ITPP=1-2 trees/plot, PLPROB=0.0,
DUPNPT=50, INGRO=1. ⇒ the AUTOES over-production has TWO distinct amount paths:
- **Disturbance tally** (stand had overstory removed): total → PLPROB (Σ per-species probs, estab.f:313/541); the
  per-tree booking prob1·300/dupnpt × capped-ITPP diverges from PLPROB = IE's ~22% stand4 residual. FIX: book
  tpaw=PLPROB/ITPP (total=PLPROB). [The earlier "book to PLPROB" recipe applies HERE.]
- **Ingrowth tally** (bare stand, itrn==0&icyc==1 — the ~182 EM real-FIA stands): PLPROB=0, so the total is NOT
  PLPROB. It = Σ_plots ITPP · per-tree-TPA. jl over-produces ~2.6× (repro 303 vs live 118) ⇒ jl's INGROWTH per-tree
  TPA (esprob·300/dupnpt with esprob=p1·newtpp/itpp) is ~2.6× live's. The live ingrowth per-tree TPA formula (a
  different scaling — possibly the SHORTY/ingrowth-time path, or a smaller constant than 300) needs tracing in
  esaddt.f / the ingrowth booking. ⇒ EM's dominant residual is the INGROWTH per-tree TPA, DISTINCT from IE's
  disturbance PLPROB residual. Both are the #143 tally-amount close-out but need SEPARATE fixes. Structural gap
  (empty stands) already RESOLVED (e1dd44e).

### AUTOES ingrowth over-production PINPOINTED (2026-08-06): ITPP-per-plot RNG realization, NOT the per-tree formula
Traced the EM ingrowth over-production (303 vs 118) to its exact source. The per-tree TPA MACHINERY is CORRECT
(matches live estab.f): ingrowth esprob = prob1·FTEMP2, FTEMP2=NEWTPP/ITPP (estab.f:945/951) = jl's
p1·newtpp/itpp; the tree TPA = ESPROB·300/DUPNPT (estab.f:1459) = jl's esprob·300/dupnpt. ⇒ the over-production
is NOT the formula. It is the ITPP (trees-per-plot) COUNT: instrumented per-plot ITPP jl=[1,1,3,1,2,…] vs live
=[1,2,1,…] — jl hits the cap (MAXING[3]=3) more often. ITPP=INT(ESTPP(DRAW)+0.5), so the divergence is the
ESTPP RNG realization (the ESRANN draw sequence entering EM's ingrowth tally) — jl's chain differs from live for
EM. prob1=0.4996 (jl); prob1 does NOT feed ITPP (ESTPP takes draw+ihab+aspect+regt only), so the count difference
is purely the RNG/ESTPP path. ⇒ THE EM AUTOES AMOUNT RESIDUAL = the ie_esrann/ie_estpp realization for EM ingrowth
(match the ESRANN sequence — the IE chain was validated bit-exact on iet01, so the EM setup enters the tally at a
different RNG state; trace the pre-tally draw count for the ingrowth path vs live). If the systematic 2.6× survives
RNG-matching it may be an accepted stochastic straddle (doctrine #3), but the consistent over-cap suggests a real
pre-tally RNG-state divergence. Structural gap RESOLVED (e1dd44e); this is the amount close-out's precise root.

### AUTOES ingrowth ITPP root — body-count nsp-scaling REFUTED (2026-08-06, tested)
Hypothesis: the per-plot RNG body (135 draws in ie_autoes_plot_seeds) + the 69 tally draws are species-dependent
(3·NOFSPE), so EM (19 sp) desyncs the seed chain → wrong ITPP. TESTED: scaled body=66+3·nsp (123 for EM) and the
tally draws 69→3·nsp (57 for EM). Result — EM repro 303→315 (essentially UNCHANGED, still ~2.6× over), IE
byte-IDENTICAL (iet01 unchanged). ⇒ the body/draw-count is NOT the ITPP root (REVERTED — inferred + unvalidated +
no effect). The ITPP-per-plot divergence (jl [1,1,3,1,2] vs live [1,2,1]) is therefore in the SEED / pre-ITPP RNG
state: the plot seed0 = ESRANN(es_stream=55329) or the pre-ITPP draws (wk6fill=50 plot-1 fill, EMSQR=2) differ from
live for EM ingrowth. NEXT: instrument live estab.f ESRANN draw-by-draw for the EM ingrowth tally (the DRAW feeding
ESTPP per plot) vs jl's ie_esrann chain — the IE chain was validated bit-exact on iet01 (disturbance path), so the
EM INGROWTH path (NTALLY=99, itrn==0&icyc==1) enters the tally at a different RNG state (possibly the ingrowth
seed0 derivation or the SHORTY/ESTIME setup draws differ). Structural gap RESOLVED (e1dd44e); this is the amount
close-out's remaining root, now narrowed to the seed/pre-ITPP RNG realization (body-count excluded).

### AUTOES ingrowth ITPP root FOUND (2026-08-06): per-plot RNG advance too large for DUPLICATE plots
Instrument-replay comparison of the ESTPP DRAW per plot (live estab.f:678 vs jl ie_estpp), EM repro (nptids=1,
idup=50 ⇒ 50 DUPLICATE plots):
```
plot:    1         2         3         4         5         6
live:  0.346302  0.350846  0.212445  0.215864  0.991703  0.932161   → ITPP [1,2,1,1,3,3]
jl:    0.346302  0.184747  0.835307  0.243140  0.462182  0.230992   → ITPP [1,1,3,1,2,…]
```
★ Plot-1 DRAW is BIT-IDENTICAL (0.346302) ⇒ seed0 + the pre-ESTPP draws (wk6fill=50 + EMSQR=2) are CORRECT. Plot-2+
diverge. ★ KEY: live's consecutive draws are CLOSE (0.346/0.351, 0.212/0.216, 0.992/0.932 — a SMALL RNG advance
between duplicate plots), while jl's are FAR APART (a full body=135 advance). ⇒ jl's ie_autoes_plot_seeds applies
the full 135-draw per-plot body between ALL plots, but for DUPLICATE plots (idup>1 from a single inventory point)
live advances the RNG by much LESS. iet01 stand4 (multiple REAL plots) validated the body=135 bit-exact, MASKING
this — the duplicate-plot advance is a different (smaller) count. ⇒ THE EM AUTOES INGROWTH AMOUNT RESIDUAL ROOT =
the DUPLICATE-plot RNG advance in ie_autoes_plot_seeds (body=135 is right for distinct plots, wrong for idup
duplicates). NEXT: instrument the ESRANN call-count between consecutive ESTPP draws in live (esrann.f counter) to
get the exact duplicate-plot advance, then make ie_autoes_plot_seeds use it for idup>1. This also affects IE
real-FIA stands with idup>1 (the ~22% IE residual may share this root). Structural gap RESOLVED (e1dd44e); amount
root now PRECISELY located (duplicate-plot RNG advance), body-count-nsp REFUTED, seed0/pre-ESTPP CONFIRMED correct.

### AUTOES ingrowth amount BIT-EXACT (2026-08-06): 303→118=live, TWO measured root causes
Full instrument-replay chain (live estab.f, EM stand 5332701010661, DB elevation NULL):
1. **Per-plot RNG body** was hardcoded 135 (=16+69+50, IE-iet01-specific). Measured live per-plot ESRANN advance
   (esrann.f call-counter): EM = 83 (constant, 137−54−…), IE iet01 = 135. Formula body = 16 + 3·nsp + 2·MAXTPP[ihab]
   fits BOTH (IE 16+69+50=135, EM 16+57+10=83). Fixed ⇒ jl ITPP [1,2,1,1,3,3]=live bit-exact. This ALONE did not
   fix the amount (still 321) — the ITPP was right but the per-tree magnitude was 2.7× high.
2. **PROB1 = ESTOCK logistic used elev=0 not 55.** Live ESPROB = PROB1·NEWTPP/ITPP, PROB1 = 1/(1+exp(−(PN+ESB−
   ESB1)))·STOADJ. Measured live: PN=−1.48962, ESB=ESB1=0, STOADJ=1 → PROB1=0.18398. jl: PN=−0.00143 → PROB1=0.4996.
   ratio 0.4996/0.18398 = 2.716 = 321/118 EXACTLY. The sole divergent ESTOCK input: elev (live 55 vs jl 0). The
   −0.027058·ELEV term = −1.488 = the entire PN gap. ROOT: em/grinit.f:190 `ELEV=55.` (variant default, hundreds-ft);
   DB overrides only when >0 (dbsstandin.f:647); this stand's ELEVATION/ELEVFT are NULL ⇒ live keeps 55. jl's
   forest_location elev fallback was SOUTHERN-GATED ⇒ western NULL-elev stands used 0. Fixed in fia_database.jl with
   the per-variant grinit ELEV default (EM 55, BM 45, IE 38, KT 35, UT 83, TT 65, CI 50). VALIDATED: 4 NULL-elev EM
   stands all jl=118=live bit-exact. NO regression: iet01 IE 536 (body inert nsp=23/ihab=10; .key skips DB reader),
   emt01 EM 1990 bit-exact. SCOPE: the 54/200 NULL-elev EM DB stands are the EMPTY establishment stands; treed growth
   stands all HAVE elevation ⇒ the cornered ~7%-BA growth-tail (treed) is UNAFFECTED (genuinely a separate phenomenon,
   verdict stands). ⇒ EM AUTOES amount [#143] CLOSED. The elev-default fix is a cluster-wide correctness fix (ESTOCK
   AND large-tree DG_EL·elev+DG_EL2·elev²) for any western NULL-elevation DB stand.

## 2026-08-06 — #137 re-confirmed SETTLED on emt01 (bit-exact early + cornered compounding tail)
Re-measured the canonical EM establishment stand emt01 (ESTAB 1992 + PLANT sp2/sp10 400 each) jl vs live FVSem
(emt01.sum), full 11-cycle TPA+BA (correct cols: BA=$4, SDI=$5):
  1990-2010: TPA + BA **BIT-EXACT** (536/77, 526/96, 517/114≈115).
  2020→2090: compounding tail — jl BA +3→+16 (246 vs live 230 = +7% by 2090), jl TPA -0→-17 (394 vs 411 = -4%).
The late TPA over-kill (-4%) is DOWNSTREAM of the +7% BA over-growth (higher BA → higher SDI → more self-thin), NOT a
separate cyc0 mortality bug — cyc0 (1990) is bit-exact. This is the accepted RDPSRT/AVHT40 BA-percentile/crown tie-break
precision compounding (the goal file's own "EM/IE growth-only ~7%-BA-by-2090 tail is CORNERED, no fix warranted").
⇒ #137 is bit-exact-or-cornered on emt01, consistent with the c7c7d2f/#152 fix (the real #137 bug was a HEIGHT-GROWTH
esgent birth-cycle omission, not mortality) and the dense-seedling real-FIA campaign ([[fvsjl-dense-seedling-sweep-
campaign]], #145 fixed). The goal-file's "#137 estab self-thin over-kill still open" is STALE — superseded by those
fixes. No mortality bug remains on the EM establishment path; the residual is the cluster-wide cornered DG-precision tail.
NOTE (breadth check, honoring "don't narrow to one variant"): this session's #140 PVREF6 fix is BM-only (fia_database.jl
BlueMountains branch) — EM emt01 is byte-unchanged by it (verified: EM reads PV_CODE via the numeric mod-1000 path, no
PVREF6), so this cornered tail is the pre-existing state, not a regression.

## 2026-08-07 — ★★ EM FIA sweep: jl AUTOES massively OVER-establishes vs live (+63-86% TPA) — #143, quantified
Fresh multi-cycle EM sweep (em40_sub.db, 3 stands, keyword-less DATABASE) vs FVSem_clean:
  31446929010690: live 1098→902 (monotone mortality) vs jl 1098→1671 (+85%); jl JUMPS at 2028 (1074→1375) + re-
    establishes every ~20yr (2048/2068/2088). 39600883010690: live 444 vs jl 828 (+86%). 42536261010690: live 409
    vs jl 667 (+63%). All jl QMD LOWER (regen dilutes).
★ ISOLATED: jl+NOAUTOES is BIT-EXACT with live (1098→902 all cycles; BA 91 vs live 99 = the accepted growth tail).
  ⇒ the ENTIRE divergence is jl's AUTOES firing where live establishes ~0. jl runs ie_autoes_establish!
  UNCONDITIONALLY for EM/IE (simulate.jl:565, no gate). Live LAUTAL defaults TRUE (esinit.f:52) so live's AUTOES is
  "on" too — but live's TALLY produces ~0 regen for these habitat/site conditions, while jl's over-produces hundreds.
  ⇒ the jl AUTOES TALLY over-produces (the #143 "tally amount / target contamination" issue), here HUGE not 22%.
★ VALIDATION-GAP META: the prior EM/IE FIA re-validation (2026-08-05) was CYC0-focused (bit-exact at cycle 0).
  AUTOES is a LATER-cycle process ⇒ cyc0 validation CANNOT see it. Multi-cycle sweeping surfaced it. LESSON: FIA
  validation must run MULTI-CYCLE (not just cyc0) to exercise establishment/mortality/self-thin, not only initial DG.
NEXT (#143): the jl AUTOES tally (ie_autoes_tally / ie_autoes_run) over-produces for EM habitat conditions where
  live yields ~0. Instrument jl vs live AUTOES tally (count/species/schedule) on 31446929010690 at the first firing
  cycle (2028). Either the tally MODEL over-counts, or a gate/threshold (seed source, habitat, stocking) that
  suppresses live's tally is missing in jl. Reuses the ie_autoes machinery. This is the stated open priority.

## 2026-08-07 — ★★ #143 MECHANISM localized: LINGRW auto-INGROWTH (not LAUTAL disturbance) + "EM AUTOES validated" was NOAUTOES+PLANT
Two findings that sharpen + partly RE-FRAME #143:
1. emt01.key (the stand where "EM AUTOES wired+bit-exact / #143 closed for EM" was claimed) uses NOAUTOES (natural
   auto-establishment OFF) + explicit PLANT (planted regen). ⇒ EM's NATURAL AUTOES/ingrowth tally was NEVER
   validated bit-exact against live — the prior "bit-exact" was the PLANT path (deterministic), not the natural tally.
   The memory "EM AUTOES bit-exact (#143 closed for EM)" is OVERSTATED for the natural path.
2. The +63-86% over-establishment fires via LINGRW automatic INGROWTH, NOT LAUTAL disturbance-tally. establishment.jl
   :1076-1081: idsdat defaults inv_year-20 (esnutr.f:113); the 40-yr-gap ingrowth rule (next_year-idsdat≥40) first
   fires at inv_year+20 = 2028 — EXACTLY the observed jl establishment cycle on 31446929010690 (then ~every 20-40yr).
   These FIA stands have NO thinning, so LAUTAL (disturbance ≥THRES removal) is NOT the trigger — it's the periodic
   ingrowth. jl's ingrowth tally produces hundreds TPA; live produces 0.
⇒ #143 for keyword-less multi-cycle FIA = jl's LINGRW auto-ingrowth tally OVER-produces vs live (~0). NEXT (measure):
  instrument live vs jl at 2028 — does live fire the ingrowth tally at all (LINGRW active on plain DATABASE stands?)
  and if so what count/BAAA/PN, vs jl's. If live's ingrowth yields 0 for these habitat/stocking conditions, jl's
  ESTOCK/BAAA ingrowth model over-produces (the known #143 baaa/tally-amount root); if live never fires ingrowth on
  keyword-less DB stands, jl's unconditional LINGRW firing (simulate.jl:565 + ie_autoes_schedule!) is a gate bug.
  Distinguish via the FVSem_g16 rebuild + an estab-tally dump. This is the stated open priority.

## 2026-08-07 (cont.) — #143 gate-vs-tally: measurement INCONCLUSIVE (wrong instrument target); mechanism source-confirmed
Source-confirmed the LINGRW rule is IDENTICAL jl↔live: esnutr.f:331-338 fires ingrowth when IY(ICYC+1)-IDSDAT≥40,
IDSDAT=IY(1)-20 ⇒ first fire inv_year+20=2028 (= jl). NTALLY=99 signals ingrowth (esnutr.f:337). ⇒ NOT a firing-rule
diff. BUT: esnutr is only reached via an establishment ACTIVITY (keyword IACTK-427, or the LAUTAL after-thinning
path fmcons.f:252) — so whether live INVOKES esnutr (thus the LINGRW branch) on a keyword-less undisturbed DATABASE
stand is the OPEN gate-vs-tally question. MEASUREMENT ATTEMPT was INCONCLUSIVE: instrumented estab.f:544
(NSTORE=ITPP) but fort.89 was EMPTY for BOTH the FIA stand AND emt01 (control, which DOES establish via ESTAB+PLANT)
⇒ estab.f:541/544 "ITPP AT INVENTORY" is NOT on the natural-ingrowth tree-booking path; the empty result is a
wrong-target artifact, NOT evidence live skips establishment. (Control emt01 correctly invalidated the reading.)
NEXT: find the ACTUAL natural-ingrowth tree-booking point (ESADDT / where ITRN is incremented for NTALLY=99), dump
its count + year on the FIA stand vs live. If live never books ingrowth on keyword-less DB stands ⇒ jl's
unconditional ie_autoes_establish! (simulate.jl:565) is a GATE bug (simple). If live books but ~0 ⇒ ESTOCK/ITPP
tally-amount over-produce (known baaa root). Still the priority; the mechanism (LINGRW inv_year+20) is source-solid.

## 2026-08-07 (cont.2) — ★★ #143 gate-vs-tally SETTLED = TALLY AMOUNT (live fires ingrowth same cycles, ITPP≈0)
Correct measurement (instrumented esnutr.f:313 LINGRW-branch inputs, FVSem_g16, FIA stand 31446929010690):
  LGCHK ICYC/IYnext/IDSDAT/ITRN — cyc1 2018/1988/24; cyc2 2028/1988/72; cyc3 2038/2008; cyc4 2048/2008; cyc5
  2058/2028; ... IDSDAT RESETS by 20 each time IYnext-IDSDAT≥40 (esnutr.f:338 IDSDAT=IYnext-20 on the NTALLY=99
  fire). The reset pattern 1988→2008→2028→2048→2068 PROVES live's LINGRW ingrowth FIRES at cyc 2/4/6/8/10 =
  2028/2048/2068/2088/2108 — EXACTLY jl's establishment cycles. So it is NOT a gate/firing bug (CALL ESNUTR is
  ungated, gradd.f:229; the rule is identical). Yet live TPA is MONOTONE-decreasing (0 net regen) ⇒ live's ingrowth
  ESTOCK/ITPP tally yields ~0 established trees where jl yields hundreds. ⇒ #143 = TALLY-AMOUNT over-production,
  DEFINITIVELY (the known baaa/ESTOCK/PLPROB root, establishment.jl:1104-1113 jl per-point BAAA ~1.5× live).
  (Caveat corrected: an earlier LINGRWFIRE dump at NTALLY=99 read empty — a whitespace-match miss on the insertion,
  NOT non-firing; the LGCHK IDSDAT-reset trace is the reliable evidence. Used the reset invariant, not the raw dump.)
NEXT: instrument the ingrowth ITPP (estab.f:589, the NTALLY=99/ingrowth path — NOT :541 which is inventory) +
PLPROB(NNID) + the ESTOCK PN vs jl on 31446929010690 cyc2(2028) — jl's PLPROB/ITPP is the over-producing term.
The fix is the ESTOCK stocking-probability / BAAA model (jl over-estimates P(stocking) for these EM habitats).

## 2026-08-07 (cont.3) — #143 jl-side quantified: ingrowth 323 TPA@2028; BAAA is RIGHT here → root is ESTOCK-PN/habitat
jl FVSJL_AUTOES_DEBUG on 31446929010690: ingrowth tally (ntally=99) total=323.5 TPA@icyc2(2028), then 163/164/175/
190 at icyc4/6/8/10. baaa=28.24@2028 ≈ the stand BA (jl BA 34) ⇒ on THIS stand the BAAA input is ~CORRECT (NOT the
documented 1.5×-low case of stand 753199439290487). So the +hundreds over-production here is NOT the BAAA input — it
is the ESTOCK stocking-probability PN itself. PN≈-1.513+1.245·ln(AGE)-0.003052·baa (establishment.jl:76): TIME=1
(ingrowth SHORTY, live-matching) ⇒ ln(AGE)≈0, baa=28 ⇒ PN≈-1.6 ⇒ logistic≈0.17 (~17% stocking/species) × ~10 estb
species × plot-expansion ≈ 323. For live's ~0, live's effective PN must be near-zero ⇒ HYPOTHESIS: EM dry-habitat
species-gating (OCURHT zeros the wet-side estb species for EM's habitats, establishment.jl:1067-1068) is not zeroing
jl's establishing species — jl establishes species live's habitat forbids. NEXT (safe jl-side first): dump jl per-
species r.tally[sp] + the ihab_code for this stand; check which species jl establishes vs the EM habitat's allowed
set (OCURHT/ie_estab_indices). Then confirm vs live (correct estab.f ingrowth-tally-count target). The BAAA-input
root (#143 title) is real on OTHER stands but NOT the driver here — the driver is the ESTOCK-PN/habitat-species set.

## 2026-08-07 (cont.4) — #143 jl per-species: over-production is SUBALPINE AF/ES/LP at ihab=870
jl AUTOES per-species tally @2028 (31446929010690, ihab=870 fcode=102): sp7=LP 15.6, sp8=ES 54.6, sp9=AF 253.4
(total 323.5). ⇒ the over-production is dominated by AF (subalpine fir, 253) + ES + LP — the SUBALPINE estb species.
So the fix question narrows to: at EM habitat 870, does live's ESTOCK establish AF/ES/LP at all (and if so how much)?
Two candidate roots: (a) SPECIES-SET — live's habitat-870 OCURHT zeros/limits AF/ES/LP but jl establishes them; or
(b) AMOUNT — species set matches but jl's ESTOCK PN over-tallies AF (~17%/species stocking → 253 TPA). DECISIVE
NEXT = live per-species ingrowth tally at 2028 on this stand (correct estab.f ingrowth-count target, e.g. where the
NTALLY=99 path books trees per species — trace from esnutr GOTO 200 → the ESTAB tally loop). If live AF≈0 → species-
gating bug (cleaner fix); if live AF>0 but «253 → ESTOCK-PN amount. BAAA confirmed ~correct on this stand (28.24).

## 2026-08-07 (cont.5) — #143 species-set-vs-amount RESOLVED = AMOUNT (subalpine species appropriate for the WB stand)
Stand 31446929010690 = PURE WHITEBARK PINE (FIA 101, 1170 TPA, avgDBH 6.0), subalpine, habitat 870, forest 2,
slope 50 aspect 110. ⇒ jl's establishing species AF(9)/ES(8)/LP(7) ARE appropriate subalpine associates of whitebark
pine — NOT a wrong species-SET. ⇒ #143 here is the AMOUNT: jl's ESTOCK stocking-probability PN over-tallies
(AF 253 TPA @2028) where live yields ~0 for this harsh subalpine WB habitat. Resolved via a DB-composition query (no
error-prone live-Fortran instrumentation). So the fix is jl's ESTOCK PN calibration over-predicting stocking for
subalpine EM habitats — jl PN≈-1.6→logistic 0.17 per species is too high; live's effective stocking prob for AF/ES/LP
at habitat 870 must be near-zero. CANDIDATE TERMS (establishment.jl:74-108 + FTEMP): the uhab(sp) habitat-adjustment
for habitat 870, and/or the ESB inventory-calibration (which the audit notes applies only at inv-year, ESB=0 for
ingrowth — so ingrowth has NO stocking down-calibration in jl; if live's ingrowth PN is effectively suppressed for
subalpine habitats, jl misses that suppression). NEXT: compare jl uhab/PN for AF@ihab-870 vs the live ESTOCK PN
(estab.f ESTOCK) — a targeted per-species-per-habitat coefficient check, jl-side first (read the EM uhab table for
habitat 870). This is the known-hard ESTOCK-model core, now bounded to AF/ES/LP @ subalpine habitat 870, AMOUNT.

## 2026-08-07 (cont.6) — #143 habitat-mapping RULED OUT (jl matches live); bound tightened to ESTOCK-PN/OCURHT@ihab16
Checked jl's habitat→ESTOCK-class mapping for the WB stand (habitat 870): jl ie_estab_indices(870,2) → ihab=16
(fallback, because 870 > max IEND 799 in the esplt2 bracket). LIVE esplt2.f:46-53 does the IDENTICAL thing —
IHTYPE=ICL5; if IHTYPE>IEND(all 33, max 799) → IHTYPE=16 (fallback, line 52). EM's live esplt2 IEND table is
byte-identical to IE's (both max 799). ⇒ jl's 870→16 MATCHES live. NOT the bug. (Note: 14/118 EM habitats, NI codes
810-999, all fall to fallback-16 in BOTH jl and live — a shared design, not a jl gap.)
★ #143 CANDIDATES RULED OUT via SAFE checks this session: (1) gate/firing (both fire same cycles, IDSDAT-reset);
(2) species-SET (AF/ES/LP appropriate for the subalpine WB stand); (3) BAAA input (jl 28.24 ≈ stand BA, correct
here); (4) habitat MAPPING (both →ihab16). ⇒ REMAINING = the ESTOCK stocking-probability PN coefficients / OCURHT
species-availability for ihab=16 (the fallback class), OR the ITPP normalization / ESRANN draw. jl produces AF 253
where live ~0 at the SAME ihab16/BAAA28/TIME1 — so a jl coefficient (uhab/OCURHT for ihab16) or the ITPP/RNG path
differs. This is the known-hard ESTOCK core; the DECISIVE step is the live per-species ingrowth tally (estab.f
NTALLY=99 booking) vs jl's r.tally for AF/ES/LP@ihab16 — needs the correct estab.f target (fresh-session task).
META: 4 candidates eliminated by DB-query + jl-eval + source-compare (all safe, no Fortran instrumentation), each a
clean negative — the diagnosis is now tightly bounded to the ESTOCK PN/OCURHT tables for the fallback habitat class.

## 2026-08-07 (cont.7) — #143 ALL ESTOCK INPUTS verified == live; bug isolated to PN-equation/ITPP or tree-persistence
Confirmed jl ELEV=87.0 (reads ELEVFT 8700×0.01; the missing-elev bug was already fixed, fia_database.jl:81-82 — this
stand has ELEVFT present so elev is correct, NOT 0). Confirmed live esplt2.f:191/267 sets IPPREP=1 & IPHYS=3 for
keyword-less DB stands (no PLOTINFO) = jl iprep=1/iphy=3. ⇒ the COMPLETE ESTOCK input set is byte-identical jl↔live:
ihab=16, iphy=3, iprep=1, ELEV=87, BAAA=28.24, species AF/ES/LP, TIME=1 (ingrowth), same firing cycles, same ESRANN
seeds (43303/61677, "ESAVE chain bit-exact"). YET jl books AF=253 TPA@2028 where live nets ~0 (.sum monotone).
★ 6 CANDIDATES RULED OUT by SAFE checks (no live-Fortran instrumentation): gate/firing, species-set, BAAA, habitat-
mapping(→16), elevation(=87), iphy/iprep(=3/1). ⇒ the divergence is NOT the inputs — it is the ESTOCK PN-EQUATION
COMPUTATION (base coeffs / uphy·uhab·upre tables for the fallback class) OR the ITPP normalization / per-species
draw, OR live-side tree-PERSISTENCE (live may book then immediately lose them). This is the tightest bound safe
checks allow. DECISIVE NEXT (needs live instrumentation): dump live's per-species ESTOCK PN + PROB1 + ITPP + the
booked tree count at the NTALLY=99 ingrowth tally (2028) vs jl's ie_estock PN / r.tally for AF/ES/LP — one clean
per-species comparison at identical inputs will expose the exact diverging term (a coefficient, the ITPP formula, or
that live's tally genuinely yields ~0 while jl's ITPP rounds up). META: 6 clean negatives via DB-query/jl-eval/
source-compare — a model of eliminating candidates cheaply before the expensive live measurement.

## 2026-08-07 (cont.8) — ★★★ #143 DECISIVE: live AF ingrowth = 0.1 TPA vs jl 253 TPA (~2500× over) = STRUCTURAL tally bug
Instrumented live estab.f tree-booking (FVSem_g16, stand 31446929010690). Live books via path p2 (estab.f:1293, the
DO 226 II=1,IBRKUP breakup loop): 80 AF(9) records + 6 LP(7) + 3 ES(8). BUT the per-record PROB = (FTEMP2*300*
XCSMAX)/DUPNPT sums to: AF total TPA=0.1, ES=0.0, LP=0.0 — essentially ZERO. jl's r.tally[AF]=253.4 TPA@2028. ⇒
live establishes AF at 0.1 TPA where jl establishes 253 TPA = ~2500× OVER-production. Two-branch question SETTLED =
branch (a) TALLY-AMOUNT (NOT persistence): live's ESTOCK stocking→TPA yields ~0.1 TPA; jl's yields 253. The ~2500×
magnitude ⇒ a STRUCTURAL error (a missing/wrong factor in jl's stocking→TPA), not a subtle coefficient. All ESTOCK
INPUTS proven identical (ihab16/iphy3/iprep1/elev87/baaa28/species/TIME1/RNG-seeds), so the divergence is jl's
PROB/ITPP→TPA COMPUTATION vs live's PROB=(FTEMP2·300·XCSMAX)/DUPNPT. CANDIDATE: jl may miss the XCSMAX (max-crown-
area normalization) and/or DUPNPT (plot-expansion) divisor, OR FTEMP2 (the realized stocking) — jl's PN→logistic≈0.17
but live's realized FTEMP2 here is ~0.0004 (0.1/80/~3). NEXT (jl-side, safe): trace jl's r.tally computation
(ie_autoes_run/ie_esnspe/the ITPP→TPA step) and compare its per-record TPA formula to live PROB=(FTEMP2·300·XCSMAX)/
DUPNPT — find the missing structural factor. This is the FIX target, now concrete (jl over-scales the tally ~2500×
for this habitat). Note iet01 stand-4 was "bit-exact" so jl's tally is right for SOME stands — the over-scale is
habitat/condition-specific (subalpine ihab16 fallback). META: pushing the g16 instrumentation (right target via the
booking-count intermediate) gave the decisive 0.1-vs-253 number that no amount of input-checking could.
