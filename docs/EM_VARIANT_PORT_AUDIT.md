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
