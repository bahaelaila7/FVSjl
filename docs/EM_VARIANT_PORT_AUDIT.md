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

## Chunk plan (mirror KT)
1. Species: `data/easternmontana/species_coefficients.csv` (19 species × ~40 cols from em/*.f DATA:
   bark bratio.f, site sitset.f, small-tree/regent, htcalc/htdbh, morts, sdimax, volume specs) +
   `species_translation.csv` (SPCTRN western crosswalk, EM target column) + `species.jl` blkdat init.
2. Site/habitat (sitset + em_habtyp, 8 habitat-type groups). 3. DG (em/dgf.f — reuse KT
   `diameter_growth!`, swap coefs + DGHAB(8)/DGFOR(6) dims). 4. Height (em/htgf.f). 5. Crown/CCF
   (em/crown.f/ccfcal.f). 6. Regent (em/regent.f). 7. Mortality (em/morts.f + varmrt.f — verify vs
   shared driver + KT Hamilton; the low-thinning behavior). 8. Volume (shared NVEL/cr_fw2_vol).
   9. Full-cycle differential vs FVSem_clean on emt01 + native stands.
