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

## Chunk plan (mirror KT)
1. Species: `data/easternmontana/species_coefficients.csv` (19 species × ~40 cols from em/*.f DATA:
   bark bratio.f, site sitset.f, small-tree/regent, htcalc/htdbh, morts, sdimax, volume specs) +
   `species_translation.csv` (SPCTRN western crosswalk, EM target column) + `species.jl` blkdat init.
2. Site/habitat (sitset + em_habtyp, 8 habitat-type groups). 3. DG (em/dgf.f — reuse KT
   `diameter_growth!`, swap coefs + DGHAB(8)/DGFOR(6) dims). 4. Height (em/htgf.f). 5. Crown/CCF
   (em/crown.f/ccfcal.f). 6. Regent (em/regent.f). 7. Mortality (em/morts.f + varmrt.f — verify vs
   shared driver + KT Hamilton; the low-thinning behavior). 8. Volume (shared NVEL/cr_fw2_vol).
   9. Full-cycle differential vs FVSem_clean on emt01 + native stands.
