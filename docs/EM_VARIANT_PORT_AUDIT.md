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

## Chunk plan (mirror KT)
1. Species: `data/easternmontana/species_coefficients.csv` (19 species × ~40 cols from em/*.f DATA:
   bark bratio.f, site sitset.f, small-tree/regent, htcalc/htdbh, morts, sdimax, volume specs) +
   `species_translation.csv` (SPCTRN western crosswalk, EM target column) + `species.jl` blkdat init.
2. Site/habitat (sitset + em_habtyp, 8 habitat-type groups). 3. DG (em/dgf.f — reuse KT
   `diameter_growth!`, swap coefs + DGHAB(8)/DGFOR(6) dims). 4. Height (em/htgf.f). 5. Crown/CCF
   (em/crown.f/ccfcal.f). 6. Regent (em/regent.f). 7. Mortality (em/morts.f + varmrt.f — verify vs
   shared driver + KT Hamilton; the low-thinning behavior). 8. Volume (shared NVEL/cr_fw2_vol).
   9. Full-cycle differential vs FVSem_clean on emt01 + native stands.
