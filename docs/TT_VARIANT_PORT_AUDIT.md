# TT (Teton) Variant Port — Audit Log

FVS **Teton** (VARACD `TT`, MAXSP 18, Region-4 Intermountain, 10-yr cycle). The next western
Rockies cluster variant after CR/KT/IE/EM. Bit-exact-or-cornered vs live `FVStt`, chunk by chunk.
Doctrine (unchanged from CR/EM): validate vs LIVE per chunk; MEASURE don't infer; per-record treelist
diff INVALID after tripling; port faithfully then validate; reuse the shared engine (TT-gate all shared
changes). Oracle + Fortran source are the SOLE ground truth (no FVSjulia oracle).

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

## Chunk 1 — Infra/species (blkdat, FIA map, grinit)  ⬜ NEXT
## Chunk 2 — Site/habitat (tt/sitset.f + tt/habtyp.f)  ⬜
## Chunk 3 — Large-tree DG (tt/dgf.f = EM clone + CR-Zeide density dispatch)  ⬜
## Chunk 4 — Height (tt/htgf.f)  ⬜   ## Chunk 5 — Crown (tt/crown+ccfcal)  ⬜
## Chunk 6 — Small-tree (tt/regent.f)  ⬜   ## Chunk 7 — Mortality (tt/morts+varmrt)  ⬜
## Chunk 8 — Volume (Region-4 FW2, likely `I00FW2W<FIAJSP>`)  ⬜   ## Chunk 9 — Full-cycle diff  ⬜

## Off-switch
`touch docs/TT_VARIANT_PORT_COMPLETE` (USER's call).
