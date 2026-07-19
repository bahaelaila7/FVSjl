# ACTIVE GOAL — Port the CR (Central Rockies) FVS variant to FVSjl

## Mission
Port the FVS **Central Rockies (CR)** variant into FVSjl, bit-exact-or-cornered vs live `FVScr`, chunk by chunk.
CR is the FIRST *western* variant (the 4 done — SN/NE/CS/LS — are all eastern). It is the recommended hub of the
western Rockies cluster (KT/IE/EM/BM/TT/UT follow at a discount). Branch: `cr-variant-port`.

## What is SHARED (already built — do NOT re-port)
The variant-agnostic engine + infra carries over verbatim: grow-cycle, **tripling/record-splitting**, the
mortality driver + **self-thinning RDPSRT** tie-break, summary/`.sum`, keyword parsing, the **RNG (seed 55329 —
CR uses the SAME)**, FFE, AND the **NVEL volume driver (NATCRS/VOLS) incl. the r9clark extreme-height fix**. The
whole validation harness + methodology is reused. CR infra: `VARACD='CR'`, **MAXSP=38**, **YR=10** (10-yr cycle,
like NE/CS/LS), seed 55329.

## What is NEW (the port surface — cr/*.f, ~35 routines)
CR is a **GENGYM** variant (Edminster's SW growth-and-yield), NOT the standard western Wykoff DDS —
CORRECTED by MEASUREMENT (doctrine #2; the DGLD/DGBAL/DGCCF coefficient-name inference was WRONG):
- **dgf.f** computes stand stats (SDI/RELSDI/DSTAG, BAL, PBAL, CR, BAUTBA, SPBA, SSITE) then delegates the DDS
  to **`CALL GEMDG(...)`** (cr/gemdg.f, line 194). `WK2 = DDS + COR + DGCON`.
- **gemdg.f (GENGYM)** — dispatch on **IMODTY** (model type 1-5: SW mixed conifer, SW ponderosa, …) then
  **SELECT CASE(species)** with per-species/per-model-type regression DDS equations (hardcoded coefficients,
  terms in DP/BA/CR/SI/slope/aspect/elev). The gem* routines (gemdg/gemht/gemcr) are the CR-specific model.
  A new dispatch on `diameter_growth!(::CentralRockies)`. This is a LARGE careful transcription-and-diff chunk.
- **htgf.f** — western height growth. **crown.f / cratet.f / ccfcal.f** — western crown + CCF.
- **regent.f** — western small-tree/regen. **sitset.f / habtyp.f** — site index + **HABITAT TYPE groups**
  (a new coefficient dimension the eastern variants lack: OCURHT is (16 habitat-type-groups, MAXSP)).
- **morts.f / varmrt.f** — mortality (check vs the shared driver; likely mostly shared). **bratio.f** bark.
  **forkod.f** forest code. Volume assignment (bfvol/cubrds/cutstk/estump) routes through the shared NVEL driver.

## DOCTRINE (hard-won — carry from the FIA campaign)
1. **Validate against LIVE FVScr, bit-exact per chunk.** Relink the oracle from `bin/FVScr_buildDir/*.o`
   (668 .o present); differential each chunk. There is NO FVSjulia oracle for CR — the live binary + Fortran
   source are the SOLE ground truth.
2. **MEASURE, don't infer.** Read the Fortran, trace the data flow, instrument if unsure; never conclude a model
   from test pass/fail.
3. **Per-record treelist comparison is INVALID after tripling** — compare per-tree DBH ONLY in the pre-split/
   pre-mortality window, or use stand-level `.sum` aggregates. (The trap that cost a retracted "17% bug".)
4. **Port FVS semantics faithfully**, then validate; a regression on a ported-faithful chunk = examine the
   oracle, don't cargo-cult to green.
5. **Reuse the shared engine** — only add CR-specific equations + data. Until a hook has a CR method, dispatching
   on `CentralRockies()` should error loudly (like LS did).
6. Document every chunk verdict in `docs/CR_VARIANT_PORT_AUDIT.md`; data lands in `data/centralrockies/`.

## Chunk plan (coarse; refine as we go)
0. Scaffold: `CentralRockies` singleton + registration (variant.jl), relink oracle, canonical CR keyfile, harness.
1. Infra/species: blkdat constants, MAXSP=38 species list + FIA map, grinit defaults.
2. Site/habitat: sitset + habtyp (habitat-type groups) — a prerequisite for the DG coefficients.
3. Large-tree diameter growth: dgf/dgdriv (the Wykoff DDS) — the core.
4. Height growth: htgf. 5. Crown: crown/cratet/ccfcal. 6. Small-tree: regent. 7. Mortality reconcile: morts/varmrt.
8. Volume assignment (VEQNNC) through the shared NVEL driver. 9. Full-cycle differential vs live on real stands.

## Off-switch
`touch docs/CR_VARIANT_PORT_COMPLETE` (USER's call). Log: `docs/CR_VARIANT_PORT_AUDIT.md`.
